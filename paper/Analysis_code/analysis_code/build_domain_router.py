"""Evaluate the frozen coverage-layer router without cross-domain score mixing."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "domain_router_01"
BOOTSTRAP_REPLICATES = 5000
BOOTSTRAP_SEED = 20260913
CUTS = (0.10, 0.25, 0.50, 0.75, 1.00)

INTERACTION_CELL = "represented_protein_unseen_chemistry"
MMSEQS_CELL = "new_protein_seen_chemistry_with_homologue"
TRANSPORT_CELLS = {
    "chemical_cold": "new_protein_unseen_chemistry_with_homologue__chemical_cold_accounting",
    "protein_cold": "new_protein_unseen_chemistry_with_homologue__protein_cold",
    "double_cold": "new_protein_unseen_chemistry_with_homologue__double_cold",
}
CONTROLS = {
    "protein_position_permutation": ROOT / "interaction_control_protein_position_permutation_01/interaction_conditional_sequence_01",
    "reaction_center_row_permutation": ROOT / "interaction_control_reaction_center_row_permutation_01/interaction_conditional_sequence_01",
}

METRICS = (
    "rr", "mean_positive_rr", "recall_at_1", "recall_at_5", "recall_at_10",
    "hit_at_1", "hit_at_5", "hit_at_10", "ndcg_at_10",
)


def read_json(path: Path):
    return json.loads(path.read_text())


def random_expectations(n: int, k: int) -> dict:
    if not 0 < k <= n:
        raise ValueError("invalid uniform denominator")
    probability = k / n
    rr = 0.0
    for rank in range(1, n - k + 2):
        rr += probability / rank
        if rank < n - k + 1:
            probability *= (n - rank - k + 1) / (n - rank)
    out = {
        "rr": rr,
        "mean_positive_rr": math.fsum(1 / rank for rank in range(1, n + 1)) / n,
        "positive_count": k,
        "covered_positives": k,
        "candidate_count": n,
        "covered_candidates": n,
        "conditional_rr": rr,
    }
    for cutoff in (1, 5, 10):
        top = min(cutoff, n)
        none = 1.0
        for offset in range(top):
            none *= max(0, n - k - offset) / (n - offset)
        out[f"recall_at_{cutoff}"] = top / n
        out[f"hit_at_{cutoff}"] = 1 - none
    discounts = [1 / math.log2(rank + 1) for rank in range(1, min(n, 10) + 1)]
    out["ndcg_at_10"] = (k / n) * math.fsum(discounts) / math.fsum(discounts[: min(k, 10)])
    return out


def rank_panel(scores: np.ndarray, panel: np.ndarray, domain: np.ndarray,
               positives: list[int], tie: np.ndarray) -> dict:
    panel = np.asarray(panel, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    truth = np.asarray(sorted(set(positives)), dtype=int)
    if not len(truth) or not np.all(panel[truth]):
        raise ValueError("positive outside routed panel")
    defined = panel & domain
    if not np.all(np.isfinite(scores[defined])):
        raise ValueError("non-finite routed score")
    selected = np.flatnonzero(defined)
    order = selected[np.lexsort((tie[selected], -scores[selected]))]
    ranks = np.full(len(scores), np.inf)
    ranks[order] = np.arange(1, len(order) + 1)
    target_ranks = ranks[truth]
    best = float(np.min(target_ranks))
    out = {
        "rr": 0.0 if not math.isfinite(best) else 1 / best,
        "mean_positive_rr": float(np.mean(1 / target_ranks)),
        "positive_count": len(truth),
        "covered_positives": int(np.sum(defined[truth])),
        "candidate_count": int(np.sum(panel)),
        "covered_candidates": len(selected),
        "conditional_rr": None if not math.isfinite(best) else 1 / best,
    }
    for cutoff in (1, 5, 10):
        out[f"recall_at_{cutoff}"] = float(np.mean(target_ranks <= cutoff))
        out[f"hit_at_{cutoff}"] = float(np.any(target_ranks <= cutoff))
    observed = target_ranks[target_ranks <= 10]
    dcg = float(np.sum(1 / np.log2(observed + 1)))
    ideal = math.fsum(1 / math.log2(rank + 1) for rank in range(1, min(len(truth), 10) + 1))
    out["ndcg_at_10"] = dcg / ideal
    return out


def seed_for(*parts: str) -> int:
    token = "|".join(parts).encode()
    return (BOOTSTRAP_SEED + int.from_bytes(hashlib.sha256(token).digest()[:8], "little")) % (2**63 - 1)


def grouped_values(rows: list[dict], keys: tuple[str, ...]) -> np.ndarray:
    grouped = defaultdict(lambda: defaultdict(list))
    for row in rows:
        grouped[row["protein_group"]][row["sequence_sha256"]].append(row)
    values = []
    for group in sorted(grouped):
        per_query = []
        for query in sorted(grouped[group]):
            panels = grouped[group][query]
            per_query.append([math.fsum(float(row[key]) for row in panels) / len(panels) for key in keys])
        values.append(np.mean(np.asarray(per_query, dtype=float), axis=0))
    return np.asarray(values, dtype=float)


def paired_summary(rows: list[dict], label: str) -> dict:
    values = grouped_values(rows, ("left_rr", "right_rr"))
    if not len(values):
        return {
            "query_panels": 0, "unique_queries": 0, "protein_groups": 0,
            "left_mrr": None, "right_mrr": None, "difference": None,
            "conditional_protein_bootstrap_95ci": None,
        }
    delta = values[:, 0] - values[:, 1]
    interval = None
    if len(values) >= 2:
        rng = np.random.default_rng(seed_for("paired", label))
        draws = rng.integers(len(values), size=(BOOTSTRAP_REPLICATES, len(values)))
        interval = [float(value) for value in np.quantile(delta[draws].mean(axis=1), (0.025, 0.975))]
    return {
        "query_panels": len(rows),
        "unique_queries": len({row["sequence_sha256"] for row in rows}),
        "protein_groups": len(values),
        "left_mrr": float(np.mean(values[:, 0])),
        "right_mrr": float(np.mean(values[:, 1])),
        "difference": float(np.mean(delta)),
        "conditional_protein_bootstrap_95ci": interval,
    }


def summarize_methods(rows: list[dict]) -> list[dict]:
    summaries = []
    for cell in sorted({row["cell"] for row in rows}):
        for method in sorted({row["method"] for row in rows if row["cell"] == cell}):
            chosen = [row for row in rows if row["cell"] == cell and row["method"] == method]
            values = grouped_values(chosen, METRICS)
            summaries.append({
                "cell": cell, "method": method,
                "query_panels": len(chosen),
                "unique_queries": len({row["sequence_sha256"] for row in chosen}),
                "protein_groups": len(values),
                "positive_instances": sum(row["positive_count"] for row in chosen),
                "candidate_count_min": min(row["candidate_count"] for row in chosen),
                "candidate_count_max": max(row["candidate_count"] for row in chosen),
                "candidate_count_mean": float(np.mean([row["candidate_count"] for row in chosen])),
                "protein_group_macro": {metric: float(np.mean(values[:, index])) for index, metric in enumerate(METRICS)},
                "micro_positive_coverage": sum(row["covered_positives"] for row in chosen) / sum(row["positive_count"] for row in chosen),
                "query_panel_mean_candidate_coverage": float(np.mean([row["covered_candidates"] / row["candidate_count"] for row in chosen])),
            })
    return summaries


def comparison(rows: list[dict], cell: str, left: str, right: str) -> tuple[dict, list[dict]]:
    index = {(row["cell"], row["source_task"], row["block_number"], row["sequence_sha256"], row["method"]): row for row in rows}
    left_rows = [row for row in rows if row["cell"] == cell and row["method"] == left]
    pairs = []
    for a in left_rows:
        key = (cell, a["source_task"], a["block_number"], a["sequence_sha256"], right)
        b = index[key]
        if a["target_reaction_indices"] != b["target_reaction_indices"] or a["candidate_count"] != b["candidate_count"]:
            raise ValueError("paired routed denominator differs")
        pairs.append({
            "cell": cell, "source_task": a["source_task"], "block_number": a["block_number"],
            "sequence_sha256": a["sequence_sha256"], "protein_group": a["protein_group"],
            "left_method": left, "right_method": right,
            "left_rr": a["rr"], "right_rr": b["rr"],
            "candidate_count": a["candidate_count"], "positive_count": a["positive_count"],
            "left_covered_positives": a["covered_positives"], "right_covered_positives": b["covered_positives"],
            "left_covered_candidates": a["covered_candidates"], "right_covered_candidates": b["covered_candidates"],
        })
    label = f"{cell}|{left}|{right}"
    result = {
        "cell": cell, "left_method": left, "right_method": right,
        "paired": paired_summary(pairs, label),
        "denominators": {
            "positive_instances": sum(row["positive_count"] for row in pairs),
            "left_covered_positive_instances": sum(row["left_covered_positives"] for row in pairs),
            "right_covered_positive_instances": sum(row["right_covered_positives"] for row in pairs),
            "mean_candidates": None if not pairs else float(np.mean([row["candidate_count"] for row in pairs])),
            "mean_left_candidate_coverage": None if not pairs else float(np.mean([row["left_covered_candidates"] / row["candidate_count"] for row in pairs])),
            "mean_right_candidate_coverage": None if not pairs else float(np.mean([row["right_covered_candidates"] / row["candidate_count"] for row in pairs])),
        },
    }
    return result, pairs


def control_summaries(records: list[dict]) -> list[dict]:
    output = []
    for control in sorted({row["control"] for row in records}):
        chosen = [row for row in records if row["control"] == control]
        gain_rows = [{**row, "left_rr": row["control_rr"], "right_rr": row["control_global_rr"]} for row in chosen]
        did_rows = [{**row, "left_rr": row["primary_rr"] - row["primary_global_rr"],
                     "right_rr": row["control_rr"] - row["control_global_rr"]} for row in chosen]
        output.append({
            "control": control,
            "control_gain_over_global": paired_summary(gain_rows, f"control_gain|{control}"),
            "primary_minus_control_gain": paired_summary(did_rows, f"control_did|{control}"),
        })
    return output


def confidence_curve(rows: list[dict], cell: str, method: str) -> dict:
    chosen = [row for row in rows if row["cell"] == cell and row["method"] == method]
    ordered = sorted(chosen, key=lambda row: (-row["confidence"], stable(f'{row["block_number"]}|{row["sequence_sha256"]}')))
    points = []
    for fraction in CUTS:
        count = max(1, math.ceil(len(ordered) * fraction))
        selected = ordered[:count]
        values = grouped_values(selected, ("rr", "hit_at_10"))
        points.append({
            "requested_coverage": fraction,
            "selected_query_panels": count,
            "realized_query_panel_coverage": count / len(ordered),
            "unique_queries": len({row["sequence_sha256"] for row in selected}),
            "protein_groups": len(values),
            "confidence_minimum": float(selected[-1]["confidence"]),
            "protein_group_macro_mrr": float(np.mean(values[:, 0])),
            "protein_group_macro_hit_at_10": float(np.mean(values[:, 1])),
            "documented_positive_retrieval_failure_at_10": float(1 - np.mean(values[:, 1])),
        })
    return {
        "cell": cell, "method": method, "confidence_definition": chosen[0]["confidence_definition"],
        "query_panels": len(chosen), "points": points,
        "not_precision_or_probability": True,
    }


def positive_lower(comparison_row: dict) -> bool:
    interval = comparison_row["paired"]["conditional_protein_bootstrap_95ci"]
    return bool(interval is not None and interval[0] > 0)


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    gates = [
        ROOT / "main_baseline_validation_01/INDEPENDENT_METRIC_QC.json",
        ROOT / "interaction_validation_sequence_01/INDEPENDENT_QC.json",
        ROOT / "interaction_paired_validation_sequence_01/INDEPENDENT_QC.json",
        ROOT / "interaction_control_protein_position_permutation_01/interaction_validation_sequence_01/INDEPENDENT_QC.json",
        ROOT / "interaction_control_reaction_center_row_permutation_01/interaction_validation_sequence_01/INDEPENDENT_QC.json",
        ROOT / "ordered_site_validation_01/INDEPENDENT_QC.json",
    ]
    for path in gates:
        if read_json(path)["status"] != "PASS":
            raise ValueError(f"unpassed input gate: {path}")

    edges = read_json(ROOT / "dataset_02/core_edges.json")
    reactions = sorted(read_json(ROOT / "dataset_02/core_reactions.json"))
    blocks = read_json(ROOT / "multiaxis_split_01/outer_inner_blocks.json")
    axes = read_json(ROOT / "multiaxis_split_01/axis_assignments.json")
    sequences = sorted({edge["sequence_sha256"] for edge in edges})
    sequence_index = {sequence: index for index, sequence in enumerate(sequences)}
    reaction_index = {reaction: index for index, reaction in enumerate(reactions)}
    tie = np.empty(len(reactions), dtype=int)
    for rank, index in enumerate(sorted(range(len(reactions)), key=lambda value: stable(reactions[value]))):
        tie[index] = rank
    with np.load(ROOT / "ordered_site_features_01/ordered_features.npz", allow_pickle=False) as saved:
        if list(saved["sequence_ids"]) != sequences:
            raise ValueError("ordered-site sequence order")
        ordered_available = saved["sequence_projected_available"].astype(bool)
    with np.load(ROOT / "main_baselines_01/fresh_representations.npz", allow_pickle=False) as saved:
        if list(saved["sequence_ids"]) != sequences or list(saved["reaction_ids"]) != reactions:
            raise ValueError("baseline representation order")
        homology = saved["homology_bits"].astype(np.float64)

    metric_rows: list[dict] = []
    decision_rows: list[dict] = []
    control_rows: list[dict] = []
    input_paths = set(gates + [
        ROOT / "DOMAIN_ROUTER_V1.md", ROOT / "build_domain_router.py", ROOT / "run_raw.py",
        ROOT / "dataset_02/core_edges.json", ROOT / "dataset_02/core_reactions.json",
        ROOT / "multiaxis_split_01/outer_inner_blocks.json", ROOT / "multiaxis_split_01/axis_assignments.json",
        ROOT / "ordered_site_features_01/ordered_features.npz", ROOT / "main_baselines_01/fresh_representations.npz",
    ])

    for block_number, block in enumerate(blocks):
        if not block["evaluable"]:
            continue
        labels = np.zeros((len(sequences), len(reactions)), dtype=bool)
        for edge_number in block["train_edge_indices"]:
            edge = edges[edge_number]
            labels[sequence_index[edge["sequence_sha256"]], reaction_index[edge["reaction_key"]]] = True
        train_sequences = np.flatnonzero(labels.any(axis=1))
        seen = labels.any(axis=0)
        unseen = ~seen
        truth = defaultdict(set)
        for edge_number in block["test_edge_indices"]:
            edge = edges[edge_number]
            truth[edge["sequence_sha256"]].add(reaction_index[edge["reaction_key"]])

        baseline_path = ROOT / f"main_baselines_01/block_{block_number:03d}_scores.npz"
        interaction_path = ROOT / f"interaction_conditional_sequence_01/block_{block_number:03d}_scores.npz"
        input_paths.update((baseline_path, interaction_path))
        baseline = np.load(baseline_path, allow_pickle=False)
        interaction = np.load(interaction_path, allow_pickle=False)
        controls = {}
        for control_name, folder in CONTROLS.items():
            path = folder / f"block_{block_number:03d}_scores.npz"
            input_paths.add(path)
            controls[control_name] = np.load(path, allow_pickle=False)
        if list(baseline["reaction_ids"]) != reactions or list(interaction["reaction_ids"]) != reactions:
            raise ValueError("score reaction order")
        baseline_query = {str(value): index for index, value in enumerate(baseline["query_ids"])}
        interaction_query = {str(value): index for index, value in enumerate(interaction["query_ids"])}
        control_query = {name: {str(value): index for index, value in enumerate(data["query_ids"])} for name, data in controls.items()}

        def score(data, query_map, method, query):
            values = data[method]
            domain = data["domain_" + method]
            if values.ndim == 2:
                values = values[query_map[query]]
            if domain.ndim == 2:
                domain = domain[query_map[query]]
            return values.astype(np.float64), domain.astype(bool)

        def add_panel(cell, source_task, query, positives, panel, methods, source, query_map,
                      confidence, confidence_definition):
            positives = sorted(positives)
            if not positives:
                return
            if not all(panel[index] for index in positives):
                raise ValueError(f"target outside {cell}")
            for method in methods:
                if method == "uniform_expectation":
                    result = random_expectations(int(np.sum(panel)), len(positives))
                else:
                    values, domain = score(source, query_map, method, query)
                    result = rank_panel(values, panel, domain, positives, tie)
                metric_rows.append({
                    "cell": cell, "source_task": source_task, "block_number": block_number,
                    "block_id": block["block_id"], "sequence_sha256": query,
                    "protein_group": axes["protein_groups"][query], "method": method,
                    "target_reaction_indices": positives,
                    "confidence": float(confidence), "confidence_definition": confidence_definition,
                    **result,
                })

        for query in sorted(truth):
            query_number = sequence_index[query]
            represented = bool(labels[query_number].any())
            available = bool(ordered_available[query_number])
            homology_scores = homology[query_number, train_sequences]
            labelled_homologue = bool(np.any(homology_scores > 0))
            maximum_homology = float(np.max(homology_scores)) if labelled_homologue else 0.0
            all_targets = truth[query]
            seen_targets = sorted(index for index in all_targets if seen[index])
            unseen_targets = sorted(index for index in all_targets if unseen[index])

            if block["task"] == "chemical_cold":
                if seen_targets:
                    raise ValueError("chemical-cold target marked seen")
                if represented and available:
                    cell = INTERACTION_CELL
                    route = "joint_global_interaction"
                    state = "candidate_route"
                    primary_score, _ = score(interaction, interaction_query, route, query)
                    ordered_scores = np.sort(primary_score[unseen])
                    confidence = float(ordered_scores[-1] - ordered_scores[-2])
                    definition = "top1_minus_top2_joint_global_interaction_within_unseen_panel"
                    add_panel(cell, block["task"], query, unseen_targets, unseen,
                              (route, "global_residual", "joint_global_ordered", "availability_chemistry_prior", "uniform_expectation"),
                              interaction, interaction_query, confidence, definition)
                    primary_global, primary_global_domain = score(interaction, interaction_query, "global_residual", query)
                    primary_metric = rank_panel(primary_score, unseen, np.ones_like(unseen), unseen_targets, tie)
                    global_metric = rank_panel(primary_global, unseen, primary_global_domain, unseen_targets, tie)
                    for control_name, control in controls.items():
                        qmap = control_query[control_name]
                        control_score, control_domain = score(control, qmap, route, query)
                        control_global, control_global_domain = score(control, qmap, "global_residual", query)
                        control_metric = rank_panel(control_score, unseen, control_domain, unseen_targets, tie)
                        control_global_metric = rank_panel(control_global, unseen, control_global_domain, unseen_targets, tie)
                        control_rows.append({
                            "control": control_name, "cell": cell, "source_task": block["task"],
                            "block_number": block_number, "sequence_sha256": query,
                            "protein_group": axes["protein_groups"][query],
                            "target_reaction_indices": unseen_targets,
                            "candidate_count": int(np.sum(unseen)),
                            "primary_rr": primary_metric["rr"], "primary_global_rr": global_metric["rr"],
                            "control_rr": control_metric["rr"], "control_global_rr": control_global_metric["rr"],
                        })
                elif not represented and labelled_homologue:
                    cell = TRANSPORT_CELLS["chemical_cold"]
                    route = "homology_chemical_transport"
                    state = "candidate_route_accounting_only"
                    add_panel(cell, block["task"], query, unseen_targets, unseen,
                              (route, "chemistry_prior", "uniform_expectation"), baseline, baseline_query,
                              maximum_homology, "maximum_qualifying_mmseqs_bitscore_to_labelled_training_protein")
                else:
                    cell = "new_protein_unseen_chemistry_without_supported_route__chemical_cold"
                    route = "abstain"
                    state = "abstain_no_homologue" if not labelled_homologue else "abstain_ordered_site_unavailable"
                decision_rows.append({
                    "source_task": block["task"], "target_subset": "unseen", "block_number": block_number,
                    "block_id": block["block_id"], "sequence_sha256": query,
                    "protein_group": axes["protein_groups"][query], "exact_protein_represented": represented,
                    "ordered_site_available": available, "labelled_homologue_available": labelled_homologue,
                    "maximum_qualifying_homology_bitscore": maximum_homology,
                    "candidate_count": int(np.sum(unseen)), "positive_count": len(unseen_targets),
                    "cell": cell, "candidate_route": route, "pre_gate_state": state,
                })

            elif block["task"] == "protein_cold":
                if represented:
                    raise ValueError("protein-cold query represented in training")
                for subset, positives, panel in (("seen", seen_targets, seen), ("unseen", unseen_targets, unseen)):
                    if not positives:
                        continue
                    if labelled_homologue and subset == "seen":
                        cell = MMSEQS_CELL
                        route = "mmseqs_weighted"
                        state = "candidate_route"
                        add_panel(cell, block["task"], query, positives, panel,
                                  (route, "mmseqs_top1", "chemistry_prior", "homology_chemical_transport", "uniform_expectation"),
                                  baseline, baseline_query, maximum_homology,
                                  "maximum_qualifying_mmseqs_bitscore_to_labelled_training_protein")
                    elif labelled_homologue:
                        cell = TRANSPORT_CELLS["protein_cold"]
                        route = "homology_chemical_transport"
                        state = "candidate_route_pending_gate"
                        add_panel(cell, block["task"], query, positives, panel,
                                  (route, "chemistry_prior", "uniform_expectation"), baseline, baseline_query,
                                  maximum_homology, "maximum_qualifying_mmseqs_bitscore_to_labelled_training_protein")
                    else:
                        cell = f"new_protein_{subset}_chemistry_without_labelled_homologue__protein_cold"
                        route = "abstain"
                        state = "abstain_no_homologue"
                    decision_rows.append({
                        "source_task": block["task"], "target_subset": subset, "block_number": block_number,
                        "block_id": block["block_id"], "sequence_sha256": query,
                        "protein_group": axes["protein_groups"][query], "exact_protein_represented": represented,
                        "ordered_site_available": available, "labelled_homologue_available": labelled_homologue,
                        "maximum_qualifying_homology_bitscore": maximum_homology,
                        "candidate_count": int(np.sum(panel)), "positive_count": len(positives),
                        "cell": cell, "candidate_route": route, "pre_gate_state": state,
                    })

            elif block["task"] == "double_cold":
                if represented or seen_targets:
                    raise ValueError("double-cold applicability invariant")
                if labelled_homologue:
                    cell = TRANSPORT_CELLS["double_cold"]
                    route = "homology_chemical_transport"
                    state = "candidate_route_pending_gate"
                    add_panel(cell, block["task"], query, unseen_targets, unseen,
                              (route, "chemistry_prior", "uniform_expectation"), baseline, baseline_query,
                              maximum_homology, "maximum_qualifying_mmseqs_bitscore_to_labelled_training_protein")
                else:
                    cell = "new_protein_unseen_chemistry_without_labelled_homologue__double_cold"
                    route = "abstain"
                    state = "abstain_no_homologue"
                decision_rows.append({
                    "source_task": block["task"], "target_subset": "unseen", "block_number": block_number,
                    "block_id": block["block_id"], "sequence_sha256": query,
                    "protein_group": axes["protein_groups"][query], "exact_protein_represented": represented,
                    "ordered_site_available": available, "labelled_homologue_available": labelled_homologue,
                    "maximum_qualifying_homology_bitscore": maximum_homology,
                    "candidate_count": int(np.sum(unseen)), "positive_count": len(unseen_targets),
                    "cell": cell, "candidate_route": route, "pre_gate_state": state,
                })
        baseline.close(); interaction.close()
        for data in controls.values():
            data.close()
        print(f"routed block {block_number + 1}/{len(blocks)}: {block['block_id']}", flush=True)

    pair_specs = [
        (INTERACTION_CELL, "joint_global_interaction", "global_residual"),
        (INTERACTION_CELL, "joint_global_interaction", "joint_global_ordered"),
        (INTERACTION_CELL, "joint_global_interaction", "availability_chemistry_prior"),
        (INTERACTION_CELL, "joint_global_interaction", "uniform_expectation"),
        (MMSEQS_CELL, "mmseqs_weighted", "chemistry_prior"),
        (MMSEQS_CELL, "mmseqs_weighted", "uniform_expectation"),
        (MMSEQS_CELL, "mmseqs_top1", "mmseqs_weighted"),
    ]
    for cell in TRANSPORT_CELLS.values():
        pair_specs.extend(((cell, "homology_chemical_transport", "chemistry_prior"),
                           (cell, "homology_chemical_transport", "uniform_expectation")))
    comparisons = []
    paired_rows = []
    for spec in pair_specs:
        result, rows = comparison(metric_rows, *spec)
        comparisons.append(result); paired_rows.extend(rows)
    control_results = control_summaries(control_rows)
    comparison_index = {(row["cell"], row["left_method"], row["right_method"]): row for row in comparisons}
    control_index = {row["control"]: row for row in control_results}

    interaction_primary = comparison_index[(INTERACTION_CELL, "joint_global_interaction", "global_residual")]
    interaction_control_not_reproduced = all(
        not positive_lower({"paired": control_index[name]["control_gain_over_global"]}) for name in CONTROLS
    )
    interaction_accepted = positive_lower(interaction_primary) and interaction_control_not_reproduced
    mmseqs_accepted = all(positive_lower(comparison_index[(MMSEQS_CELL, "mmseqs_weighted", right)])
                          for right in ("chemistry_prior", "uniform_expectation"))
    transport_acceptance = {}
    for task in ("protein_cold", "double_cold"):
        cell = TRANSPORT_CELLS[task]
        transport_acceptance[cell] = all(
            positive_lower(comparison_index[(cell, "homology_chemical_transport", right)])
            for right in ("chemistry_prior", "uniform_expectation")
        )
    transport_acceptance[TRANSPORT_CELLS["chemical_cold"]] = False

    accepted_cells = {INTERACTION_CELL: interaction_accepted, MMSEQS_CELL: mmseqs_accepted, **transport_acceptance}
    for row in decision_rows:
        if row["candidate_route"] == "abstain":
            row["final_state"] = row["pre_gate_state"]
        elif row["cell"] == TRANSPORT_CELLS["chemical_cold"]:
            row["final_state"] = "abstain_accounting_panel_not_a_strict_protein_novelty_gate"
        elif accepted_cells.get(row["cell"], False):
            row["final_state"] = "accepted_route"
        else:
            row["final_state"] = "abstain_failed_development_gate"

    coverage = []
    for task in sorted({row["source_task"] for row in decision_rows}):
        chosen = [row for row in decision_rows if row["source_task"] == task]
        accepted = [row for row in chosen if row["final_state"] == "accepted_route"]
        coverage.append({
            "source_task": task,
            "target_subpanels": len(chosen),
            "unique_queries": len({row["sequence_sha256"] for row in chosen}),
            "protein_groups": len({row["protein_group"] for row in chosen}),
            "positive_instances": sum(row["positive_count"] for row in chosen),
            "accepted_target_subpanels": len(accepted),
            "accepted_positive_instances": sum(row["positive_count"] for row in accepted),
            "accepted_target_subpanel_fraction": len(accepted) / len(chosen),
            "accepted_positive_instance_fraction": sum(row["positive_count"] for row in accepted) / sum(row["positive_count"] for row in chosen),
            "states": {state: sum(row["final_state"] == state for row in chosen) for state in sorted({row["final_state"] for row in chosen})},
            "task_denominators_are_not_additive_across_tasks": True,
        })

    curves = []
    if interaction_accepted:
        curves.append(confidence_curve(metric_rows, INTERACTION_CELL, "joint_global_interaction"))
    if mmseqs_accepted:
        curves.append(confidence_curve(metric_rows, MMSEQS_CELL, "mmseqs_weighted"))
    for cell, accepted in transport_acceptance.items():
        if accepted:
            curves.append(confidence_curve(metric_rows, cell, "homology_chemical_transport"))

    protein_did = control_index["protein_position_permutation"]["primary_minus_control_gain"]
    reaction_did = control_index["reaction_center_row_permutation"]["primary_minus_control_gain"]
    acceptance = {
        "status": "COMPUTED_QC_PENDING",
        "interaction_route": {
            "accepted": interaction_accepted,
            "positive_vs_global": positive_lower(interaction_primary),
            "permutation_gains_do_not_have_positive_intervals": interaction_control_not_reproduced,
            "position_specificity_supported": bool(protein_did["conditional_protein_bootstrap_95ci"] and protein_did["conditional_protein_bootstrap_95ci"][0] > 0),
            "reaction_identity_specificity_supported": bool(reaction_did["conditional_protein_bootstrap_95ci"] and reaction_did["conditional_protein_bootstrap_95ci"][0] > 0),
        },
        "mmseqs_seen_chemistry_route": {
            "accepted": mmseqs_accepted,
            "positive_vs_chemistry_prior": positive_lower(comparison_index[(MMSEQS_CELL, "mmseqs_weighted", "chemistry_prior")]),
            "positive_vs_uniform_expectation": positive_lower(comparison_index[(MMSEQS_CELL, "mmseqs_weighted", "uniform_expectation")]),
        },
        "homology_transport_unseen_chemistry": {
            cell: {"accepted": accepted, "chemical_cold_accounting_only": cell.endswith("accounting")}
            for cell, accepted in transport_acceptance.items()
        },
        "cross_domain_score_fusion_performed": False,
        "structure_route_enabled": False,
        "independent_biological_validation": False,
        "whole_tool_complete": False,
    }

    OUTPUT.mkdir()
    for name, rows in (("query_panel_metrics.jsonl", metric_rows), ("route_decisions.jsonl", decision_rows),
                       ("paired_query_records.jsonl", paired_rows), ("control_query_records.jsonl", control_rows)):
        with (OUTPUT / name).open("w") as handle:
            for row in rows:
                handle.write(json.dumps(row, allow_nan=False) + "\n")
    write_json(OUTPUT / "method_summary.json", {
        "status": "COMPUTED_QC_PENDING", "created_utc": now(), "summaries": summarize_methods(metric_rows),
        "metric_rows": len(metric_rows), "candidate_catalogue": len(reactions),
    })
    write_json(OUTPUT / "comparison_summary.json", {
        "status": "COMPUTED_QC_PENDING", "created_utc": now(), "comparisons": comparisons,
        "control_comparisons": control_results, "paired_query_rows": len(paired_rows),
        "control_query_rows": len(control_rows), "bootstrap_seed": BOOTSTRAP_SEED,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "interval_scope": "Conditional fixed-fit protein-group bootstrap; no split, chemical-component, publication, refit, external or multiplicity uncertainty",
    })
    write_json(OUTPUT / "routing_coverage.json", {
        "status": "COMPUTED_QC_PENDING", "created_utc": now(), "by_source_task": coverage,
        "target_subpanels": len(decision_rows), "task_denominators_are_not_additive": True,
    })
    write_json(OUTPUT / "selective_risk_curves.json", {
        "status": "COMPUTED_QC_PENDING", "created_utc": now(), "curves": curves,
        "coverage_cuts": list(CUTS), "outer_threshold_selected": False,
        "risk_definition": "one minus protein-group-macro hit@10 among documented-positive query panels",
        "not_precision_probability_or_calibration": True,
    })
    write_json(OUTPUT / "acceptance.json", acceptance)
    manifest = {path.relative_to(ROOT).as_posix(): digest_file(path) for path in sorted(input_paths)}
    write_json(OUTPUT / "input_manifest.json", manifest)
    print(json.dumps({
        "status": "COMPUTED_QC_PENDING", "metric_rows": len(metric_rows),
        "decision_rows": len(decision_rows), "paired_rows": len(paired_rows),
        "control_rows": len(control_rows), "acceptance": acceptance,
    }, indent=2))


if __name__ == "__main__":
    main()
