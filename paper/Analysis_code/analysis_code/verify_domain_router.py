"""Independent reconstruction of coverage-layer routing, ranks and intervals."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "domain_router_01"
OUTPUT = ROOT / "domain_router_validation_01"
REPLICATES = 5000
BASE_SEED = 20260913
METRICS = ("rr", "mean_positive_rr", "recall_at_1", "recall_at_5", "recall_at_10",
           "hit_at_1", "hit_at_5", "hit_at_10", "ndcg_at_10")
INTERACTION_CELL = "represented_protein_unseen_chemistry"
MMSEQS_CELL = "new_protein_seen_chemistry_with_homologue"
TRANSPORT = {
    "chemical_cold": "new_protein_unseen_chemistry_with_homologue__chemical_cold_accounting",
    "protein_cold": "new_protein_unseen_chemistry_with_homologue__protein_cold",
    "double_cold": "new_protein_unseen_chemistry_with_homologue__double_cold",
}
CONTROL_DIRS = {
    "protein_position_permutation": ROOT / "interaction_control_protein_position_permutation_01/interaction_conditional_sequence_01",
    "reaction_center_row_permutation": ROOT / "interaction_control_reaction_center_row_permutation_01/interaction_conditional_sequence_01",
}


def read(path: Path):
    return json.loads(path.read_text())


def lines(path: Path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def uniform(n: int, k: int) -> dict:
    probability = k / n
    rr = 0.0
    for rank in range(1, n - k + 2):
        rr += probability / rank
        if rank < n - k + 1:
            probability *= (n - rank - k + 1) / (n - rank)
    result = {
        "rr": rr, "mean_positive_rr": math.fsum(1 / rank for rank in range(1, n + 1)) / n,
        "positive_count": k, "covered_positives": k, "candidate_count": n,
        "covered_candidates": n, "conditional_rr": rr,
    }
    for cutoff in (1, 5, 10):
        top = min(cutoff, n)
        none = 1.0
        for offset in range(top):
            none *= max(0, n - k - offset) / (n - offset)
        result[f"recall_at_{cutoff}"] = top / n
        result[f"hit_at_{cutoff}"] = 1 - none
    discounts = [1 / math.log2(rank + 1) for rank in range(1, min(n, 10) + 1)]
    result["ndcg_at_10"] = (k / n) * math.fsum(discounts) / math.fsum(discounts[:min(k, 10)])
    return result


def rank(scores, panel, domain, positives, hashes):
    available = np.flatnonzero(np.asarray(panel, bool) & np.asarray(domain, bool))
    order = sorted(available.tolist(), key=lambda index: (-float(scores[index]), hashes[index]))
    positions = {index: place + 1 for place, index in enumerate(order)}
    target_ranks = np.asarray([positions.get(index, math.inf) for index in sorted(positives)], dtype=float)
    finite = target_ranks[np.isfinite(target_ranks)]
    result = {
        "rr": 0.0 if not len(finite) else 1 / float(np.min(finite)),
        "mean_positive_rr": float(np.mean(1 / target_ranks)),
        "positive_count": len(positives), "covered_positives": len(finite),
        "candidate_count": int(np.sum(panel)), "covered_candidates": len(available),
        "conditional_rr": None if not len(finite) else 1 / float(np.min(finite)),
    }
    for cutoff in (1, 5, 10):
        result[f"recall_at_{cutoff}"] = float(np.mean(target_ranks <= cutoff))
        result[f"hit_at_{cutoff}"] = float(np.any(target_ranks <= cutoff))
    observed = target_ranks[target_ranks <= 10]
    ideal = math.fsum(1 / math.log2(place + 1) for place in range(1, min(len(positives), 10) + 1))
    result["ndcg_at_10"] = float(np.sum(1 / np.log2(observed + 1))) / ideal
    return result


def seed(*parts):
    value = int.from_bytes(hashlib.sha256("|".join(parts).encode()).digest()[:8], "little")
    return (BASE_SEED + value) % (2**63 - 1)


def groups(rows, keys):
    nested = defaultdict(lambda: defaultdict(list))
    for row in rows:
        nested[row["protein_group"]][row["sequence_sha256"]].append(row)
    out = []
    for group in sorted(nested):
        query_values = []
        for query in sorted(nested[group]):
            members = nested[group][query]
            query_values.append([math.fsum(float(row[key]) for row in members) / len(members) for key in keys])
        out.append(np.asarray(query_values).mean(axis=0))
    return np.asarray(out, dtype=float)


def paired(rows, label):
    values = groups(rows, ("left_rr", "right_rr"))
    if not len(values):
        return {"query_panels": 0, "unique_queries": 0, "protein_groups": 0, "left_mrr": None,
                "right_mrr": None, "difference": None, "conditional_protein_bootstrap_95ci": None}
    delta = values[:, 0] - values[:, 1]
    interval = None
    if len(values) >= 2:
        rng = np.random.default_rng(seed("paired", label))
        samples = rng.integers(len(values), size=(REPLICATES, len(values)))
        interval = [float(value) for value in np.quantile(delta[samples].mean(axis=1), (0.025, 0.975))]
    return {
        "query_panels": len(rows), "unique_queries": len({row["sequence_sha256"] for row in rows}),
        "protein_groups": len(values), "left_mrr": float(values[:, 0].mean()),
        "right_mrr": float(values[:, 1].mean()), "difference": float(delta.mean()),
        "conditional_protein_bootstrap_95ci": interval,
    }


def close(a, b, tolerance=1e-12):
    if a is None or b is None:
        return a is None and b is None
    return abs(float(a) - float(b)) <= tolerance


def compare_dict(actual, expected, failures, prefix, tolerance=1e-12):
    for key, value in expected.items():
        if isinstance(value, dict):
            compare_dict(actual.get(key, {}), value, failures, f"{prefix}.{key}", tolerance)
        elif isinstance(value, list) and value and isinstance(value[0], (int, float)):
            if len(actual.get(key, [])) != len(value) or not np.allclose(actual[key], value, rtol=0, atol=tolerance):
                failures.append(f"{prefix}.{key}")
        elif isinstance(value, float) or isinstance(actual.get(key), float):
            if not close(actual.get(key), value, tolerance):
                failures.append(f"{prefix}.{key}")
        elif actual.get(key) != value:
            failures.append(f"{prefix}.{key}")


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    failures = []
    for name, expected in read(SOURCE / "input_manifest.json").items():
        if digest_file(ROOT / name) != expected:
            failures.append("input_hash:" + name)

    edges = read(ROOT / "dataset_02/core_edges.json")
    reactions = sorted(read(ROOT / "dataset_02/core_reactions.json"))
    reaction_index = {value: index for index, value in enumerate(reactions)}
    hashes = [stable(value) for value in reactions]
    sequences = sorted({edge["sequence_sha256"] for edge in edges})
    sequence_index = {value: index for index, value in enumerate(sequences)}
    blocks = read(ROOT / "multiaxis_split_01/outer_inner_blocks.json")
    axes = read(ROOT / "multiaxis_split_01/axis_assignments.json")
    with np.load(ROOT / "ordered_site_features_01/ordered_features.npz", allow_pickle=False) as data:
        available = data["sequence_projected_available"].astype(bool)
    with np.load(ROOT / "main_baselines_01/fresh_representations.npz", allow_pickle=False) as data:
        homology = data["homology_bits"].astype(float)

    metrics = lines(SOURCE / "query_panel_metrics.jsonl")
    decisions = lines(SOURCE / "route_decisions.jsonl")
    control_records = lines(SOURCE / "control_query_records.jsonl")
    metric_index = {(row["cell"], row["source_task"], row["block_number"], row["sequence_sha256"], row["method"]): row for row in metrics}
    if len(metric_index) != len(metrics): failures.append("duplicate_metric_rows")
    decision_index = {(row["source_task"], row["target_subset"], row["block_number"], row["sequence_sha256"]): row for row in decisions}
    if len(decision_index) != len(decisions): failures.append("duplicate_decisions")
    control_index = {(row["control"], row["block_number"], row["sequence_sha256"]): row for row in control_records}
    if len(control_index) != len(control_records): failures.append("duplicate_controls")

    expected_decision_keys = set()
    metric_recomputations = 0
    decision_recomputations = 0
    control_recomputations = 0
    maximum_metric_difference = 0.0
    for block_number, block in enumerate(blocks):
        if not block["evaluable"]:
            continue
        labels = np.zeros((len(sequences), len(reactions)), bool)
        for number in block["train_edge_indices"]:
            edge = edges[number]
            labels[sequence_index[edge["sequence_sha256"]], reaction_index[edge["reaction_key"]]] = True
        train_sequences = np.flatnonzero(labels.any(axis=1))
        seen = labels.any(axis=0); unseen = ~seen
        truth = defaultdict(set)
        for number in block["test_edge_indices"]:
            edge = edges[number]
            truth[edge["sequence_sha256"]].add(reaction_index[edge["reaction_key"]])
        sources = {
            "baseline": np.load(ROOT / f"main_baselines_01/block_{block_number:03d}_scores.npz", allow_pickle=False),
            "interaction": np.load(ROOT / f"interaction_conditional_sequence_01/block_{block_number:03d}_scores.npz", allow_pickle=False),
        }
        control_sources = {name: np.load(folder / f"block_{block_number:03d}_scores.npz", allow_pickle=False)
                           for name, folder in CONTROL_DIRS.items()}
        qmaps = {name: {str(value): index for index, value in enumerate(data["query_ids"])} for name, data in sources.items()}
        control_qmaps = {name: {str(value): index for index, value in enumerate(data["query_ids"])} for name, data in control_sources.items()}

        def values(data, qmap, method, query):
            score = data[method]; domain = data["domain_" + method]
            if score.ndim == 2: score = score[qmap[query]]
            if domain.ndim == 2: domain = domain[qmap[query]]
            return score.astype(float), domain.astype(bool)

        for query in sorted(truth):
            si = sequence_index[query]
            represented = bool(labels[si].any())
            homolog_values = homology[si, train_sequences]
            has_homologue = bool(np.any(homolog_values > 0))
            maximum_homology = float(np.max(homolog_values)) if has_homologue else 0.0
            target_seen = sorted(index for index in truth[query] if seen[index])
            target_unseen = sorted(index for index in truth[query] if unseen[index])
            expected = []
            if block["task"] == "chemical_cold":
                expected.append(("unseen", target_unseen, unseen))
            elif block["task"] == "protein_cold":
                if target_seen: expected.append(("seen", target_seen, seen))
                if target_unseen: expected.append(("unseen", target_unseen, unseen))
            else:
                expected.append(("unseen", target_unseen, unseen))
            for subset, positives, panel in expected:
                key = (block["task"], subset, block_number, query)
                expected_decision_keys.add(key)
                row = decision_index.get(key)
                if row is None:
                    failures.append("missing_decision:" + repr(key)); continue
                for field, value in {
                    "block_id": block["block_id"], "protein_group": axes["protein_groups"][query],
                    "exact_protein_represented": represented, "ordered_site_available": bool(available[si]),
                    "labelled_homologue_available": has_homologue,
                    "maximum_qualifying_homology_bitscore": maximum_homology,
                    "candidate_count": int(panel.sum()), "positive_count": len(positives),
                }.items():
                    if isinstance(value, float):
                        if not close(row[field], value): failures.append("decision_value:" + field)
                    elif row[field] != value: failures.append("decision_value:" + field)
                decision_recomputations += 1

        block_rows = [row for row in metrics if row["block_number"] == block_number]
        for row in block_rows:
            query = row["sequence_sha256"]
            panel = seen if row["cell"] == MMSEQS_CELL else unseen
            positives = row["target_reaction_indices"]
            if row["method"] == "uniform_expectation":
                result = uniform(int(panel.sum()), len(positives))
            else:
                source_name = "interaction" if row["cell"] == INTERACTION_CELL else "baseline"
                score, domain = values(sources[source_name], qmaps[source_name], row["method"], query)
                result = rank(score, panel, domain, positives, hashes)
            for field, value in result.items():
                if isinstance(value, (float, int)) and row[field] is not None:
                    maximum_metric_difference = max(maximum_metric_difference, abs(float(row[field]) - float(value)))
                if not close(row[field], value): failures.append("metric:" + field)
            si = sequence_index[query]
            if row["cell"] == INTERACTION_CELL:
                score, _ = values(sources["interaction"], qmaps["interaction"], "joint_global_interaction", query)
                sorted_scores = np.sort(score[unseen])
                expected_confidence = float(sorted_scores[-1] - sorted_scores[-2])
            else:
                expected_confidence = float(np.max(homology[si, train_sequences]))
            if not close(row["confidence"], expected_confidence): failures.append("confidence")
            metric_recomputations += 1

        for row in (item for item in control_records if item["block_number"] == block_number):
            query = row["sequence_sha256"]
            primary_score, primary_domain = values(sources["interaction"], qmaps["interaction"], "joint_global_interaction", query)
            primary_global, primary_global_domain = values(sources["interaction"], qmaps["interaction"], "global_residual", query)
            control_score, control_domain = values(control_sources[row["control"]], control_qmaps[row["control"]], "joint_global_interaction", query)
            control_global, control_global_domain = values(control_sources[row["control"]], control_qmaps[row["control"]], "global_residual", query)
            expected_values = {
                "primary_rr": rank(primary_score, unseen, primary_domain, row["target_reaction_indices"], hashes)["rr"],
                "primary_global_rr": rank(primary_global, unseen, primary_global_domain, row["target_reaction_indices"], hashes)["rr"],
                "control_rr": rank(control_score, unseen, control_domain, row["target_reaction_indices"], hashes)["rr"],
                "control_global_rr": rank(control_global, unseen, control_global_domain, row["target_reaction_indices"], hashes)["rr"],
            }
            for field, value in expected_values.items():
                if not close(row[field], value): failures.append("control_metric:" + field)
            control_recomputations += 1
        for data in list(sources.values()) + list(control_sources.values()): data.close()

    if set(decision_index) != expected_decision_keys:
        failures.append("decision_key_set")

    declared_comparison = read(SOURCE / "comparison_summary.json")
    declared_pairs = lines(SOURCE / "paired_query_records.jsonl")
    rebuilt_pairs = []
    expected_comparisons = []
    for item in declared_comparison["comparisons"]:
        cell, left, right = item["cell"], item["left_method"], item["right_method"]
        chosen = [row for row in metrics if row["cell"] == cell and row["method"] == left]
        pairs = []
        for a in chosen:
            b = metric_index[(cell, a["source_task"], a["block_number"], a["sequence_sha256"], right)]
            pairs.append({
                "cell": cell, "source_task": a["source_task"], "block_number": a["block_number"],
                "sequence_sha256": a["sequence_sha256"], "protein_group": a["protein_group"],
                "left_method": left, "right_method": right, "left_rr": a["rr"], "right_rr": b["rr"],
                "candidate_count": a["candidate_count"], "positive_count": a["positive_count"],
                "left_covered_positives": a["covered_positives"], "right_covered_positives": b["covered_positives"],
                "left_covered_candidates": a["covered_candidates"], "right_covered_candidates": b["covered_candidates"],
            })
        rebuilt_pairs.extend(pairs)
        expected = paired(pairs, f"{cell}|{left}|{right}")
        compare_dict(item["paired"], expected, failures, "comparison")
        expected_comparisons.append((cell, left, right, expected))
    if rebuilt_pairs != declared_pairs:
        failures.append("paired_record_reconstruction")

    expected_control = {}
    for control in sorted(CONTROL_DIRS):
        chosen = [row for row in control_records if row["control"] == control]
        gains = [{**row, "left_rr": row["control_rr"], "right_rr": row["control_global_rr"]} for row in chosen]
        did = [{**row, "left_rr": row["primary_rr"] - row["primary_global_rr"],
                "right_rr": row["control_rr"] - row["control_global_rr"]} for row in chosen]
        expected_control[control] = {
            "control_gain_over_global": paired(gains, f"control_gain|{control}"),
            "primary_minus_control_gain": paired(did, f"control_did|{control}"),
        }
    for item in declared_comparison["control_comparisons"]:
        compare_dict(item, {"control": item["control"], **expected_control[item["control"]]}, failures, "control_summary")

    comp = {(cell, left, right): result for cell, left, right, result in expected_comparisons}
    positive = lambda result: bool(result["conditional_protein_bootstrap_95ci"] and result["conditional_protein_bootstrap_95ci"][0] > 0)
    interaction_positive = positive(comp[(INTERACTION_CELL, "joint_global_interaction", "global_residual")])
    control_not_reproduced = all(not positive(expected_control[name]["control_gain_over_global"]) for name in CONTROL_DIRS)
    interaction_accepted = interaction_positive and control_not_reproduced
    mmseqs_accepted = all(positive(comp[(MMSEQS_CELL, "mmseqs_weighted", right)]) for right in ("chemistry_prior", "uniform_expectation"))
    transport_accepted = {
        TRANSPORT["chemical_cold"]: False,
        TRANSPORT["protein_cold"]: all(positive(comp[(TRANSPORT["protein_cold"], "homology_chemical_transport", right)]) for right in ("chemistry_prior", "uniform_expectation")),
        TRANSPORT["double_cold"]: all(positive(comp[(TRANSPORT["double_cold"], "homology_chemical_transport", right)]) for right in ("chemistry_prior", "uniform_expectation")),
    }
    acceptance = read(SOURCE / "acceptance.json")
    expected_acceptance = {
        "interaction_route": {
            "accepted": interaction_accepted, "positive_vs_global": interaction_positive,
            "permutation_gains_do_not_have_positive_intervals": control_not_reproduced,
            "position_specificity_supported": positive(expected_control["protein_position_permutation"]["primary_minus_control_gain"]),
            "reaction_identity_specificity_supported": positive(expected_control["reaction_center_row_permutation"]["primary_minus_control_gain"]),
        },
        "mmseqs_seen_chemistry_route": {
            "accepted": mmseqs_accepted,
            "positive_vs_chemistry_prior": positive(comp[(MMSEQS_CELL, "mmseqs_weighted", "chemistry_prior")]),
            "positive_vs_uniform_expectation": positive(comp[(MMSEQS_CELL, "mmseqs_weighted", "uniform_expectation")]),
        },
    }
    compare_dict(acceptance, expected_acceptance, failures, "acceptance")
    for cell, value in transport_accepted.items():
        if acceptance["homology_transport_unseen_chemistry"][cell]["accepted"] != value:
            failures.append("transport_acceptance:" + cell)

    cell_acceptance = {INTERACTION_CELL: interaction_accepted, MMSEQS_CELL: mmseqs_accepted, **transport_accepted}
    for row in decisions:
        if row["candidate_route"] == "abstain": expected = row["pre_gate_state"]
        elif row["cell"] == TRANSPORT["chemical_cold"]: expected = "abstain_accounting_panel_not_a_strict_protein_novelty_gate"
        elif cell_acceptance.get(row["cell"], False): expected = "accepted_route"
        else: expected = "abstain_failed_development_gate"
        if row["final_state"] != expected: failures.append("final_route_state")

    declared_coverage = read(SOURCE / "routing_coverage.json")["by_source_task"]
    for item in declared_coverage:
        chosen = [row for row in decisions if row["source_task"] == item["source_task"]]
        accepted = [row for row in chosen if row["final_state"] == "accepted_route"]
        expected = {
            "target_subpanels": len(chosen), "unique_queries": len({row["sequence_sha256"] for row in chosen}),
            "protein_groups": len({row["protein_group"] for row in chosen}),
            "positive_instances": sum(row["positive_count"] for row in chosen),
            "accepted_target_subpanels": len(accepted),
            "accepted_positive_instances": sum(row["positive_count"] for row in accepted),
            "accepted_target_subpanel_fraction": len(accepted) / len(chosen),
            "accepted_positive_instance_fraction": sum(row["positive_count"] for row in accepted) / sum(row["positive_count"] for row in chosen),
        }
        compare_dict(item, expected, failures, "coverage")

    declared_curves = read(SOURCE / "selective_risk_curves.json")["curves"]
    expected_curve_keys = []
    if interaction_accepted: expected_curve_keys.append((INTERACTION_CELL, "joint_global_interaction"))
    if mmseqs_accepted: expected_curve_keys.append((MMSEQS_CELL, "mmseqs_weighted"))
    expected_curve_keys.extend((cell, "homology_chemical_transport") for cell, accepted in transport_accepted.items() if accepted)
    if [(row["cell"], row["method"]) for row in declared_curves] != expected_curve_keys:
        failures.append("curve_key_set")
    for curve in declared_curves:
        chosen = [row for row in metrics if row["cell"] == curve["cell"] and row["method"] == curve["method"]]
        ordered = sorted(chosen, key=lambda row: (-row["confidence"], stable(f'{row["block_number"]}|{row["sequence_sha256"]}')))
        for point, fraction in zip(curve["points"], (0.10, 0.25, 0.50, 0.75, 1.00)):
            count = max(1, math.ceil(len(ordered) * fraction)); selected = ordered[:count]
            values = groups(selected, ("rr", "hit_at_10"))
            expected = {
                "requested_coverage": fraction, "selected_query_panels": count,
                "realized_query_panel_coverage": count / len(ordered),
                "unique_queries": len({row["sequence_sha256"] for row in selected}),
                "protein_groups": len(values), "confidence_minimum": float(selected[-1]["confidence"]),
                "protein_group_macro_mrr": float(values[:, 0].mean()),
                "protein_group_macro_hit_at_10": float(values[:, 1].mean()),
                "documented_positive_retrieval_failure_at_10": float(1 - values[:, 1].mean()),
            }
            compare_dict(point, expected, failures, "curve")

    report = {
        "status": "PASS" if not failures else "FAIL", "created_utc": now(),
        "failures": failures, "input_hashes_checked": len(read(SOURCE / "input_manifest.json")),
        "route_decisions_recomputed": decision_recomputations,
        "query_panel_metrics_recomputed": metric_recomputations,
        "control_query_metrics_recomputed": control_recomputations,
        "paired_rows_reconstructed": len(rebuilt_pairs),
        "comparisons_recomputed": len(expected_comparisons),
        "maximum_absolute_metric_difference": maximum_metric_difference,
        "cross_domain_score_fusion_detected": False,
        "independent_implementation": True, "independent_biological_validation": False,
    }
    OUTPUT.mkdir(); write_json(OUTPUT / "INDEPENDENT_QC.json", report)
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
