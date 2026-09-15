"""Nested CLEAN candidate pruning over the frozen coverage-layer experts."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "clean_pruning_01"
MODES = ("none", "open_L2", "open_L3", "open_L4", "strict_L2", "strict_L3", "strict_L4")
INTERACTION_CELL = "represented_protein_unseen_chemistry"
MMSEQS_CELL = "new_protein_seen_chemistry_with_homologue"
TRANSPORT_CELLS = {
    "chemical_cold": "new_protein_unseen_chemistry_with_homologue__chemical_cold_accounting",
    "protein_cold": "new_protein_unseen_chemistry_with_homologue__protein_cold",
    "double_cold": "new_protein_unseen_chemistry_with_homologue__double_cold",
}
BOOTSTRAP_REPLICATES = 5000
BOOTSTRAP_SEED = 20260914


def read(path):
    return json.loads(path.read_text())


def construct(split, edges, sequence_index, reaction_index, shape):
    labels = np.zeros(shape, dtype=np.float32)
    truth = defaultdict(set)
    for edge_number in split["train_edge_indices"]:
        edge = edges[edge_number]
        labels[sequence_index[edge["sequence_sha256"]], reaction_index[edge["reaction_key"]]] = 1
    for edge_number in split["test_edge_indices"]:
        edge = edges[edge_number]
        truth[edge["sequence_sha256"]].add(reaction_index[edge["reaction_key"]])
    return labels, truth


def train_scores(labels, chemical_kernel, homology, query_indices, training_indices):
    frequency = labels.sum(axis=0)
    balanced = labels / np.maximum(labels.sum(axis=1, keepdims=True), 1)
    query_homology = homology[np.asarray(query_indices)]
    weighted = query_homology @ labels
    relevant = query_homology[:, training_indices]
    transport = (relevant @ balanced[training_indices] / np.maximum(relevant.sum(axis=1, keepdims=True), 1)) @ chemical_kernel
    return {
        "mmseqs_weighted": weighted,
        "homology_chemical_transport": transport,
    }, {
        "mmseqs_weighted": weighted > 0,
        "homology_chemical_transport": np.broadcast_to(relevant.sum(axis=1, keepdims=True) > 0, transport.shape).copy(),
    }


def reaction_mappings(reaction_ids):
    payload = read(ROOT / "clean_ec_bridge_01/reaction_ec_support.json")
    output = []
    for reaction in reaction_ids:
        output.append({ec: set(value["supporting_sequences"])
                       for ec, value in payload[reaction].items() if value["in_clean_catalog"]})
    return output


def compatible(left, right, level):
    return tuple(left.split(".")[:level]) == tuple(right.split(".")[:level])


def clean_mask(mode, query, calls, mappings):
    if mode == "none" or not calls:
        return np.ones(len(mappings), dtype=bool)
    rule, level_name = mode.split("_")
    level = int(level_name[1:])
    mask = np.zeros(len(mappings), dtype=bool)
    for index, mapping in enumerate(mappings):
        available = [ec for ec, supporters in mapping.items() if any(value != query for value in supporters)]
        match = any(compatible(ec, call, level) for ec in available for call in calls)
        mask[index] = (not available or match) if rule == "open" else bool(available and match)
    return mask


def rank(scores, panel, domain, pruning, positives, tie):
    truth = np.asarray(sorted(positives), dtype=int)
    if not len(truth) or not np.all(np.asarray(panel, dtype=bool)[truth]):
        raise ValueError("positive outside candidate panel")
    defined = np.asarray(panel, dtype=bool) & np.asarray(domain, dtype=bool) & np.asarray(pruning, dtype=bool)
    if not np.all(np.isfinite(scores[defined])):
        raise ValueError("non-finite expert score")
    selected = np.flatnonzero(defined)
    order = selected[np.lexsort((tie[selected], -scores[selected]))]
    ranks = np.full(len(scores), np.inf)
    ranks[order] = np.arange(1, len(order) + 1)
    target = ranks[truth]
    best = float(np.min(target))
    result = {
        "rr": 0.0 if not math.isfinite(best) else 1 / best,
        "positive_count": len(truth), "covered_positives": int(np.sum(defined[truth])),
        "candidate_count": int(np.sum(panel)), "covered_candidates": len(selected),
    }
    for cutoff in (1, 5, 10):
        result[f"recall_at_{cutoff}"] = float(np.mean(target <= cutoff))
        result[f"hit_at_{cutoff}"] = float(np.any(target <= cutoff))
    return result


def route_panels(task, query, labels, truth, homology, ordered_available, sequence_index):
    query_number = sequence_index[query]
    represented = bool(labels[query_number].any())
    training_sequences = np.flatnonzero(labels.any(axis=1))
    seen = labels.any(axis=0)
    homolog = bool(np.any(homology[query_number, training_sequences] > 0))
    targets = truth[query]
    seen_targets = set(index for index in targets if seen[index])
    unseen_targets = set(targets) - seen_targets
    output = []
    if task == "chemical_cold":
        if represented and ordered_available[query_number] and unseen_targets:
            output.append((INTERACTION_CELL, "joint_global_interaction", ~seen, unseen_targets))
        elif not represented and homolog and unseen_targets:
            output.append((TRANSPORT_CELLS[task], "homology_chemical_transport", ~seen, unseen_targets))
    elif task == "protein_cold":
        if not represented and homolog and seen_targets:
            output.append((MMSEQS_CELL, "mmseqs_weighted", seen, seen_targets))
        if not represented and homolog and unseen_targets:
            output.append((TRANSPORT_CELLS[task], "homology_chemical_transport", ~seen, unseen_targets))
    elif task == "double_cold":
        if not represented and homolog and unseen_targets:
            output.append((TRANSPORT_CELLS[task], "homology_chemical_transport", ~seen, unseen_targets))
    return output


def group_macro(rows, key="rr"):
    by_group = defaultdict(lambda: defaultdict(list))
    for row in rows:
        by_group[row["protein_group"]][row["sequence_sha256"]].append(float(row[key]))
    values = []
    for group in sorted(by_group):
        per_query = [float(np.mean(values)) for values in by_group[group].values()]
        values.append(float(np.mean(per_query)))
    return values


def choose_mode(rows):
    if not rows:
        return "none", {mode: None for mode in MODES}
    values = {}
    for mode in MODES:
        chosen = [row for row in rows if row["mode"] == mode]
        groups = group_macro(chosen)
        values[mode] = None if not groups else float(np.mean(groups))
    best = max(value for value in values.values() if value is not None)
    selected = next(mode for mode in MODES if values[mode] is not None and abs(values[mode] - best) <= 1e-15)
    return selected, values


def seed_for(value):
    return (BOOTSTRAP_SEED + int.from_bytes(hashlib.sha256(value.encode()).digest()[:8], "little")) % (2**63 - 1)


def paired_summary(rows, label):
    by_group = defaultdict(lambda: defaultdict(list))
    for row in rows:
        by_group[row["protein_group"]][row["sequence_sha256"]].append((row["selected_rr"], row["none_rr"]))
    pairs = []
    for group in sorted(by_group):
        query_pairs = [np.mean(values, axis=0) for values in by_group[group].values()]
        pairs.append(np.mean(query_pairs, axis=0))
    array = np.asarray(pairs, dtype=float)
    if not len(array):
        return {"query_panels": 0, "unique_queries": 0, "protein_groups": 0,
                "selected_mrr": None, "unpruned_mrr": None, "difference": None,
                "conditional_protein_bootstrap_95ci": None}
    delta = array[:, 0] - array[:, 1]
    interval = None
    if len(array) >= 2:
        rng = np.random.default_rng(seed_for(label))
        draws = rng.integers(len(array), size=(BOOTSTRAP_REPLICATES, len(array)))
        interval = [float(value) for value in np.quantile(delta[draws].mean(axis=1), (0.025, 0.975))]
    return {
        "query_panels": len(rows), "unique_queries": len({row["sequence_sha256"] for row in rows}),
        "protein_groups": len(array), "selected_mrr": float(np.mean(array[:, 0])),
        "unpruned_mrr": float(np.mean(array[:, 1])), "difference": float(np.mean(delta)),
        "conditional_protein_bootstrap_95ci": interval,
    }


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    for gate in (ROOT / "clean_ec_bridge_validation_01/INDEPENDENT_QC.json",
                 ROOT / "clean_core_features_validation_01/INDEPENDENT_QC.json",
                 ROOT / "domain_router_validation_01/INDEPENDENT_QC.json"):
        if read(gate)["status"] != "PASS": raise ValueError("unpassed gate: " + str(gate))
    edges = read(ROOT / "dataset_02/core_edges.json")
    reactions = sorted(read(ROOT / "dataset_02/core_reactions.json"))
    sequences = sorted({edge["sequence_sha256"] for edge in edges})
    sequence_index = {value: index for index, value in enumerate(sequences)}
    reaction_index = {value: index for index, value in enumerate(reactions)}
    axes = read(ROOT / "multiaxis_split_01/axis_assignments.json")
    blocks = read(ROOT / "multiaxis_split_01/outer_inner_blocks.json")
    mappings = reaction_mappings(reactions)
    calls_payload = read(ROOT / "clean_core_features_01/maximum_separation_calls.json")
    calls = {query: [row["ec"] for row in values] for query, values in calls_payload.items()}
    feature_manifest = {row["sequence_sha256"]: row for row in
                        [json.loads(line) for line in (ROOT / "clean_core_features_01/sequence_manifest.jsonl").read_text().splitlines()]}
    with np.load(ROOT / "main_baselines_01/fresh_representations.npz", allow_pickle=False) as saved:
        if list(saved["sequence_ids"]) != sequences or list(saved["reaction_ids"]) != reactions: raise ValueError("representation order")
        chemical_kernel = saved["chemical_kernel"].astype(np.float64)
        homology = saved["homology_bits"].astype(np.float64)
    with np.load(ROOT / "ordered_site_features_01/ordered_features.npz", allow_pickle=False) as saved:
        if list(saved["sequence_ids"]) != sequences: raise ValueError("ordered-site order")
        ordered_available = saved["sequence_projected_available"].astype(bool)
    tie = np.empty(len(reactions), dtype=int)
    for rank_number, index in enumerate(sorted(range(len(reactions)), key=lambda value: stable(reactions[value]))): tie[index] = rank_number

    inner_rows = []
    selections = []
    outer_rows = []
    input_paths = {
        ROOT / "CLEAN_CANDIDATE_PRUNING_V1.md", Path(__file__), ROOT / "run_raw.py",
        ROOT / "dataset_02/core_edges.json", ROOT / "dataset_02/core_reactions.json",
        ROOT / "multiaxis_split_01/outer_inner_blocks.json", ROOT / "multiaxis_split_01/axis_assignments.json",
        ROOT / "main_baselines_01/fresh_representations.npz", ROOT / "ordered_site_features_01/ordered_features.npz",
        ROOT / "clean_ec_bridge_01/reaction_ec_support.json", ROOT / "clean_ec_bridge_validation_01/INDEPENDENT_QC.json",
        ROOT / "clean_core_features_01/maximum_separation_calls.json", ROOT / "clean_core_features_01/sequence_manifest.jsonl",
        ROOT / "clean_core_features_01/input_manifest.json", ROOT / "clean_core_features_validation_01/INDEPENDENT_QC.json",
        ROOT / "domain_router_01/route_decisions.jsonl", ROOT / "domain_router_validation_01/INDEPENDENT_QC.json",
    }
    router_candidates = {(row["source_task"], row["block_number"], row["sequence_sha256"], row["cell"])
                         for row in [json.loads(line) for line in (ROOT / "domain_router_01/route_decisions.jsonl").read_text().splitlines()]
                         if row["candidate_route"] != "abstain"}

    for block_number, block in enumerate(blocks):
        if not block["evaluable"]: continue
        interaction_manifest_path = ROOT / f"interaction_conditional_sequence_01/block_{block_number:03d}_inner_manifest.json"
        interaction_inner_path = ROOT / f"interaction_conditional_sequence_01/block_{block_number:03d}_inner_predictions.npz"
        input_paths.update((interaction_manifest_path, interaction_inner_path))
        interaction_manifest = read(interaction_manifest_path)
        with np.load(interaction_inner_path, allow_pickle=False) as saved:
            interaction_inner = {key: saved[key] for key in saved.files}
        joint_lambdas = interaction_manifest["joint_global_interaction_selection"]["lambdas"]
        interaction_scores = (interaction_inner["log_prior"] + joint_lambdas[0] * interaction_inner["delta_global"] +
                              joint_lambdas[1] * interaction_inner["delta_interaction"])
        interaction_row = {(row["inner_number"], row["query"]): index for index, row in enumerate(interaction_manifest["rows"])}

        block_inner = []
        for inner_number, inner in enumerate(block["inner_blocks"]):
            if not inner["evaluable"]: continue
            labels, truth = construct(inner, edges, sequence_index, reaction_index, (len(sequences), len(reactions)))
            queries = sorted(truth)
            training = np.flatnonzero(labels.any(axis=1))
            baseline_scores, baseline_domains = train_scores(labels, chemical_kernel, homology,
                                                               [sequence_index[query] for query in queries], training)
            query_position = {query: index for index, query in enumerate(queries)}
            seen = labels.any(axis=0)
            for query in queries:
                for cell, method, panel, positives in route_panels(inner["task"], query, labels, truth, homology,
                                                                   ordered_available, sequence_index):
                    if method == "joint_global_interaction":
                        row_number = interaction_row[(inner_number, query)]
                        score = interaction_scores[row_number]
                        domain = np.ones(len(reactions), dtype=bool)
                    else:
                        row_number = query_position[query]
                        score = baseline_scores[method][row_number]
                        domain = baseline_domains[method][row_number]
                    for mode in MODES:
                        metric = rank(score, panel, domain, clean_mask(mode, query, calls.get(query, []), mappings), positives, tie)
                        block_inner.append({
                            "block_number": block_number, "block_id": block["block_id"], "inner_number": inner_number,
                            "cell": cell, "method": method, "mode": mode, "sequence_sha256": query,
                            "protein_group": axes["protein_groups"][query], "target_reaction_indices": sorted(positives),
                            "clean_available": query in calls, "clean_exact_exposed": feature_manifest[query]["clean_exact_exposed"],
                            **metric,
                        })
        inner_rows.extend(block_inner)
        expected_cells = {
            "chemical_cold": (INTERACTION_CELL, TRANSPORT_CELLS["chemical_cold"]),
            "protein_cold": (MMSEQS_CELL, TRANSPORT_CELLS["protein_cold"]),
            "double_cold": (TRANSPORT_CELLS["double_cold"],),
        }[block["task"]]
        chosen_by_cell = {}
        for cell in expected_cells:
            selected, mode_scores = choose_mode([row for row in block_inner if row["cell"] == cell])
            chosen_by_cell[cell] = selected
            selections.append({
                "block_number": block_number, "block_id": block["block_id"], "cell": cell,
                "selected_mode": selected, "inner_protein_group_macro_mrr": mode_scores,
                "inner_query_panels": sum(row["mode"] == "none" for row in block_inner if row["cell"] == cell),
            })

        labels, truth = construct(block, edges, sequence_index, reaction_index, (len(sequences), len(reactions)))
        baseline_path = ROOT / f"main_baselines_01/block_{block_number:03d}_scores.npz"
        interaction_path = ROOT / f"interaction_conditional_sequence_01/block_{block_number:03d}_scores.npz"
        input_paths.update((baseline_path, interaction_path))
        with np.load(baseline_path, allow_pickle=False) as saved:
            baseline = {key: saved[key] for key in saved.files}
        with np.load(interaction_path, allow_pickle=False) as saved:
            interaction = {key: saved[key] for key in saved.files}
        baseline_query = {str(query): index for index, query in enumerate(baseline["query_ids"])}
        interaction_query = {str(query): index for index, query in enumerate(interaction["query_ids"])}
        for query in sorted(truth):
            for cell, method, panel, positives in route_panels(block["task"], query, labels, truth, homology,
                                                               ordered_available, sequence_index):
                if (block["task"], block_number, query, cell) not in router_candidates:
                    raise ValueError("router reconstruction mismatch")
                if method == "joint_global_interaction": source, query_map = interaction, interaction_query
                else: source, query_map = baseline, baseline_query
                score = source[method]
                domain = source["domain_" + method]
                if score.ndim == 2: score = score[query_map[query]]
                if domain.ndim == 2: domain = domain[query_map[query]]
                mode = chosen_by_cell.get(cell, "none")
                none_metric = rank(score, panel, domain, np.ones(len(reactions), dtype=bool), positives, tie)
                selected_metric = rank(score, panel, domain, clean_mask(mode, query, calls.get(query, []), mappings), positives, tie)
                outer_rows.append({
                    "block_number": block_number, "block_id": block["block_id"], "source_task": block["task"],
                    "cell": cell, "method": method, "selected_mode": mode, "sequence_sha256": query,
                    "protein_group": axes["protein_groups"][query], "target_reaction_indices": sorted(positives),
                    "clean_available": query in calls, "clean_exact_exposed": feature_manifest[query]["clean_exact_exposed"],
                    "candidate_count": none_metric["candidate_count"], "positive_count": none_metric["positive_count"],
                    "none_rr": none_metric["rr"], "selected_rr": selected_metric["rr"],
                    "none_covered_candidates": none_metric["covered_candidates"],
                    "selected_covered_candidates": selected_metric["covered_candidates"],
                    "none_covered_positives": none_metric["covered_positives"],
                    "selected_covered_positives": selected_metric["covered_positives"],
                    "none_hit_at_10": none_metric["hit_at_10"], "selected_hit_at_10": selected_metric["hit_at_10"],
                })
        print(f"CLEAN pruning block {block_number + 1}/{len(blocks)} {block['block_id']}", flush=True)

    comparisons = []
    acceptance = {}
    for cell in sorted({row["cell"] for row in outer_rows}):
        selected = [row for row in outer_rows if row["cell"] == cell]
        all_summary = paired_summary(selected, cell + "|all")
        exposed = paired_summary([row for row in selected if row["clean_exact_exposed"]], cell + "|exposed")
        unexposed = paired_summary([row for row in selected if not row["clean_exact_exposed"]], cell + "|unexposed")
        comparisons.append({
            "cell": cell, "all": all_summary, "clean_exact_exposed": exposed, "clean_exact_unexposed": unexposed,
            "pre_pruning_positive_instances": sum(row["positive_count"] for row in selected),
            "selected_covered_positive_instances": sum(row["selected_covered_positives"] for row in selected),
            "mean_selected_candidate_coverage": float(np.mean([row["selected_covered_candidates"] / row["candidate_count"] for row in selected])),
            "selected_mode_frequencies": {mode: sum(row["selected_mode"] == mode for row in selected) for mode in MODES},
        })
        interval = all_summary["conditional_protein_bootstrap_95ci"]
        useful = bool(interval and interval[0] > 0 and unexposed["protein_groups"] >= 2 and unexposed["difference"] > 0)
        acceptance[cell] = {
            "clean_pruning_useful": useful,
            "overall_conditional_interval_lower_above_zero": bool(interval and interval[0] > 0),
            "exact_unexposed_point_improvement": bool(unexposed["difference"] is not None and unexposed["difference"] > 0),
            "does_not_by_itself_accept_underlying_router": True,
        }

    OUTPUT.mkdir()
    for filename, rows in (("inner_mode_metrics.jsonl", inner_rows), ("mode_selections.jsonl", selections),
                           ("outer_paired_records.jsonl", outer_rows)):
        with (OUTPUT / filename).open("w") as handle:
            for row in rows: handle.write(json.dumps(row, allow_nan=False) + "\n")
    write_json(OUTPUT / "comparison_summary.json", {
        "status": "COMPUTED_QC_PENDING", "created_utc": now(), "comparisons": comparisons,
        "bootstrap_replicates": BOOTSTRAP_REPLICATES, "bootstrap_seed": BOOTSTRAP_SEED,
        "interval_scope": "Conditional fixed-fit protein-group bootstrap; no split, chemistry, mapping, model-refit, external or multiplicity uncertainty",
    })
    write_json(OUTPUT / "acceptance.json", {
        "status": "COMPUTED_QC_PENDING", "by_cell": acceptance,
        "old_integrated_label_table_used": False, "outer_information_selected_pruning_mode": False,
        "uncovered_positive_removed_from_denominator": False, "cross_domain_score_fusion_performed": False,
        "underlying_router_routes_previously_accepted": 0, "independent_biological_validation": False,
    })
    write_json(OUTPUT / "input_manifest.json", {str(path.resolve()): digest_file(path) for path in sorted(input_paths)})
    print(json.dumps({"status": "COMPUTED_QC_PENDING", "inner_mode_rows": len(inner_rows),
                      "outer_paired_rows": len(outer_rows), "selections": len(selections),
                      "acceptance": acceptance}, indent=2))


if __name__ == "__main__": main()
