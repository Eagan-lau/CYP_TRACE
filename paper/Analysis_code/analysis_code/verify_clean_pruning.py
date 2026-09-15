"""Independent denominator, selection and outer-ranking checks for CLEAN pruning."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "clean_pruning_01"
OUTPUT = ROOT / "clean_pruning_validation_01"
MODES = ("none", "open_L2", "open_L3", "open_L4", "strict_L2", "strict_L3", "strict_L4")
BOOTSTRAP_SEED = 20260914
BOOTSTRAP_REPLICATES = 5000


def load(path): return json.loads(path.read_text())
def rows(path): return [json.loads(line) for line in path.read_text().splitlines()]


def prefix_equal(a, b, width): return a.split(".")[:width] == b.split(".")[:width]


def mask_for(mode, query, predictions, mappings):
    if mode == "none" or query not in predictions: return np.ones(len(mappings), dtype=bool)
    policy, depth = mode.split("_"); width = int(depth[1:]); calls = predictions[query]
    keep = []
    for mapping in mappings:
        available = [ec for ec, supporters in mapping.items() if supporters - {query}]
        agrees = any(prefix_equal(ec, call, width) for ec in available for call in calls)
        keep.append((not available or agrees) if policy == "open" else bool(available and agrees))
    return np.asarray(keep, dtype=bool)


def first_rr(score, allowed, positives, tie):
    selected = np.flatnonzero(allowed)
    order = selected[np.lexsort((tie[selected], -score[selected]))]
    ranks = np.full(len(score), np.inf); ranks[order] = np.arange(1, len(order) + 1)
    target = ranks[np.asarray(positives, dtype=int)]
    return (0.0 if not np.isfinite(target).any() else float(1 / np.min(target))), int(np.isfinite(target).sum()), len(selected), float(np.any(target <= 10))


def group_values(records, left, right=None):
    grouped = defaultdict(lambda: defaultdict(list))
    for row in records:
        value = [row[left]] if right is None else [row[left], row[right]]
        grouped[row["protein_group"]][row["sequence_sha256"]].append(value)
    result = []
    for group in sorted(grouped):
        per_query = [np.mean(value, axis=0) for value in grouped[group].values()]
        result.append(np.mean(per_query, axis=0))
    return np.asarray(result, dtype=float)


def choose(records):
    if not records: return "none"
    scores = {}
    for mode in MODES:
        array = group_values([row for row in records if row["mode"] == mode], "rr")
        scores[mode] = float(array.mean()) if len(array) else None
    best = max(value for value in scores.values() if value is not None)
    return next(mode for mode in MODES if scores[mode] is not None and abs(scores[mode] - best) <= 1e-15)


def seed(label):
    return (BOOTSTRAP_SEED + int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "little")) % (2**63 - 1)


def summary(records, label):
    pairs = group_values(records, "selected_rr", "none_rr")
    if not len(pairs): return {"protein_groups": 0, "difference": None, "ci": None}
    difference = pairs[:, 0] - pairs[:, 1]
    interval = None
    if len(pairs) >= 2:
        rng = np.random.default_rng(seed(label))
        draws = rng.integers(len(pairs), size=(BOOTSTRAP_REPLICATES, len(pairs)))
        interval = [float(x) for x in np.quantile(difference[draws].mean(axis=1), (0.025, 0.975))]
    return {"protein_groups": len(pairs), "difference": float(difference.mean()), "ci": interval,
            "selected": float(pairs[:, 0].mean()), "none": float(pairs[:, 1].mean())}


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    failures = []
    for filename, expected in load(SOURCE / "input_manifest.json").items():
        if digest_file(Path(filename)) != expected: failures.append("input_hash:" + filename)
    inner = rows(SOURCE / "inner_mode_metrics.jsonl")
    selections = rows(SOURCE / "mode_selections.jsonl")
    outer = rows(SOURCE / "outer_paired_records.jsonl")
    selection_index = {(row["block_number"], row["cell"]): row for row in selections}
    if len(selection_index) != len(selections): failures.append("duplicate_selection")

    inner_panels = defaultdict(list)
    for row in inner:
        inner_panels[(row["block_number"], row["inner_number"], row["cell"], row["sequence_sha256"])].append(row)
    for panel in inner_panels.values():
        if tuple(row["mode"] for row in panel) != MODES: failures.append("inner_mode_order")
        if len({tuple(row["target_reaction_indices"]) for row in panel}) != 1 or len({row["candidate_count"] for row in panel}) != 1: failures.append("inner_denominator")
        if not panel[0]["clean_available"] and len({(row["rr"], row["covered_candidates"], row["covered_positives"]) for row in panel}) != 1: failures.append("unsupported_clean_changed_inner")
    for key, declared in selection_index.items():
        chosen = choose([row for row in inner if row["block_number"] == key[0] and row["cell"] == key[1]])
        if declared["selected_mode"] != chosen: failures.append("nested_selection")

    reactions = sorted(load(ROOT / "dataset_02/core_reactions.json"))
    reaction_index = {value: index for index, value in enumerate(reactions)}
    edges = load(ROOT / "dataset_02/core_edges.json")
    blocks = load(ROOT / "multiaxis_split_01/outer_inner_blocks.json")
    declared_mapping = load(ROOT / "clean_ec_bridge_01/reaction_ec_support.json")
    mappings = [{ec: set(value["supporting_sequences"]) for ec, value in declared_mapping[reaction].items() if value["in_clean_catalog"]} for reaction in reactions]
    predictions = {query: [item["ec"] for item in values] for query, values in load(ROOT / "clean_core_features_01/maximum_separation_calls.json").items()}
    tie = np.empty(len(reactions), dtype=int)
    for rank, index in enumerate(sorted(range(len(reactions)), key=lambda value: stable(reactions[value]))): tie[index] = rank
    outer_keys = set()
    score_cache = {}
    for record in outer:
        block_number = record["block_number"]
        key = (block_number, record["cell"])
        if record["selected_mode"] != selection_index[key]["selected_mode"]: failures.append("outer_mode")
        if record["candidate_count"] <= 0 or record["positive_count"] != len(record["target_reaction_indices"]): failures.append("outer_denominator")
        if record["selected_covered_candidates"] > record["none_covered_candidates"] or record["selected_covered_positives"] > record["none_covered_positives"]: failures.append("coverage_increase")
        if (record["selected_mode"] == "none" or not record["clean_available"]) and (record["selected_rr"] != record["none_rr"] or record["selected_covered_candidates"] != record["none_covered_candidates"]): failures.append("inactive_clean_changed_outer")

        block = blocks[block_number]
        train_reactions = {reaction_index[edges[index]["reaction_key"]] for index in block["train_edge_indices"]}
        panel = np.zeros(len(reactions), dtype=bool)
        if record["cell"] == "new_protein_seen_chemistry_with_homologue": panel[list(train_reactions)] = True
        else: panel[:] = True; panel[list(train_reactions)] = False
        source_name = "interaction_conditional_sequence_01" if record["method"] == "joint_global_interaction" else "main_baselines_01"
        cache_key = (source_name, block_number)
        if cache_key not in score_cache:
            with np.load(ROOT / source_name / f"block_{block_number:03d}_scores.npz", allow_pickle=False) as saved:
                score_cache[cache_key] = {name: saved[name] for name in saved.files}
        source = score_cache[cache_key]
        query_map = {str(value): index for index, value in enumerate(source["query_ids"])}
        score = source[record["method"]]; domain = source["domain_" + record["method"]]
        if score.ndim == 2: score = score[query_map[record["sequence_sha256"]]]
        if domain.ndim == 2: domain = domain[query_map[record["sequence_sha256"]]]
        base = panel & domain.astype(bool)
        selected = base & mask_for(record["selected_mode"], record["sequence_sha256"], predictions, mappings)
        none_rr, none_positive, none_candidates, none_hit10 = first_rr(score, base, record["target_reaction_indices"], tie)
        selected_rr, selected_positive, selected_candidates, selected_hit10 = first_rr(score, selected, record["target_reaction_indices"], tie)
        observed = (record["none_rr"], record["none_covered_positives"], record["none_covered_candidates"], record["none_hit_at_10"], record["selected_rr"], record["selected_covered_positives"], record["selected_covered_candidates"], record["selected_hit_at_10"])
        expected = (none_rr, none_positive, none_candidates, none_hit10, selected_rr, selected_positive, selected_candidates, selected_hit10)
        if observed != expected: failures.append("outer_ranking_reconstruction")
        outer_keys.add((record["source_task"], block_number, record["sequence_sha256"], record["cell"]))
    router_keys = {(row["source_task"], row["block_number"], row["sequence_sha256"], row["cell"]) for row in rows(ROOT / "domain_router_01/route_decisions.jsonl") if row["candidate_route"] != "abstain"}
    if outer_keys != router_keys: failures.append("router_panel_set")

    declared_comparisons = {row["cell"]: row for row in load(SOURCE / "comparison_summary.json")["comparisons"]}
    declared_acceptance = load(SOURCE / "acceptance.json")["by_cell"]
    for cell, declared in declared_comparisons.items():
        chosen = [row for row in outer if row["cell"] == cell]
        all_result = summary(chosen, cell + "|all")
        unexposed = summary([row for row in chosen if not row["clean_exact_exposed"]], cell + "|unexposed")
        if abs(declared["all"]["selected_mrr"] - all_result["selected"]) > 1e-15 or abs(declared["all"]["unpruned_mrr"] - all_result["none"]) > 1e-15 or declared["all"]["conditional_protein_bootstrap_95ci"] != all_result["ci"]: failures.append("comparison_reconstruction")
        useful = bool(all_result["ci"] and all_result["ci"][0] > 0 and unexposed["protein_groups"] >= 2 and unexposed["difference"] > 0)
        if declared_acceptance[cell]["clean_pruning_useful"] != useful: failures.append("acceptance_reconstruction")

    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
              "input_hashes_checked": len(load(SOURCE / "input_manifest.json")),
              "inner_query_panels_checked": len(inner_panels), "nested_selections_checked": len(selections),
              "outer_query_panels_reconstructed": len(outer), "comparison_cells_reconstructed": len(declared_comparisons),
              "same_pre_pruning_denominator_enforced": True, "independent_implementation": True,
              "independent_biological_validation": False}
    OUTPUT.mkdir(); write_json(OUTPUT / "INDEPENDENT_QC.json", report)
    print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": main()
