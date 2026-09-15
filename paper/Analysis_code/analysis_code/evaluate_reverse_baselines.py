"""Evaluate species-matched reverse retrieval on frozen source-observed panels."""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "reverse_baselines_01"
PERMUTATIONS = 99
BOOTSTRAPS = 5000
SEED = 20260915
SEEN_METHODS = ("mmseqs_nearest_positive", "mmseqs_weighted_positive", "blast_nearest_positive", "blast_weighted_positive")
UNSEEN_METHODS = ("mmseqs_homology_chemical_transport", "blast_homology_chemical_transport")
BASELINES = ("training_label_count_prior", "uniform_expectation")


def read(path): return json.loads(path.read_text())
def jsonl(path): return [json.loads(line) for line in path.read_text().splitlines()]


def alignment_matrix(path, candidates, core, spec):
    qi = {value: index for index, value in enumerate(candidates)}
    ti = {value: index for index, value in enumerate(core)}
    matrix = np.zeros((len(candidates), len(core)), dtype=np.float32)
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            query, target = row[:2]
            identity, qcov, tcov, length, evalue, bits = map(float, row[2:])
            if evalue <= spec["homology_evalue_maximum"] and min(qcov, tcov) >= spec["homology_minimum_bilateral_coverage"]:
                matrix[qi[query], ti[target]] = max(matrix[qi[query], ti[target]], bits)
    return matrix


def rank(score, domain, candidates, positives, tie):
    selected = np.flatnonzero(domain)
    order = selected[np.lexsort((tie[selected], -score[selected]))]
    positions = np.full(len(candidates), np.inf); positions[order] = np.arange(1, len(order) + 1)
    target = positions[np.asarray(positives, dtype=int)]
    best = float(np.min(target))
    result = {"rr": 0.0 if not math.isfinite(best) else 1 / best,
              "covered_candidates": len(selected), "covered_positives": int(np.isfinite(target).sum())}
    for cutoff in (1, 5, 10):
        result[f"recall_at_{cutoff}"] = float(np.mean(target <= cutoff))
        result[f"hit_at_{cutoff}"] = float(np.any(target <= cutoff))
    return result


def random_metrics(n, k):
    probability = k / n; rr = 0.0
    for position in range(1, n - k + 2):
        rr += probability / position
        if position < n - k + 1: probability *= (n - position - k + 1) / (n - position)
    result = {"rr": rr, "covered_candidates": n, "covered_positives": k}
    for cutoff in (1, 5, 10):
        top = min(cutoff, n); none = 1.0
        for offset in range(top): none *= max(0, n - k - offset) / (n - offset)
        result[f"recall_at_{cutoff}"] = top / n
        result[f"hit_at_{cutoff}"] = 1 - none
    return result


def taxid_values(records, key):
    grouped = defaultdict(list)
    for row in records: grouped[row["taxid"]].append(float(row[key]))
    return {taxid: float(np.mean(values)) for taxid, values in grouped.items()}


def paired(records, left, right, label):
    index = {(row["block_number"], row["taxid"], row["reaction_key"], row["positive_state"], row["method"]): row for row in records}
    left_rows = [row for row in records if row["method"] == left]
    by_taxid = defaultdict(list); pairs = []
    for row in left_rows:
        key = (row["block_number"], row["taxid"], row["reaction_key"], row["positive_state"], right)
        other = index[key]
        if row["candidate_count"] != other["candidate_count"] or row["positive_sequence_sha256"] != other["positive_sequence_sha256"]: raise ValueError("paired denominator")
        delta = row["rr"] - other["rr"]
        by_taxid[row["taxid"]].append((row["rr"], other["rr"], delta))
        pairs.append({"block_number": row["block_number"], "taxid": row["taxid"], "reaction_key": row["reaction_key"],
                      "positive_state": row["positive_state"], "left_method": left, "right_method": right,
                      "left_rr": row["rr"], "right_rr": other["rr"], "difference": delta,
                      "candidate_count": row["candidate_count"], "positive_count": row["positive_count"]})
    values = np.asarray([np.mean(items, axis=0) for taxid, items in sorted(by_taxid.items())])
    interval = None
    if len(values) >= 2:
        token = int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "little")
        rng = np.random.default_rng((SEED + token) % (2**63 - 1))
        draws = rng.integers(len(values), size=(BOOTSTRAPS, len(values)))
        interval = [float(value) for value in np.quantile(values[draws, 2].mean(axis=1), (0.025, 0.975))]
    return {"query_panels": len(left_rows), "taxids": len(values), "left_mrr": float(values[:, 0].mean()),
            "right_mrr": float(values[:, 1].mean()), "difference": float(values[:, 2].mean()),
            "taxid_bootstrap_95ci": interval}, pairs


def permutation(candidates, taxid, replicate):
    token = int.from_bytes(hashlib.sha256(f"{SEED}|{taxid}|{replicate}".encode()).digest()[:8], "little")
    return np.random.default_rng(token).permutation(len(candidates))


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    for gate in (ROOT / "reverse_candidate_census_validation_01/INDEPENDENT_QC.json",
                 ROOT / "reverse_alignment_validation_01/INDEPENDENT_QC.json",
                 ROOT / "main_baseline_validation_01/INDEPENDENT_METRIC_QC.json"):
        if read(gate)["status"] != "PASS": raise ValueError("unpassed gate: " + str(gate))
    edges = read(ROOT / "dataset_02/core_edges.json")
    reactions = sorted(read(ROOT / "dataset_02/core_reactions.json"))
    core = sorted({edge["sequence_sha256"] for edge in edges})
    candidates = sorted(row["sequence_sha256"] for row in jsonl(ROOT / "reverse_candidate_census_01/candidate_records.jsonl"))
    panels = read(ROOT / "reverse_candidate_census_01/taxid_candidate_panels.json")
    queries = jsonl(ROOT / "reverse_candidate_census_01/query_panels.jsonl")
    ci = {value: index for index, value in enumerate(candidates)}
    si = {value: index for index, value in enumerate(core)}
    ri = {value: index for index, value in enumerate(reactions)}
    spec = read(ROOT / "main_baseline_spec_v1.json")
    mmseqs = alignment_matrix(ROOT / "reverse_alignment_01/mmseqs_normalized.tsv", candidates, core, spec)
    blast = alignment_matrix(ROOT / "reverse_alignment_01/blast_normalized.tsv", candidates, core, spec)
    with np.load(ROOT / "main_baselines_01/fresh_representations.npz", allow_pickle=False) as saved:
        if list(saved["reaction_ids"]) != reactions: raise ValueError("chemical kernel order")
        chemical_kernel = saved["chemical_kernel"].astype(np.float64)
    blocks = read(ROOT / "multiaxis_split_01/outer_inner_blocks.json")
    tie = np.empty(len(candidates), dtype=int)
    for position, index in enumerate(sorted(range(len(candidates)), key=lambda value: stable(candidates[value]))): tie[index] = position

    metric_rows = []; permutation_values = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for block_number, block in enumerate(blocks):
        if not block["evaluable"]: continue
        labels = np.zeros((len(core), len(reactions)), dtype=np.float32)
        for edge_number in block["train_edge_indices"]:
            edge = edges[edge_number]; labels[si[edge["sequence_sha256"]], ri[edge["reaction_key"]]] = 1
        training_indices = np.flatnonzero(labels.any(axis=1))
        balanced = labels / np.maximum(labels.sum(axis=1, keepdims=True), 1)
        represented = set(core[index] for index in training_indices)
        label_count_all = np.asarray([labels[si[value]].sum() if value in si else 0 for value in candidates], dtype=np.float64)
        for query in (row for row in queries if row["block_number"] == block_number):
            taxid_candidates = panels[query["taxid"]]
            candidate_indices = np.asarray([ci[value] for value in taxid_candidates])
            reaction_number = ri[query["reaction_key"]]
            seen = bool(labels[:, reaction_number].any())
            subsets = {
                "represented_positive": [value for value in query["positive_sequence_sha256"] if value in represented],
                "new_positive": [value for value in query["positive_sequence_sha256"] if value not in represented],
            }
            matrix_by_name = {"mmseqs": mmseqs[candidate_indices], "blast": blast[candidate_indices]}
            scores = {"training_label_count_prior": label_count_all[candidate_indices]}
            domains = {"training_label_count_prior": np.ones(len(candidate_indices), dtype=bool)}
            if seen:
                positive_training = np.flatnonzero(labels[:, reaction_number] > 0)
                for name, matrix in matrix_by_name.items():
                    values = matrix[:, positive_training]
                    scores[name + "_nearest_positive"] = values.max(axis=1)
                    scores[name + "_weighted_positive"] = values.sum(axis=1)
                    domains[name + "_nearest_positive"] = scores[name + "_nearest_positive"] > 0
                    domains[name + "_weighted_positive"] = scores[name + "_weighted_positive"] > 0
                methods = SEEN_METHODS
            else:
                for name, matrix in matrix_by_name.items():
                    relevant = matrix[:, training_indices]
                    weights = relevant.sum(axis=1, keepdims=True)
                    scores[name + "_homology_chemical_transport"] = ((relevant @ balanced[training_indices]) / np.maximum(weights, 1)) @ chemical_kernel[:, reaction_number]
                    domains[name + "_homology_chemical_transport"] = np.broadcast_to(weights[:, 0] > 0, (len(candidate_indices),)).copy()
                methods = UNSEEN_METHODS
            for positive_state, positive_ids in subsets.items():
                if not positive_ids: continue
                positive_positions = [taxid_candidates.index(value) for value in positive_ids]
                cell = ("seen" if seen else "unseen") + "_reaction__" + positive_state
                method_names = (*BASELINES, *methods)
                for method in method_names:
                    if method == "uniform_expectation": result = random_metrics(len(taxid_candidates), len(positive_ids))
                    else: result = rank(scores[method], domains[method], taxid_candidates, positive_positions, tie[candidate_indices])
                    metric_rows.append({
                        "cell": cell, "reaction_seen_in_training": seen, "positive_state": positive_state,
                        "block_number": block_number, "block_id": block["block_id"], "taxid": query["taxid"],
                        "reaction_key": query["reaction_key"], "method": method,
                        "positive_sequence_sha256": sorted(positive_ids), "positive_count": len(positive_ids),
                        "candidate_count": len(taxid_candidates), **result,
                    })
                for method in methods:
                    for replicate in range(PERMUTATIONS):
                        order = permutation(taxid_candidates, query["taxid"], replicate)
                        result = rank(scores[method][order], domains[method][order], taxid_candidates, positive_positions, tie[candidate_indices])
                        permutation_values[(cell, method)][replicate][query["taxid"]].append(result["rr"])
        print(f"reverse block {block_number + 1}/{len(blocks)} {block['block_id']}", flush=True)

    summaries = []
    for cell in sorted({row["cell"] for row in metric_rows}):
        chosen = [row for row in metric_rows if row["cell"] == cell]
        for method in sorted({row["method"] for row in chosen}):
            rows_method = [row for row in chosen if row["method"] == method]
            values = taxid_values(rows_method, "rr")
            summaries.append({
                "cell": cell, "method": method, "query_panels": len(rows_method), "taxids": len(values),
                "taxid_macro_mrr": float(np.mean(list(values.values()))),
                "documented_positive_coverage": sum(row["covered_positives"] for row in rows_method) / sum(row["positive_count"] for row in rows_method),
                "mean_candidate_coverage": float(np.mean([row["covered_candidates"] / row["candidate_count"] for row in rows_method])),
            })
    comparisons = []; pair_rows = []
    comparison_index = {}
    for cell in sorted({row["cell"] for row in metric_rows}):
        methods = SEEN_METHODS if cell.startswith("seen_") else UNSEEN_METHODS
        for left in methods:
            for right in BASELINES:
                result, pairs = paired([row for row in metric_rows if row["cell"] == cell], left, right, f"{cell}|{left}|{right}")
                item = {"cell": cell, "left_method": left, "right_method": right, **result}
                comparisons.append(item); pair_rows.extend(pairs); comparison_index[(cell, left, right)] = item

    controls = []
    for (cell, method), by_replicate in sorted(permutation_values.items()):
        observed = next(row["taxid_macro_mrr"] for row in summaries if row["cell"] == cell and row["method"] == method)
        replicate_mrr = []
        for replicate in range(PERMUTATIONS):
            replicate_mrr.append(float(np.mean([np.mean(values) for values in by_replicate[replicate].values()])))
        pvalue = (1 + sum(value >= observed for value in replicate_mrr)) / (PERMUTATIONS + 1)
        controls.append({"cell": cell, "method": method, "observed_taxid_macro_mrr": observed,
                         "permutation_taxid_macro_mrr": replicate_mrr, "exact_upper_tail_p": pvalue})

    acceptance = {}
    summary_index = {(row["cell"], row["method"]): row for row in summaries}
    control_index = {(row["cell"], row["method"]): row for row in controls}
    for cell in sorted({row["cell"] for row in metric_rows}):
        primary = "mmseqs_weighted_positive" if cell.startswith("seen_") else "mmseqs_homology_chemical_transport"
        versus = [comparison_index[(cell, primary, baseline)] for baseline in BASELINES]
        coverage = summary_index[(cell, primary)]["documented_positive_coverage"]
        taxids = summary_index[(cell, primary)]["taxids"]
        gates = {
            "positive_vs_training_label_count_prior": bool(versus[0]["taxid_bootstrap_95ci"] and versus[0]["taxid_bootstrap_95ci"][0] > 0),
            "positive_vs_uniform_expectation": bool(versus[1]["taxid_bootstrap_95ci"] and versus[1]["taxid_bootstrap_95ci"][0] > 0),
            "at_least_20_taxids": taxids >= 20,
            "documented_positive_coverage_at_least_0_8": coverage >= 0.8,
            "protein_permutation_exact_p_at_most_0_05": control_index[(cell, primary)]["exact_upper_tail_p"] <= 0.05,
        }
        acceptance[cell] = {"primary_method": primary, "accepted": all(gates.values()), "gates": gates,
                            "taxids": taxids, "documented_positive_coverage": coverage,
                            "permutation_exact_p": control_index[(cell, primary)]["exact_upper_tail_p"]}

    OUTPUT.mkdir()
    for name, values in (("per_query_metrics.jsonl", metric_rows), ("paired_query_records.jsonl", pair_rows)):
        with (OUTPUT / name).open("w") as handle:
            for row in values: handle.write(json.dumps(row, allow_nan=False) + "\n")
    write_json(OUTPUT / "method_summary.json", {"status": "COMPUTED_QC_PENDING", "created_utc": now(), "summaries": summaries})
    write_json(OUTPUT / "comparison_summary.json", {"status": "COMPUTED_QC_PENDING", "created_utc": now(), "comparisons": comparisons,
                                                     "bootstrap_replicates": BOOTSTRAPS, "bootstrap_unit": "taxid"})
    write_json(OUTPUT / "permutation_controls.json", {"status": "COMPUTED_QC_PENDING", "created_utc": now(), "controls": controls,
                                                       "permutations": PERMUTATIONS, "permutation_unit": "candidate_scores_within_taxid"})
    write_json(OUTPUT / "acceptance.json", {"status": "COMPUTED_QC_PENDING", "by_cell": acceptance,
                                             "source_observed_candidate_panel_not_complete_proteome": True,
                                             "independent_biological_validation": False})
    inputs = [ROOT / "REVERSE_RETRIEVAL_V1.md", ROOT / "REVERSE_RETRIEVAL_ACCEPTANCE_V1.md", Path(__file__), ROOT / "run_raw.py",
              ROOT / "reverse_candidate_census_01/candidate_records.jsonl", ROOT / "reverse_candidate_census_01/query_panels.jsonl",
              ROOT / "reverse_candidate_census_01/taxid_candidate_panels.json", ROOT / "reverse_candidate_census_validation_01/INDEPENDENT_QC.json",
              ROOT / "reverse_alignment_01/mmseqs_normalized.tsv", ROOT / "reverse_alignment_01/blast_normalized.tsv",
              ROOT / "reverse_alignment_validation_01/INDEPENDENT_QC.json", ROOT / "dataset_02/core_edges.json",
              ROOT / "dataset_02/core_reactions.json", ROOT / "multiaxis_split_01/outer_inner_blocks.json",
              ROOT / "main_baseline_spec_v1.json", ROOT / "main_baselines_01/fresh_representations.npz",
              ROOT / "main_baseline_validation_01/INDEPENDENT_METRIC_QC.json"]
    write_json(OUTPUT / "input_manifest.json", {str(path.resolve()): digest_file(path) for path in inputs})
    print(json.dumps({"status": "COMPUTED_QC_PENDING", "metric_rows": len(metric_rows), "paired_rows": len(pair_rows),
                      "acceptance": acceptance}, indent=2))


if __name__ == "__main__": main()
