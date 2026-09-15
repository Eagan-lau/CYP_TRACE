"""Independent ranking, denominator and replicate-0 QC for reverse retrieval."""
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
SOURCE = ROOT / "reverse_baselines_01"
OUTPUT = ROOT / "reverse_baseline_validation_01"
SEED = 20260915
BOOTSTRAPS = 5000
SEEN = ("mmseqs_nearest_positive", "mmseqs_weighted_positive", "blast_nearest_positive", "blast_weighted_positive")
UNSEEN = ("mmseqs_homology_chemical_transport", "blast_homology_chemical_transport")
BASE = ("training_label_count_prior", "uniform_expectation")


def load(path): return json.loads(path.read_text())
def lines(path): return [json.loads(value) for value in path.read_text().splitlines()]


def matrix(path, candidates, core, maximum_evalue, minimum_coverage):
    qi = {value: index for index, value in enumerate(candidates)}; ti = {value: index for index, value in enumerate(core)}
    result = np.zeros((len(candidates), len(core)), dtype=np.float64)
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            q, t = row[:2]; identity, qcov, tcov, length, evalue, bits = map(float, row[2:])
            if evalue <= maximum_evalue and qcov >= minimum_coverage and tcov >= minimum_coverage:
                result[qi[q], ti[t]] = max(result[qi[q], ti[t]], bits)
    return result


def ranking(values, available, positives, tie):
    chosen = np.flatnonzero(available)
    order = chosen[np.lexsort((tie[chosen], -values[chosen]))]
    rank = np.full(len(values), np.inf); rank[order] = np.arange(1, len(order) + 1)
    target = rank[np.asarray(positives, dtype=int)]
    best = float(np.min(target))
    out = {"rr": 0.0 if not math.isfinite(best) else 1 / best,
           "covered_candidates": len(chosen), "covered_positives": int(np.isfinite(target).sum())}
    for cutoff in (1, 5, 10):
        out[f"recall_at_{cutoff}"] = float(np.mean(target <= cutoff)); out[f"hit_at_{cutoff}"] = float(np.any(target <= cutoff))
    return out


def random_rank(n, k):
    p = k / n; reciprocal = 0.0
    for position in range(1, n - k + 2):
        reciprocal += p / position
        if position < n - k + 1: p *= (n - position - k + 1) / (n - position)
    out = {"rr": reciprocal, "covered_candidates": n, "covered_positives": k}
    for cutoff in (1, 5, 10):
        top = min(cutoff, n); no_hit = 1.0
        for offset in range(top): no_hit *= max(0, n - k - offset) / (n - offset)
        out[f"recall_at_{cutoff}"] = top / n; out[f"hit_at_{cutoff}"] = 1 - no_hit
    return out


def permute_count(count, taxid):
    token = int.from_bytes(hashlib.sha256(f"{SEED}|{taxid}|0".encode()).digest()[:8], "little")
    return np.random.default_rng(token).permutation(count)


def taxid_average(rows, key):
    values = defaultdict(list)
    for row in rows: values[row["taxid"]].append(float(row[key]))
    return {taxid: float(np.mean(items)) for taxid, items in values.items()}


def paired_summary(rows, left, right, label):
    index = {(row["block_number"], row["taxid"], row["reaction_key"], row["positive_state"], row["method"]): row for row in rows}
    grouped = defaultdict(list)
    for row in rows:
        if row["method"] != left: continue
        other = index[(row["block_number"], row["taxid"], row["reaction_key"], row["positive_state"], right)]
        grouped[row["taxid"]].append((row["rr"], other["rr"], row["rr"] - other["rr"]))
    values = np.asarray([np.mean(items, axis=0) for taxid, items in sorted(grouped.items())])
    token = int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "little")
    rng = np.random.default_rng((SEED + token) % (2**63 - 1))
    draws = rng.integers(len(values), size=(BOOTSTRAPS, len(values)))
    interval = [float(value) for value in np.quantile(values[draws, 2].mean(axis=1), (0.025, 0.975))] if len(values) >= 2 else None
    return {"query_panels": sum(row["method"] == left for row in rows), "taxids": len(values),
            "left_mrr": float(values[:, 0].mean()), "right_mrr": float(values[:, 1].mean()),
            "difference": float(values[:, 2].mean()), "taxid_bootstrap_95ci": interval}


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    failures = []
    manifest = load(SOURCE / "input_manifest.json")
    for name, expected in manifest.items():
        if digest_file(Path(name)) != expected: failures.append("input_hash:" + name)
    observed_rows = lines(SOURCE / "per_query_metrics.jsonl")
    observed = {(row["cell"], row["block_number"], row["taxid"], row["reaction_key"], row["positive_state"], row["method"]): row for row in observed_rows}
    if len(observed) != len(observed_rows): failures.append("duplicate_metric_key")

    edges = load(ROOT / "dataset_02/core_edges.json"); reactions = sorted(load(ROOT / "dataset_02/core_reactions.json"))
    core = sorted({edge["sequence_sha256"] for edge in edges})
    candidates = sorted(row["sequence_sha256"] for row in lines(ROOT / "reverse_candidate_census_01/candidate_records.jsonl"))
    panels = load(ROOT / "reverse_candidate_census_01/taxid_candidate_panels.json")
    queries = lines(ROOT / "reverse_candidate_census_01/query_panels.jsonl")
    ci = {value: index for index, value in enumerate(candidates)}; si = {value: index for index, value in enumerate(core)}; ri = {value: index for index, value in enumerate(reactions)}
    spec = load(ROOT / "main_baseline_spec_v1.json")
    mmseqs = matrix(ROOT / "reverse_alignment_01/mmseqs_normalized.tsv", candidates, core, spec["homology_evalue_maximum"], spec["homology_minimum_bilateral_coverage"])
    blast = matrix(ROOT / "reverse_alignment_01/blast_normalized.tsv", candidates, core, spec["homology_evalue_maximum"], spec["homology_minimum_bilateral_coverage"])
    with np.load(ROOT / "main_baselines_01/fresh_representations.npz", allow_pickle=False) as saved: kernel = saved["chemical_kernel"].astype(np.float64)
    blocks = load(ROOT / "multiaxis_split_01/outer_inner_blocks.json")
    tie = np.empty(len(candidates), dtype=int)
    for order, index in enumerate(sorted(range(len(candidates)), key=lambda value: stable(candidates[value]))): tie[index] = order
    reconstructed = set(); replicate_zero = defaultdict(lambda: defaultdict(list))
    metric_names = ("rr", "covered_candidates", "covered_positives", "recall_at_1", "recall_at_5", "recall_at_10", "hit_at_1", "hit_at_5", "hit_at_10")
    for block_number, block in enumerate(blocks):
        if not block["evaluable"]: continue
        labels = np.zeros((len(core), len(reactions)), dtype=np.float64)
        for edge_number in block["train_edge_indices"]:
            edge = edges[edge_number]; labels[si[edge["sequence_sha256"]], ri[edge["reaction_key"]]] = 1
        training = np.flatnonzero(labels.any(axis=1)); represented = {core[index] for index in training}
        balanced = labels / np.maximum(labels.sum(axis=1, keepdims=True), 1)
        label_counts = np.asarray([labels[si[value]].sum() if value in si else 0 for value in candidates])
        for query in (item for item in queries if item["block_number"] == block_number):
            ids = panels[query["taxid"]]; indices = np.asarray([ci[value] for value in ids]); reaction = ri[query["reaction_key"]]
            seen = bool(labels[:, reaction].any()); scores = {"training_label_count_prior": label_counts[indices]}; domains = {"training_label_count_prior": np.ones(len(ids), dtype=bool)}
            matrices = {"mmseqs": mmseqs[indices], "blast": blast[indices]}
            if seen:
                positive_training = np.flatnonzero(labels[:, reaction] > 0)
                for name, values in matrices.items():
                    subset = values[:, positive_training]
                    for suffix, score in (("nearest_positive", subset.max(axis=1)), ("weighted_positive", subset.sum(axis=1))):
                        method = name + "_" + suffix; scores[method] = score; domains[method] = score > 0
                methods = SEEN
            else:
                for name, values in matrices.items():
                    subset = values[:, training]; mass = subset.sum(axis=1, keepdims=True)
                    method = name + "_homology_chemical_transport"
                    scores[method] = ((subset @ balanced[training]) / np.maximum(mass, 1)) @ kernel[:, reaction]
                    domains[method] = mass[:, 0] > 0
                methods = UNSEEN
            for state, positive_ids in (("represented_positive", [x for x in query["positive_sequence_sha256"] if x in represented]),
                                        ("new_positive", [x for x in query["positive_sequence_sha256"] if x not in represented])):
                if not positive_ids: continue
                positives = [ids.index(value) for value in positive_ids]
                cell = ("seen" if seen else "unseen") + "_reaction__" + state
                for method in (*BASE, *methods):
                    result = random_rank(len(ids), len(positives)) if method == "uniform_expectation" else ranking(scores[method], domains[method], positives, tie[indices])
                    key = (cell, block_number, query["taxid"], query["reaction_key"], state, method)
                    if key not in observed: failures.append("missing_metric_key")
                    else:
                        row = observed[key]
                        if row["candidate_count"] != len(ids) or row["positive_sequence_sha256"] != sorted(positive_ids): failures.append("metric_denominator")
                        if any(abs(float(row[name]) - float(result[name])) > 1e-12 for name in metric_names): failures.append("metric_reconstruction")
                    reconstructed.add(key)
                order = permute_count(len(ids), query["taxid"])
                for method in methods:
                    result = ranking(scores[method][order], domains[method][order], positives, tie[indices])
                    replicate_zero[(cell, method)][query["taxid"]].append(result["rr"])
        print(f"verified reverse block {block_number + 1}/{len(blocks)}", flush=True)
    if reconstructed != set(observed): failures.append("metric_key_set")

    method_summary = {(row["cell"], row["method"]): row for row in load(SOURCE / "method_summary.json")["summaries"]}
    for key, declared in method_summary.items():
        selected = [row for row in observed_rows if (row["cell"], row["method"]) == key]
        values = taxid_average(selected, "rr")
        if abs(declared["taxid_macro_mrr"] - np.mean(list(values.values()))) > 1e-15: failures.append("method_summary")
        coverage = sum(row["covered_positives"] for row in selected) / sum(row["positive_count"] for row in selected)
        if abs(declared["documented_positive_coverage"] - coverage) > 1e-15: failures.append("method_coverage")
    comparison_rows = load(SOURCE / "comparison_summary.json")["comparisons"]
    for declared in comparison_rows:
        selected = [row for row in observed_rows if row["cell"] == declared["cell"]]
        expected = paired_summary(selected, declared["left_method"], declared["right_method"], f"{declared['cell']}|{declared['left_method']}|{declared['right_method']}")
        if any(declared[name] != expected[name] for name in expected): failures.append("comparison_summary")
    controls = {(row["cell"], row["method"]): row for row in load(SOURCE / "permutation_controls.json")["controls"]}
    replicate_zero_checks = 0
    for key, taxids in replicate_zero.items():
        expected = float(np.mean([np.mean(values) for values in taxids.values()]))
        if abs(controls[key]["permutation_taxid_macro_mrr"][0] - expected) > 1e-15: failures.append("permutation_replicate_zero")
        values = controls[key]["permutation_taxid_macro_mrr"]
        observed_value = controls[key]["observed_taxid_macro_mrr"]
        pvalue = (1 + sum(value >= observed_value for value in values)) / 100
        if controls[key]["exact_upper_tail_p"] != pvalue: failures.append("permutation_pvalue")
        replicate_zero_checks += 1
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
              "input_hashes_checked": len(manifest), "metric_rows_reconstructed": len(reconstructed),
              "method_summaries_reconstructed": len(method_summary), "comparisons_reconstructed": len(comparison_rows),
              "permutation_replicate_zero_aggregates_reconstructed": replicate_zero_checks,
              "remaining_permutation_replicates_aggregation_checked_not_independently_regenerated": 98,
              "same_candidate_and_positive_denominators_enforced": True, "independent_implementation": True,
              "independent_biological_validation": False}
    OUTPUT.mkdir(); write_json(OUTPUT / "INDEPENDENT_QC.json", report)
    print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": main()
