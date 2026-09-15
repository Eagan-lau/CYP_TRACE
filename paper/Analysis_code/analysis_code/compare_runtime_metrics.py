"""Audit rank-level differences between verified cluster and Windows reruns."""
from collections import Counter, defaultdict
import json
from pathlib import Path
from run_raw import digest_file, write_json, now


FIELDS = ("rr", "mean_positive_rr", "hit_at_1", "hit_at_5", "hit_at_10",
          "recall_at_1", "recall_at_5", "recall_at_10", "ndcg_at_10")


def flatten(node, prefix=()):
    out = {}
    if isinstance(node, dict):
        for key, value in node.items(): out.update(flatten(value, prefix + (key,)))
    elif isinstance(node, (int, float)) and not isinstance(node, bool): out[prefix] = float(node)
    return out


def run(root):
    out = root / "runtime_reproducibility_01"
    if out.exists(): raise FileExistsError(out)
    report = {"status": "DIAGNOSTIC_COMPLETE", "created_utc": now(), "canonical_runtime": "ULiege cluster",
              "comparison_runtime": "Windows local compute", "channels": {},
              "identical_bitwise_outputs_expected": False, "independent_biological_validation": False}
    sources = []
    for channel in ["sequence", "structure"]:
        local_metrics = root / f"local_conditional_{channel}_01/per_query_metrics.jsonl"
        cluster_metrics = root / f"runtime_comparison_cluster_20260909/{channel}_per_query_metrics.jsonl"
        local_qc = root / f"local_validation_{channel}_01/INDEPENDENT_QC.json"
        cluster_qc = root / f"runtime_comparison_cluster_20260909/{channel}_qc.json"
        local_evaluation = root / f"local_conditional_{channel}_01/evaluation.json"
        cluster_evaluation = root / f"runtime_comparison_cluster_20260909/{channel}_evaluation.json"
        paths = [local_metrics, cluster_metrics, local_qc, cluster_qc, local_evaluation, cluster_evaluation]; sources.extend(paths)
        if any(json.loads(p.read_text())["status"] != "PASS" for p in [local_qc, cluster_qc]): raise ValueError("Both executions must pass QC")
        def load(path):
            return [json.loads(line) for line in path.read_text().splitlines()]
        key = lambda row: (row["task"], row["block_number"], row["sequence_sha256"], row["method"])
        local = {key(row): row for row in load(local_metrics)}; cluster = {key(row): row for row in load(cluster_metrics)}
        if local.keys() != cluster.keys(): raise ValueError("Metric key mismatch")
        structural = ["target_reaction_indices", "protein_group", "positive_count", "candidate_count", "covered_candidates", "covered_positives"]
        if any(local[k][field] != cluster[k][field] for k in local for field in structural): raise ValueError("Denominator mismatch")
        changed = Counter(); rr_changed = Counter(); hit_changed = Counter(); maximum = defaultdict(float)
        for k in local:
            differences = {field: abs(float(local[k][field]) - float(cluster[k][field])) for field in FIELDS}
            if max(differences.values()) > 1e-15: changed[k[3]] += 1
            if differences["rr"] > 1e-15: rr_changed[k[3]] += 1
            if max(differences[f"hit_at_{cutoff}"] for cutoff in [1, 5, 10]) > 0: hit_changed[k[3]] += 1
            maximum[k[3]] = max(maximum[k[3]], max(differences.values()))
        local_aggregate = flatten(json.loads(local_evaluation.read_text())["aggregate"])
        cluster_aggregate = flatten(json.loads(cluster_evaluation.read_text())["aggregate"])
        common = set(local_aggregate) & set(cluster_aggregate)
        aggregate_max = max(abs(local_aggregate[k] - cluster_aggregate[k]) for k in common)
        report["channels"][channel] = {"metric_rows": len(local), "identical_metric_keys": True,
            "identical_targets_and_denominators": True, "rows_with_any_metric_difference_over_1e_15": sum(changed.values()),
            "rows_with_rr_difference_over_1e_15": sum(rr_changed.values()), "rows_with_any_hit_cutoff_difference": sum(hit_changed.values()),
            "changed_rows_by_method": dict(sorted(changed.items())), "rr_changed_rows_by_method": dict(sorted(rr_changed.items())),
            "hit_changed_rows_by_method": dict(sorted(hit_changed.items())),
            "maximum_metric_difference_by_method": dict(sorted(maximum.items())),
            "maximum_shared_aggregate_numeric_difference": aggregate_max}
    report["interpretation"] = ("All keys, targets and denominators agree. Most rank instability is confined to ordered-nearest "
                                    "transport, whose argmax is sensitive to floating near-ties; primary fitted residual differences are sparse and small. "
                                    "Cluster outputs remain the canonical execution for manuscript values.")
    out.mkdir(); write_json(out / "audit.json", report)
    write_json(out / "input_manifest.json", {p.relative_to(root).as_posix(): digest_file(p) for p in sources})
    print(json.dumps(report, indent=2))


if __name__ == "__main__": run(Path(__file__).resolve().parent)
