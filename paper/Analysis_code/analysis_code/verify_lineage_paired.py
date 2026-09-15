"""Independent ranks and protein-group bootstrap for lineage contrasts."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import numpy as np
from run_raw import digest_file, stable, write_json, now


def run(root, channel):
    result = root / f"lineage_paired_{channel}_01"; source_dir = root / f"lineage_conditional_{channel}_01"
    out = root / f"lineage_paired_validation_{channel}_01"
    if out.exists(): raise FileExistsError(out)
    read = lambda path: json.loads(path.read_text()); declared = read(result / "comparison_summary.json"); failures = []
    for name, expected in read(result / "input_manifest.json").items():
        if digest_file(root / name) != expected: failures.append("source:" + name)
    rows = [json.loads(line) for line in (result / "paired_query_records.jsonl").read_text().splitlines()]
    original = [json.loads(line) for line in (source_dir / "per_query_metrics.jsonl").read_text().splitlines()]
    source = {(row["method"], row["task"], row["block_number"], row["sequence_sha256"]): row for row in original}
    rids = sorted(read(root / "dataset_02/core_reactions.json")); hashes = [stable(r) for r in rids]; byblock = defaultdict(list)
    for row in rows: byblock[row["block_number"]].append(row)
    keys = {(r["left_method"], r["right_method"], r["task"], r["block_number"], r["sequence_sha256"]) for r in rows}
    if len(keys) != len(rows): failures.append("duplicate_paired_rows")
    checked = intervals = 0
    for bn, members in sorted(byblock.items()):
        with np.load(source_dir / f"block_{bn:03d}_scores.npz", allow_pickle=False) as saved: data = {name: saved[name] for name in saved.files}
        if list(data["reaction_ids"]) != rids: failures.append("reaction_order")
        for row in members:
            values = []; masks = []; targets = []
            for side in ["left", "right"]:
                name = row[side + "_method"]; original_row = source[(name, row["task"], bn, row["sequence_sha256"])]
                if original_row["rr"] != row[side + "_rr"] or original_row["covered_positives"] != row[side + "_covered_positives"]: failures.append("source_metric")
                if original_row["protein_group"] != row["protein_group"] or original_row["positive_count"] != row["target_positives"]: failures.append("source_denominator")
                qi = original_row["query_index_in_block"]
                if str(data["query_ids"][qi]) != row["sequence_sha256"]: failures.append("query_order")
                value = data[name]; mask = data["domain_" + name]
                values.append(value[qi] if value.ndim == 2 else value); masks.append(mask[qi] if mask.ndim == 2 else mask)
                targets.append(original_row["target_reaction_indices"])
                if int(masks[-1].sum()) != row[side + "_candidate_count"]: failures.append("native_domain_count")
            if targets[0] != targets[1]: failures.append("target_disagreement")
            common = [i for i in range(len(rids)) if masks[0][i] and masks[1][i]]; commonset = set(common)
            truth = [i for i in targets[0] if i in commonset]
            if len(common) != row["common_candidate_count"] or len(truth) != row["common_positive_count"] or len(rids) != row["full_candidate_count"]: failures.append("intersection_denominator")
            if bool(truth) != row["common_reranking_defined"]: failures.append("undefined_flag")
            for side, value in zip(["left", "right"], values):
                if not truth:
                    if row["common_" + side + "_rr"] is not None: failures.append("undefined_value")
                    continue
                order = sorted(common, key=lambda i: (-float(value[i]), hashes[i])); truthset = set(truth)
                expected = 1 / next(j + 1 for j, i in enumerate(order) if i in truthset)
                if abs(expected - row["common_" + side + "_rr"]) > 1e-12: failures.append("rank")
            checked += 1
        print(f"{channel} lineage independently ranked block {bn}", flush=True)
    rng = np.random.Generator(np.random.PCG64(declared["bootstrap_seed"])); scopes = 0
    for comparison in declared["comparisons"]:
        selected = [r for r in rows if (r["task"], r["left_method"], r["right_method"]) ==
                    (comparison["task"], comparison["left_method"], comparison["right_method"])]
        coverage = [("original_target_positive_instances", "target_positives"), ("left_covered_positive_instances", "left_covered_positives"),
                    ("right_covered_positive_instances", "right_covered_positives"), ("common_positive_instances", "common_positive_count")]
        for field, key in coverage:
            if comparison["coverage"][field] != sum(r[key] for r in selected): failures.append("coverage_summary")
        if comparison["coverage"]["common_undefined_query_panels"] != sum(not r["common_reranking_defined"] for r in selected): failures.append("undefined_summary")
        for scope, left, right, chosen in [("end_to_end", "left_rr", "right_rr", selected),
                                           ("common_domain_reranking", "common_left_rr", "common_right_rr", [r for r in selected if r["common_reranking_defined"]])]:
            groups = defaultdict(lambda: defaultdict(list))
            for row in chosen: groups[row["protein_group"]][row["sequence_sha256"]].append(row)
            values = []
            for group, queries in sorted(groups.items()):
                values.append([math.fsum(math.fsum(r[key] for r in panels) / len(panels) for panels in queries.values()) / len(queries) for key in [left, right]])
            n = len(values); expected = comparison[scope]
            if n != expected["protein_groups"] or len(chosen) != expected["query_panels"] or len({r["sequence_sha256"] for r in chosen}) != expected["unique_queries"]: failures.append("summary_denominator")
            if n:
                matrix = np.array(values); delta = matrix[:, 0] - matrix[:, 1]
                for key, value in [("left_mrr", math.fsum(matrix[:, 0]) / n), ("right_mrr", math.fsum(matrix[:, 1]) / n), ("difference", math.fsum(delta) / n)]:
                    if abs(value - expected[key]) > 1e-12: failures.append("summary_mean")
            elif any(expected[key] is not None for key in ["left_mrr", "right_mrr", "difference"]): failures.append("undefined_mean")
            if n >= 2:
                draws = rng.integers(n, size=(declared["bootstrap_replicates"], n)); counts = np.zeros((len(draws), n), dtype=int)
                np.add.at(counts, (np.repeat(np.arange(len(draws)), n), draws.ravel()), 1); boot = np.sort(counts @ delta / n); interval = []
                for probability in [.025, .975]:
                    index = probability * (len(boot) - 1); low = math.floor(index); high = math.ceil(index)
                    interval.append(boot[low] + (index - low) * (boot[high] - boot[low]))
                if not np.allclose(interval, expected["conditional_protein_bootstrap_95ci"], rtol=0, atol=1e-12): failures.append("interval")
                intervals += 1
            elif expected["conditional_protein_bootstrap_95ci"] is not None: failures.append("undefined_interval")
            scopes += 1
    if checked != declared["paired_query_rows"] or len(declared["comparisons"]) != 40: failures.append("total_count")
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "cohort": channel, "failures": failures,
              "paired_rows_recomputed": checked, "summary_scopes_recomputed": scopes, "conditional_intervals_recomputed": intervals,
              "verifier_sha256": digest_file(Path(__file__)), "independent_biological_validation": False, "multiplicity_adjusted": False}
    out.mkdir(); write_json(out / "INDEPENDENT_QC.json", report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    run(Path(__file__).resolve().parent, parser.parse_args().channel)
