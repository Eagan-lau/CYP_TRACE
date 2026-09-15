"""Independent ranking and bootstrap audit for interaction contrasts."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

import numpy as np

from run_raw import digest_file, stable, write_json, now


def run(root: Path, channel: str):
    result = root / f"interaction_paired_{channel}_01"
    source_dir = root / f"interaction_conditional_{channel}_01"
    output = root / f"interaction_paired_validation_{channel}_01"
    if output.exists(): raise FileExistsError(output)
    read = lambda path: json.loads(Path(path).read_text())
    declared = read(result / "comparison_summary.json"); failures = []
    for name, expected in read(result / "input_manifest.json").items():
        if digest_file(root / name) != expected: failures.append("input_hash:" + name)
    rows = [json.loads(line) for line in (result / "paired_query_records.jsonl").read_text().splitlines()]
    original = [json.loads(line) for line in (source_dir / "per_query_metrics.jsonl").read_text().splitlines()]
    source = {(row["method"], row["task"], row["block_number"], row["sequence_sha256"]): row for row in original}
    reaction_ids = sorted(read(root / "dataset_02/core_reactions.json")); hashes = [stable(value) for value in reaction_ids]
    warning = np.load(root / "interaction_transform_01/reaction_transform.npz", allow_pickle=False)["mapping_warning"].astype(bool)
    warning_candidates = np.flatnonzero(~warning); warning_set = set(warning_candidates.tolist())
    byblock = defaultdict(list)
    for row in rows: byblock[row["block_number"]].append(row)
    unique = {(row["left_method"], row["right_method"], row["task"], row["block_number"], row["sequence_sha256"]) for row in rows}
    if len(unique) != len(rows): failures.append("duplicate_rows")
    ranked = 0
    for block_number, members in sorted(byblock.items()):
        with np.load(source_dir / f"block_{block_number:03d}_scores.npz", allow_pickle=False) as saved:
            data = {name: saved[name] for name in saved.files}
        if list(data["reaction_ids"]) != reaction_ids: failures.append("reaction_order")
        for row in members:
            scores = []
            targets = []
            for side in ("left", "right"):
                method = row[side + "_method"]
                source_row = source[(method, row["task"], block_number, row["sequence_sha256"])]
                if source_row["rr"] != row[side + "_rr"] or source_row["positive_count"] != row["target_positives"]: failures.append("source_metric")
                query_index = source_row["query_index_in_block"]
                if str(data["query_ids"][query_index]) != row["sequence_sha256"]: failures.append("query_order")
                score = data[method]; scores.append(score[query_index] if score.ndim == 2 else score)
                targets.append(source_row["target_reaction_indices"])
            if targets[0] != targets[1]: failures.append("target_disagreement")
            positives = [index for index in targets[0] if index in warning_set]
            if len(warning_candidates) != row["mapping_warning_excluded_candidate_count"] or len(positives) != row["mapping_warning_excluded_positive_count"]: failures.append("warning_denominator")
            if bool(positives) != row["mapping_warning_excluded_defined"]: failures.append("warning_defined")
            for side, score in zip(("left", "right"), scores):
                saved_value = row["mapping_warning_excluded_" + side + "_rr"]
                if not positives:
                    if saved_value is not None: failures.append("warning_undefined_value")
                    continue
                order = sorted(warning_candidates, key=lambda index: (-float(score[index]), hashes[index]))
                positive_set = set(positives)
                expected = 1 / next(rank + 1 for rank, index in enumerate(order) if index in positive_set)
                if abs(expected - saved_value) > 1e-12: failures.append("warning_rank")
            ranked += 1
        print(f"{channel} interaction paired QA block {block_number}", flush=True)
    rng = np.random.Generator(np.random.PCG64(declared["bootstrap_seed"])); scopes = intervals = 0
    for comparison in declared["comparisons"]:
        selected = [row for row in rows if (row["task"], row["left_method"], row["right_method"]) == (comparison["task"], comparison["left_method"], comparison["right_method"])]
        sensitivity = [row for row in selected if row["mapping_warning_excluded_defined"]]
        expected_coverage = comparison["coverage"]
        if expected_coverage["original_positive_instances"] != sum(row["target_positives"] for row in selected): failures.append("coverage_original")
        if expected_coverage["warning_excluded_positive_instances"] != sum(row["mapping_warning_excluded_positive_count"] for row in selected): failures.append("coverage_warning")
        if expected_coverage["warning_excluded_undefined_query_panels"] != len(selected) - len(sensitivity): failures.append("coverage_undefined")
        for scope, left, right, chosen in [
            ("end_to_end", "left_rr", "right_rr", selected),
            ("mapping_warning_excluded", "mapping_warning_excluded_left_rr", "mapping_warning_excluded_right_rr", sensitivity),
        ]:
            grouped = defaultdict(lambda: defaultdict(list))
            for row in chosen: grouped[row["protein_group"]][row["sequence_sha256"]].append(row)
            values = []
            for group, queries in sorted(grouped.items()):
                values.append([math.fsum(math.fsum(row[key] for row in panels) / len(panels) for panels in queries.values()) / len(queries) for key in (left, right)])
            expected = comparison[scope]; n = len(values)
            if n != expected["protein_groups"] or len(chosen) != expected["query_panels"] or len({row["sequence_sha256"] for row in chosen}) != expected["unique_queries"]: failures.append("scope_denominator")
            if n:
                matrix = np.asarray(values); delta = matrix[:, 0] - matrix[:, 1]
                for key, value in (("left_mrr", matrix[:, 0].mean()), ("right_mrr", matrix[:, 1].mean()), ("difference", delta.mean())):
                    if abs(float(value) - expected[key]) > 1e-12: failures.append("scope_mean")
            elif any(expected[key] is not None for key in ("left_mrr", "right_mrr", "difference")): failures.append("scope_undefined")
            if n >= 2:
                draws = rng.integers(n, size=(declared["bootstrap_replicates"], n)); bootstrap = np.sort(delta[draws].mean(axis=1))
                interval = []
                for probability in (0.025, 0.975):
                    position = probability * (len(bootstrap) - 1); low = math.floor(position); high = math.ceil(position)
                    interval.append(float(bootstrap[low] + (position - low) * (bootstrap[high] - bootstrap[low])))
                if not np.allclose(interval, expected["conditional_protein_bootstrap_95ci"], rtol=0, atol=1e-12): failures.append("scope_interval")
                intervals += 1
            elif expected["conditional_protein_bootstrap_95ci"] is not None: failures.append("scope_interval_undefined")
            scopes += 1
    if ranked != declared["paired_query_rows"]: failures.append("row_count")
    report = {
        "status": "PASS" if not failures else "FAIL", "created_utc": now(), "cohort": channel,
        "failures": failures, "paired_rows_recomputed": ranked, "summary_scopes_recomputed": scopes,
        "conditional_intervals_recomputed": intervals, "independent_implementation": True,
        "independent_biological_validation": False, "multiplicity_adjusted": False,
    }
    output.mkdir(); write_json(output / "INDEPENDENT_QC.json", report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    args = parser.parse_args(); run(Path(__file__).resolve().parent, args.channel)
