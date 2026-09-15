"""Independent arithmetic and bootstrap audit of control gain differences."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

import numpy as np

from run_raw import digest_file, now, write_json


def run(root: Path, channel: str):
    directory = root / "interaction_control_comparison_01"
    summary = json.loads((directory / f"{channel}_summary.json").read_text()); failures = []
    for name, expected in json.loads((directory / f"{channel}_input_manifest.json").read_text()).items():
        if digest_file(root / name) != expected: failures.append("input_hash:" + name)
    rows = [json.loads(line) for line in (directory / f"{channel}_rows.jsonl").read_text().splitlines()]
    for row in rows:
        if abs((row["primary_rr"] - row["global_rr"]) - row["primary_gain"]) > 1e-15: failures.append("primary_gain")
        if abs((row["control_rr"] - row["global_rr"]) - row["control_gain"]) > 1e-15: failures.append("control_gain")
        if abs((row["primary_gain"] - row["control_gain"]) - row["gain_difference"]) > 1e-15: failures.append("gain_difference")
    rng = np.random.Generator(np.random.PCG64(summary["bootstrap_seed"])); intervals = 0
    for declared in summary["summaries"]:
        selected = [row for row in rows if (row["control"], row["method"], row["task"]) == (declared["control"], declared["method"], declared["task"])]
        grouped = defaultdict(lambda: defaultdict(list))
        for row in selected: grouped[row["protein_group"]][row["sequence_sha256"]].append(row)
        values = []
        for group, queries in sorted(grouped.items()):
            values.append([math.fsum(math.fsum(row[key] for row in panels) / len(panels) for panels in queries.values()) / len(queries) for key in ("primary_gain", "control_gain")])
        expected = declared["paired_gain_comparison"]; n = len(values)
        if n != expected["protein_groups"] or len(selected) != expected["query_panels"] or len({row["sequence_sha256"] for row in selected}) != expected["unique_queries"]: failures.append("denominator")
        matrix = np.asarray(values); difference = matrix[:, 0] - matrix[:, 1]
        for key, value in (("left_mrr", matrix[:, 0].mean()), ("right_mrr", matrix[:, 1].mean()), ("difference", difference.mean())):
            if abs(float(value) - expected[key]) > 1e-12: failures.append("mean")
        if n >= 2:
            draws = rng.integers(n, size=(summary["bootstrap_replicates"], n)); bootstrap = np.sort(difference[draws].mean(axis=1)); interval = []
            for probability in (0.025, 0.975):
                position = probability * (len(bootstrap) - 1); low = math.floor(position); high = math.ceil(position)
                interval.append(float(bootstrap[low] + (position - low) * (bootstrap[high] - bootstrap[low])))
            if not np.allclose(interval, expected["conditional_protein_bootstrap_95ci"], rtol=0, atol=1e-12): failures.append("interval")
            intervals += 1
    if len(rows) != summary["row_count"]: failures.append("row_count")
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "cohort": channel,
              "failures": failures, "rows_recomputed": len(rows), "summaries_recomputed": len(summary["summaries"]),
              "intervals_recomputed": intervals, "independent_implementation": True}
    output = directory / f"{channel}_validation.json"; write_json(output, report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    args = parser.parse_args(); run(Path(__file__).resolve().parent, args.channel)
