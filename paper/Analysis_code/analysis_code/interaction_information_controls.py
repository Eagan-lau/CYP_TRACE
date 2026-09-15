"""Protocol-required paired contrasts to non-interaction information controls."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from paired_expert_comparisons import group_summary
from run_raw import digest_file, now, write_json


PAIRS = [
    ("interaction_residual", "missing_residual"),
    ("interaction_residual", "composition_residual"),
    ("interaction_residual", "lineage_residual"),
    ("joint_global_interaction", "missing_residual"),
    ("joint_global_interaction", "composition_residual"),
    ("joint_global_interaction", "lineage_residual"),
]


def run(root: Path, channel: str):
    source_dir = root / f"interaction_conditional_{channel}_01"
    output = root / "interaction_information_controls_01"; output.mkdir(exist_ok=True)
    summary_path = output / f"{channel}_summary.json"; rows_path = output / f"{channel}_rows.jsonl"
    if summary_path.exists() or rows_path.exists(): raise FileExistsError(summary_path)
    qc_path = root / f"interaction_validation_{channel}_01/INDEPENDENT_QC.json"
    if json.loads(qc_path.read_text())["status"] != "PASS": raise ValueError("interaction QC gate")
    metrics = [json.loads(line) for line in (source_dir / "per_query_metrics.jsonl").read_text().splitlines()]
    lookup = {(row["method"], row["task"], row["block_number"], row["sequence_sha256"]): row for row in metrics}
    rows = []
    for left, right in PAIRS:
        for key, left_row in lookup.items():
            if key[0] != left: continue
            _, task, block_number, sequence = key; right_row = lookup[(right, task, block_number, sequence)]
            if left_row["target_reaction_indices"] != right_row["target_reaction_indices"] or left_row["protein_group"] != right_row["protein_group"]: raise ValueError("denominator mismatch")
            rows.append({"cohort": channel, "task": task, "block_number": block_number, "sequence_sha256": sequence,
                         "protein_group": left_row["protein_group"], "left_method": left, "right_method": right,
                         "left_rr": left_row["rr"], "right_rr": right_row["rr"], "positive_count": left_row["positive_count"]})
    rng = np.random.default_rng(20260909); comparisons = []
    for task in sorted({row["task"] for row in rows}):
        for left, right in PAIRS:
            selected = [row for row in rows if (row["task"], row["left_method"], row["right_method"]) == (task, left, right)]
            comparisons.append({"task": task, "left_method": left, "right_method": right,
                                "end_to_end": group_summary(selected, "left_rr", "right_rr", rng)})
    with rows_path.open("w") as handle:
        for row in rows: handle.write(json.dumps(row, allow_nan=False) + "\n")
    write_json(summary_path, {"status": "COMPUTED_QC_PENDING", "created_utc": now(), "cohort": channel,
        "pairs": [list(pair) for pair in PAIRS], "comparisons": comparisons, "row_count": len(rows),
        "bootstrap_seed": 20260909, "bootstrap_replicates": 5000,
        "comparison_list_frozen_before_primary_readout": False,
        "protocol_required_information_contrasts_computed_post_readout": True,
        "confirmatory_hypothesis_test": False, "multiplicity_adjusted": False, "independent_biological_validation": False})
    write_json(output / f"{channel}_input_manifest.json", {
        str((source_dir / "per_query_metrics.jsonl").relative_to(root)).replace("\\", "/"): digest_file(source_dir / "per_query_metrics.jsonl"),
        str(qc_path.relative_to(root)).replace("\\", "/"): digest_file(qc_path),
        "interaction_information_controls.py": digest_file(Path(__file__)),
        "paired_expert_comparisons.py": digest_file(root / "paired_expert_comparisons.py"),
    })
    print(json.dumps({"cohort": channel, "rows": len(rows), "comparisons": len(comparisons)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    args = parser.parse_args(); run(Path(__file__).resolve().parent, args.channel)
