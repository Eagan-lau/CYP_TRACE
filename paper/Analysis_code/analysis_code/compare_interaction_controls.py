"""Paired difference-in-gain comparisons against both negative controls."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from paired_expert_comparisons import group_summary
from run_raw import digest_file, now, write_json


CONTROLS = {
    "protein_position_permutation": "interaction_control_protein_position_permutation_01",
    "reaction_center_row_permutation": "interaction_control_reaction_center_row_permutation_01",
}
METHODS = ["joint_global_interaction", "interaction_residual"]


def run(root: Path, channel: str):
    output = root / "interaction_control_comparison_01"
    output.mkdir(exist_ok=True)
    summary_path = output / f"{channel}_summary.json"
    rows_path = output / f"{channel}_rows.jsonl"
    if summary_path.exists() or rows_path.exists(): raise FileExistsError(summary_path)
    primary_dir = root / f"interaction_conditional_{channel}_01"
    if json.loads((root / f"interaction_validation_{channel}_01/INDEPENDENT_QC.json").read_text())["status"] != "PASS": raise ValueError("primary QC gate")
    primary_rows = [json.loads(line) for line in (primary_dir / "per_query_metrics.jsonl").read_text().splitlines()]
    primary = {(row["method"], row["task"], row["block_number"], row["sequence_sha256"]): row for row in primary_rows}
    rows = []; input_paths = [primary_dir / "per_query_metrics.jsonl", root / f"interaction_validation_{channel}_01/INDEPENDENT_QC.json"]
    for control_name, workspace_name in CONTROLS.items():
        control_dir = root / workspace_name / f"interaction_conditional_{channel}_01"
        control_qc = root / workspace_name / f"interaction_validation_{channel}_01/INDEPENDENT_QC.json"
        if json.loads(control_qc.read_text())["status"] != "PASS": raise ValueError("control QC gate")
        control_rows = [json.loads(line) for line in (control_dir / "per_query_metrics.jsonl").read_text().splitlines()]
        control = {(row["method"], row["task"], row["block_number"], row["sequence_sha256"]): row for row in control_rows}
        input_paths.extend([control_dir / "per_query_metrics.jsonl", control_qc])
        for method in METHODS:
            for key, primary_method in primary.items():
                if key[0] != method: continue
                _, task, block_number, sequence = key
                control_method = control[key]
                global_key = ("global_residual", task, block_number, sequence)
                primary_global = primary[global_key]; control_global = control[global_key]
                if primary_method["target_reaction_indices"] != control_method["target_reaction_indices"] or primary_method["protein_group"] != control_method["protein_group"]:
                    raise ValueError("primary/control denominator mismatch")
                if primary_global["rr"] != control_global["rr"] or primary_global["target_reaction_indices"] != control_global["target_reaction_indices"]:
                    raise ValueError("shared global branch mismatch")
                primary_gain = primary_method["rr"] - primary_global["rr"]
                control_gain = control_method["rr"] - control_global["rr"]
                rows.append({
                    "cohort": channel, "control": control_name, "method": method, "task": task,
                    "block_number": block_number, "sequence_sha256": sequence,
                    "protein_group": primary_method["protein_group"],
                    "primary_rr": primary_method["rr"], "control_rr": control_method["rr"],
                    "global_rr": primary_global["rr"], "primary_gain": primary_gain,
                    "control_gain": control_gain, "gain_difference": primary_gain - control_gain,
                    "positive_count": primary_method["positive_count"],
                })
    rng = np.random.default_rng(20260909); summaries = []
    for control_name in CONTROLS:
        for method in METHODS:
            for task in sorted({row["task"] for row in rows}):
                selected = [row for row in rows if (row["control"], row["method"], row["task"]) == (control_name, method, task)]
                summaries.append({
                    "control": control_name, "method": method, "task": task,
                    "paired_gain_comparison": group_summary(selected, "primary_gain", "control_gain", rng),
                })
    with rows_path.open("w") as handle:
        for row in rows: handle.write(json.dumps(row, allow_nan=False) + "\n")
    summary = {
        "status": "COMPUTED_QC_PENDING", "created_utc": now(), "cohort": channel,
        "summaries": summaries, "row_count": len(rows), "bootstrap_seed": 20260909,
        "bootstrap_replicates": 5000, "comparison": "(primary method - shared global) - (permuted method - shared global)",
        "comparison_list_frozen_before_control_readout": False,
        "comparison_is_protocol_implied_but_computed_after_control_readout": True,
        "confirmatory_hypothesis_test": False,
        "multiplicity_adjusted": False, "independent_biological_validation": False,
    }
    write_json(summary_path, summary)
    manifest = {str(path.relative_to(root)).replace("\\", "/"): digest_file(path) for path in input_paths}
    manifest.update({"compare_interaction_controls.py": digest_file(Path(__file__)), "paired_expert_comparisons.py": digest_file(root / "paired_expert_comparisons.py")})
    write_json(output / f"{channel}_input_manifest.json", manifest)
    print(json.dumps({"cohort": channel, "rows": len(rows), "summaries": len(summaries)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    args = parser.parse_args(); run(Path(__file__).resolve().parent, args.channel)
