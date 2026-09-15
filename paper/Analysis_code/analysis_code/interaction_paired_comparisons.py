"""Frozen paired contrasts for the candidate-specific interaction model."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from paired_expert_comparisons import group_summary
from run_main_baselines import rank_result
from run_raw import write_json, digest_file, stable, now


PAIRS = [
    ("joint_global_interaction", "global_residual"),
    ("joint_global_interaction", "joint_global_ordered"),
    ("joint_global_interaction", "interaction_residual"),
    ("joint_lineage_interaction", "lineage_residual"),
    ("joint_lineage_interaction", "joint_global_interaction"),
    ("interaction_residual", "ordered_residual"),
    ("interaction_residual", "global_residual"),
    ("interaction_residual", "availability_chemistry_prior"),
]


def run(root: Path, channel: str):
    result = root / f"interaction_conditional_{channel}_01"
    output = root / f"interaction_paired_{channel}_01"
    qc = root / f"interaction_validation_{channel}_01/INDEPENDENT_QC.json"
    if output.exists(): raise FileExistsError(output)
    if json.loads(qc.read_text())["status"] != "PASS": raise ValueError("interaction QC gate")
    original = [json.loads(line) for line in (result / "per_query_metrics.jsonl").read_text().splitlines()]
    source = {(row["method"], row["task"], row["block_number"], row["sequence_sha256"]): row for row in original}
    if len(source) != len(original): raise ValueError("duplicate source metrics")
    byblock = defaultdict(list)
    for row in original: byblock[row["block_number"]].append(row)
    reaction_ids = sorted(json.loads((root / "dataset_02/core_reactions.json").read_text()))
    tie = np.empty(len(reaction_ids), dtype=int)
    for order, index in enumerate(sorted(range(len(reaction_ids)), key=lambda i: stable(reaction_ids[i]))): tie[index] = order
    transform = np.load(root / "interaction_transform_01/reaction_transform.npz", allow_pickle=False)
    mapping_warning = transform["mapping_warning"].astype(bool)
    if len(mapping_warning) != len(reaction_ids): raise ValueError("warning mask order")
    rows = []; score_files = []
    for block_number, block_rows in sorted(byblock.items()):
        score_path = result / f"block_{block_number:03d}_scores.npz"; score_files.append(score_path.relative_to(root).as_posix())
        with np.load(score_path, allow_pickle=False) as saved: data = {name: saved[name] for name in saved.files}
        if list(data["reaction_ids"]) != reaction_ids: raise ValueError("reaction order")
        for left, right in PAIRS:
            for left_row in (row for row in block_rows if row["method"] == left):
                right_row = source[(right, left_row["task"], block_number, left_row["sequence_sha256"])]
                if left_row["target_reaction_indices"] != right_row["target_reaction_indices"] or left_row["protein_group"] != right_row["protein_group"]:
                    raise ValueError("paired denominator")
                query_index = left_row["query_index_in_block"]
                if str(data["query_ids"][query_index]) != left_row["sequence_sha256"]: raise ValueError("query order")
                left_score = data[left]; right_score = data[right]
                left_score = left_score[query_index] if left_score.ndim == 2 else left_score
                right_score = right_score[query_index] if right_score.ndim == 2 else right_score
                warning_mask = ~mapping_warning
                warning_positives = [index for index in left_row["target_reaction_indices"] if warning_mask[index]]
                warning_left = rank_result(left_score, warning_mask, warning_positives, tie)["rr"] if warning_positives else None
                warning_right = rank_result(right_score, warning_mask, warning_positives, tie)["rr"] if warning_positives else None
                rows.append({
                    "cohort": channel, "task": left_row["task"], "block_number": block_number,
                    "sequence_sha256": left_row["sequence_sha256"], "protein_group": left_row["protein_group"],
                    "left_method": left, "right_method": right, "left_rr": left_row["rr"], "right_rr": right_row["rr"],
                    "full_candidate_count": len(reaction_ids), "target_positives": left_row["positive_count"],
                    "mapping_warning_excluded_candidate_count": int(warning_mask.sum()),
                    "mapping_warning_excluded_positive_count": len(warning_positives),
                    "mapping_warning_excluded_left_rr": warning_left,
                    "mapping_warning_excluded_right_rr": warning_right,
                    "mapping_warning_excluded_defined": bool(warning_positives),
                })
        print(f"{channel} interaction paired block {block_number}", flush=True)
    rng = np.random.default_rng(20260909); comparisons = []
    for task in sorted({row["task"] for row in rows}):
        for left, right in PAIRS:
            selected = [row for row in rows if (row["task"], row["left_method"], row["right_method"]) == (task, left, right)]
            sensitivity = [row for row in selected if row["mapping_warning_excluded_defined"]]
            comparisons.append({
                "task": task, "left_method": left, "right_method": right,
                "end_to_end": group_summary(selected, "left_rr", "right_rr", rng),
                "mapping_warning_excluded": group_summary(sensitivity, "mapping_warning_excluded_left_rr", "mapping_warning_excluded_right_rr", rng),
                "coverage": {
                    "original_positive_instances": sum(row["target_positives"] for row in selected),
                    "warning_excluded_positive_instances": sum(row["mapping_warning_excluded_positive_count"] for row in selected),
                    "warning_excluded_undefined_query_panels": len(selected) - len(sensitivity),
                    "warning_excluded_candidate_count": int((~mapping_warning).sum()),
                },
            })
    output.mkdir()
    with (output / "paired_query_records.jsonl").open("w") as handle:
        for row in rows: handle.write(json.dumps(row, allow_nan=False) + "\n")
    summary = {
        "status": "COMPUTED_QC_PENDING", "created_utc": now(), "cohort": channel,
        "pairs": [list(pair) for pair in PAIRS], "comparisons": comparisons,
        "paired_query_rows": len(rows), "bootstrap_seed": 20260909, "bootstrap_replicates": 5000,
        "comparison_list_frozen_before_aggregate_readout": True,
        "confirmatory_hypothesis_test": False, "independent_biological_validation": False,
        "multiplicity_adjusted": False,
        "interval_scope": "Conditional fixed-fit protein-group bootstrap; no chemical, publication, refit, external or multiplicity uncertainty",
    }
    write_json(output / "comparison_summary.json", summary)
    names = [
        "interaction_paired_comparisons.py", "paired_expert_comparisons.py", "CANDIDATE_INTERACTION_MODEL_V1.md",
        "run_main_baselines.py", "run_raw.py",
        result.relative_to(root).as_posix() + "/per_query_metrics.jsonl", qc.relative_to(root).as_posix(),
        "dataset_02/core_reactions.json", "interaction_transform_01/reaction_transform.npz",
    ] + score_files
    write_json(output / "input_manifest.json", {name: digest_file(root / name) for name in names})
    print(json.dumps({"cohort": channel, "paired_rows": len(rows), "comparisons": len(comparisons)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    args = parser.parse_args(); run(Path(__file__).resolve().parent, args.channel)
