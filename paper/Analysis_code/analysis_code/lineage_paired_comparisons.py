"""Matched-cohort taxonomy versus ordered-site exploratory contrasts."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import numpy as np
from paired_expert_comparisons import group_summary
from run_main_baselines import rank_result
from run_raw import write_json, digest_file, stable, now


PAIRS = [
    ("ordered_residual", "lineage_residual"),
    ("ordered_residual", "missing_residual"),
    ("ordered_residual", "global_residual"),
    ("lineage_residual", "global_residual"),
    ("joint_lineage_ordered", "lineage_residual"),
    ("joint_lineage_ordered", "ordered_residual"),
    ("joint_lineage_ordered", "global_residual"),
    ("joint_lineage_ordered", "availability_chemistry_prior"),
]


def run(root, channel):
    result = root / f"lineage_conditional_{channel}_01"; out = root / f"lineage_paired_{channel}_01"
    if out.exists(): raise FileExistsError(out)
    read = lambda path: json.loads(path.read_text())
    qc = root / f"lineage_validation_{channel}_01/INDEPENDENT_QC.json"
    if read(qc)["status"] != "PASS": raise ValueError("Lineage model QC gate")
    original = [json.loads(line) for line in (result / "per_query_metrics.jsonl").read_text().splitlines()]
    records = {(row["method"], row["task"], row["block_number"], row["sequence_sha256"]): row for row in original}
    if len(records) != len(original): raise ValueError("Duplicate records")
    byblock = defaultdict(list)
    for row in original: byblock[row["block_number"]].append(row)
    rids = sorted(read(root / "dataset_02/core_reactions.json")); tie = np.empty(len(rids), int)
    for j, i in enumerate(sorted(range(len(rids)), key=lambda i: stable(rids[i]))): tie[i] = j
    rows = []; scorefiles = []
    for bn, blockrows in sorted(byblock.items()):
        file = result / f"block_{bn:03d}_scores.npz"; scorefiles.append(file.relative_to(root).as_posix())
        with np.load(file, allow_pickle=False) as saved: data = {name: saved[name] for name in saved.files}
        if list(data["reaction_ids"]) != rids: raise ValueError("Reaction order")
        for left, right in PAIRS:
            for a in (row for row in blockrows if row["method"] == left):
                b = records[(right, a["task"], bn, a["sequence_sha256"])]
                if a["target_reaction_indices"] != b["target_reaction_indices"] or a["protein_group"] != b["protein_group"]: raise ValueError("Pair denominator")
                qi = a["query_index_in_block"]
                if str(data["query_ids"][qi]) != a["sequence_sha256"]: raise ValueError("Query order")
                left_values = data[left]; right_values = data[right]
                left_values = left_values[qi] if left_values.ndim == 2 else left_values
                right_values = right_values[qi] if right_values.ndim == 2 else right_values
                left_mask = data["domain_" + left]; right_mask = data["domain_" + right]
                left_mask = left_mask[qi] if left_mask.ndim == 2 else left_mask
                right_mask = right_mask[qi] if right_mask.ndim == 2 else right_mask
                common = left_mask & right_mask; positive = [i for i in a["target_reaction_indices"] if common[i]]
                common_left = rank_result(left_values, common, positive, tie)["rr"] if positive else None
                common_right = rank_result(right_values, common, positive, tie)["rr"] if positive else None
                rows.append({"cohort": channel, "task": a["task"], "block_number": bn, "sequence_sha256": a["sequence_sha256"],
                             "protein_group": a["protein_group"], "left_method": left, "right_method": right,
                             "left_rr": a["rr"], "right_rr": b["rr"], "full_candidate_count": len(rids),
                             "target_positives": a["positive_count"], "left_covered_positives": a["covered_positives"],
                             "right_covered_positives": b["covered_positives"], "left_candidate_count": int(left_mask.sum()),
                             "right_candidate_count": int(right_mask.sum()), "common_candidate_count": int(common.sum()),
                             "common_positive_count": len(positive), "common_left_rr": common_left,
                             "common_right_rr": common_right, "common_reranking_defined": bool(positive)})
        print(f"{channel} lineage paired block {bn}", flush=True)
    rng = np.random.default_rng(20260909); comparisons = []
    for task in sorted({row["task"] for row in rows}):
        for left, right in PAIRS:
            chosen = [row for row in rows if (row["task"], row["left_method"], row["right_method"]) == (task, left, right)]
            common = [row for row in chosen if row["common_reranking_defined"]]
            comparisons.append({"task": task, "left_method": left, "right_method": right,
                                "end_to_end": group_summary(chosen, "left_rr", "right_rr", rng),
                                "common_domain_reranking": group_summary(common, "common_left_rr", "common_right_rr", rng),
                                "coverage": {"original_target_positive_instances": sum(r["target_positives"] for r in chosen),
                                             "left_covered_positive_instances": sum(r["left_covered_positives"] for r in chosen),
                                             "right_covered_positive_instances": sum(r["right_covered_positives"] for r in chosen),
                                             "common_positive_instances": sum(r["common_positive_count"] for r in chosen),
                                             "common_undefined_query_panels": len(chosen) - len(common)}})
    out.mkdir()
    with (out / "paired_query_records.jsonl").open("w") as stream:
        for row in rows: stream.write(json.dumps(row, allow_nan=False) + "\n")
    summary = {"status": "COMPUTED_QC_PENDING", "created_utc": now(), "cohort": channel,
               "comparisons": comparisons, "paired_query_rows": len(rows), "bootstrap_seed": 20260909,
               "bootstrap_replicates": 5000, "comparison_list_frozen_before_lineage_results": True,
               "confirmatory_hypothesis_test": False, "independent_validation": False, "multiplicity_adjusted": False,
               "interval_scope": "Conditional fixed-fit protein-group bootstrap; no chemical/publication/refit/external uncertainty"}
    write_json(out / "comparison_summary.json", summary)
    names = ["lineage_paired_comparisons.py", "paired_expert_comparisons.py", "TAXONOMIC_LINEAGE_CONTROL_V1.md",
             result.relative_to(root).as_posix() + "/per_query_metrics.jsonl",
             qc.relative_to(root).as_posix(), "dataset_02/core_reactions.json"] + scorefiles
    write_json(out / "input_manifest.json", {name: digest_file(root / name) for name in names})
    print(json.dumps({"cohort": channel, "paired_rows": len(rows), "comparisons": len(comparisons)}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    run(Path(__file__).resolve().parent, parser.parse_args().channel)
