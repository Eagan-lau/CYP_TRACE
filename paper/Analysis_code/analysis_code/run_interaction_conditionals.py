"""Nested candidate-specific ordered-site x reaction-centre models."""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import time

import numpy as np

from conditional_core import fit_residual, predict_delta, normalized_features
from interaction_core import fit_interaction, predict_interaction, interaction_receipt
from local_conditional_core import local_descriptions, select_joint
from conditional_core import select_lambda
from run_main_baselines import rank_result, random_expectations, summarize
from run_raw import write_json, digest_file, stable, now


COPIED_LINEAGE_METHODS = [
    "availability_chemistry_prior", "global_residual", "ordered_residual",
    "missing_residual", "lineage_residual",
]
COPIED_LOCAL_METHODS = ["composition_residual", "joint_global_ordered"]


def run(root: Path, channel: str) -> None:
    output = root / f"interaction_conditional_{channel}_01"
    lineage_result = root / f"lineage_conditional_{channel}_01"
    local_result = root / f"local_conditional_{channel}_01"
    prefix = channel + "_projected"
    if output.exists():
        raise FileExistsError(output)
    read = lambda path: json.loads(Path(path).read_text())
    gates = [
        root / "reaction_center_environment_validation_01/audit.json",
        root / "interaction_transform_validation_01/audit.json",
        root / f"lineage_validation_{channel}_01/INDEPENDENT_QC.json",
        root / f"local_validation_{channel}_01/INDEPENDENT_QC.json",
    ]
    if any(read(path)["status"] != "PASS" for path in gates):
        raise ValueError("interaction input QC gate")
    for name, expected in read(root / "multiaxis_split_01/output_checksums.json").items():
        if digest_file(root / "multiaxis_split_01" / name) != expected:
            raise ValueError("frozen split changed")

    edges = read(root / "dataset_02/core_edges.json")
    reaction_ids = sorted(read(root / "dataset_02/core_reactions.json"))
    sequences = sorted({edge["sequence_sha256"] for edge in edges})
    sequence_index = {value: index for index, value in enumerate(sequences)}
    reaction_index = {value: index for index, value in enumerate(reaction_ids)}
    blocks = read(root / "multiaxis_split_01/outer_inner_blocks.json")
    axes = read(root / "multiaxis_split_01/axis_assignments.json")
    site = np.load(root / "ordered_site_features_01/ordered_features.npz", allow_pickle=False)
    esm = np.load(root / "esm_global_01/global_features.npz", allow_pickle=False)
    lineage = np.load(root / "taxonomy_lineages_01/lineage_features.npz", allow_pickle=False)
    representation = np.load(root / "main_baselines_01/fresh_representations.npz", allow_pickle=False)
    reaction_archive = np.load(root / "interaction_transform_01/reaction_transform.npz", allow_pickle=False)
    reaction_transform = {name: reaction_archive[name] for name in reaction_archive.files}
    lineage_receipts = read(lineage_result / "block_receipts.json")
    if any(list(data["sequence_ids"]) != sequences for data in (site, esm, lineage)):
        raise ValueError("protein feature order")
    if list(representation["reaction_ids"]) != reaction_ids or list(reaction_transform["reaction_keys"]) != reaction_ids:
        raise ValueError("reaction feature order")
    availability = site[prefix + "_available"]
    available = {sequence for sequence, flag in zip(sequences, availability) if flag}
    descriptions = local_descriptions(site[prefix + "_onehot"], site[prefix + "_present_mask"], esm["features"])
    descriptions["lineage"] = normalized_features(lineage["features"])
    protein_interaction = descriptions["ordered"]
    chemical_kernel = representation["chemical_kernel"].astype(np.float64)
    tie = np.empty(len(reaction_ids), dtype=int)
    for order, index in enumerate(sorted(range(len(reaction_ids)), key=lambda i: stable(reaction_ids[i]))):
        tie[index] = order

    def construct(block):
        labels = np.zeros((len(sequences), len(reaction_ids)), dtype=np.float64)
        truth = defaultdict(set)
        full_seen = set()
        for edge_index in block["train_edge_indices"]:
            edge = edges[edge_index]
            full_seen.add(reaction_index[edge["reaction_key"]])
            if edge["sequence_sha256"] in available:
                labels[sequence_index[edge["sequence_sha256"]], reaction_index[edge["reaction_key"]]] = 1.0
        for edge_index in block["test_edge_indices"]:
            edge = edges[edge_index]
            if edge["sequence_sha256"] in available:
                truth[edge["sequence_sha256"]].add(reaction_index[edge["reaction_key"]])
        return labels, truth, full_seen

    output.mkdir()
    all_metrics = []
    receipts = []
    maximum_reference_difference = 0.0
    for block_number, block in enumerate(blocks):
        started = time.time()
        labels, truth, seen = construct(block)
        queries = sorted(truth)
        training = np.flatnonzero(labels.sum(axis=1))
        receipt = {
            "block_number": block_number,
            "block_id": block["block_id"],
            "query_ids": queries,
            "training_sequence_ids": [sequences[index] for index in training],
            "retained_train_edges": int(labels.sum()),
            "retained_test_edges": sum(map(len, truth.values())),
            "full_training_seen_labels": len(seen),
            "available_training_seen_labels": int(np.sum(labels.sum(axis=0) > 0)),
        }
        if not len(training) or not queries:
            receipts.append({**receipt, "status": "NOT_EVALUABLE", "reason": "empty_available_training_or_query"})
            continue

        reference_manifest = read(lineage_result / f"block_{block_number:03d}_inner_manifest.json")
        with np.load(lineage_result / f"block_{block_number:03d}_inner_predictions.npz", allow_pickle=False) as saved:
            reference_inner = {name: saved[name] for name in saved.files}
        records = reference_manifest["rows"]
        targets = reference_manifest["targets"]
        interaction_inner = np.zeros_like(reference_inner["log_prior"])
        indices_by_inner = defaultdict(list)
        for index, record in enumerate(records):
            indices_by_inner[record["inner_number"]].append(index)
        inner_receipts = []
        for inner_number, inner in enumerate(block["inner_blocks"]):
            inner_labels, inner_truth, _ = construct(inner)
            inner_queries = sorted(inner_truth)
            selected_rows = indices_by_inner[inner_number]
            if not np.any(inner_labels) or not inner_queries:
                if selected_rows:
                    raise ValueError("reference contains rows for unavailable inner split")
                inner_receipts.append({"inner_number": inner_number, "status": "NOT_EVALUABLE"})
                continue
            if [records[index]["query"] for index in selected_rows] != inner_queries:
                raise ValueError("reference inner row order")
            interaction_model = fit_interaction(protein_interaction, inner_labels, chemical_kernel, reaction_transform)
            predicted = predict_interaction(interaction_model, protein_interaction[[sequence_index[q] for q in inner_queries]], reaction_transform)
            interaction_inner[selected_rows] = predicted
            if not np.allclose(interaction_model["prior_log"], reference_inner["log_prior"][selected_rows], rtol=0, atol=1e-12):
                raise ValueError("interaction and reference chemical prior differ")
            inner_receipts.append({
                "inner_number": inner_number,
                "status": "COMPLETE",
                "training_sequence_ids": [sequences[index] for index in np.flatnonzero(inner_labels.sum(axis=1))],
                "query_ids": inner_queries,
                "fit": interaction_receipt(interaction_model),
            })

        log_prior = reference_inner["log_prior"]
        if records:
            interaction_selection = select_lambda(log_prior, interaction_inner, targets, records)
            joint_global = select_joint(log_prior, reference_inner["delta_global"], interaction_inner, targets, records)
            joint_lineage = select_joint(log_prior, reference_inner["delta_lineage"], interaction_inner, targets, records)
        else:
            interaction_selection = {"lambda": 0.0, "reason": "no_inner_predictions"}
            joint_global = {"lambdas": [0.0, 0.0], "reason": "no_inner_predictions", "joint_optimizer_accepted": False}
            joint_lineage = {"lambdas": [0.0, 0.0], "reason": "no_inner_predictions", "joint_optimizer_accepted": False}
        np.savez_compressed(
            output / f"block_{block_number:03d}_inner_predictions.npz",
            log_prior=log_prior,
            delta_global=reference_inner["delta_global"],
            delta_lineage=reference_inner["delta_lineage"],
            delta_ordered=reference_inner["delta_ordered"],
            delta_interaction=interaction_inner,
        )
        write_json(output / f"block_{block_number:03d}_inner_manifest.json", {
            "rows": records,
            "targets": targets,
            "fits": inner_receipts,
            "interaction_selection": interaction_selection,
            "joint_global_interaction_selection": joint_global,
            "joint_lineage_interaction_selection": joint_lineage,
            "reference_lineage_manifest_sha256": digest_file(lineage_result / f"block_{block_number:03d}_inner_manifest.json"),
        })

        interaction_model = fit_interaction(protein_interaction, labels, chemical_kernel, reaction_transform)
        query_indices = [sequence_index[q] for q in queries]
        delta_interaction = predict_interaction(interaction_model, protein_interaction[query_indices], reaction_transform)
        global_model = fit_residual(descriptions["global"], labels, chemical_kernel)
        lineage_model = fit_residual(descriptions["lineage"], labels, chemical_kernel)
        delta_global = predict_delta(global_model, descriptions["global"][query_indices])
        delta_lineage = predict_delta(lineage_model, descriptions["lineage"][query_indices])
        prior = interaction_model["prior_log"]
        if not np.allclose(prior, global_model["prior_log"], rtol=0, atol=1e-12):
            raise ValueError("outer chemical prior differs")

        scores = {}
        masks = {}
        with np.load(lineage_result / f"block_{block_number:03d}_scores.npz", allow_pickle=False) as reference:
            if list(reference["query_ids"]) != queries or list(reference["reaction_ids"]) != reaction_ids:
                raise ValueError("lineage reference order")
            for name in COPIED_LINEAGE_METHODS:
                scores[name] = reference[name]
                masks[name] = reference["domain_" + name]
        with np.load(local_result / f"block_{block_number:03d}_scores.npz", allow_pickle=False) as reference:
            if list(reference["query_ids"]) != queries or list(reference["reaction_ids"]) != reaction_ids:
                raise ValueError("local reference order")
            for name in COPIED_LOCAL_METHODS:
                scores[name] = reference[name]
                masks[name] = reference["domain_" + name]
        reference_receipt = lineage_receipts[block_number]
        reconstructed_global = prior[None, :] + reference_receipt["selections"]["global"]["lambda"] * delta_global
        reconstructed_lineage = prior[None, :] + reference_receipt["selections"]["lineage"]["lambda"] * delta_lineage
        for name, expected in (("availability_chemistry_prior", prior), ("global_residual", reconstructed_global), ("lineage_residual", reconstructed_lineage)):
            difference = float(np.max(np.abs(scores[name] - expected)))
            maximum_reference_difference = max(maximum_reference_difference, difference)
            if difference > 1e-5:
                raise ValueError(f"reference score reconstruction failed: {name}; max_abs_difference={difference}")

        scores["interaction_residual"] = prior[None, :] + interaction_selection["lambda"] * delta_interaction
        scores["joint_global_interaction"] = prior[None, :] + joint_global["lambdas"][0] * delta_global + joint_global["lambdas"][1] * delta_interaction
        scores["joint_lineage_interaction"] = prior[None, :] + joint_lineage["lambdas"][0] * delta_lineage + joint_lineage["lambdas"][1] * delta_interaction
        for name in ("interaction_residual", "joint_global_interaction", "joint_lineage_interaction"):
            masks[name] = np.ones_like(scores[name], dtype=bool)
        np.savez_compressed(
            output / f"block_{block_number:03d}_scores.npz",
            query_ids=np.asarray(queries), reaction_ids=np.asarray(reaction_ids),
            **scores, **{"domain_" + name: mask for name, mask in masks.items()},
        )
        model_arrays = {
            "training_sequence_ids": np.asarray([sequences[index] for index in training]),
            "prior_log": prior,
            "interaction_lambda": np.asarray(interaction_selection["lambda"]),
            "joint_global_interaction_lambdas": np.asarray(joint_global["lambdas"]),
            "joint_lineage_interaction_lambdas": np.asarray(joint_lineage["lambdas"]),
        }
        for key in ("x_center", "x_scale", "x_singular_values", "x_vt", "core_coefficients", "residual_rms"):
            model_arrays["interaction__" + key] = np.asarray(interaction_model[key])
        model_arrays["interaction__alpha"] = np.asarray(np.nan if interaction_model["alpha"] is None else interaction_model["alpha"])
        model_arrays["interaction__zero"] = np.asarray(interaction_model["zero_interaction"])
        np.savez_compressed(output / f"block_{block_number:03d}_model.npz", **model_arrays)

        for query_number, query in enumerate(queries):
            tasks = {block["task"] + "_all": truth[query]}
            if block["task"] == "protein_cold":
                tasks.update({"protein_cold_seen": truth[query] & seen, "protein_cold_unseen": truth[query] - seen})
            for task, positives in tasks.items():
                if not positives:
                    continue
                for method in list(scores) + ["uniform_expectation"]:
                    if method == "uniform_expectation":
                        metric = dict(random_expectations(len(reaction_ids), len(positives)))
                    else:
                        score = scores[method]
                        domain = masks[method]
                        metric = rank_result(
                            score[query_number] if score.ndim == 2 else score,
                            domain[query_number] if domain.ndim == 2 else domain,
                            positives, tie,
                        )
                    all_metrics.append({
                        "cohort": channel, "task": task, "block_number": block_number,
                        "block_id": block["block_id"], "sequence_sha256": query,
                        "protein_group": axes["protein_groups"][query], "method": method,
                        "query_index_in_block": query_number,
                        "target_reaction_indices": sorted(positives), **metric,
                    })
        receipts.append({
            **receipt, "status": "COMPLETE", "interaction_fit": interaction_receipt(interaction_model),
            "interaction_selection": interaction_selection,
            "joint_global_interaction_selection": joint_global,
            "joint_lineage_interaction_selection": joint_lineage,
            "inner_completed_blocks": sum(item["status"] == "COMPLETE" for item in inner_receipts),
            "inner_unavailable_blocks": sum(item["status"] != "COMPLETE" for item in inner_receipts),
            "inner_rows": len(records), "elapsed_seconds": time.time() - started,
        })
        write_json(output / "block_receipts.json", receipts)
        write_json(output / "progress.json", {"status": "RUNNING", "completed_outer_entries": len(receipts), "updated_utc": now()})
        print(f"{channel} interaction block={block_number + 1}/35 queries={len(queries)} alpha={interaction_model['alpha']} lambda={interaction_selection['lambda']:.6g} seconds={time.time() - started:.1f}", flush=True)

    with (output / "per_query_metrics.jsonl").open("w") as handle:
        for row in all_metrics:
            handle.write(json.dumps(row, allow_nan=False) + "\n")
    write_json(output / "block_receipts.json", receipts)
    evaluation = {
        "status": "COMPLETE_INDEPENDENT_QC_PENDING", "created_utc": now(), "cohort": channel,
        "available_sequences": len(available),
        "available_protein_groups": len({axes["protein_groups"][q] for q in available}),
        "core_sequences": len(sequences), "candidate_reactions": len(reaction_ids),
        "reaction_centres_available": int(np.sum(reaction_transform["center_available"])),
        "mapping_warning_candidates": int(np.sum(reaction_transform["mapping_warning"])),
        "reaction_numerical_rank": int(len(reaction_transform["singular_values"])),
        "metric_rows": len(all_metrics), "eligible_outer_blocks": sum(row["status"] == "COMPLETE" for row in receipts),
        "aggregate": summarize(all_metrics),
        "interaction_zero_outer_fits": sum(row["interaction_fit"]["zero_interaction"] for row in receipts if row["status"] == "COMPLETE"),
        "interaction_zero_selected_coefficients": sum(row["interaction_selection"]["lambda"] == 0 for row in receipts if row["status"] == "COMPLETE"),
        "joint_global_zero_interaction_coefficients": sum(row["joint_global_interaction_selection"]["lambdas"][1] == 0 for row in receipts if row["status"] == "COMPLETE"),
        "joint_lineage_zero_interaction_coefficients": sum(row["joint_lineage_interaction_selection"]["lambdas"][1] == 0 for row in receipts if row["status"] == "COMPLETE"),
        "maximum_existing_score_difference": maximum_reference_difference,
        "historical_exposure": "DEVELOPMENT_EXPOSED_OR_UNVERIFIED",
        "independent_biological_validation": False, "final_router_or_tool_complete": False,
    }
    write_json(output / "evaluation.json", evaluation)
    manifest_names = [
        "run_interaction_conditionals.py", "interaction_core.py", "CANDIDATE_INTERACTION_MODEL_V1.md",
        "conditional_core.py", "local_conditional_core.py", "run_main_baselines.py", "run_raw.py",
        "dataset_02/core_edges.json", "dataset_02/core_reactions.json",
        "ordered_site_features_01/ordered_features.npz", "esm_global_01/global_features.npz",
        "taxonomy_lineages_01/lineage_features.npz", "main_baselines_01/fresh_representations.npz",
        "interaction_transform_01/reaction_transform.npz", "interaction_transform_validation_01/audit.json",
        "multiaxis_split_01/outer_inner_blocks.json", "multiaxis_split_01/axis_assignments.json",
        f"lineage_validation_{channel}_01/INDEPENDENT_QC.json", f"local_validation_{channel}_01/INDEPENDENT_QC.json",
        f"lineage_conditional_{channel}_01/block_receipts.json",
    ]
    for block_number in range(len(blocks)):
        candidates = [
            f"lineage_conditional_{channel}_01/block_{block_number:03d}_inner_manifest.json",
            f"lineage_conditional_{channel}_01/block_{block_number:03d}_inner_predictions.npz",
            f"lineage_conditional_{channel}_01/block_{block_number:03d}_scores.npz",
            f"local_conditional_{channel}_01/block_{block_number:03d}_scores.npz",
        ]
        manifest_names.extend(name for name in candidates if (root / name).exists())
    write_json(output / "input_manifest.json", {name: digest_file(root / name) for name in manifest_names})
    write_json(output / "progress.json", {"status": evaluation["status"], "outer_entries": len(receipts), "updated_utc": now()})
    print(json.dumps({key: value for key, value in evaluation.items() if key != "aggregate"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    arguments = parser.parse_args()
    run(Path(__file__).resolve().parent, arguments.channel)
