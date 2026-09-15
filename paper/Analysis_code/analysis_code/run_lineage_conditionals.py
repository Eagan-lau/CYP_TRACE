"""Availability-matched taxonomy, local-site and global conditional models."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import time
import numpy as np
from conditional_core import fit_residual, predict_delta, select_lambda, model_receipt, normalized_features
from local_conditional_core import local_descriptions, select_joint
from run_main_baselines import rank_result, random_expectations, summarize
from run_raw import write_json, digest_file, stable, now


FEATURES = ("global", "ordered", "missing", "lineage")


def run(root, channel):
    out = root / f"lineage_conditional_{channel}_01"
    reference = root / f"local_conditional_{channel}_01"
    prefix = channel + "_projected"
    if out.exists(): raise FileExistsError(out)
    read = lambda name: json.loads((root / name).read_text())
    if read("taxonomy_lineages_validation_01/INDEPENDENT_QC.json")["status"] != "PASS": raise ValueError("Taxonomy QC gate")
    if read(f"local_validation_{channel}_01/INDEPENDENT_QC.json")["status"] != "PASS": raise ValueError("Matched local QC gate")
    for name, expected in read("multiaxis_split_01/output_checksums.json").items():
        if digest_file(root / "multiaxis_split_01" / name) != expected: raise ValueError("Frozen split changed")
    edges = read("dataset_02/core_edges.json"); rids = sorted(read("dataset_02/core_reactions.json"))
    seqs = sorted({edge["sequence_sha256"] for edge in edges}); si = {q: i for i, q in enumerate(seqs)}
    ri = {r: i for i, r in enumerate(rids)}; blocks = read("multiaxis_split_01/outer_inner_blocks.json")
    axes = read("multiaxis_split_01/axis_assignments.json")
    site = np.load(root / "ordered_site_features_01/ordered_features.npz", allow_pickle=False)
    esm = np.load(root / "esm_global_01/global_features.npz", allow_pickle=False)
    rep = np.load(root / "main_baselines_01/fresh_representations.npz", allow_pickle=False)
    lineage = np.load(root / "taxonomy_lineages_01/lineage_features.npz", allow_pickle=False)
    if any(list(data["sequence_ids"]) != seqs for data in [site, esm, rep, lineage]) or list(rep["reaction_ids"]) != rids:
        raise ValueError("Feature order")
    available = site[prefix + "_available"]; availset = {q for q, flag in zip(seqs, available) if flag}
    base = local_descriptions(site[prefix + "_onehot"], site[prefix + "_present_mask"], esm["features"])
    features = {name: base[name] for name in ["global", "ordered", "missing"]}
    features["lineage"] = normalized_features(lineage["features"])
    K = rep["chemical_kernel"].astype(float); tie = np.empty(len(rids), dtype=int)
    for j, i in enumerate(sorted(range(len(rids)), key=lambda i: stable(rids[i]))): tie[i] = j

    def construct(block):
        y = np.zeros((len(seqs), len(rids))); truth = defaultdict(set); full_seen = set()
        for index in block["train_edge_indices"]:
            edge = edges[index]; full_seen.add(ri[edge["reaction_key"]])
            if edge["sequence_sha256"] in availset: y[si[edge["sequence_sha256"]], ri[edge["reaction_key"]]] = 1
        for index in block["test_edge_indices"]:
            edge = edges[index]
            if edge["sequence_sha256"] in availset: truth[edge["sequence_sha256"]].add(ri[edge["reaction_key"]])
        return y, truth, full_seen

    out.mkdir(); allrows = []; receipts = []; maximum_reference_difference = 0.
    for bn, block in enumerate(blocks):
        started = time.time(); y, truth, seen = construct(block); qs = sorted(truth); training = np.flatnonzero(y.sum(1))
        receipt = {"block_number": bn, "block_id": block["block_id"], "query_ids": qs,
                   "training_sequence_ids": [seqs[i] for i in training], "retained_train_edges": int(y.sum()),
                   "retained_test_edges": sum(map(len, truth.values())), "full_training_seen_labels": len(seen),
                   "available_training_seen_labels": int(np.sum(y.sum(0) > 0))}
        if not len(training) or not qs:
            receipts.append({**receipt, "status": "NOT_EVALUABLE", "reason": "empty_available_training_or_query"}); continue
        logs = []; deltas = {name: [] for name in FEATURES}; records = []; targets = []; innerfits = []
        for j, inner in enumerate(block["inner_blocks"]):
            if not set(inner["train_edge_indices"]).issubset(block["train_edge_indices"]) or not set(inner["test_edge_indices"]).issubset(block["train_edge_indices"]):
                raise ValueError("Inner leakage")
            iy, it, _ = construct(inner); iq = sorted(it)
            if not np.any(iy) or not iq:
                innerfits.append({"inner_number": j, "status": "NOT_EVALUABLE", "reason": "empty_available_training_or_query"}); continue
            models = {name: fit_residual(features[name], iy, K) for name in FEATURES}; prior = models["global"]["prior_log"]
            for name, model in models.items():
                if not np.allclose(prior, model["prior_log"], rtol=0, atol=1e-12): raise ValueError("Feature-specific chemical backbone")
                deltas[name].extend(predict_delta(model, features[name][[si[q] for q in iq]]))
            logs.extend([prior] * len(iq)); targets.extend([sorted(it[q]) for q in iq])
            records.extend([{"group": axes["protein_groups"][q], "query": q, "inner_number": j} for q in iq])
            innerfits.append({"inner_number": j, "status": "COMPLETE", "training_sequence_ids": [seqs[i] for i in np.flatnonzero(iy.sum(1))],
                              "query_ids": iq, "fits": {name: model_receipt(model) for name, model in models.items()}})
        logmat = np.array(logs).reshape(-1, len(rids)); mats = {name: np.array(values).reshape(-1, len(rids)) for name, values in deltas.items()}
        choices = {name: select_lambda(logmat, mats[name], targets, records) if records else {"lambda": 0., "reason": "no_inner_predictions"} for name in FEATURES}
        joint = select_joint(logmat, mats["lineage"], mats["ordered"], targets, records) if records else {"lambdas": [0., 0.], "reason": "no_inner_predictions", "joint_optimizer_accepted": False}
        np.savez_compressed(out / f"block_{bn:03d}_inner_predictions.npz", log_prior=logmat, **{"delta_" + name: values for name, values in mats.items()})
        write_json(out / f"block_{bn:03d}_inner_manifest.json", {"rows": records, "targets": targets, "fits": innerfits,
                                                                          "selections": choices, "joint_selection": joint})
        models = {name: fit_residual(features[name], y, K) for name in FEATURES}; prior = models["global"]["prior_log"]
        qidx = [si[q] for q in qs]; ds = {name: predict_delta(models[name], features[name][qidx]) for name in FEATURES}
        scores = {"availability_chemistry_prior": prior}
        for name in FEATURES: scores[name + "_residual"] = prior[None, :] + choices[name]["lambda"] * ds[name]
        scores["joint_lineage_ordered"] = prior[None, :] + joint["lambdas"][0] * ds["lineage"] + joint["lambdas"][1] * ds["ordered"]
        masks = {name: np.ones_like(values, dtype=bool) for name, values in scores.items()}
        with np.load(reference / f"block_{bn:03d}_scores.npz", allow_pickle=False) as ref:
            if list(ref["query_ids"]) != qs or list(ref["reaction_ids"]) != rids: raise ValueError("Reference order")
            for name in ["availability_chemistry_prior", "global_residual", "ordered_residual", "missing_residual"]:
                difference = float(np.max(np.abs(scores[name] - ref[name]))); maximum_reference_difference = max(maximum_reference_difference, difference)
                if difference > 1e-10 or not np.array_equal(masks[name], ref["domain_" + name]): raise ValueError("Existing-score reconstruction: " + name)
        np.savez_compressed(out / f"block_{bn:03d}_scores.npz", query_ids=np.array(qs), reaction_ids=np.array(rids),
                            **scores, **{"domain_" + name: mask for name, mask in masks.items()})
        saved = {"training_sequence_ids": np.array([seqs[i] for i in training]), "prior_log": prior,
                 "joint_lambdas": np.array(joint["lambdas"])}
        for name, model in models.items():
            for key in ["center", "training_features_centered", "kernel_scale", "residual_rms", "coefficients"]:
                saved[name + "__" + key] = model[key]
            saved[name + "__lambda"] = choices[name]["lambda"]
        np.savez_compressed(out / f"block_{bn:03d}_models.npz", **saved)
        for qi, query in enumerate(qs):
            tasks = {block["task"] + "_all": truth[query]}
            if block["task"] == "protein_cold": tasks.update({"protein_cold_seen": truth[query] & seen, "protein_cold_unseen": truth[query] - seen})
            for task, positive in tasks.items():
                if not positive: continue
                for method in list(scores) + ["uniform_expectation"]:
                    metric = dict(random_expectations(len(rids), len(positive))) if method == "uniform_expectation" else rank_result(
                        scores[method][qi] if scores[method].ndim == 2 else scores[method],
                        masks[method][qi] if masks[method].ndim == 2 else masks[method], positive, tie)
                    allrows.append({"cohort": channel, "task": task, "block_number": bn, "block_id": block["block_id"],
                                    "sequence_sha256": query, "protein_group": axes["protein_groups"][query], "method": method,
                                    "query_index_in_block": qi, "target_reaction_indices": sorted(positive), **metric})
        receipts.append({**receipt, "status": "COMPLETE", "fits": {name: model_receipt(model) for name, model in models.items()},
                         "selections": choices, "joint_selection": joint,
                         "inner_completed_blocks": sum(r["status"] == "COMPLETE" for r in innerfits),
                         "inner_unavailable_blocks": sum(r["status"] != "COMPLETE" for r in innerfits), "inner_rows": len(records),
                         "maximum_existing_score_difference": maximum_reference_difference, "elapsed_seconds": time.time() - started})
        write_json(out / "block_receipts.json", receipts)
        write_json(out / "progress.json", {"status": "RUNNING", "completed_outer_entries": len(receipts), "updated_utc": now()})
        print(f"{channel} lineage block={bn + 1}/35 queries={len(qs)} joint={joint['lambdas']}", flush=True)
    with (out / "per_query_metrics.jsonl").open("w") as stream:
        for row in allrows: stream.write(json.dumps(row, allow_nan=False) + "\n")
    write_json(out / "block_receipts.json", receipts)
    evaluation = {"status": "COMPLETE_INDEPENDENT_QC_PENDING", "created_utc": now(), "cohort": channel,
                  "available_sequences": len(availset), "available_protein_groups": len({axes["protein_groups"][q] for q in availset}),
                  "core_sequences": len(seqs), "candidate_reactions": len(rids), "metric_rows": len(allrows),
                  "eligible_outer_blocks": sum(r["status"] == "COMPLETE" for r in receipts), "aggregate": summarize(allrows),
                  "feature_dimensions": {name: features[name].shape[1] for name in FEATURES},
                  "joint_optimizer_fallback_blocks": sum(not r["joint_selection"]["joint_optimizer_accepted"] for r in receipts if r["status"] == "COMPLETE"),
                  "maximum_existing_score_difference": maximum_reference_difference,
                  "confirmatory_hypothesis_test": False, "independent_biological_validation": False, "final_router_or_tool_complete": False}
    write_json(out / "evaluation.json", evaluation)
    names = ["run_lineage_conditionals.py", "TAXONOMIC_LINEAGE_CONTROL_V1.md", "conditional_core.py", "local_conditional_core.py",
             "run_main_baselines.py", "dataset_02/core_edges.json", "dataset_02/core_reactions.json",
             "ordered_site_features_01/ordered_features.npz", "ordered_site_validation_01/INDEPENDENT_QC.json",
             "esm_global_01/global_features.npz", "main_baselines_01/fresh_representations.npz",
             "taxonomy_lineages_01/lineage_features.npz", "taxonomy_lineages_validation_01/INDEPENDENT_QC.json",
             "multiaxis_split_01/outer_inner_blocks.json", "multiaxis_split_01/axis_assignments.json",
             f"local_validation_{channel}_01/INDEPENDENT_QC.json"]
    write_json(out / "input_manifest.json", {name: digest_file(root / name) for name in names})
    write_json(out / "progress.json", {"status": evaluation["status"], "outer_entries": len(receipts), "updated_utc": now()})
    print(json.dumps({k: v for k, v in evaluation.items() if k != "aggregate"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    run(Path(__file__).resolve().parent, parser.parse_args().channel)
