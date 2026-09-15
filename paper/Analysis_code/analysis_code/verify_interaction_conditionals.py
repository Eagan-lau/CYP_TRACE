"""Independent reconstruction of nested candidate-specific interaction fits."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.sparse import csr_matrix
from scipy.special import logsumexp

from run_raw import digest_file, stable, write_json, now
from verify_main_baselines import FIELDS


METHODS = [
    "availability_chemistry_prior", "global_residual", "ordered_residual", "missing_residual",
    "lineage_residual", "composition_residual", "joint_global_ordered", "interaction_residual",
    "joint_global_interaction", "joint_lineage_interaction",
]


def feature_matrices(root, channel, sequences):
    site = np.load(root / "ordered_site_features_01/ordered_features.npz", allow_pickle=False)
    onehot = site[channel + "_projected_onehot"].astype(float)
    absent = (~site[channel + "_projected_present_mask"]).astype(float)
    ordered = np.concatenate([onehot, absent[..., None]], axis=2).reshape(len(sequences), -1)
    global_features = np.load(root / "esm_global_01/global_features.npz", allow_pickle=False)["features"].astype(float)
    lineage = np.load(root / "taxonomy_lineages_01/lineage_features.npz", allow_pickle=False)["features"].astype(float)
    for matrix in (ordered, global_features, lineage):
        matrix /= np.maximum(np.linalg.norm(matrix, axis=1, keepdims=True), 1e-12)
    return site, ordered, global_features, lineage


def conditional_delta(features, labels, kernel, query_features, declared):
    training = np.flatnonzero(labels.sum(axis=1))
    balanced = labels[training] / labels[training].sum(axis=1, keepdims=True)
    target = csr_matrix(balanced) @ kernel
    prior = np.asarray(target).mean(axis=0)
    probability = np.maximum(prior / prior.sum(), 1e-12); probability /= probability.sum()
    residual = np.asarray(target) - prior
    rms = float(np.sqrt(np.mean(residual * residual)))
    x = features[training]; center = x.mean(axis=0); centered = x - center
    gram = centered @ centered.T; scale = float(np.trace(gram) / len(training))
    if declared["zero_residual"]:
        delta = np.zeros((len(query_features), kernel.shape[1]))
    else:
        coefficient = np.linalg.solve(gram / scale + declared["alpha"] * np.eye(len(training)), residual)
        delta = ((query_features - center) @ centered.T / scale) @ coefficient / rms
    return np.log(probability), delta


def interaction_delta(protein, labels, kernel, transform, query_features, declared):
    training = np.flatnonzero(labels.sum(axis=1))
    balanced = labels[training] / labels[training].sum(axis=1, keepdims=True)
    target = np.asarray(csr_matrix(balanced) @ kernel)
    prior = target.mean(axis=0)
    probability = np.maximum(prior / prior.sum(), 1e-12); probability /= probability.sum()
    residual = target - prior
    rms = float(np.sqrt(np.mean(residual * residual)))
    x = protein[training]; center = x.mean(axis=0); centered = x - center
    scale = float(np.trace(centered @ centered.T) / len(training))
    checks = {
        "residual_rms": abs(rms - declared["residual_rms"]),
        "x_scale": abs(scale - declared["x_scale"]),
    }
    if declared["zero_interaction"] and declared["zero_reason"] == "insufficient_training_or_feature_variation":
        return np.log(probability), np.zeros((len(query_features), kernel.shape[1])), checks, None
    xs = centered / np.sqrt(scale)
    ux, sx, vtx = np.linalg.svd(xs, full_matrices=False)
    tolerance = float(np.finfo(float).eps * max(xs.shape) * sx[0])
    keep = sx > tolerance; ux, sx, vtx = ux[:, keep], sx[keep], vtx[keep]
    checks.update({"x_rank_tolerance": abs(tolerance - declared["x_rank_tolerance"]),
                   "x_rank": abs(len(sx) - declared["x_rank"])})
    zu = transform["u"].astype(float); sz = transform["singular_values"].astype(float)
    target_u = np.asarray(csr_matrix(balanced) @ transform["chemical_kernel_times_u"].astype(float))
    projected = ux.T @ (target_u - target_u.mean(axis=0))
    spectrum = (sx * sx)[:, None] * (sz * sz)[None, :]
    energy = projected * projected; observations = residual.size
    def evaluate(value):
        shrink = np.ones_like(spectrum) if value == 0 else spectrum / (spectrum + value)
        rss = float(np.sum(residual * residual)) - 2 * float(np.sum(energy * shrink)) + float(np.sum(energy * shrink * shrink))
        denominator = max(1 - 1 / len(training) - float(np.sum(shrink)) / observations, 1e-12)
        return max(rss, 0) / observations / denominator ** 2
    zero = float(np.mean(residual * residual)) / max(1 - 1 / len(training), 1e-12) ** 2
    checks["gcv_zero"] = abs(zero - declared["gcv_zero"])
    alphas = np.asarray(declared["gcv_curve_alpha"], dtype=float)
    expected_curve = np.asarray([evaluate(value) if np.isfinite(value) else zero for value in alphas])
    checks["gcv_curve"] = float(np.max(np.abs(expected_curve - np.asarray(declared["gcv_curve_value"]))))
    optimum = minimize_scalar(lambda loga: evaluate(float(np.exp(loga))), bounds=(-14, 14), method="bounded", options={"xatol": 1e-6})
    best_value, best_alpha = min([(evaluate(0.0), 0.0), (float(optimum.fun), float(np.exp(optimum.x)))])
    checks["gcv_optimum_value"] = abs(best_value - declared["gcv_selected"])
    if declared["zero_interaction"]:
        checks["gcv_zero_decision"] = max(0.0, zero * (1.0 - 1e-10) - best_value)
        return np.log(probability), np.zeros((len(query_features), kernel.shape[1])), checks, None
    alpha = float(declared["alpha"])
    core = projected * (sx[:, None] * sz[None, :]) / (spectrum + alpha)
    xlatent = ((query_features - center) / np.sqrt(scale)) @ vtx.T
    delta = ((xlatent @ core) * sz[None, :]) @ zu.T / rms
    checks["gcv_selected"] = abs(evaluate(alpha) - declared["gcv_selected"])
    checks["gcv_optimum_alpha"] = abs(best_alpha - alpha) / max(1.0, abs(alpha))
    return np.log(probability), delta, checks, core


def weighted_objective(prior, deltas, targets, records, coefficients):
    panels = Counter((row["group"], row["query"]) for row in records)
    groups = Counter(group for group, query in panels)
    weights = np.asarray([1 / len(groups) / groups[row["group"]] / panels[(row["group"], row["query"])] for row in records])
    score = prior.copy()
    for coefficient, delta in zip(coefficients, deltas): score += coefficient * delta
    normalizer = logsumexp(score, axis=1)
    probability = np.exp(score - normalizer[:, None])
    truth = np.asarray([score[index, target].mean() for index, target in enumerate(targets)])
    gradients = []
    for delta in deltas:
        expected = np.sum(probability * delta, axis=1)
        observed = np.asarray([delta[index, target].mean() for index, target in enumerate(targets)])
        gradients.append(float(weights @ (expected - observed)))
    return float(weights @ (normalizer - truth)), np.asarray(gradients)


def uniform_metrics(candidates, positives):
    harmonic = math.fsum(1 / rank for rank in range(1, candidates + 1)) / candidates
    denominator = math.comb(candidates, positives)
    reciprocal_rank = math.fsum(
        math.comb(candidates - rank, positives - 1) / denominator / rank
        for rank in range(1, candidates - positives + 2)
    )
    result = {"rr": reciprocal_rank, "mean_positive_rr": harmonic}
    for cutoff in (1, 5, 10):
        expected_hits = positives * min(cutoff, candidates) / candidates
        result["hit_at_" + str(cutoff)] = 1 - math.comb(candidates - min(cutoff, candidates), positives) / math.comb(candidates, positives) if positives <= candidates - min(cutoff, candidates) else 1.0
        result["recall_at_" + str(cutoff)] = expected_hits / positives
    result["ndcg_at_10"] = sum((positives / candidates) / math.log2(rank + 1) for rank in range(1, min(10, candidates) + 1)) / sum(1 / math.log2(rank + 1) for rank in range(1, min(positives, 10) + 1))
    return result


def run(root: Path, channel: str):
    result = root / f"interaction_conditional_{channel}_01"
    output = root / f"interaction_validation_{channel}_01"
    if output.exists(): raise FileExistsError(output)
    read = lambda path: json.loads(Path(path).read_text())
    failures = []
    for name, expected in read(result / "input_manifest.json").items():
        if digest_file(root / name) != expected: failures.append("input_hash:" + name)
    edges = read(root / "dataset_02/core_edges.json"); reaction_ids = sorted(read(root / "dataset_02/core_reactions.json"))
    sequences = sorted({edge["sequence_sha256"] for edge in edges}); si = {q: i for i, q in enumerate(sequences)}
    ri = {r: i for i, r in enumerate(reaction_ids)}; blocks = read(root / "multiaxis_split_01/outer_inner_blocks.json")
    axes = read(root / "multiaxis_split_01/axis_assignments.json")
    site, ordered, global_features, lineage = feature_matrices(root, channel, sequences)
    available = {q for q, flag in zip(sequences, site[channel + "_projected_available"]) if flag}
    representation = np.load(root / "main_baselines_01/fresh_representations.npz", allow_pickle=False); kernel = representation["chemical_kernel"].astype(float)
    loaded = np.load(root / "interaction_transform_01/reaction_transform.npz", allow_pickle=False)
    transform = {name: loaded[name] for name in loaded.files}
    receipts = read(result / "block_receipts.json"); evaluation = read(result / "evaluation.json")
    metrics = [json.loads(line) for line in (result / "per_query_metrics.jsonl").read_text().splitlines()]
    byblock = defaultdict(list)
    for row in metrics: byblock[row["block_number"]].append(row)
    ties = [stable(value) for value in reaction_ids]

    def construct(block):
        labels = np.zeros((len(sequences), len(reaction_ids))); truth = defaultdict(set); seen = set()
        for index in block["train_edge_indices"]:
            edge = edges[index]; seen.add(ri[edge["reaction_key"]])
            if edge["sequence_sha256"] in available: labels[si[edge["sequence_sha256"]], ri[edge["reaction_key"]]] = 1
        for index in block["test_edge_indices"]:
            edge = edges[index]
            if edge["sequence_sha256"] in available: truth[edge["sequence_sha256"]].add(ri[edge["reaction_key"]])
        return labels, truth, seen

    maximum_difference = 0.0; interaction_fits = score_vectors = metric_checks = coefficient_checks = 0
    aggregates = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for block_number, block in enumerate(blocks):
        labels, truth, seen = construct(block); queries = sorted(truth); training = np.flatnonzero(labels.sum(axis=1)); receipt = receipts[block_number]
        if receipt["query_ids"] != queries or receipt["training_sequence_ids"] != [sequences[index] for index in training]: failures.append("outer_pool")
        if not queries or not len(training):
            if receipt["status"] != "NOT_EVALUABLE" or block_number in byblock: failures.append("outer_unavailable")
            continue
        manifest = read(result / f"block_{block_number:03d}_inner_manifest.json")
        with np.load(result / f"block_{block_number:03d}_inner_predictions.npz", allow_pickle=False) as saved:
            inner_saved = {name: saved[name] for name in saved.files}
        indices = defaultdict(list)
        for index, row in enumerate(manifest["rows"]): indices[row["inner_number"]].append(index)
        for inner_number, inner in enumerate(block["inner_blocks"]):
            iy, it, _ = construct(inner); iq = sorted(it); declared = manifest["fits"][inner_number]
            selected = indices[inner_number]
            if not iq or not np.any(iy):
                if declared["status"] != "NOT_EVALUABLE" or selected: failures.append("inner_unavailable")
                continue
            if [manifest["rows"][index]["query"] for index in selected] != iq: failures.append("inner_order")
            for index, query in zip(selected, iq):
                if manifest["targets"][index] != sorted(it[query]): failures.append("inner_target")
            prior, delta, checks, _ = interaction_delta(ordered, iy, kernel, transform, ordered[[si[q] for q in iq]], declared["fit"])
            interaction_fits += 1
            difference = float(np.max(np.abs(delta - inner_saved["delta_interaction"][selected]))); maximum_difference = max(maximum_difference, difference)
            if difference > 1e-8 or not np.allclose(prior, inner_saved["log_prior"][selected], rtol=0, atol=1e-12): failures.append(f"inner_prediction:{block_number}:{inner_number}")
            if any(value > 1e-8 for value in checks.values()): failures.append(f"inner_fit_receipt:{block_number}:{inner_number}")
        if manifest["rows"]:
            prior = inner_saved["log_prior"]; target = manifest["targets"]; records = manifest["rows"]
            base, _ = weighted_objective(prior, [inner_saved["delta_interaction"]], target, records, [0])
            scalar = manifest["interaction_selection"]; value, gradient = weighted_objective(prior, [inner_saved["delta_interaction"]], target, records, [scalar["lambda"]]); coefficient_checks += 1
            if scalar["lambda"] < 0 or value > base + 1e-8: failures.append("scalar_selection")
            if scalar["reason"] == "convex_inner_optimum" and abs(gradient[0]) > 1e-6: failures.append("scalar_stationarity")
            for key, delta_name in (("joint_global_interaction_selection", "delta_global"), ("joint_lineage_interaction_selection", "delta_lineage")):
                choice = manifest[key]; coefficients = np.asarray(choice["lambdas"])
                value, gradient = weighted_objective(prior, [inner_saved[delta_name], inner_saved["delta_interaction"]], target, records, coefficients); coefficient_checks += 1
                alternative = min(base,
                    weighted_objective(prior, [inner_saved[delta_name]], target, records, [0])[0],
                    weighted_objective(prior, [inner_saved["delta_interaction"]], target, records, [scalar["lambda"]])[0])
                if np.any(coefficients < 0) or value > alternative + 1e-8: failures.append("joint_selection")
                projected = np.where(coefficients <= 1e-8, np.minimum(gradient, 0), gradient)
                if choice["joint_optimizer_accepted"] and np.max(np.abs(projected)) > 1e-6: failures.append("joint_stationarity")

        qfeatures = ordered[[si[q] for q in queries]]
        prior, interaction, checks, reconstructed_core = interaction_delta(ordered, labels, kernel, transform, qfeatures, receipt["interaction_fit"])
        interaction_fits += 1
        if any(value > 1e-8 for value in checks.values()): failures.append("outer_fit_receipt")
        global_prior, global_delta = conditional_delta(global_features, labels, kernel, global_features[[si[q] for q in queries]], read(root / f"lineage_conditional_{channel}_01/block_receipts.json")[block_number]["fits"]["global"])
        lineage_prior, lineage_delta = conditional_delta(lineage, labels, kernel, lineage[[si[q] for q in queries]], read(root / f"lineage_conditional_{channel}_01/block_receipts.json")[block_number]["fits"]["lineage"])
        if not np.allclose(prior, global_prior, rtol=0, atol=1e-12) or not np.allclose(prior, lineage_prior, rtol=0, atol=1e-12): failures.append("outer_prior")
        with np.load(result / f"block_{block_number:03d}_scores.npz", allow_pickle=False) as saved:
            data = {name: saved[name] for name in saved.files}
        expected = {
            "interaction_residual": prior[None, :] + receipt["interaction_selection"]["lambda"] * interaction,
            "joint_global_interaction": prior[None, :] + receipt["joint_global_interaction_selection"]["lambdas"][0] * global_delta + receipt["joint_global_interaction_selection"]["lambdas"][1] * interaction,
            "joint_lineage_interaction": prior[None, :] + receipt["joint_lineage_interaction_selection"]["lambdas"][0] * lineage_delta + receipt["joint_lineage_interaction_selection"]["lambdas"][1] * interaction,
        }
        with np.load(result / f"block_{block_number:03d}_model.npz", allow_pickle=False) as model_saved:
            if reconstructed_core is not None:
                difference = float(np.max(np.abs(reconstructed_core - model_saved["interaction__core_coefficients"]))); maximum_difference = max(maximum_difference, difference)
                if difference > 1e-8: failures.append("outer_core_coefficients")
        for name, value in expected.items():
            difference = float(np.max(np.abs(value - data[name]))); maximum_difference = max(maximum_difference, difference)
            if difference > 1e-5 or not np.all(data["domain_" + name]): failures.append("outer_score:" + name)
            score_vectors += len(queries)
        if list(data["query_ids"]) != queries or list(data["reaction_ids"]) != reaction_ids: failures.append("score_order")
        rank_cache = {}
        for row in byblock[block_number]:
            query = row["sequence_sha256"]; qi = row["query_index_in_block"]; positives = truth[query]
            if row["task"] == "protein_cold_seen": positives = positives & seen
            elif row["task"] == "protein_cold_unseen": positives = positives - seen
            if row["target_reaction_indices"] != sorted(positives): failures.append("metric_target")
            if row["method"] == "uniform_expectation": values = uniform_metrics(len(reaction_ids), len(positives))
            else:
                key = (row["method"], qi)
                if key not in rank_cache:
                    score = data[row["method"]]; score = score[qi] if score.ndim == 2 else score
                    order = sorted(range(len(reaction_ids)), key=lambda index: (-float(score[index]), ties[index]))
                    rank_cache[key] = {index: rank + 1 for rank, index in enumerate(order)}
                ranks = [rank_cache[key][index] for index in positives]
                values = {"rr": 1 / min(ranks), "mean_positive_rr": math.fsum(1 / rank for rank in ranks) / len(ranks)}
                for cutoff in (1, 5, 10):
                    hits = sum(rank <= cutoff for rank in ranks); values["hit_at_" + str(cutoff)] = float(hits > 0); values["recall_at_" + str(cutoff)] = hits / len(ranks)
                values["ndcg_at_10"] = math.fsum(1 / math.log2(rank + 1) for rank in ranks if rank <= 10) / math.fsum(1 / math.log2(rank + 1) for rank in range(1, min(len(ranks), 10) + 1))
            for field, expected_value in values.items():
                if abs(expected_value - row[field]) > 1e-10: failures.append("metric:" + field)
            if row["positive_count"] != len(positives) or row["covered_positives"] != len(positives) or row["candidate_count"] != len(reaction_ids) or row["covered_candidates"] != len(reaction_ids): failures.append("metric_denominator")
            metric_checks += 1
            aggregates[(row["task"], row["method"])][row["protein_group"]][query].append(values)
        print(f"{channel} interaction independent QA {block_number + 1}/35", flush=True)
    for (task, method), groups in aggregates.items():
        for field in FIELDS:
            value = math.fsum(math.fsum(math.fsum(panel[field] for panel in panels) / len(panels) for panels in queries.values()) / len(queries) for queries in groups.values()) / len(groups)
            if abs(value - evaluation["aggregate"][task][method]["protein_group_macro_" + field]) > 1e-10: failures.append("aggregate:" + field)
    if len(metrics) != evaluation["metric_rows"]: failures.append("metric_row_count")
    report = {
        "status": "PASS" if not failures else "FAIL", "created_utc": now(), "cohort": channel,
        "failures": failures, "interaction_fits_reconstructed": interaction_fits,
        "coefficient_selections_checked": coefficient_checks, "score_vectors_reconstructed": score_vectors,
        "metric_rows_recomputed": metric_checks, "maximum_prediction_or_coefficient_difference": maximum_difference,
        "independent_implementation": True,
        "scope": "all evaluable nested fits, GCV receipts, coefficients, outer scores, ranks, denominators and macro aggregates",
    }
    output.mkdir(); write_json(output / "INDEPENDENT_QC.json", report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    args = parser.parse_args(); run(Path(__file__).resolve().parent, args.channel)
