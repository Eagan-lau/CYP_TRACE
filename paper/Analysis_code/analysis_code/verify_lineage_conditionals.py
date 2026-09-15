"""Independent reconstruction of taxonomy/local/global conditional models."""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from scipy.special import logsumexp
from verify_main_baselines import uniform, FIELDS
from run_raw import digest_file, stable, write_json, now


FEATURES = ("global", "ordered", "missing", "lineage")


def run(root, channel):
    result = root / f"lineage_conditional_{channel}_01"; out = root / f"lineage_validation_{channel}_01"
    reference = root / f"local_conditional_{channel}_01"; prefix = channel + "_projected"
    if out.exists(): raise FileExistsError(out)
    read = lambda path: json.loads(path.read_text()); failures = []
    for name, expected in read(result / "input_manifest.json").items():
        if digest_file(root / name) != expected: failures.append("source:" + name)
    edges = read(root / "dataset_02/core_edges.json"); rids = sorted(read(root / "dataset_02/core_reactions.json"))
    seqs = sorted({edge["sequence_sha256"] for edge in edges}); si = {q: i for i, q in enumerate(seqs)}
    ri = {r: i for i, r in enumerate(rids)}; blocks = read(root / "multiaxis_split_01/outer_inner_blocks.json")
    receipts = read(result / "block_receipts.json"); evaluation = read(result / "evaluation.json")
    axes = read(root / "multiaxis_split_01/axis_assignments.json")
    site = np.load(root / "ordered_site_features_01/ordered_features.npz", allow_pickle=False)
    available = {q for q, flag in zip(seqs, site[prefix + "_available"]) if flag}
    onehot = site[prefix + "_onehot"].astype(float); absent = (~site[prefix + "_present_mask"]).astype(float)
    ordered = np.empty((len(seqs), 35, 22)); ordered[..., :21] = onehot; ordered[..., 21] = absent
    missing = np.empty((len(seqs), 35, 2)); missing[..., 0] = onehot[..., 20]; missing[..., 1] = absent
    features = {"global": np.load(root / "esm_global_01/global_features.npz", allow_pickle=False)["features"].astype(float),
                "ordered": ordered.reshape(len(seqs), -1), "missing": missing.reshape(len(seqs), -1),
                "lineage": np.load(root / "taxonomy_lineages_01/lineage_features.npz", allow_pickle=False)["features"].astype(float)}
    for name in FEATURES:
        features[name] /= np.maximum(np.sqrt(np.sum(features[name] * features[name], axis=1, keepdims=True)), 1e-12)
    rep = np.load(root / "main_baselines_01/fresh_representations.npz", allow_pickle=False); K = rep["chemical_kernel"].astype(float)
    ties = [stable(r) for r in rids]

    def construct(block):
        y = np.zeros((len(seqs), len(rids))); truth = defaultdict(set); seen = set()
        for index in block["train_edge_indices"]:
            edge = edges[index]; seen.add(ri[edge["reaction_key"]])
            if edge["sequence_sha256"] in available: y[si[edge["sequence_sha256"]], ri[edge["reaction_key"]]] = 1
        for index in block["test_edge_indices"]:
            edge = edges[index]
            if edge["sequence_sha256"] in available: truth[edge["sequence_sha256"]].add(ri[edge["reaction_key"]])
        return y, truth, seen

    def reconstruct(y, queries, name, declared):
        training = np.flatnonzero(y.sum(1)); target = (y[training] / y[training].sum(1, keepdims=True)) @ K
        prior = target.mean(0); probability = np.maximum(prior / prior.sum(), 1e-12); probability /= probability.sum()
        residual = target - prior; rms = float(np.sqrt(np.mean(residual ** 2))); x = features[name][training]
        center = x.mean(0); centered = x - center; gram = centered @ centered.T; scale = np.trace(gram) / len(training)
        if declared["zero_residual"]: delta = np.zeros((len(queries), len(rids)))
        else:
            coefficient = solve(gram / scale + declared["alpha"] * np.eye(len(training)), residual, assume_a="pos", check_finite=False)
            delta = ((features[name][[si[q] for q in queries]] - center) @ centered.T / scale) @ coefficient / rms
        return np.log(probability), delta

    def objective(prior, deltas, targets, records, coefficients):
        panels = Counter((r["group"], r["query"]) for r in records); queries = Counter(g for g, q in panels)
        weights = np.array([1 / len(queries) / queries[r["group"]] / panels[(r["group"], r["query"])] for r in records])
        score = prior.copy()
        for coefficient, delta in zip(coefficients, deltas): score += coefficient * delta
        normalizer = logsumexp(score, axis=1); probability = np.exp(score - normalizer[:, None])
        true_score = np.array([score[i, target].mean() for i, target in enumerate(targets)])
        gradients = [float(weights @ (np.sum(probability * delta, axis=1) - np.array([delta[i, target].mean() for i, target in enumerate(targets)]))) for delta in deltas]
        return float(weights @ (normalizer - true_score)), np.array(gradients)

    metrics = [json.loads(line) for line in (result / "per_query_metrics.jsonl").read_text().splitlines()]
    byblock = defaultdict(list); aggregate = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in metrics: byblock[row["block_number"]].append(row)
    inner_count = outer_count = scalar_checks = joint_checks = vectors = 0; maximum_difference = 0.; reference_checks = 0
    for bn, block in enumerate(blocks):
        y, truth, seen = construct(block); queries = sorted(truth); training = np.flatnonzero(y.sum(1)); receipt = receipts[bn]
        if receipt["query_ids"] != queries or receipt["training_sequence_ids"] != [seqs[i] for i in training]: failures.append("outer_pool")
        if not queries or not len(training):
            if receipt["status"] != "NOT_EVALUABLE" or bn in byblock: failures.append("unevaluable")
            continue
        info = read(result / f"block_{bn:03d}_inner_manifest.json")
        with np.load(result / f"block_{bn:03d}_inner_predictions.npz", allow_pickle=False) as saved: inner = {name: saved[name] for name in saved.files}
        indices = defaultdict(list)
        for i, row in enumerate(info["rows"]): indices[row["inner_number"]].append(i)
        for j, innerblock in enumerate(block["inner_blocks"]):
            iy, it, _ = construct(innerblock); iq = sorted(it); declared = info["fits"][j]
            if not iq or not np.any(iy):
                if declared["status"] != "NOT_EVALUABLE" or j in indices: failures.append("inner_unevaluable")
                continue
            selected = indices[j]
            if declared["query_ids"] != iq or declared["training_sequence_ids"] != [seqs[i] for i in np.flatnonzero(iy.sum(1))]: failures.append("inner_pool")
            if [info["rows"][i]["query"] for i in selected] != iq: failures.append("inner_order")
            for i, query in zip(selected, iq):
                if info["targets"][i] != sorted(it[query]) or info["rows"][i]["group"] != axes["protein_groups"][query]: failures.append("inner_target")
            for name in FEATURES:
                prior, delta = reconstruct(iy, iq, name, declared["fits"][name]); inner_count += 1
                difference = float(np.max(np.abs(delta - inner["delta_" + name][selected]))); maximum_difference = max(maximum_difference, difference)
                if difference > 1e-5 or not np.allclose(prior, inner["log_prior"][selected], atol=1e-10, rtol=0): failures.append(f"inner_fit:{bn}:{j}:{name}")
        if info["rows"]:
            prior = inner["log_prior"]; targets = info["targets"]; records = info["rows"]
            deltas = {name: inner["delta_" + name] for name in FEATURES}; base, _ = objective(prior, [deltas["global"]], targets, records, [0.])
            for name, choice in info["selections"].items():
                value, gradient = objective(prior, [deltas[name]], targets, records, [choice["lambda"]]); scalar_checks += 1
                if choice["lambda"] < 0 or value > base + 1e-8: failures.append("scalar_loss")
                if choice["reason"] == "convex_inner_optimum" and abs(gradient[0]) > 1e-6: failures.append("scalar_stationarity")
                if choice["reason"] == "zero_boundary_optimum" and gradient[0] < -1e-8: failures.append("scalar_boundary")
            choice = info["joint_selection"]; coefficients = np.array(choice["lambdas"])
            value, gradient = objective(prior, [deltas["lineage"], deltas["ordered"]], targets, records, coefficients); joint_checks += 1
            alternatives = [base] + [objective(prior, [deltas[name]], targets, records, [info["selections"][name]["lambda"]])[0] for name in ["lineage", "ordered"]]
            if np.any(coefficients < 0) or value > min(alternatives) + 1e-8: failures.append("joint_loss")
            projected = np.where(coefficients <= 1e-8, np.minimum(gradient, 0), gradient)
            if choice["joint_optimizer_accepted"] and np.max(np.abs(projected)) > 1e-6: failures.append("joint_stationarity")
        with np.load(result / f"block_{bn:03d}_scores.npz", allow_pickle=False) as saved: data = {name: saved[name] for name in saved.files}
        expected = {}; deltas = {}
        for name in FEATURES:
            prior, delta = reconstruct(y, queries, name, receipt["fits"][name]); deltas[name] = delta; outer_count += 1
            expected[name + "_residual"] = prior[None, :] + receipt["selections"][name]["lambda"] * delta
        expected["availability_chemistry_prior"] = prior
        coefficients = receipt["joint_selection"]["lambdas"]
        expected["joint_lineage_ordered"] = prior[None, :] + coefficients[0] * deltas["lineage"] + coefficients[1] * deltas["ordered"]
        if list(data["query_ids"]) != queries or list(data["reaction_ids"]) != rids: failures.append("candidate_order")
        with np.load(reference / f"block_{bn:03d}_scores.npz", allow_pickle=False) as ref:
            for name in ["availability_chemistry_prior", "global_residual", "ordered_residual", "missing_residual"]:
                difference = float(np.max(np.abs(data[name] - ref[name]))); maximum_difference = max(maximum_difference, difference); reference_checks += 1
                if difference > 1e-10: failures.append("reference_score:" + name)
        masks = {name: np.ones_like(value, dtype=bool) for name, value in expected.items()}
        for name, value in expected.items():
            difference = float(np.max(np.abs(value - data[name]))); maximum_difference = max(maximum_difference, difference)
            if not np.allclose(value, data[name], rtol=1e-7, atol=1e-5): failures.append("outer_score:" + name)
            if not np.array_equal(masks[name], data["domain_" + name]): failures.append("outer_domain:" + name)
            vectors += len(queries) if value.ndim == 2 else 1
        cache = {}
        for row in byblock[bn]:
            query = row["sequence_sha256"]; qi = row["query_index_in_block"]; name = row["method"]; positive = truth[query]
            if row["task"] == "protein_cold_seen": positive = positive & seen
            elif row["task"] == "protein_cold_unseen": positive = positive - seen
            if row["target_reaction_indices"] != sorted(positive) or row["protein_group"] != axes["protein_groups"][query]: failures.append("outer_target")
            if name == "uniform_expectation": value = uniform(len(rids), len(positive)); covered = len(positive); conditional = value["rr"]
            else:
                if (name, qi) not in cache:
                    score = data[name]; score = score[qi] if score.ndim == 2 else score
                    order = sorted(range(len(rids)), key=lambda i: (-float(score[i]), ties[i])); cache[(name, qi)] = {i: j + 1 for j, i in enumerate(order)}
                ranks = cache[(name, qi)]; positions = [ranks[i] for i in positive]; covered = len(positive); best = min(positions)
                value = {"rr": 1 / best, "mean_positive_rr": math.fsum(1 / i for i in positions) / len(positions)}; conditional = 1 / best
                for k in [1, 5, 10]:
                    hits = sum(i <= k for i in positions); value["hit_at_" + str(k)] = float(hits > 0); value["recall_at_" + str(k)] = hits / len(positions)
                value["ndcg_at_10"] = math.fsum(1 / math.log2(i + 1) for i in positions if i <= 10) / math.fsum(1 / math.log2(i + 1) for i in range(1, min(len(positions), 10) + 1))
            for field, expected_value in value.items():
                if abs(expected_value - row[field]) > 1e-10: failures.append("rank:" + field)
            if row["positive_count"] != len(positive) or row["covered_positives"] != covered or row["covered_candidates"] != len(rids) or row["candidate_count"] != len(rids): failures.append("coverage")
            if row["conditional_rr"] is None or abs(row["conditional_rr"] - conditional) > 1e-10: failures.append("conditional_rr")
            aggregate[(row["task"], name)][row["protein_group"]][query].append(value)
        print(f"{channel} lineage independent QA {bn + 1}/35", flush=True)
    for (task, name), groups in aggregate.items():
        for field in FIELDS:
            value = math.fsum(math.fsum(math.fsum(row[field] for row in panels) / len(panels) for panels in queries.values()) / len(queries) for queries in groups.values()) / len(groups)
            if abs(value - evaluation["aggregate"][task][name]["protein_group_macro_" + field]) > 1e-10: failures.append("aggregate:" + field)
    if len(metrics) != evaluation["metric_rows"] or len(byblock) != evaluation["eligible_outer_blocks"]: failures.append("overall_counts")
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "cohort": channel, "failures": failures,
              "inner_feature_models_reconstructed": inner_count, "outer_feature_models_reconstructed": outer_count,
              "scalar_coefficient_checks": scalar_checks, "joint_coefficient_checks": joint_checks,
              "existing_score_identity_checks": reference_checks, "score_vectors_reconstructed": vectors,
              "metric_rows_recomputed": len(metrics), "aggregates_recomputed": len(aggregate),
              "maximum_prediction_difference": maximum_difference, "verifier_sha256": digest_file(Path(__file__)),
              "independent_biological_validation": False,
              "scope": "All evaluable lineage/local/global pools, fits, coefficients, scores, full-catalogue ranks and macro metrics"}
    out.mkdir(); write_json(out / "INDEPENDENT_QC.json", report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--channel", choices=["sequence", "structure"], required=True)
    run(Path(__file__).resolve().parent, parser.parse_args().channel)
