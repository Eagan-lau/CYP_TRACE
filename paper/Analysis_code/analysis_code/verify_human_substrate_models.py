"""Independent reconstruction of human-substrate features, fits and inference."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "human_substrate_models_01"
OUTPUT = ROOT / "human_substrate_model_validation_01"
KS = (1, 3, 5, 11, 25, 51, 101)
WS = tuple(round(x / 10, 1) for x in range(11))
SEED = 20260916
BOOTSTRAPS = 5000
PERMUTATIONS = 99
METHODS = ("pooled_chemical_knn", "isoform_specific_knn", "data_driven_shrinkage")


def load(path): return json.loads(path.read_text(encoding="utf-8"))
def close(a, b, tolerance=2e-12): return abs(float(a) - float(b)) <= tolerance


def ap(y, score, weight=None):
    y = np.asarray(y, dtype=np.int8); score = np.asarray(score, dtype=np.float64)
    weight = np.ones(len(y)) if weight is None else np.asarray(weight, dtype=np.float64)
    retained = weight > 0; y, score, weight = y[retained], score[retained], weight[retained]
    total = float(np.dot(weight, y))
    if total == 0: return math.nan
    order = np.argsort(-score, kind="mergesort"); s = score[order]; yy = y[order]; ww = weight[order]
    starts = np.r_[0, np.flatnonzero(np.diff(s) != 0) + 1]
    mass = np.add.reduceat(ww, starts); positive = np.add.reduceat(ww * yy, starts)
    return float(np.dot(np.cumsum(positive) / np.cumsum(mass), positive) / total)


def auc(y, score, weight=None):
    y = np.asarray(y, dtype=np.int8); score = np.asarray(score, dtype=np.float64)
    weight = np.ones(len(y)) if weight is None else np.asarray(weight, dtype=np.float64)
    retained = weight > 0; y, score, weight = y[retained], score[retained], weight[retained]
    pos = float(np.dot(weight, y)); neg = float(np.dot(weight, 1 - y))
    if pos == 0 or neg == 0: return math.nan
    order = np.argsort(score, kind="mergesort"); s = score[order]; yy = y[order]; ww = weight[order]
    starts = np.r_[0, np.flatnonzero(np.diff(s) != 0) + 1]
    gp = np.add.reduceat(ww * yy, starts); gn = np.add.reduceat(ww * (1 - yy), starts)
    below = np.r_[0.0, np.cumsum(gn)[:-1]]
    return float(np.dot(gp, below + 0.5 * gn) / (pos * neg))


def macro(y, score, protein, proteins, weight=None):
    aps, aucs, briers = [], [], []
    for p in range(proteins):
        mask = protein == p; w = None if weight is None else weight[mask]
        aps.append(ap(y[mask], score[mask], w)); aucs.append(auc(y[mask], score[mask], w))
        if w is None: briers.append(float(np.mean((score[mask] - y[mask]) ** 2)))
        else: briers.append(float(np.average((score[mask] - y[mask]) ** 2, weights=w)) if w.sum() else math.nan)
    return float(np.nanmean(aps)), float(np.nanmean(aucs)), float(np.nanmean(briers))


def make_features(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = []; bits = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for i, value in enumerate(smiles):
        molecule = Chem.MolFromSmiles(value)
        if molecule is None: raise ValueError("invalid normalized molecule")
        fp = generator.GetFingerprint(molecule); fps.append(fp); DataStructs.ConvertToNumpyArray(fp, bits[i])
    similarity = np.asarray([DataStructs.BulkTanimotoSimilarity(fp, fps) for fp in fps], dtype=np.float32)
    return bits, similarity


def predict(query, train, labels, similarity, tie, ks):
    train = np.asarray(train, dtype=int)
    pool_fallback = float(np.nanmean(np.nanmean(labels[train], axis=1)))
    iso_fallback = np.nanmean(labels[train], axis=0)
    pooled = {k: np.empty(len(query)) for k in ks}; specific = {k: np.empty((len(query), labels.shape[1])) for k in ks}
    for qn, q in enumerate(query):
        ordered = train[np.lexsort((tie[train], -similarity[q, train]))]
        for k in ks:
            chosen = ordered[:min(k, len(ordered))]; w = similarity[q, chosen].astype(np.float64)
            pooled[k][qn] = np.dot(w, np.nanmean(labels[chosen], axis=1)) / w.sum() if w.sum() else pool_fallback
        for p in range(labels.shape[1]):
            valid = ordered[np.isfinite(labels[ordered, p])]
            for k in ks:
                chosen = valid[:min(k, len(valid))]; w = similarity[q, chosen].astype(np.float64)
                specific[k][qn, p] = np.dot(w, labels[chosen, p]) / w.sum() if w.sum() else iso_fallback[p]
    return pooled, specific


def inner_arrays(outer, folds, labels, similarity, tie):
    ys, ps = [], []; pools = defaultdict(list); specific = defaultdict(list)
    for validation in range(5):
        if validation == outer: continue
        train = np.flatnonzero((folds != outer) & (folds != validation)); query = np.flatnonzero(folds == validation)
        a, b = predict(query, train, labels, similarity, tie, KS)
        for qn, q in enumerate(query):
            for p in range(labels.shape[1]):
                if not np.isfinite(labels[q, p]): continue
                ys.append(labels[q, p]); ps.append(p)
                for k in KS: pools[k].append(a[k][qn]); specific[k].append(b[k][qn, p])
    return np.asarray(ys, dtype=np.int8), np.asarray(ps, dtype=int), {k: np.asarray(v) for k, v in pools.items()}, {k: np.asarray(v) for k, v in specific.items()}


def choose(y, protein, pooled, specific, proteins):
    pm = {k: macro(y, pooled[k], protein, proteins) for k in KS}; im = {k: macro(y, specific[k], protein, proteins) for k in KS}
    pk = min(KS, key=lambda k: (-pm[k][0], -pm[k][1], k)); ik = min(KS, key=lambda k: (-im[k][0], -im[k][1], k))
    best = None
    for kp in KS:
        for ki in KS:
            for w in WS:
                metrics = macro(y, (1 - w) * pooled[kp] + w * specific[ki], protein, proteins)
                key = (-metrics[0], -metrics[1], kp + ki, w, kp, ki)
                if best is None or key < best[0]: best = (key, kp, ki, w, metrics)
    return {"pooled_k": pk, "pooled": pm[pk], "isoform_k": ik, "isoform": im[ik],
            "shrinkage_k_pool": best[1], "shrinkage_k_iso": best[2], "shrinkage_w": best[3], "shrinkage": best[4], "rows": len(y)}


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    failures = []
    manifest = load(SOURCE / "input_manifest.json")
    for name, expected in manifest.items():
        if digest_file(Path(name)) != expected: failures.append("input_hash:" + name)
    rows = load(ROOT / "human_substrate_01/normalized_labels.json")
    isoforms = sorted({r["isoform"] for r in rows}); compounds = sorted({r["compound_inchikey"] for r in rows})
    ci = {v: i for i, v in enumerate(compounds)}; pi = {v: i for i, v in enumerate(isoforms)}
    labels = np.full((len(compounds), len(isoforms)), np.nan); smiles = [None] * len(compounds); scaffolds = [None] * len(compounds); folds = np.full(len(compounds), -1, dtype=int)
    for r in rows:
        c, p = ci[r["compound_inchikey"]], pi[r["isoform"]]
        if np.isfinite(labels[c, p]): failures.append("duplicate_label")
        labels[c, p] = r["label"]
        for target, value in ((smiles, r["canonical_smiles"]), (scaffolds, r["scaffold_group"])):
            if target[c] not in (None, value): failures.append("compound_metadata")
            target[c] = value
        if folds[c] not in (-1, r["development_fold"]): failures.append("compound_fold")
        folds[c] = r["development_fold"]
    if any(len({folds[i] for i, value in enumerate(scaffolds) if value == scaffold}) != 1 for scaffold in set(scaffolds)):
        failures.append("scaffold_fold_leakage")
    bits, similarity = make_features(smiles)
    with np.load(SOURCE / "fingerprint_similarity.npz", allow_pickle=False) as saved:
        feature_bit_difference = int(np.count_nonzero(bits != saved["fingerprint_bits"]))
        similarity_difference = float(np.max(np.abs(similarity - saved["similarity"])))
        if list(saved["compound_ids"]) != compounds or not np.array_equal(saved["folds"], folds): failures.append("feature_order")
        if feature_bit_difference: failures.append("fingerprint_reconstruction")
        if similarity_difference: failures.append("similarity_reconstruction")
    tie = np.empty(len(compounds), dtype=int)
    for n, i in enumerate(sorted(range(len(compounds)), key=lambda value: stable(compounds[value]))): tie[i] = n

    declared_selections = {r["outer_fold"]: r for r in load(SOURCE / "nested_selections.json")}
    reconstructed_selections = {}; reconstructed_predictions = []
    for outer in range(5):
        y, p, pooled, specific = inner_arrays(outer, folds, labels, similarity, tie)
        chosen = choose(y, p, pooled, specific, len(isoforms)); reconstructed_selections[outer] = chosen
        declared = declared_selections[outer]
        for key in ("pooled_k", "isoform_k", "shrinkage_k_pool", "shrinkage_k_iso", "shrinkage_w", "inner_rows"):
            expected = chosen["rows"] if key == "inner_rows" else chosen[key]
            if declared[key] != expected: failures.append(f"selection:{outer}:{key}")
        metric_pairs = (("pooled_inner_metrics", "pooled"), ("isoform_inner_metrics", "isoform"), ("shrinkage_inner_metrics", "shrinkage"))
        for declared_name, reconstructed_name in metric_pairs:
            for field, position in (("macro_average_precision", 0), ("macro_roc_auc", 1), ("macro_brier", 2)):
                if not close(declared[declared_name][field], chosen[reconstructed_name][position]): failures.append(f"selection_metric:{outer}:{field}")
        train = np.flatnonzero(folds != outer); query = np.flatnonzero(folds == outer)
        requested = sorted({chosen["pooled_k"], chosen["isoform_k"], chosen["shrinkage_k_pool"], chosen["shrinkage_k_iso"]})
        a, b = predict(query, train, labels, similarity, tie, requested)
        for qn, q in enumerate(query):
            for protein in range(len(isoforms)):
                if not np.isfinite(labels[q, protein]): continue
                shrink = (1 - chosen["shrinkage_w"]) * a[chosen["shrinkage_k_pool"]][qn] + chosen["shrinkage_w"] * b[chosen["shrinkage_k_iso"]][qn, protein]
                reconstructed_predictions.append({"compound_number": int(q), "compound_inchikey": compounds[q], "scaffold_group": scaffolds[q],
                    "isoform_number": protein, "isoform": isoforms[protein], "label": int(labels[q, protein]), "outer_fold": outer,
                    "pooled_chemical_knn": float(a[chosen["pooled_k"]][qn]), "isoform_specific_knn": float(b[chosen["isoform_k"]][qn, protein]),
                    "data_driven_shrinkage": float(shrink)})
        print(f"verified human outer fold {outer + 1}/5", flush=True)
    reconstructed_predictions.sort(key=lambda r: (r["compound_inchikey"], r["isoform"]))
    observed = [json.loads(line) for line in (SOURCE / "outer_predictions.jsonl").read_text().splitlines()]
    if len(observed) != len(reconstructed_predictions): failures.append("prediction_count")
    for left, right in zip(observed, reconstructed_predictions):
        for field in ("compound_number", "compound_inchikey", "scaffold_group", "isoform_number", "isoform", "label", "outer_fold"):
            if left[field] != right[field]: failures.append("prediction_key")
        for method in METHODS:
            if not close(left[method], right[method]): failures.append("prediction_score")

    y = np.asarray([r["label"] for r in reconstructed_predictions], dtype=np.int8)
    protein = np.asarray([r["isoform_number"] for r in reconstructed_predictions], dtype=int)
    scores = {m: np.asarray([r[m] for r in reconstructed_predictions]) for m in METHODS}
    summary = {(r["method"], r["scope"]): r for r in load(SOURCE / "method_summary.json")["summaries"]}
    point = {}
    for method in METHODS:
        point[method] = macro(y, scores[method], protein, len(isoforms))
        declared = summary[(method, "isoform_macro")]
        for field, position in (("macro_average_precision", 0), ("macro_roc_auc", 1), ("macro_brier", 2)):
            if not close(declared[field], point[method][position]): failures.append("macro_summary")
        for p, isoform in enumerate(isoforms):
            mask = protein == p; local = summary[(method, isoform)]
            metrics = (ap(y[mask], scores[method][mask]), auc(y[mask], scores[method][mask]), float(np.mean((scores[method][mask] - y[mask]) ** 2)))
            for field, value in zip(("average_precision", "roc_auc", "brier"), metrics):
                if not close(local[field], value): failures.append("isoform_summary")

    scaffold_names = sorted(set(scaffolds)); gi = {v: i for i, v in enumerate(scaffold_names)}
    row_group = np.asarray([gi[r["scaffold_group"]] for r in reconstructed_predictions], dtype=int)
    rng = np.random.default_rng(SEED); boot = defaultdict(list)
    for replicate in range(BOOTSTRAPS):
        count = np.bincount(rng.integers(len(scaffold_names), size=len(scaffold_names)), minlength=len(scaffold_names)); weight = count[row_group]
        metrics = {m: macro(y, scores[m], protein, len(isoforms), weight) for m in METHODS}
        for right in ("pooled_chemical_knn", "isoform_specific_knn"):
            boot[(right, "average_precision")].append(metrics["data_driven_shrinkage"][0] - metrics[right][0])
            boot[(right, "roc_auc")].append(metrics["data_driven_shrinkage"][1] - metrics[right][1])
        if (replicate + 1) % 500 == 0: print(f"verified human bootstrap {replicate + 1}/{BOOTSTRAPS}", flush=True)
    declared_comparisons = {r["right_method"]: r for r in load(SOURCE / "comparison_summary.json")["comparisons"]}
    for right, declared in declared_comparisons.items():
        for metric, position in (("average_precision", 0), ("roc_auc", 1)):
            expected = {"left": point["data_driven_shrinkage"][position], "right": point[right][position],
                        "difference": point["data_driven_shrinkage"][position] - point[right][position],
                        "scaffold_bootstrap_95ci": [float(x) for x in np.quantile(boot[(right, metric)], (0.025, 0.975))]}
            for field in ("left", "right", "difference"):
                if not close(declared[metric][field], expected[field]): failures.append("comparison_point")
            if any(not close(a, b) for a, b in zip(declared[metric]["scaffold_bootstrap_95ci"], expected["scaffold_bootstrap_95ci"])):
                failures.append("comparison_interval")

    declared_control = load(SOURCE / "score_permutation_control.json"); permutation_values = []
    by_compound = defaultdict(list)
    for i, r in enumerate(reconstructed_predictions): by_compound[r["compound_number"]].append(i)
    for replicate in range(PERMUTATIONS):
        permuted = scores["data_driven_shrinkage"].copy()
        for compound, indices in by_compound.items():
            token = int.from_bytes(hashlib.sha256(f"{SEED}|{compound}|{replicate}".encode()).digest()[:8], "little")
            order = np.random.default_rng(token).permutation(len(indices)); loc = np.asarray(indices)
            permuted[loc] = scores["data_driven_shrinkage"][loc[order]]
        permutation_values.append(macro(y, permuted, protein, len(isoforms))[0] - point["pooled_chemical_knn"][0])
    observed_difference = point["data_driven_shrinkage"][0] - point["pooled_chemical_knn"][0]
    pvalue = (1 + sum(v >= observed_difference for v in permutation_values)) / 100
    if any(not close(a, b) for a, b in zip(permutation_values, declared_control["permuted_differences"])): failures.append("permutation_values")
    if declared_control["plus_one_upper_tail_monte_carlo_p"] != pvalue: failures.append("permutation_p")

    comparison_by_right = {r["right_method"]: r for r in declared_comparisons.values()}
    class_ok = all(r["positives"] >= 20 and r["negatives"] >= 20 for r in load(SOURCE / "audit.json")["class_counts"])
    signal_gates = {"positive_ap_interval_vs_pooled": comparison_by_right["pooled_chemical_knn"]["average_precision"]["scaffold_bootstrap_95ci"][0] > 0,
                    "score_permutation_p_at_most_0_05": pvalue <= 0.05,
                    "all_normalized_rows_scored": len(reconstructed_predictions) == len(rows) == 14955,
                    "at_least_20_each_class_each_isoform": class_ok}
    both_gates = {**signal_gates, "positive_ap_interval_vs_isoform_specific": comparison_by_right["isoform_specific_knn"]["average_precision"]["scaffold_bootstrap_95ci"][0] > 0}
    declared_acceptance = load(SOURCE / "acceptance.json")
    expected_acceptance = {"isoform_condition_signal": {"accepted": all(signal_gates.values()), "gates": signal_gates},
                           "shrinkage_improves_both_components": {"accepted": all(both_gates.values()), "gates": both_gates}}
    for name, value in expected_acceptance.items():
        if declared_acceptance[name] != value: failures.append("acceptance:" + name)

    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": sorted(set(failures)),
              "input_hashes_checked": len(manifest), "normalized_rows_reconstructed": len(rows),
              "fingerprint_bits_reconstructed": int(bits.size), "fingerprint_bit_differences": feature_bit_difference,
              "similarity_cells_reconstructed": int(similarity.size), "maximum_similarity_difference": similarity_difference,
              "nested_outer_selections_reconstructed": len(reconstructed_selections), "outer_scores_reconstructed": len(reconstructed_predictions) * len(METHODS),
              "method_scopes_reconstructed": len(summary), "bootstrap_replicates_reconstructed": BOOTSTRAPS,
              "score_permutations_reconstructed": PERMUTATIONS, "scaffold_fold_separation_enforced": True,
              "independent_implementation": True, "independent_biological_validation": False}
    OUTPUT.mkdir(); write_json(OUTPUT / "INDEPENDENT_QC.json", report)
    print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": main()
