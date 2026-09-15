"""Nested scaffold-held-out evaluation of the fixed nine-isoform endpoint."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "human_substrate_01"
OUTPUT = ROOT / "human_substrate_models_01"
KS = (1, 3, 5, 11, 25, 51, 101)
WEIGHTS = tuple(round(value / 10, 1) for value in range(11))
BOOTSTRAPS = 5000
PERMUTATIONS = 99
SEED = 20260916
METHODS = ("pooled_chemical_knn", "isoform_specific_knn", "data_driven_shrinkage")


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def weighted_average_precision(y, score, weight=None):
    y = np.asarray(y, dtype=np.int8); score = np.asarray(score, dtype=np.float64)
    weight = np.ones(len(y), dtype=np.float64) if weight is None else np.asarray(weight, dtype=np.float64)
    retained = weight > 0
    y, score, weight = y[retained], score[retained], weight[retained]
    positive = float(np.sum(weight * y))
    if positive <= 0: return math.nan
    order = np.argsort(-score, kind="mergesort")
    s = score[order]; yy = y[order]; ww = weight[order]
    starts = np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1]
    group_weight = np.add.reduceat(ww, starts)
    group_positive = np.add.reduceat(ww * yy, starts)
    precision = np.cumsum(group_positive) / np.cumsum(group_weight)
    return float(np.sum(precision * group_positive) / positive)


def weighted_roc_auc(y, score, weight=None):
    y = np.asarray(y, dtype=np.int8); score = np.asarray(score, dtype=np.float64)
    weight = np.ones(len(y), dtype=np.float64) if weight is None else np.asarray(weight, dtype=np.float64)
    retained = weight > 0
    y, score, weight = y[retained], score[retained], weight[retained]
    positive = float(np.sum(weight * y)); negative = float(np.sum(weight * (1 - y)))
    if positive <= 0 or negative <= 0: return math.nan
    order = np.argsort(score, kind="mergesort")
    s = score[order]; yy = y[order]; ww = weight[order]
    starts = np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1]
    group_positive = np.add.reduceat(ww * yy, starts)
    group_negative = np.add.reduceat(ww * (1 - yy), starts)
    negatives_below = np.r_[0.0, np.cumsum(group_negative)[:-1]]
    return float(np.sum(group_positive * (negatives_below + 0.5 * group_negative)) / (positive * negative))


def macro_metrics(y, score, protein, isoform_count, weight=None):
    aps, aucs, briers = [], [], []
    per_isoform = []
    for p in range(isoform_count):
        mask = protein == p
        local_weight = None if weight is None else weight[mask]
        ap = weighted_average_precision(y[mask], score[mask], local_weight)
        auc = weighted_roc_auc(y[mask], score[mask], local_weight)
        if local_weight is None:
            brier = float(np.mean((score[mask] - y[mask]) ** 2))
        else:
            brier = float(np.average((score[mask] - y[mask]) ** 2, weights=local_weight)) if local_weight.sum() else math.nan
        aps.append(ap); aucs.append(auc); briers.append(brier)
        per_isoform.append((ap, auc, brier))
    return {"macro_average_precision": float(np.nanmean(aps)),
            "macro_roc_auc": float(np.nanmean(aucs)),
            "macro_brier": float(np.nanmean(briers)),
            "per_isoform": per_isoform}


def fingerprint_and_similarity(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = []
    bits = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for number, value in enumerate(smiles):
        molecule = Chem.MolFromSmiles(value)
        if molecule is None: raise ValueError(f"normalized structure failed: {number}")
        fingerprint = generator.GetFingerprint(molecule)
        fingerprints.append(fingerprint)
        DataStructs.ConvertToNumpyArray(fingerprint, bits[number])
    similarity = np.empty((len(smiles), len(smiles)), dtype=np.float32)
    for number, fingerprint in enumerate(fingerprints):
        similarity[number] = DataStructs.BulkTanimotoSimilarity(fingerprint, fingerprints)
    return bits, similarity


def component_predictions(query, train, labels, similarity, tie, ks):
    """Return pooled and isoform-specific predictions for every requested k."""
    train = np.asarray(train, dtype=int)
    pooled_label = np.nanmean(labels[train], axis=1)
    pooled_fallback = float(np.nanmean(pooled_label))
    isoform_fallback = np.nanmean(labels[train], axis=0)
    pool = {k: np.empty(len(query), dtype=np.float64) for k in ks}
    iso = {k: np.empty((len(query), labels.shape[1]), dtype=np.float64) for k in ks}
    for q_position, compound in enumerate(query):
        order = train[np.lexsort((tie[train], -similarity[compound, train]))]
        pooled_values = np.nanmean(labels[order], axis=1)
        for k in ks:
            chosen = order[:min(k, len(order))]
            weights = similarity[compound, chosen].astype(np.float64)
            denominator = float(weights.sum())
            pool[k][q_position] = float(np.dot(weights, np.nanmean(labels[chosen], axis=1)) / denominator) if denominator > 0 else pooled_fallback
        for p in range(labels.shape[1]):
            valid_order = order[~np.isnan(labels[order, p])]
            for k in ks:
                chosen = valid_order[:min(k, len(valid_order))]
                weights = similarity[compound, chosen].astype(np.float64)
                denominator = float(weights.sum())
                iso[k][q_position, p] = float(np.dot(weights, labels[chosen, p]) / denominator) if denominator > 0 else float(isoform_fallback[p])
    return pool, iso


def rows_from_components(query, labels, pool, iso, compounds, isoforms, fold):
    records = []
    for q_position, compound in enumerate(query):
        for p in range(labels.shape[1]):
            if np.isnan(labels[compound, p]): continue
            records.append({"compound_number": int(compound), "compound_inchikey": compounds[compound],
                            "isoform_number": p, "isoform": isoforms[p], "label": int(labels[compound, p]),
                            "fold": int(fold), "pool": {str(k): float(pool[k][q_position]) for k in KS},
                            "iso": {str(k): float(iso[k][q_position, p]) for k in KS}})
    return records


def arrays(records, field, k):
    y = np.asarray([row["label"] for row in records], dtype=np.int8)
    p = np.asarray([row["isoform_number"] for row in records], dtype=int)
    score = np.asarray([row[field][str(k)] for row in records], dtype=np.float64)
    return y, p, score


def select_parameters(records, isoform_count):
    y = np.asarray([row["label"] for row in records], dtype=np.int8)
    protein = np.asarray([row["isoform_number"] for row in records], dtype=int)
    pool_scores = {k: np.asarray([row["pool"][str(k)] for row in records]) for k in KS}
    iso_scores = {k: np.asarray([row["iso"][str(k)] for row in records]) for k in KS}
    pool_results = {k: macro_metrics(y, score, protein, isoform_count) for k, score in pool_scores.items()}
    iso_results = {k: macro_metrics(y, score, protein, isoform_count) for k, score in iso_scores.items()}
    pool_k = min(KS, key=lambda k: (-pool_results[k]["macro_average_precision"], -pool_results[k]["macro_roc_auc"], k))
    iso_k = min(KS, key=lambda k: (-iso_results[k]["macro_average_precision"], -iso_results[k]["macro_roc_auc"], k))
    best = None
    for kp in KS:
        for ki in KS:
            for w in WEIGHTS:
                score = (1 - w) * pool_scores[kp] + w * iso_scores[ki]
                result = macro_metrics(y, score, protein, isoform_count)
                key = (-result["macro_average_precision"], -result["macro_roc_auc"], kp + ki, w, kp, ki)
                if best is None or key < best[0]: best = (key, kp, ki, w, result)
    return {"pooled_k": pool_k, "pooled_inner_metrics": pool_results[pool_k],
            "isoform_k": iso_k, "isoform_inner_metrics": iso_results[iso_k],
            "shrinkage_k_pool": best[1], "shrinkage_k_iso": best[2], "shrinkage_w": best[3],
            "shrinkage_inner_metrics": best[4], "inner_rows": len(records)}


def summarize(predictions, isoforms):
    y = np.asarray([row["label"] for row in predictions], dtype=np.int8)
    protein = np.asarray([row["isoform_number"] for row in predictions], dtype=int)
    summaries = []
    for method in METHODS:
        score = np.asarray([row[method] for row in predictions], dtype=np.float64)
        result = macro_metrics(y, score, protein, len(isoforms))
        summaries.append({"method": method, "scope": "isoform_macro", "rows": len(y),
                          **{key: result[key] for key in ("macro_average_precision", "macro_roc_auc", "macro_brier")}})
        for p, isoform in enumerate(isoforms):
            mask = protein == p; ap, auc, brier = result["per_isoform"][p]
            summaries.append({"method": method, "scope": isoform, "rows": int(mask.sum()),
                              "positives": int(y[mask].sum()), "negatives": int(mask.sum() - y[mask].sum()),
                              "average_precision": ap, "roc_auc": auc, "brier": brier})
    return summaries


def bootstrap_comparisons(predictions, scaffold_by_compound, isoform_count):
    y = np.asarray([row["label"] for row in predictions], dtype=np.int8)
    protein = np.asarray([row["isoform_number"] for row in predictions], dtype=int)
    scores = {method: np.asarray([row[method] for row in predictions]) for method in METHODS}
    scaffold_names = sorted(set(scaffold_by_compound))
    scaffold_index = {value: number for number, value in enumerate(scaffold_names)}
    row_scaffold = np.asarray([scaffold_index[scaffold_by_compound[row["compound_number"]]] for row in predictions], dtype=int)
    rng = np.random.default_rng(SEED)
    replicates = defaultdict(list)
    for replicate in range(BOOTSTRAPS):
        counts = np.bincount(rng.integers(len(scaffold_names), size=len(scaffold_names)), minlength=len(scaffold_names))
        row_weight = counts[row_scaffold].astype(np.float64)
        metrics = {method: macro_metrics(y, scores[method], protein, isoform_count, row_weight) for method in METHODS}
        for right in ("pooled_chemical_knn", "isoform_specific_knn"):
            replicates[(right, "average_precision")].append(metrics["data_driven_shrinkage"]["macro_average_precision"] - metrics[right]["macro_average_precision"])
            replicates[(right, "roc_auc")].append(metrics["data_driven_shrinkage"]["macro_roc_auc"] - metrics[right]["macro_roc_auc"])
        if (replicate + 1) % 500 == 0: print(f"human bootstrap {replicate + 1}/{BOOTSTRAPS}", flush=True)
    comparisons = []
    summary_index = {(row["method"], row["scope"]): row for row in summarize(predictions, [str(i) for i in range(isoform_count)])}
    for right in ("pooled_chemical_knn", "isoform_specific_knn"):
        left = summary_index[("data_driven_shrinkage", "isoform_macro")]
        other = summary_index[(right, "isoform_macro")]
        item = {"left_method": "data_driven_shrinkage", "right_method": right}
        for metric, field in (("average_precision", "macro_average_precision"), ("roc_auc", "macro_roc_auc")):
            values = np.asarray(replicates[(right, metric)])
            item[metric] = {"left": left[field], "right": other[field], "difference": left[field] - other[field],
                            "scaffold_bootstrap_95ci": [float(x) for x in np.quantile(values, (0.025, 0.975))]}
        comparisons.append(item)
    return comparisons, len(scaffold_names)


def permutation_control(predictions, isoform_count, observed_difference, pooled_macro_ap):
    y = np.asarray([row["label"] for row in predictions], dtype=np.int8)
    protein = np.asarray([row["isoform_number"] for row in predictions], dtype=int)
    original = np.asarray([row["data_driven_shrinkage"] for row in predictions], dtype=np.float64)
    by_compound = defaultdict(list)
    for number, row in enumerate(predictions): by_compound[row["compound_number"]].append(number)
    values = []
    for replicate in range(PERMUTATIONS):
        permuted = original.copy()
        for compound, indices in by_compound.items():
            token = int.from_bytes(hashlib.sha256(f"{SEED}|{compound}|{replicate}".encode()).digest()[:8], "little")
            order = np.random.default_rng(token).permutation(len(indices))
            permuted[np.asarray(indices)] = original[np.asarray(indices)[order]]
        macro = macro_metrics(y, permuted, protein, isoform_count)["macro_average_precision"]
        values.append(macro - pooled_macro_ap)
    pvalue = (1 + sum(value >= observed_difference for value in values)) / (PERMUTATIONS + 1)
    return {"permutations": PERMUTATIONS, "unit": "conditional_scores_within_compound_available_isoforms",
            "refitted": False, "observed_macro_average_precision_difference_vs_pooled": observed_difference,
            "permuted_differences": values, "plus_one_upper_tail_monte_carlo_p": pvalue}


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    rows = read(SOURCE / "normalized_labels.json")
    source_audit = read(SOURCE / "human_substrate_audit.json")
    isoforms = sorted({row["isoform"] for row in rows})
    compounds = sorted({row["compound_inchikey"] for row in rows})
    ci = {value: number for number, value in enumerate(compounds)}; pi = {value: number for number, value in enumerate(isoforms)}
    labels = np.full((len(compounds), len(isoforms)), np.nan, dtype=np.float64)
    smiles = [None] * len(compounds); scaffolds = [None] * len(compounds); folds = np.full(len(compounds), -1, dtype=int)
    for row in rows:
        c = ci[row["compound_inchikey"]]; p = pi[row["isoform"]]
        if not np.isnan(labels[c, p]): raise ValueError("duplicate normalized label")
        labels[c, p] = row["label"]
        for target, value, label in ((smiles, row["canonical_smiles"], "smiles"), (scaffolds, row["scaffold_group"], "scaffold")):
            if target[c] is not None and target[c] != value: raise ValueError("compound " + label)
            target[c] = value
        if folds[c] >= 0 and folds[c] != row["development_fold"]: raise ValueError("compound fold")
        folds[c] = row["development_fold"]
    bits, similarity = fingerprint_and_similarity(smiles)
    if not np.allclose(np.diag(similarity), 1): raise ValueError("similarity diagonal")
    tie = np.empty(len(compounds), dtype=int)
    for position, index in enumerate(sorted(range(len(compounds)), key=lambda value: stable(compounds[value]))): tie[index] = position

    predictions, selections = [], []
    for outer in range(5):
        inner_records = []
        for validation in range(5):
            if validation == outer: continue
            train = np.flatnonzero((folds != outer) & (folds != validation)); query = np.flatnonzero(folds == validation)
            pool, iso = component_predictions(query, train, labels, similarity, tie, KS)
            inner_records.extend(rows_from_components(query, labels, pool, iso, compounds, isoforms, validation))
        selection = select_parameters(inner_records, len(isoforms)); selection["outer_fold"] = outer
        selections.append(selection)
        train = np.flatnonzero(folds != outer); query = np.flatnonzero(folds == outer)
        requested = sorted(set((selection["pooled_k"], selection["isoform_k"], selection["shrinkage_k_pool"], selection["shrinkage_k_iso"])))
        pool, iso = component_predictions(query, train, labels, similarity, tie, requested)
        for q_position, compound in enumerate(query):
            for p in range(len(isoforms)):
                if np.isnan(labels[compound, p]): continue
                pooled = float(pool[selection["pooled_k"]][q_position])
                specific = float(iso[selection["isoform_k"]][q_position, p])
                shrink_pool = float(pool[selection["shrinkage_k_pool"]][q_position])
                shrink_iso = float(iso[selection["shrinkage_k_iso"]][q_position, p])
                predictions.append({"compound_number": int(compound), "compound_inchikey": compounds[compound],
                    "scaffold_group": scaffolds[compound], "isoform_number": p, "isoform": isoforms[p],
                    "label": int(labels[compound, p]), "outer_fold": outer,
                    "pooled_chemical_knn": pooled, "isoform_specific_knn": specific,
                    "data_driven_shrinkage": (1 - selection["shrinkage_w"]) * shrink_pool + selection["shrinkage_w"] * shrink_iso})
        print(f"human outer fold {outer + 1}/5", flush=True)

    predictions.sort(key=lambda row: (row["compound_inchikey"], row["isoform"]))
    summaries = summarize(predictions, isoforms)
    comparisons, scaffold_count = bootstrap_comparisons(predictions, scaffolds, len(isoforms))
    comparison_index = {row["right_method"]: row for row in comparisons}
    observed = comparison_index["pooled_chemical_knn"]["average_precision"]["difference"]
    pooled_ap = comparison_index["pooled_chemical_knn"]["average_precision"]["right"]
    control = permutation_control(predictions, len(isoforms), observed, pooled_ap)
    class_counts = []
    for p, isoform in enumerate(isoforms):
        values = labels[:, p]; valid = ~np.isnan(values)
        class_counts.append({"isoform": isoform, "rows": int(valid.sum()), "positives": int(np.nansum(values)),
                             "negatives": int(valid.sum() - np.nansum(values))})
    gates_signal = {
        "positive_ap_interval_vs_pooled": comparison_index["pooled_chemical_knn"]["average_precision"]["scaffold_bootstrap_95ci"][0] > 0,
        "score_permutation_p_at_most_0_05": control["plus_one_upper_tail_monte_carlo_p"] <= 0.05,
        "all_normalized_rows_scored": len(predictions) == len(rows) == 14955,
        "at_least_20_each_class_each_isoform": all(row["positives"] >= 20 and row["negatives"] >= 20 for row in class_counts),
    }
    gates_both = {**gates_signal,
                  "positive_ap_interval_vs_isoform_specific": comparison_index["isoform_specific_knn"]["average_precision"]["scaffold_bootstrap_95ci"][0] > 0}
    acceptance = {"isoform_condition_signal": {"accepted": all(gates_signal.values()), "gates": gates_signal},
                  "shrinkage_improves_both_components": {"accepted": all(gates_both.values()), "gates": gates_both},
                  "development_endpoint_not_independent_validation": True}

    OUTPUT.mkdir()
    with (OUTPUT / "outer_predictions.jsonl").open("w", encoding="utf-8") as handle:
        for row in predictions: handle.write(json.dumps(row, allow_nan=False) + "\n")
    write_json(OUTPUT / "nested_selections.json", selections)
    write_json(OUTPUT / "method_summary.json", {"status": "COMPUTED_QC_PENDING", "summaries": summaries})
    write_json(OUTPUT / "comparison_summary.json", {"status": "COMPUTED_QC_PENDING", "bootstrap_replicates": BOOTSTRAPS,
                                                       "bootstrap_unit": "scaffold_component", "comparisons": comparisons})
    write_json(OUTPUT / "score_permutation_control.json", {"status": "COMPUTED_QC_PENDING", **control})
    write_json(OUTPUT / "acceptance.json", {"status": "COMPUTED_QC_PENDING", **acceptance})
    np.savez_compressed(OUTPUT / "fingerprint_similarity.npz", compound_ids=np.asarray(compounds),
                        fingerprint_bits=bits, similarity=similarity, folds=folds)
    audit = {"status": "COMPUTED_QC_PENDING", "created_utc": now(), "normalized_rows": len(rows),
             "compounds": len(compounds), "scaffold_components": scaffold_count, "isoforms": isoforms,
             "class_counts": class_counts, "outer_predictions": len(predictions), "inner_rows_by_outer": [r["inner_rows"] for r in selections],
             "rdkit_version": rdBase.rdkitVersion, "morgan_radius": 2, "morgan_bits": 2048,
             "source_unresolved_rows_not_scored": source_audit["counts"]["structure_unresolved_rows"],
             "source_labels_not_assay_harmonized": True, "independent_validation": False}
    write_json(OUTPUT / "audit.json", audit)
    raw_source = ROOT / "raw_additions/molecules-26-04678-s001.zip"
    inputs = [ROOT / "HUMAN_SUBSTRATE_ENDPOINT_V1.md", Path(__file__), ROOT / "run_raw.py",
              SOURCE / "normalized_labels.json", SOURCE / "human_substrate_audit.json"]
    if raw_source.exists(): inputs.append(raw_source)
    write_json(OUTPUT / "input_manifest.json", {str(path.resolve()): digest_file(path) for path in inputs})
    print(json.dumps({"audit": audit, "selections": selections, "comparisons": comparisons,
                      "permutation_p": control["plus_one_upper_tail_monte_carlo_p"], "acceptance": acceptance}, indent=2))


if __name__ == "__main__": main()
