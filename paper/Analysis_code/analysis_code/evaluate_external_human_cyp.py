"""One-time evaluation of the frozen external six-isoform human CYP endpoint."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import rdFingerprintGenerator

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "human_substrate_01" / "normalized_labels.json"
EXTERNAL = ROOT / "external_human_cyp_01" / "normalized_external_labels.json"
EXTERNAL_AUDIT = ROOT / "external_human_cyp_01" / "overlap_audit.json"
ACQUISITION = ROOT / "external_human_cyp_01" / "acquisition_manifest.json"
PROTOCOL = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1.md"
AMENDMENT = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1_AMENDMENT_1.md"
AMENDMENT_2 = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1_AMENDMENT_2.md"
AMENDMENT_3 = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1_AMENDMENT_3.md"
OUTPUT = ROOT / "external_human_cyp_models_01"
ISOFORM_K = 25
POOLED_K = 51
BOOTSTRAPS = 5000
PERMUTATIONS = 99
SEED = 20260917
METHODS = ("pooled_chemical_knn", "isoform_specific_knn")
STRATA = (
    "all_parseable_unique", "exact_novel", "scaffold_novel",
    "lineage_eligible", "scaffold_and_lineage_eligible",
)


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def weighted_average_precision(y, score, weight=None):
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=np.float64)
    weight = np.ones(len(y), dtype=np.float64) if weight is None else np.asarray(weight, dtype=np.float64)
    keep = weight > 0
    y, score, weight = y[keep], score[keep], weight[keep]
    positive = float(np.sum(weight * y))
    if positive <= 0:
        return math.nan
    order = np.argsort(-score, kind="mergesort")
    score, y, weight = score[order], y[order], weight[order]
    starts = np.r_[0, np.flatnonzero(score[1:] != score[:-1]) + 1]
    group_weight = np.add.reduceat(weight, starts)
    group_positive = np.add.reduceat(weight * y, starts)
    precision = np.cumsum(group_positive) / np.cumsum(group_weight)
    return float(np.sum(precision * group_positive) / positive)


def weighted_roc_auc(y, score, weight=None):
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=np.float64)
    weight = np.ones(len(y), dtype=np.float64) if weight is None else np.asarray(weight, dtype=np.float64)
    keep = weight > 0
    y, score, weight = y[keep], score[keep], weight[keep]
    positive = float(np.sum(weight * y))
    negative = float(np.sum(weight * (1 - y)))
    if positive <= 0 or negative <= 0:
        return math.nan
    order = np.argsort(score, kind="mergesort")
    score, y, weight = score[order], y[order], weight[order]
    starts = np.r_[0, np.flatnonzero(score[1:] != score[:-1]) + 1]
    group_positive = np.add.reduceat(weight * y, starts)
    group_negative = np.add.reduceat(weight * (1 - y), starts)
    below = np.r_[0.0, np.cumsum(group_negative)[:-1]]
    return float(np.sum(group_positive * (below + 0.5 * group_negative)) / (positive * negative))


def metrics(rows, method, isoforms, weights=None):
    per_isoform = []
    for isoform in isoforms:
        indices = np.asarray([n for n, row in enumerate(rows) if row["isoform"] == isoform], dtype=int)
        y = np.asarray([rows[n]["label"] for n in indices], dtype=np.int8)
        score = np.asarray([rows[n][method] for n in indices], dtype=np.float64)
        local_weight = None if weights is None else np.asarray(weights, dtype=np.float64)[indices]
        positives = int(y.sum())
        negatives = int(len(y) - positives)
        if positives < 20 or negatives < 20:
            per_isoform.append({"isoform": isoform, "rows": len(y), "positives": positives,
                                "negatives": negatives, "evaluable": False,
                                "average_precision": None, "roc_auc": None, "brier": None})
            continue
        ap = weighted_average_precision(y, score, local_weight)
        auc = weighted_roc_auc(y, score, local_weight)
        if local_weight is None:
            brier = float(np.mean((score - y) ** 2))
        else:
            brier = float(np.average((score - y) ** 2, weights=local_weight)) if local_weight.sum() else math.nan
        per_isoform.append({"isoform": isoform, "rows": len(y), "positives": positives,
                            "negatives": negatives, "evaluable": True,
                            "average_precision": ap, "roc_auc": auc, "brier": brier})
    evaluable = [row for row in per_isoform if row["evaluable"]]
    return {
        "method": method,
        "rows": len(rows),
        "evaluable_isoforms": len(evaluable),
        "macro_average_precision": float(np.nanmean([row["average_precision"] for row in evaluable])) if evaluable else None,
        "macro_roc_auc": float(np.nanmean([row["roc_auc"] for row in evaluable])) if evaluable else None,
        "macro_brier": float(np.nanmean([row["brier"] for row in evaluable])) if evaluable else None,
        "per_isoform": per_isoform,
    }


def fingerprints(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = []
    bits = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for number, value in enumerate(smiles):
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            raise ValueError(f"normalized structure failed: {number}")
        fp = generator.GetFingerprint(molecule)
        fps.append(fp)
        DataStructs.ConvertToNumpyArray(fp, bits[number])
    return fps, bits


def score_queries(train_rows, external_rows):
    all_isoforms = sorted({row["isoform"] for row in train_rows})
    shared_isoforms = sorted({row["isoform"] for row in external_rows})
    train_ids = sorted({row["compound_inchikey"] for row in train_rows})
    external_ids = sorted({row["compound_inchikey"] for row in external_rows})
    ti = {value: number for number, value in enumerate(train_ids)}
    pi = {value: number for number, value in enumerate(all_isoforms)}
    train_smiles = [None] * len(train_ids)
    labels = np.full((len(train_ids), len(all_isoforms)), np.nan, dtype=np.float64)
    for row in train_rows:
        c, p = ti[row["compound_inchikey"]], pi[row["isoform"]]
        labels[c, p] = row["label"]
        if train_smiles[c] is not None and train_smiles[c] != row["canonical_smiles"]:
            raise ValueError("inconsistent source structure")
        train_smiles[c] = row["canonical_smiles"]
    external_smiles_by_id = {}
    for row in external_rows:
        previous = external_smiles_by_id.setdefault(row["compound_inchikey"], row["canonical_smiles"])
        if previous != row["canonical_smiles"]:
            raise ValueError("inconsistent external structure")
    external_smiles = [external_smiles_by_id[value] for value in external_ids]

    train_fps, train_bits = fingerprints(train_smiles)
    external_fps, external_bits = fingerprints(external_smiles)
    similarity = np.empty((len(external_ids), len(train_ids)), dtype=np.float32)
    for number, fp in enumerate(external_fps):
        similarity[number] = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
    tie = np.empty(len(train_ids), dtype=int)
    for position, index in enumerate(sorted(range(len(train_ids)), key=lambda value: stable(train_ids[value]))):
        tie[index] = position

    pooled_labels = np.nanmean(labels, axis=1)
    pooled_fallback = float(np.nanmean(pooled_labels))
    iso_fallback = np.nanmean(labels, axis=0)
    pooled_scores = np.empty(len(external_ids), dtype=np.float64)
    specific_scores = np.empty((len(external_ids), len(shared_isoforms)), dtype=np.float64)
    for q in range(len(external_ids)):
        order = np.lexsort((tie, -similarity[q]))
        chosen = order[:min(POOLED_K, len(order))]
        weight = similarity[q, chosen].astype(np.float64)
        pooled_scores[q] = float(np.dot(weight, pooled_labels[chosen]) / weight.sum()) if weight.sum() else pooled_fallback
        for p, isoform in enumerate(shared_isoforms):
            source_p = pi[isoform]
            valid = order[~np.isnan(labels[order, source_p])]
            chosen = valid[:min(ISOFORM_K, len(valid))]
            weight = similarity[q, chosen].astype(np.float64)
            specific_scores[q, p] = float(np.dot(weight, labels[chosen, source_p]) / weight.sum()) if weight.sum() else float(iso_fallback[source_p])
        if (q + 1) % 500 == 0:
            print(f"external query {q + 1}/{len(external_ids)}", flush=True)
    return {
        "train_ids": train_ids, "external_ids": external_ids,
        "train_bits": train_bits, "external_bits": external_bits,
        "similarity": similarity, "pooled_scores": pooled_scores,
        "specific_scores": specific_scores, "shared_isoforms": shared_isoforms,
        "all_isoforms": all_isoforms,
    }


def bootstrap(rows, isoforms, stratum):
    scaffolds = sorted({row["scaffold_group"] for row in rows})
    si = {value: number for number, value in enumerate(scaffolds)}
    row_scaffold = np.asarray([si[row["scaffold_group"]] for row in rows], dtype=int)
    token = int.from_bytes(hashlib.sha256(f"{SEED}|{stratum}|bootstrap".encode()).digest()[:8], "little")
    rng = np.random.default_rng(token)
    curve_cache = {}
    for method in METHODS:
        for isoform in isoforms:
            indices = np.asarray([n for n, row in enumerate(rows) if row["isoform"] == isoform], dtype=int)
            y = np.asarray([rows[n]["label"] for n in indices], dtype=np.int8)
            score = np.asarray([rows[n][method] for n in indices], dtype=np.float64)
            descending = np.argsort(-score, kind="mergesort")
            ascending = np.argsort(score, kind="mergesort")
            curve_cache[(method, isoform)] = {
                "indices": indices,
                "y_desc": y[descending],
                "desc": descending,
                "desc_starts": np.r_[0, np.flatnonzero(score[descending][1:] != score[descending][:-1]) + 1],
                "y_asc": y[ascending],
                "asc": ascending,
                "asc_starts": np.r_[0, np.flatnonzero(score[ascending][1:] != score[ascending][:-1]) + 1],
            }

    def cached_curves(method, weights):
        aps, aucs = [], []
        for isoform in isoforms:
            cache = curve_cache[(method, isoform)]
            local = weights[cache["indices"]]
            desc_weight = local[cache["desc"]]
            desc_y = cache["y_desc"]
            total_positive = float(np.dot(desc_weight, desc_y))
            if total_positive <= 0:
                aps.append(math.nan)
            else:
                level_weight = np.add.reduceat(desc_weight, cache["desc_starts"])
                level_positive = np.add.reduceat(desc_weight * desc_y, cache["desc_starts"])
                retained = level_weight > 0
                level_weight = level_weight[retained]
                level_positive = level_positive[retained]
                precision = np.cumsum(level_positive) / np.cumsum(level_weight)
                aps.append(float(np.dot(precision, level_positive) / total_positive))
            asc_weight = local[cache["asc"]]
            asc_y = cache["y_asc"]
            positive = float(np.dot(asc_weight, asc_y))
            negative = float(np.dot(asc_weight, 1 - asc_y))
            if positive <= 0 or negative <= 0:
                aucs.append(math.nan)
            else:
                level_positive = np.add.reduceat(asc_weight * asc_y, cache["asc_starts"])
                level_negative = np.add.reduceat(asc_weight * (1 - asc_y), cache["asc_starts"])
                below = np.r_[0.0, np.cumsum(level_negative)[:-1]]
                aucs.append(float(np.dot(level_positive, below + 0.5 * level_negative) / (positive * negative)))
        return float(np.nanmean(aps)), float(np.nanmean(aucs))

    ap_differences, auc_differences = [], []
    for replicate in range(BOOTSTRAPS):
        counts = np.bincount(rng.integers(len(scaffolds), size=len(scaffolds)), minlength=len(scaffolds))
        weights = counts[row_scaffold].astype(np.float64)
        pooled_ap, pooled_auc = cached_curves("pooled_chemical_knn", weights)
        specific_ap, specific_auc = cached_curves("isoform_specific_knn", weights)
        ap_differences.append(specific_ap - pooled_ap)
        auc_differences.append(specific_auc - pooled_auc)
    pooled = metrics(rows, "pooled_chemical_knn", isoforms)
    specific = metrics(rows, "isoform_specific_knn", isoforms)
    return {
        "stratum": stratum,
        "bootstrap_replicates": BOOTSTRAPS,
        "bootstrap_unit": "external_bemis_murcko_scaffold",
        "scaffold_groups": len(scaffolds),
        "average_precision": {
            "isoform_specific": specific["macro_average_precision"],
            "pooled": pooled["macro_average_precision"],
            "difference": specific["macro_average_precision"] - pooled["macro_average_precision"],
            "scaffold_bootstrap_95ci": [float(x) for x in np.quantile(ap_differences, (0.025, 0.975))],
        },
        "roc_auc": {
            "isoform_specific": specific["macro_roc_auc"],
            "pooled": pooled["macro_roc_auc"],
            "difference": specific["macro_roc_auc"] - pooled["macro_roc_auc"],
            "scaffold_bootstrap_95ci": [float(x) for x in np.quantile(auc_differences, (0.025, 0.975))],
        },
    }


def permutation(rows, isoforms, observed_difference, pooled_ap, stratum):
    by_compound = defaultdict(list)
    for number, row in enumerate(rows):
        by_compound[row["compound_inchikey"]].append(number)
    original = np.asarray([row["isoform_specific_knn"] for row in rows], dtype=np.float64)
    values = []
    for replicate in range(PERMUTATIONS):
        permuted = original.copy()
        for compound, indices in by_compound.items():
            token = int.from_bytes(hashlib.sha256(f"{SEED}|{stratum}|{compound}|{replicate}".encode()).digest()[:8], "little")
            order = np.random.default_rng(token).permutation(len(indices))
            permuted[np.asarray(indices)] = original[np.asarray(indices)[order]]
        copied = [{**row, "permuted": float(permuted[n])} for n, row in enumerate(rows)]
        macro = metrics(copied, "permuted", isoforms)["macro_average_precision"]
        values.append(macro - pooled_ap)
    pvalue = (1 + sum(value >= observed_difference for value in values)) / (PERMUTATIONS + 1)
    return {
        "stratum": stratum,
        "permutations": PERMUTATIONS,
        "unit": "isoform_specific_scores_within_compound_available_isoforms",
        "refitted": False,
        "observed_macro_average_precision_difference_vs_pooled": observed_difference,
        "permuted_differences": values,
        "plus_one_upper_tail_monte_carlo_p": pvalue,
    }


def main():
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to replace {OUTPUT}")
    RDLogger.DisableLog("rdApp.*")
    source_rows = read(SOURCE)
    external_rows = [row for row in read(EXTERNAL) if row["strata"]["all_parseable_unique"]]
    audit = read(EXTERNAL_AUDIT)
    acquisition = read(ACQUISITION)
    if audit["status"] != "PASS" or not acquisition["complete"]:
        raise RuntimeError("external inputs did not pass audit")
    scored = score_queries(source_rows, external_rows)
    ei = {value: number for number, value in enumerate(scored["external_ids"])}
    pi = {value: number for number, value in enumerate(scored["shared_isoforms"])}
    predictions = []
    for row in external_rows:
        c, p = ei[row["compound_inchikey"]], pi[row["isoform"]]
        predictions.append({
            "isoform": row["isoform"], "compound_inchikey": row["compound_inchikey"],
            "scaffold_group": row["scaffold_group"], "label": row["label"],
            "sources": row["sources"], "strata": row["strata"],
            "pooled_chemical_knn": float(scored["pooled_scores"][c]),
            "isoform_specific_knn": float(scored["specific_scores"][c, p]),
        })
    predictions.sort(key=lambda row: (row["compound_inchikey"], row["isoform"]))

    summaries, comparisons, controls = {}, {}, {}
    for stratum in STRATA:
        rows = [row for row in predictions if row["strata"][stratum]]
        summaries[stratum] = {method: metrics(rows, method, scored["shared_isoforms"]) for method in METHODS}
        comparisons[stratum] = bootstrap(rows, scored["shared_isoforms"], stratum)
        observed = comparisons[stratum]["average_precision"]["difference"]
        pooled_ap = comparisons[stratum]["average_precision"]["pooled"]
        controls[stratum] = permutation(rows, scored["shared_isoforms"], observed, pooled_ap, stratum)
        print(f"external stratum complete: {stratum}", flush=True)

    def gate(stratum):
        info = audit["strata"][stratum]
        comparison = comparisons[stratum]
        control = controls[stratum]
        class_ok = sum(row["positives"] >= 20 and row["negatives"] >= 20 for row in info["by_isoform"]) >= 4
        gates = {
            "at_least_four_isoforms_with_20_each_class": class_ok,
            "at_least_100_unique_compounds": info["unique_compounds"] >= 100,
            "at_least_20_scaffold_groups": info["scaffold_groups"] >= 20,
            "positive_ap_interval_vs_pooled": comparison["average_precision"]["scaffold_bootstrap_95ci"][0] > 0,
            "score_permutation_p_at_most_0_05": control["plus_one_upper_tail_monte_carlo_p"] <= 0.05,
        }
        return {"accepted": all(gates.values()), "gates": gates}

    independent = gate("lineage_eligible")
    scaffold = gate("scaffold_and_lineage_eligible")
    if not independent["accepted"]:
        scaffold["accepted"] = False
        scaffold["gates"]["independent_external_transfer_accepted"] = False
    else:
        scaffold["gates"]["independent_external_transfer_accepted"] = True
    acceptance = {
        "independent_external_transfer": independent,
        "scaffold_external_transfer": scaffold,
        "source_labels_not_assay_harmonized": True,
        "fixed_human_isoforms_only": True,
        "no_new_protein_or_exact_reaction_claim": True,
    }

    OUTPUT.mkdir()
    with (OUTPUT / "external_predictions.jsonl").open("w", encoding="utf-8") as stream:
        for row in predictions:
            stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + "\n")
    write_json(OUTPUT / "method_summary.json", {"status": "COMPUTED_QC_PENDING", "strata": summaries})
    write_json(OUTPUT / "comparison_summary.json", {"status": "COMPUTED_QC_PENDING", "strata": comparisons})
    write_json(OUTPUT / "score_permutation_control.json", {"status": "COMPUTED_QC_PENDING", "strata": controls})
    write_json(OUTPUT / "acceptance.json", {"status": "COMPUTED_QC_PENDING", **acceptance})
    np.savez_compressed(
        OUTPUT / "fingerprint_similarity.npz",
        train_ids=np.asarray(scored["train_ids"]), external_ids=np.asarray(scored["external_ids"]),
        train_bits=scored["train_bits"], external_bits=scored["external_bits"],
        similarity=scored["similarity"],
    )
    run_audit = {
        "status": "COMPUTED_QC_PENDING", "created_utc": now(),
        "source_rows": len(source_rows), "external_scored_rows": len(predictions),
        "source_compounds": len(scored["train_ids"]), "external_compounds": len(scored["external_ids"]),
        "source_isoforms": scored["all_isoforms"], "external_isoforms": scored["shared_isoforms"],
        "isoform_k": ISOFORM_K, "pooled_k": POOLED_K,
        "morgan_radius": 2, "morgan_bits": 2048, "rdkit_version": rdBase.rdkitVersion,
        "source_labels_not_assay_harmonized": True,
    }
    write_json(OUTPUT / "audit.json", run_audit)
    inputs = [SOURCE, EXTERNAL, EXTERNAL_AUDIT, ACQUISITION, PROTOCOL, AMENDMENT, AMENDMENT_2, AMENDMENT_3,
              Path(__file__), ROOT / "run_raw.py"]
    write_json(OUTPUT / "input_manifest.json", {str(path.resolve()): digest_file(path) for path in inputs})
    print(json.dumps({
        "audit": run_audit,
        "comparisons": comparisons,
        "permutation_p": {key: value["plus_one_upper_tail_monte_carlo_p"] for key, value in controls.items()},
        "acceptance": acceptance,
    }, indent=2))


if __name__ == "__main__":
    main()
