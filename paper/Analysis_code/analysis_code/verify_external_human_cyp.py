"""Independent reconstruction of the frozen external human CYP evaluation."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import AllChem

from run_raw import digest_file, now, stable, write_json


ROOT = Path(__file__).resolve().parent
TRAIN_PATH = ROOT / "human_substrate_01" / "normalized_labels.json"
TEST_PATH = ROOT / "external_human_cyp_01" / "normalized_external_labels.json"
AUDIT_PATH = ROOT / "external_human_cyp_01" / "overlap_audit.json"
RESULTS = ROOT / "external_human_cyp_models_01"
VALIDATION = ROOT / "external_human_cyp_validation_01"
STRATA = (
    "all_parseable_unique", "exact_novel", "scaffold_novel",
    "lineage_eligible", "scaffold_and_lineage_eligible",
)
METHODS = ("pooled_chemical_knn", "isoform_specific_knn")
SEED = 20260917
BOOTSTRAPS = 5000
PERMUTATIONS = 99


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def assert_close(left, right, path="root", tolerance=1e-11):
    if isinstance(left, dict) and isinstance(right, dict):
        if set(left) != set(right):
            raise AssertionError(f"keys differ at {path}: {set(left) ^ set(right)}")
        for key in left:
            assert_close(left[key], right[key], f"{path}.{key}", tolerance)
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            raise AssertionError(f"length differs at {path}")
        for number, (a, b) in enumerate(zip(left, right)):
            assert_close(a, b, f"{path}[{number}]", tolerance)
    elif isinstance(left, (int, float)) and not isinstance(left, bool) and isinstance(right, (int, float)) and not isinstance(right, bool):
        if not math.isclose(float(left), float(right), rel_tol=tolerance, abs_tol=tolerance):
            raise AssertionError(f"numeric mismatch at {path}: {left} != {right}")
    elif left != right:
        raise AssertionError(f"mismatch at {path}: {left!r} != {right!r}")


def average_precision(labels, scores, weights=None):
    labels = np.asarray(labels, dtype=np.int8)
    scores = np.asarray(scores, dtype=float)
    weights = np.ones(labels.size, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    mask = weights > 0
    labels, scores, weights = labels[mask], scores[mask], weights[mask]
    total_positive = np.dot(weights, labels)
    if total_positive == 0:
        return math.nan
    order = np.argsort(-scores, kind="mergesort")
    labels, scores, weights = labels[order], scores[order], weights[order]
    boundaries = np.r_[0, np.nonzero(np.diff(scores))[0] + 1]
    positive_at_level = np.add.reduceat(weights * labels, boundaries)
    weight_at_level = np.add.reduceat(weights, boundaries)
    precision = np.cumsum(positive_at_level) / np.cumsum(weight_at_level)
    return float(np.dot(precision, positive_at_level) / total_positive)


def roc_auc(labels, scores, weights=None):
    labels = np.asarray(labels, dtype=np.int8)
    scores = np.asarray(scores, dtype=float)
    weights = np.ones(labels.size, dtype=float) if weights is None else np.asarray(weights, dtype=float)
    mask = weights > 0
    labels, scores, weights = labels[mask], scores[mask], weights[mask]
    positive = np.dot(weights, labels)
    negative = np.dot(weights, 1 - labels)
    if positive == 0 or negative == 0:
        return math.nan
    order = np.argsort(scores, kind="mergesort")
    labels, scores, weights = labels[order], scores[order], weights[order]
    boundaries = np.r_[0, np.nonzero(np.diff(scores))[0] + 1]
    positives = np.add.reduceat(weights * labels, boundaries)
    negatives = np.add.reduceat(weights * (1 - labels), boundaries)
    lower_negatives = np.r_[0.0, np.cumsum(negatives)[:-1]]
    return float(np.dot(positives, lower_negatives + 0.5 * negatives) / (positive * negative))


def method_metrics(rows, method, isoforms, weights=None):
    details = []
    for isoform in isoforms:
        positions = [number for number, row in enumerate(rows) if row["isoform"] == isoform]
        y = np.asarray([rows[number]["label"] for number in positions], dtype=np.int8)
        score = np.asarray([rows[number][method] for number in positions], dtype=float)
        local_weights = None if weights is None else np.asarray(weights, dtype=float)[positions]
        positive, negative = int(y.sum()), int(y.size - y.sum())
        if positive < 20 or negative < 20:
            details.append({"isoform": isoform, "rows": len(y), "positives": positive,
                            "negatives": negative, "evaluable": False,
                            "average_precision": None, "roc_auc": None, "brier": None})
            continue
        ap = average_precision(y, score, local_weights)
        auc = roc_auc(y, score, local_weights)
        brier = float(np.mean((score - y) ** 2)) if local_weights is None else (
            float(np.average((score - y) ** 2, weights=local_weights)) if local_weights.sum() else math.nan
        )
        details.append({"isoform": isoform, "rows": len(y), "positives": positive,
                        "negatives": negative, "evaluable": True,
                        "average_precision": ap, "roc_auc": auc, "brier": brier})
    evaluable = [row for row in details if row["evaluable"]]
    return {
        "method": method, "rows": len(rows), "evaluable_isoforms": len(evaluable),
        "macro_average_precision": float(np.mean([row["average_precision"] for row in evaluable])) if evaluable else None,
        "macro_roc_auc": float(np.mean([row["roc_auc"] for row in evaluable])) if evaluable else None,
        "macro_brier": float(np.mean([row["brier"] for row in evaluable])) if evaluable else None,
        "per_isoform": details,
    }


def make_fingerprints(smiles):
    fingerprints = []
    bits = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for number, value in enumerate(smiles):
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            raise AssertionError(f"invalid normalized structure {number}")
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, 2, nBits=2048)
        fingerprints.append(fingerprint)
        DataStructs.ConvertToNumpyArray(fingerprint, bits[number])
    return fingerprints, bits


def reconstruct_predictions(train_rows, external_rows, saved_arrays):
    train_ids = sorted({row["compound_inchikey"] for row in train_rows})
    query_ids = sorted({row["compound_inchikey"] for row in external_rows})
    all_isoforms = sorted({row["isoform"] for row in train_rows})
    shared_isoforms = sorted({row["isoform"] for row in external_rows})
    if train_ids != saved_arrays["train_ids"].tolist() or query_ids != saved_arrays["external_ids"].tolist():
        raise AssertionError("saved array identifiers differ")
    ti = {value: number for number, value in enumerate(train_ids)}
    pi = {value: number for number, value in enumerate(all_isoforms)}
    qi = {value: number for number, value in enumerate(query_ids)}
    training_smiles = [None] * len(train_ids)
    labels = np.full((len(train_ids), len(all_isoforms)), np.nan)
    for row in train_rows:
        c, p = ti[row["compound_inchikey"]], pi[row["isoform"]]
        labels[c, p] = row["label"]
        training_smiles[c] = row["canonical_smiles"]
    query_smiles = [None] * len(query_ids)
    for row in external_rows:
        c = qi[row["compound_inchikey"]]
        if query_smiles[c] is not None and query_smiles[c] != row["canonical_smiles"]:
            raise AssertionError("external identity has multiple representations")
        query_smiles[c] = row["canonical_smiles"]

    train_fps, train_bits = make_fingerprints(training_smiles)
    query_fps, query_bits = make_fingerprints(query_smiles)
    if not np.array_equal(train_bits, saved_arrays["train_bits"]):
        raise AssertionError("training fingerprints differ")
    if not np.array_equal(query_bits, saved_arrays["external_bits"]):
        raise AssertionError("external fingerprints differ")
    similarities = np.vstack([
        np.asarray(DataStructs.BulkTanimotoSimilarity(fp, train_fps), dtype=np.float32)
        for fp in query_fps
    ])
    maximum_similarity_difference = float(np.max(np.abs(similarities - saved_arrays["similarity"])))
    if maximum_similarity_difference != 0:
        raise AssertionError("similarity matrix differs")

    tiebreak = np.empty(len(train_ids), dtype=int)
    for rank, position in enumerate(sorted(range(len(train_ids)), key=lambda index: stable(train_ids[index]))):
        tiebreak[position] = rank
    pooled_target = np.nanmean(labels, axis=1)
    pooled_default = float(np.nanmean(pooled_target))
    isoform_default = np.nanmean(labels, axis=0)
    scores = {}
    for c, compound in enumerate(query_ids):
        order = np.lexsort((tiebreak, -similarities[c]))
        selected = order[:51]
        weight = similarities[c, selected].astype(float)
        pooled = float(np.dot(weight, pooled_target[selected]) / weight.sum()) if weight.sum() else pooled_default
        for isoform in shared_isoforms:
            p = pi[isoform]
            eligible = order[~np.isnan(labels[order, p])][:25]
            weight = similarities[c, eligible].astype(float)
            specific = float(np.dot(weight, labels[eligible, p]) / weight.sum()) if weight.sum() else float(isoform_default[p])
            scores[(compound, isoform)] = (pooled, specific)

    predictions = []
    for row in external_rows:
        pooled, specific = scores[(row["compound_inchikey"], row["isoform"])]
        predictions.append({
            "isoform": row["isoform"], "compound_inchikey": row["compound_inchikey"],
            "scaffold_group": row["scaffold_group"], "label": row["label"],
            "sources": row["sources"], "strata": row["strata"],
            "pooled_chemical_knn": pooled, "isoform_specific_knn": specific,
        })
    predictions.sort(key=lambda row: (row["compound_inchikey"], row["isoform"]))
    return predictions, shared_isoforms, maximum_similarity_difference, int(train_bits.sum()), int(query_bits.sum())


def reconstruct_comparison(rows, isoforms, stratum):
    scaffold_names = sorted({row["scaffold_group"] for row in rows})
    scaffold_number = {name: number for number, name in enumerate(scaffold_names)}
    row_group = np.asarray([scaffold_number[row["scaffold_group"]] for row in rows], dtype=int)
    token = int.from_bytes(hashlib.sha256(f"{SEED}|{stratum}|bootstrap".encode()).digest()[:8], "little")
    generator = np.random.default_rng(token)
    cache = {}
    for method in METHODS:
        for isoform in isoforms:
            positions = np.asarray([n for n, row in enumerate(rows) if row["isoform"] == isoform], dtype=int)
            labels = np.asarray([rows[n]["label"] for n in positions], dtype=np.int8)
            scores = np.asarray([rows[n][method] for n in positions], dtype=float)
            down = np.argsort(-scores, kind="mergesort")
            up = np.argsort(scores, kind="mergesort")
            cache[(method, isoform)] = {
                "positions": positions, "labels_down": labels[down], "down": down,
                "down_starts": np.r_[0, np.nonzero(np.diff(scores[down]))[0] + 1],
                "labels_up": labels[up], "up": up,
                "up_starts": np.r_[0, np.nonzero(np.diff(scores[up]))[0] + 1],
            }

    def cached_metrics(method, row_weights):
        ap_values, auc_values = [], []
        for isoform in isoforms:
            item = cache[(method, isoform)]
            local = row_weights[item["positions"]]
            down_weight = local[item["down"]]
            down_labels = item["labels_down"]
            total_positive = float(np.dot(down_weight, down_labels))
            if total_positive <= 0:
                ap_values.append(math.nan)
            else:
                group_weight = np.add.reduceat(down_weight, item["down_starts"])
                group_positive = np.add.reduceat(down_weight * down_labels, item["down_starts"])
                keep = group_weight > 0
                group_weight, group_positive = group_weight[keep], group_positive[keep]
                precision = np.cumsum(group_positive) / np.cumsum(group_weight)
                ap_values.append(float(np.dot(precision, group_positive) / total_positive))
            up_weight = local[item["up"]]
            up_labels = item["labels_up"]
            positive = float(np.dot(up_weight, up_labels))
            negative = float(np.dot(up_weight, 1 - up_labels))
            if positive <= 0 or negative <= 0:
                auc_values.append(math.nan)
            else:
                group_positive = np.add.reduceat(up_weight * up_labels, item["up_starts"])
                group_negative = np.add.reduceat(up_weight * (1 - up_labels), item["up_starts"])
                below = np.r_[0.0, np.cumsum(group_negative)[:-1]]
                auc_values.append(float(np.dot(group_positive, below + 0.5 * group_negative) / (positive * negative)))
        return float(np.nanmean(ap_values)), float(np.nanmean(auc_values))

    differences_ap, differences_auc = [], []
    for _ in range(BOOTSTRAPS):
        sampled = generator.integers(len(scaffold_names), size=len(scaffold_names))
        group_weights = np.bincount(sampled, minlength=len(scaffold_names))
        row_weights = group_weights[row_group].astype(float)
        left_ap, left_auc = cached_metrics("isoform_specific_knn", row_weights)
        right_ap, right_auc = cached_metrics("pooled_chemical_knn", row_weights)
        differences_ap.append(left_ap - right_ap)
        differences_auc.append(left_auc - right_auc)
    left = method_metrics(rows, "isoform_specific_knn", isoforms)
    right = method_metrics(rows, "pooled_chemical_knn", isoforms)
    return {
        "stratum": stratum, "bootstrap_replicates": BOOTSTRAPS,
        "bootstrap_unit": "external_bemis_murcko_scaffold", "scaffold_groups": len(scaffold_names),
        "average_precision": {
            "isoform_specific": left["macro_average_precision"], "pooled": right["macro_average_precision"],
            "difference": left["macro_average_precision"] - right["macro_average_precision"],
            "scaffold_bootstrap_95ci": [float(x) for x in np.quantile(differences_ap, (0.025, 0.975))],
        },
        "roc_auc": {
            "isoform_specific": left["macro_roc_auc"], "pooled": right["macro_roc_auc"],
            "difference": left["macro_roc_auc"] - right["macro_roc_auc"],
            "scaffold_bootstrap_95ci": [float(x) for x in np.quantile(differences_auc, (0.025, 0.975))],
        },
    }


def reconstruct_permutation(rows, isoforms, stratum, comparison):
    positions_by_compound = defaultdict(list)
    for number, row in enumerate(rows):
        positions_by_compound[row["compound_inchikey"]].append(number)
    original = np.asarray([row["isoform_specific_knn"] for row in rows], dtype=float)
    null_differences = []
    pooled_ap = comparison["average_precision"]["pooled"]
    observed = comparison["average_precision"]["difference"]
    for replicate in range(PERMUTATIONS):
        shuffled = original.copy()
        for compound, positions in positions_by_compound.items():
            token = int.from_bytes(hashlib.sha256(f"{SEED}|{stratum}|{compound}|{replicate}".encode()).digest()[:8], "little")
            permutation = np.random.default_rng(token).permutation(len(positions))
            shuffled[np.asarray(positions)] = original[np.asarray(positions)[permutation]]
        temporary = [{**row, "shuffled": float(shuffled[number])} for number, row in enumerate(rows)]
        null_differences.append(method_metrics(temporary, "shuffled", isoforms)["macro_average_precision"] - pooled_ap)
    pvalue = (1 + sum(value >= observed for value in null_differences)) / (PERMUTATIONS + 1)
    return {
        "stratum": stratum, "permutations": PERMUTATIONS,
        "unit": "isoform_specific_scores_within_compound_available_isoforms", "refitted": False,
        "observed_macro_average_precision_difference_vs_pooled": observed,
        "permuted_differences": null_differences,
        "plus_one_upper_tail_monte_carlo_p": pvalue,
    }


def main():
    if VALIDATION.exists():
        raise FileExistsError(f"refusing to replace {VALIDATION}")
    RDLogger.DisableLog("rdApp.*")
    manifest = load(RESULTS / "input_manifest.json")
    for path, expected in manifest.items():
        if digest_file(path) != expected:
            raise AssertionError(f"input hash mismatch: {path}")
    train_rows = load(TRAIN_PATH)
    all_external = load(TEST_PATH)
    external_rows = [row for row in all_external if row["strata"]["all_parseable_unique"]]
    saved = np.load(RESULTS / "fingerprint_similarity.npz")
    rebuilt, isoforms, max_difference, train_on_bits, external_on_bits = reconstruct_predictions(train_rows, external_rows, saved)
    stored_predictions = [json.loads(line) for line in (RESULTS / "external_predictions.jsonl").read_text(encoding="utf-8").splitlines()]
    assert_close(rebuilt, stored_predictions, "predictions")

    expected_summary = load(RESULTS / "method_summary.json")["strata"]
    expected_comparisons = load(RESULTS / "comparison_summary.json")["strata"]
    expected_controls = load(RESULTS / "score_permutation_control.json")["strata"]
    rebuilt_summaries, rebuilt_comparisons, rebuilt_controls = {}, {}, {}
    for stratum in STRATA:
        rows = [row for row in rebuilt if row["strata"][stratum]]
        rebuilt_summaries[stratum] = {method: method_metrics(rows, method, isoforms) for method in METHODS}
        rebuilt_comparisons[stratum] = reconstruct_comparison(rows, isoforms, stratum)
        rebuilt_controls[stratum] = reconstruct_permutation(rows, isoforms, stratum, rebuilt_comparisons[stratum])
        assert_close(rebuilt_summaries[stratum], expected_summary[stratum], f"summary.{stratum}")
        assert_close(rebuilt_comparisons[stratum], expected_comparisons[stratum], f"comparison.{stratum}")
        assert_close(rebuilt_controls[stratum], expected_controls[stratum], f"permutation.{stratum}")
        print(f"verified external stratum: {stratum}", flush=True)

    audit = load(AUDIT_PATH)
    acceptance = load(RESULTS / "acceptance.json")
    for stratum, claim in (("lineage_eligible", "independent_external_transfer"),
                           ("scaffold_and_lineage_eligible", "scaffold_external_transfer")):
        information = audit["strata"][stratum]
        comparison = rebuilt_comparisons[stratum]
        control = rebuilt_controls[stratum]
        base_gates = {
            "at_least_four_isoforms_with_20_each_class": sum(
                row["positives"] >= 20 and row["negatives"] >= 20 for row in information["by_isoform"]
            ) >= 4,
            "at_least_100_unique_compounds": information["unique_compounds"] >= 100,
            "at_least_20_scaffold_groups": information["scaffold_groups"] >= 20,
            "positive_ap_interval_vs_pooled": comparison["average_precision"]["scaffold_bootstrap_95ci"][0] > 0,
            "score_permutation_p_at_most_0_05": control["plus_one_upper_tail_monte_carlo_p"] <= 0.05,
        }
        for key, value in base_gates.items():
            if acceptance[claim]["gates"].get(key) != value:
                raise AssertionError(f"acceptance gate mismatch: {claim}.{key}")
    independent = acceptance["independent_external_transfer"]["accepted"]
    if independent != all(acceptance["independent_external_transfer"]["gates"].values()):
        raise AssertionError("independent acceptance aggregation mismatch")
    scaffold_base = all(value for key, value in acceptance["scaffold_external_transfer"]["gates"].items()
                        if key != "independent_external_transfer_accepted")
    if acceptance["scaffold_external_transfer"]["accepted"] != (independent and scaffold_base):
        raise AssertionError("scaffold acceptance aggregation mismatch")

    result_files = [
        RESULTS / "external_predictions.jsonl", RESULTS / "fingerprint_similarity.npz",
        RESULTS / "method_summary.json", RESULTS / "comparison_summary.json",
        RESULTS / "score_permutation_control.json", RESULTS / "acceptance.json",
        RESULTS / "audit.json", RESULTS / "input_manifest.json",
    ]
    report = {
        "status": "PASS", "created_utc": now(),
        "input_hashes_verified": len(manifest), "predictions_verified": len(rebuilt),
        "train_fingerprint_bits_verified": int(saved["train_bits"].size),
        "external_fingerprint_bits_verified": int(saved["external_bits"].size),
        "train_on_bits": train_on_bits, "external_on_bits": external_on_bits,
        "similarity_cells_verified": int(saved["similarity"].size),
        "maximum_similarity_difference": max_difference,
        "strata_verified": len(STRATA), "bootstrap_replicates_per_stratum": BOOTSTRAPS,
        "permutations_per_stratum": PERMUTATIONS,
        "rdkit_version": rdBase.rdkitVersion,
        "acceptance": {
            "independent_external_transfer": acceptance["independent_external_transfer"]["accepted"],
            "scaffold_external_transfer": acceptance["scaffold_external_transfer"]["accepted"],
        },
        "result_hashes": {path.name: digest_file(path) for path in result_files},
        "validator_sha256": digest_file(Path(__file__)),
    }
    VALIDATION.mkdir()
    write_json(VALIDATION / "validation.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
