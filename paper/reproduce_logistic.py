"""Reproduce exploratory logistic selection and external evaluation separately.

`tune` reads development labels only; `evaluate` consumes its frozen selection.
The algorithm reproduces logistic_baseline in reanalyse_frozen_outputs.py:
six C/weight combinations, five stored development folds, mean fold AP,
then smaller C and lexical str(class_weight) as deterministic tie breakers.
This is a reproduction of an exploratory comparison, not a new validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import warnings

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
SEED = 20260914
ISOFORMS = ("CYP1A2", "CYP2C19", "CYP2C9", "CYP2D6", "CYP2E1", "CYP3A4")
GRID = tuple((c, weight) for c in (0.01, 0.1, 1.0) for weight in (None, "balanced"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def table(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def close(a, b, tolerance=1e-12):
    if not math.isclose(float(a), float(b), abs_tol=tolerance, rel_tol=0):
        raise AssertionError(f"Numerical mismatch: {a} != {b}")


def average_precision(y, score):
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=np.float64)
    if not y.sum():
        return math.nan
    order = np.argsort(-score, kind="mergesort")
    s, yy = score[order], y[order]
    starts = np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1]
    mass = np.diff(np.r_[starts, len(y)])
    positive = np.add.reduceat(yy, starts)
    return float(np.dot(np.cumsum(positive) / np.cumsum(mass), positive) / y.sum())


def fingerprints(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    bits = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for i, text in enumerate(smiles):
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            raise ValueError(f"Invalid normalized SMILES at index {i}")
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(mol), bits[i])
    return bits


def development_arrays(path):
    rows = read(path)
    compounds = sorted({r["compound_inchikey"] for r in rows})
    isoforms = sorted({r["isoform"] for r in rows})
    ci = {key: i for i, key in enumerate(compounds)}
    pi = {key: i for i, key in enumerate(isoforms)}
    labels = np.full((len(compounds), len(isoforms)), np.nan)
    compounds_meta = {}
    scaffold_folds = {}
    for row in rows:
        key, iso = row["compound_inchikey"], row["isoform"]
        c, p = ci[key], pi[iso]
        if not np.isnan(labels[c, p]):
            raise ValueError("Duplicate development pair")
        if row["label"] not in (0, 1):
            raise ValueError("Non-binary label")
        labels[c, p] = row["label"]
        meta = (row["canonical_smiles"], row["scaffold_group"], row["development_fold"])
        if key in compounds_meta and compounds_meta[key] != meta:
            raise ValueError("Conflicting compound representation or fold")
        compounds_meta[key] = meta
        scaffold, fold = meta[1:]
        if fold not in range(5):
            raise ValueError("Development fold outside 0..4")
        if scaffold in scaffold_folds and scaffold_folds[scaffold] != fold:
            raise ValueError("Scaffold crosses development folds")
        scaffold_folds[scaffold] = fold
    bits = fingerprints([compounds_meta[key][0] for key in compounds])
    folds = np.array([compounds_meta[key][2] for key in compounds], dtype=np.int64)
    return compounds, isoforms, bits, folds, labels


def choose(candidates):
    return sorted(candidates, key=lambda r: (-r["mean_fold_AP"], r["C"], str(r["class_weight"])))[0]


def fit(bits, labels, c, weight):
    model = LogisticRegression(C=c, solver="liblinear", class_weight=weight,
                               max_iter=2000, random_state=SEED)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        model.fit(bits.astype(np.float64), labels.astype(np.int8))
    return model


def tune(development, output):
    compounds, isoforms, bits, folds, labels = development_arrays(development)
    candidates, selections = [], {}
    for iso in ISOFORMS:
        known = ~np.isnan(labels[:, isoforms.index(iso)])
        x = bits[known].astype(np.float64)
        y = labels[known, isoforms.index(iso)].astype(np.int8)
        local_folds = folds[known]
        trials = []
        for c, weight in GRID:
            fold_aps = []
            for fold in range(5):
                train, test = local_folds != fold, local_folds == fold
                if not test.any() or len(np.unique(y[train])) != 2 or not y[test].sum():
                    raise ValueError(f"Invalid CV support for {iso}, fold {fold}")
                model = fit(x[train], y[train], c, weight)
                fold_aps.append(average_precision(y[test], model.predict_proba(x[test])[:, 1]))
            trial = dict(isoform=iso, C=c, class_weight=weight,
                         fold_AP=fold_aps, mean_fold_AP=float(np.mean(fold_aps)))
            trials.append(trial)
            candidates.append(trial)
        selections[iso] = dict(choose(trials), development_rows=int(known.sum()),
                              development_positives=int(y.sum()))
        print(f"Selected {iso}: C={selections[iso]['C']}, weight={selections[iso]['class_weight']}", flush=True)
    result = dict(analysis_status="exploratory_after_external_inspection", seed=SEED,
                  selection_inputs=["development_labels.json"], external_data_read=False,
                  development_sha256=sha(development), compound_count=len(compounds),
                  development_label_count=int(np.isfinite(labels).sum()),
                  fingerprint_bits_sha256=hashlib.sha256(bits.tobytes()).hexdigest(),
                  folds_int64_le_sha256=hashlib.sha256(folds.astype("<i8").tobytes()).hexdigest(),
                  grid_fits=6*6*5, selected=selections, trials=candidates,
                  versions={p: importlib.metadata.version(p) for p in ("numpy", "rdkit", "scikit-learn")})
    write(output, result)
    return result


def ap_template(y, scores):
    order = np.argsort(-scores, kind="mergesort")
    s = scores[order]
    return y[order], order, np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1]


def weighted_ap(template, weights):
    y, order, starts = template
    w = weights[order]
    positive = np.add.reduceat(w*y, starts)
    mass = np.add.reduceat(w, starts)
    if positive.sum() == 0:
        return math.nan
    denominator = np.cumsum(mass)
    precision = np.divide(np.cumsum(positive), denominator,
                          out=np.zeros_like(denominator), where=denominator > 0)
    return float(np.dot(precision, positive) / positive.sum())


def evaluate(root, selection_path, output, bootstrap=5000):
    selected = read(selection_path)
    development = root / "data/development_labels.json"
    if selected["development_sha256"] != sha(development):
        raise ValueError("Selection and fitting use different development data")
    ids, isoforms, bits, _, labels = development_arrays(development)
    rows = [json.loads(line) for line in (root / "data/external_predictions.jsonl").read_text(encoding="utf-8").splitlines() if line]
    strict = [r for r in rows if r["strata"]["scaffold_and_lineage_eligible"]]
    external = read(root / "data/external_labels.json")
    smiles = {}
    for row in external:
        key, value = row["compound_inchikey"], row["canonical_smiles"]
        if key in smiles and smiles[key] != value:
            raise ValueError("Conflicting external structure")
        smiles[key] = value
    external_ids = sorted(smiles)
    external_bits = fingerprints([smiles[key] for key in external_ids])
    ei = {key: i for i, key in enumerate(external_ids)}
    expected = {r["isoform"]: r for r in table(root / "paper/Tables/Table_S14b_logistic_regression_supplement.tsv")}
    previous = read(root / "paper/Source_data/analysis_tables/comparison_revision_summary.json")
    provenance = read(root / "paper/Provenance/LOGISTIC_INPUT_PROVENANCE.json")
    close(len(strict), 3035)
    if hashlib.sha256(bits.tobytes()).hexdigest() != provenance["development_fingerprint_bits_sha256"]:
        raise AssertionError("Development fingerprints differ from the original frozen matrix")
    if hashlib.sha256(external_bits.tobytes()).hexdigest() != provenance["external_fingerprint_bits_sha256"]:
        raise AssertionError("External fingerprints differ from the original frozen matrix")
    if selected["folds_int64_le_sha256"] != provenance["development_folds_int64_le_sha256"]:
        raise AssertionError("Development folds differ from the original matrix")
    scaffolds = sorted({r["scaffold_group"] for r in strict})
    si = {key: i for i, key in enumerate(scaffolds)}
    close(len(scaffolds), 844)
    predictions, comparison, templates = [], [], {}
    for iso in ISOFORMS:
        choice = selected["selected"][iso]
        exp = expected[iso]
        close(choice["C"], exp["selected_C"])
        if (choice["class_weight"] or "none") != exp["selected_class_weight"]:
            raise AssertionError("Selected class weight changed")
        close(choice["mean_fold_AP"], exp["mean_five_fold_development_AP"])
        p = isoforms.index(iso)
        known = ~np.isnan(labels[:, p])
        model = fit(bits[known], labels[known, p], choice["C"], choice["class_weight"])
        local = [r for r in strict if r["isoform"] == iso]
        scores = model.predict_proba(external_bits[[ei[r["compound_inchikey"]] for r in local]].astype(np.float64))[:, 1]
        y = np.array([r["label"] for r in local], dtype=np.int8)
        knn = np.array([r["isoform_specific_knn"] for r in local])
        ap = average_precision(y, scores)
        close(ap, exp["strict_external_AP"])
        comparison.append(dict(isoform=iso, labels=len(local), logistic_AP=ap,
                               knn_AP=average_precision(y, knn)))
        templates[iso] = ap_template(y, knn), ap_template(y, scores), np.array([si[r["scaffold_group"]] for r in local])
        for r, score in zip(local, scores):
            predictions.append({**{k: r[k] for k in ("isoform", "compound_inchikey", "scaffold_group", "label")}, "logistic_score":float(score)})
    original_scores = {(r["isoform"], r["compound_inchikey"]): float(r["isoform_logistic"]) for r in table(root / "paper/Source_data/analysis_tables/strict_logistic_paired_predictions.tsv")}
    max_error = max(abs(r["logistic_score"]-original_scores[(r["isoform"], r["compound_inchikey"])]) for r in predictions)
    close(max_error, 0, 1e-10)
    macro = float(np.mean([r["logistic_AP"] for r in comparison]))
    delta = float(np.mean([r["knn_AP"]-r["logistic_AP"] for r in comparison]))
    close(macro, previous["macro"]["logistic_AP"])
    close(delta, previous["knn_minus_logistic"])
    rng = np.random.default_rng(SEED)
    deltas = []
    for _ in range(bootstrap):
        counts = np.bincount(rng.integers(0, len(scaffolds), size=len(scaffolds)), minlength=len(scaffolds)).astype(float)
        pair = [(weighted_ap(a, counts[ix]), weighted_ap(b, counts[ix])) for a, b, ix in templates.values()]
        if np.isfinite(pair).all():
            deltas.append(float(np.mean([a-b for a, b in pair])))
    interval = np.quantile(deltas, [.025, .975]).tolist()
    if bootstrap == 5000:
        for a, b in zip(interval, previous["paired_fixed_fit_95_interval"]):
            close(a, b)
    result = dict(status="PASS", analysis_status="exploratory_after_external_inspection",
                  selection_matches_all_six_isoforms=True, selection_grid_fits=selected["grid_fits"],
                  source_fingerprints_and_folds_match=True, strict_external_rows=len(strict),
                  external_score_max_abs_error=max_error, logistic_macro_AP=macro,
                  knn_minus_logistic=delta, paired_fixed_fit_95_interval=interval,
                  bootstrap_requested=bootstrap, bootstrap_valid=len(deltas),
                  refitted_within_bootstrap=False, per_isoform=comparison,
                  selection_sha256=sha(selection_path), raw_source_normalization_repeated=False)
    write(output / "logistic_evaluation.json", result)
    write_table(output / "logistic_predictions.tsv", predictions)
    print(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("tune", "evaluate"))
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, required=True, help="New directory; existing paths are never overwritten")
    parser.add_argument("--selection", type=Path)
    args = parser.parse_args()
    if args.mode == "evaluate" and args.selection is None:
        parser.error("evaluate requires --selection")
    args.output.mkdir(parents=True, exist_ok=False)
    if args.mode == "tune":
        tune(args.root / "data/development_labels.json", args.output / "logistic_selection.json")
    else:
        evaluate(args.root, args.selection, args.output)


if __name__ == "__main__":
    main()
