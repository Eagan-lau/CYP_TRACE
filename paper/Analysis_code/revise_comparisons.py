#!/usr/bin/env python3
"""Exploratory context from frozen data; never changes the locked model or split.

Reproduce the six already selected logistic fits, retain individual predictions,
and compare them with frozen kNN scores using paired, fixed-fit scaffold resampling.
The model-class comparison and its interval were added after external inspection.
"""
import argparse
import csv
import hashlib
import itertools
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parent
SEED = 20260914


def read_jsonl(path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def random_ap(n, m):
    """Uniform random strict ordering, conditional on n and m > 0."""
    assert 0 < m <= n
    if n == 1:
        return 1.0
    harmonic = np.sum(1.0 / np.arange(1, n + 1))
    return float(harmonic / n + (m - 1) * (n - harmonic) / (n * (n - 1)))


def ap_template(y, scores):
    order = np.argsort(-scores, kind="mergesort")
    s = scores[order]
    starts = np.r_[0, np.flatnonzero(s[1:] != s[:-1]) + 1]
    return y[order], order, starts


def weighted_ap(template, weights):
    y, order, starts = template
    w = weights[order]
    positive = np.add.reduceat(w * y, starts)
    mass = np.add.reduceat(w, starts)
    if positive.sum() == 0:
        return np.nan
    cumulative_mass = np.cumsum(mass)
    precision = np.divide(np.cumsum(positive), cumulative_mass,
                          out=np.zeros_like(cumulative_mass), where=cumulative_mass > 0)
    return float(np.dot(precision, positive) / positive.sum())


def tests():
    # Exhaustive label-position enumeration checks the finite-sample expression.
    for n in range(2, 8):
        for m in range(1, n + 1):
            aps = []
            for positive_positions in itertools.combinations(range(n), m):
                y = np.zeros(n)
                y[list(positive_positions)] = 1
                aps.append(average_precision_score(y, -np.arange(n)))
            assert abs(np.mean(aps) - random_ap(n, m)) < 1e-12
    y = np.array([0, 1, 1, 0, 1])
    scores = np.array([.2, .2, .7, .1, .1])
    for weights in (np.ones(5), np.array([0., 2., 3., 1., 0.])):
        assert abs(weighted_ap(ap_template(y, scores), weights) -
                   average_precision_score(y, scores, sample_weight=weights)) < 1e-12


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=ROOT.parents[1] / "paper_rebuild_20260909")
    parser.add_argument("--bootstrap", type=int, default=5000)
    args = parser.parse_args()
    source = args.source.resolve()
    out = ROOT / "tables"
    out.mkdir(exist_ok=True)
    tests()
    source_files = [
        "human_substrate_01/normalized_labels.json",
        "human_substrate_models_01/fingerprint_similarity.npz",
        "external_human_cyp_models_01/fingerprint_similarity.npz",
        "external_human_cyp_models_01/external_predictions.jsonl",
        "reverse_candidate_census_01/query_panels.jsonl",
        "reverse_baselines_01/input_manifest.json",
        "REVERSE_RETRIEVAL_V1.md", "REVERSE_RETRIEVAL_ACCEPTANCE_V1.md",
    ]
    rows = json.loads((source / source_files[0]).read_text(encoding="utf-8"))
    with np.load(source / source_files[1], allow_pickle=False) as saved:
        train_ids = list(saved["compound_ids"])
        bits = saved["fingerprint_bits"].astype(np.float64)
    with np.load(source / source_files[2], allow_pickle=False) as saved:
        assert list(saved["train_ids"]) == train_ids
        external_ids = list(saved["external_ids"])
        external_bits = saved["external_bits"].astype(np.float64)
    train_index = {key: i for i, key in enumerate(train_ids)}
    external_index = {key: i for i, key in enumerate(external_ids)}
    labels = {iso: np.full(len(train_ids), np.nan) for iso in sorted({r["isoform"] for r in rows})}
    for row in rows:
        labels[row["isoform"]][train_index[row["compound_inchikey"]]] = row["label"]
    strict = [r for r in read_jsonl(source / source_files[3]) if r["strata"]["scaffold_and_lineage_eligible"]]
    assert len(strict) == 3035
    with (ROOT / "inputs/logistic_regression_supplement.tsv").open(encoding="utf-8") as handle:
        settings = {r["isoform"]: r for r in csv.DictReader(handle, delimiter="\t") if r["isoform"] != "MACRO"}
    comparison = []
    predictions = []
    templates = {}
    scaffolds = sorted({r["scaffold_group"] for r in strict})
    scaffold_index = {key: i for i, key in enumerate(scaffolds)}
    assert len(scaffolds) == 844
    for iso, parameters in settings.items():
        known = ~np.isnan(labels[iso])
        weight = parameters["selected_class_weight"]
        model = LogisticRegression(C=float(parameters["selected_C"]), solver="liblinear",
                                   class_weight=None if weight == "none" else weight,
                                   max_iter=2000, random_state=SEED)
        model.fit(bits[known], labels[iso][known].astype(np.int8))
        local = [r for r in strict if r["isoform"] == iso]
        ix = [external_index[r["compound_inchikey"]] for r in local]
        scores = model.predict_proba(external_bits[ix])[:, 1]
        y = np.array([r["label"] for r in local], dtype=np.int8)
        knn = np.array([r["isoform_specific_knn"] for r in local])
        pooled = np.array([r["pooled_chemical_knn"] for r in local])
        logistic_ap = average_precision_score(y, scores)
        assert abs(logistic_ap - float(parameters["strict_external_AP"])) < 1e-12, iso
        assert abs(average_precision_score(y, knn) - float(parameters["production_isoform_specific_kNN_AP"])) < 1e-12
        n, m = len(y), int(y.sum())
        comparison.append(dict(isoform=iso, labels=n, positives=m, prevalence=m/n,
                               exact_random_order_AP=random_ap(n, m), pooled_AP=average_precision_score(y, pooled),
                               specific_AP=average_precision_score(y, knn), logistic_AP=logistic_ap))
        templates[iso] = (ap_template(y, knn), ap_template(y, scores),
                          np.array([scaffold_index[r["scaffold_group"]] for r in local]))
        for row, score in zip(local, scores):
            predictions.append({**{k: row[k] for k in ("isoform", "compound_inchikey", "scaffold_group", "label",
                                                     "isoform_specific_knn", "pooled_chemical_knn")},
                                "isoform_logistic": float(score)})
    macro = {key: float(np.mean([r[key] for r in comparison])) for key in
             ("prevalence", "exact_random_order_AP", "pooled_AP", "specific_AP", "logistic_AP")}
    comparison.append(dict(isoform="MACRO", labels=len(strict), positives=sum(r["positives"] for r in comparison), **macro))
    write_tsv(out / "external_AP_context.tsv", comparison)
    write_tsv(out / "strict_logistic_paired_predictions.tsv", predictions)
    rng = np.random.default_rng(SEED)
    deltas = []
    for _ in range(args.bootstrap):
        counts = np.bincount(rng.integers(0, len(scaffolds), size=len(scaffolds)), minlength=len(scaffolds)).astype(float)
        pair = [(weighted_ap(a, counts[ix]), weighted_ap(b, counts[ix])) for a, b, ix in templates.values()]
        if np.isfinite(pair).all():
            deltas.append(float(np.mean([a-b for a, b in pair])))
    interval = np.quantile(deltas, [.025, .975]).tolist()
    panels = read_jsonl(source / source_files[4])
    assert len(panels) == 4880
    panel_stats = {}
    for field in ("candidate_count", "positive_count"):
        values = np.array([r[field] for r in panels])
        panel_stats[field] = dict(min=int(values.min()), q1=float(np.quantile(values, .25)),
                                 median=float(np.median(values)), q3=float(np.quantile(values, .75)), max=int(values.max()))
    frozen = json.loads((source / source_files[5]).read_text(encoding="utf-8"))
    protocols = {}
    for name in source_files[-2:]:
        payload = (source / name).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        expected = next(v for k, v in frozen.items() if k.endswith("/" + name))
        assert digest == expected, f"Protocol hash mismatch: {name}"
        protocols[name] = dict(sha256=digest, matches_frozen_reverse_input_manifest=True,
                               chronology="Sequence described in the contemporaneous protocol; hash verifies content, not independent timestamp.")
    result = dict(analysis_status="exploratory_after_external_inspection", seed=SEED,
                  analysis_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  versions={p: importlib.metadata.version(p) for p in ("numpy", "scikit-learn")},
                  strict_rows=len(strict), scaffold_clusters=len(scaffolds), macro=macro,
                  knn_minus_logistic=macro["specific_AP"]-macro["logistic_AP"],
                  paired_fixed_fit_95_interval=interval, bootstrap_requested=args.bootstrap,
                  bootstrap_valid=len(deltas), refit_within_bootstrap=False, retuned_on_external=False,
                  interval_scope="Fixed predictions; shared scaffold multiplicities preserve isoform pairing; no refitting uncertainty.",
                  reverse_panel_statistics=panel_stats, reverse_protocols=protocols,
                  verification="Six logistic AP values reproduce the prior table to 1e-12; finite-sample AP and tied weighted AP unit checks passed.",
                  input_sha256={name: hashlib.sha256((source/name).read_bytes()).hexdigest() for name in source_files},
                  parameter_table_sha256=hashlib.sha256((ROOT/"inputs/logistic_regression_supplement.tsv").read_bytes()).hexdigest())
    (out / "comparison_revision_summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
