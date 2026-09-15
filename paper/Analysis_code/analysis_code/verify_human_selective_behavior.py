"""Independent reconstruction of fixed-human selective-behavior outputs."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
DEVELOPMENT = ROOT / "human_substrate_models_01" / "outer_predictions.jsonl"
EXTERNAL = ROOT / "external_human_cyp_models_01" / "external_predictions.jsonl"
PROTOCOL = ROOT / "HUMAN_SELECTIVE_BEHAVIOR_V1.md"
ANALYSIS_DIR = Path(os.environ.get("CYPTRACE_SELECTIVE_ANALYSIS_DIR", ROOT / "human_selective_behavior_01")).resolve()
ANALYSIS = ANALYSIS_DIR / "analysis.json"
TABLE = ANALYSIS_DIR / "coverage_risk_curves.tsv"
OUTPUT = Path(os.environ.get("CYPTRACE_SELECTIVE_VALIDATION_OUTPUT", ROOT / "human_selective_behavior_validation_01")).resolve()
METHODS = ("pooled_chemical_knn", "isoform_specific_knn")
COVERAGES = (1.0, 0.75, 0.5, 0.25, 0.1)
STRATA = ("lineage_eligible", "scaffold_and_lineage_eligible")
BOOTSTRAPS = 5000
SEED = 20260918


def read_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def tie(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def metrics(rows: list[dict], accepted: set[int], isoforms: list[str]) -> dict:
    values = []
    for isoform in isoforms:
        all_rows = [index for index, row in enumerate(rows) if row["isoform"] == isoform]
        selected = [index for index in all_rows if index in accepted]
        positive = sum(rows[index]["label"] for index in all_rows)
        selected_positive = sum(rows[index]["label"] for index in selected)
        values.append((len(selected) / len(all_rows), selected_positive / len(selected), selected_positive / positive))
    return {
        "selected_rows": len(accepted),
        "macro_coverage": float(np.mean([value[0] for value in values])),
        "macro_precision": float(np.mean([value[1] for value in values])),
        "macro_selective_risk": float(np.mean([1 - value[1] for value in values])),
        "macro_positive_recall": float(np.mean([value[2] for value in values])),
    }


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to replace {OUTPUT}")
    result = json.loads(ANALYSIS.read_text(encoding="utf-8"))
    development = read_jsonl(DEVELOPMENT)
    external = read_jsonl(EXTERNAL)
    isoforms = sorted({row["isoform"] for row in external})
    checks = {
        "analysis_status_pending_qc": result["status"] == "COMPUTED_QC_PENDING",
        "protocol_hash": result["input_sha256"]["protocol"] == digest(PROTOCOL),
        "development_hash": result["input_sha256"]["development"] == digest(DEVELOPMENT),
        "external_hash": result["input_sha256"]["external"] == digest(EXTERNAL),
        "external_six_isoforms": result["isoforms"] == isoforms and len(isoforms) == 6,
    }

    threshold_map = {(row["method"], row["isoform"], row["target_development_coverage"]): row for row in result["thresholds"]}
    threshold_differences = []
    for method in METHODS:
        for isoform in isoforms:
            rows = [row for row in development if row["isoform"] == isoform]
            ordered = sorted(rows, key=lambda row: (-row[method], tie(row["compound_inchikey"])))
            for coverage in COVERAGES:
                count = min(len(rows), max(1, math.ceil(coverage * len(rows))))
                threshold_differences.append(abs(float(ordered[count - 1][method]) - threshold_map[(method, isoform, coverage)]["threshold"]))
    checks["all_thresholds_reconstructed_from_development"] = max(threshold_differences) == 0

    max_metric_difference = 0.0
    max_bootstrap_difference = 0.0
    label_blind_masks_unchanged = True
    for stratum_number, stratum in enumerate(STRATA):
        rows = [row for row in external if row["strata"][stratum]]
        stored = result["strata"][stratum]
        matched_lookup = {(row["method"], row["target_coverage"]): row for row in stored["label_blind_matched_coverage_curve"]}
        frozen_lookup = {(row["method"], row["target_coverage"]): row for row in stored["development_frozen_threshold_curve"]}
        masks = {}
        for method in METHODS:
            for coverage in COVERAGES:
                accepted = set()
                flipped_accepted = set()
                for isoform in isoforms:
                    indices = [index for index, row in enumerate(rows) if row["isoform"] == isoform]
                    order = sorted(indices, key=lambda index: (-rows[index][method], tie(rows[index]["compound_inchikey"])))
                    count = min(len(order), max(1, math.ceil(coverage * len(order))))
                    accepted.update(order[:count])
                    # Selection is rebuilt after changing labels; ranks must remain identical.
                    flipped = [{**rows[index], "label": 1 - rows[index]["label"]} for index in indices]
                    local_order = sorted(range(len(flipped)), key=lambda index: (-flipped[index][method], tie(flipped[index]["compound_inchikey"])))
                    flipped_accepted.update(indices[index] for index in local_order[:count])
                label_blind_masks_unchanged &= accepted == flipped_accepted
                masks[(method, coverage)] = accepted
                actual = metrics(rows, accepted, isoforms)
                expected = matched_lookup[(method, coverage)]
                for field in ("macro_coverage", "macro_precision", "macro_selective_risk", "macro_positive_recall"):
                    max_metric_difference = max(max_metric_difference, abs(actual[field] - expected[field]))
                max_metric_difference = max(max_metric_difference, abs(actual["selected_rows"] - expected["selected_rows"]))

                frozen = {
                    index for index, row in enumerate(rows)
                    if row[method] >= threshold_map[(method, row["isoform"], coverage)]["threshold"]
                }
                actual_frozen = metrics(rows, frozen, isoforms)
                expected_frozen = frozen_lookup[(method, coverage)]
                for field in ("macro_coverage", "macro_precision", "macro_selective_risk", "macro_positive_recall"):
                    max_metric_difference = max(max_metric_difference, abs(actual_frozen[field] - expected_frozen[field]))
                max_metric_difference = max(max_metric_difference, abs(actual_frozen["selected_rows"] - expected_frozen["selected_rows"]))

        scaffolds = sorted({row["scaffold_group"] for row in rows})
        scaffold_index = {value: index for index, value in enumerate(scaffolds)}
        rng = np.random.default_rng(SEED + stratum_number)
        weights = rng.multinomial(len(scaffolds), np.full(len(scaffolds), 1 / len(scaffolds)), size=BOOTSTRAPS)
        for coverage in COVERAGES:
            bootstrap_precision = {}
            for method in METHODS:
                isoform_values = []
                for isoform in isoforms:
                    selected_count = np.zeros(len(scaffolds), dtype=np.int64)
                    selected_positive = np.zeros(len(scaffolds), dtype=np.int64)
                    for index, row in enumerate(rows):
                        if row["isoform"] == isoform and index in masks[(method, coverage)]:
                            position = scaffold_index[row["scaffold_group"]]
                            selected_count[position] += 1
                            selected_positive[position] += row["label"]
                    denominator = weights @ selected_count
                    numerator = weights @ selected_positive
                    isoform_values.append(np.divide(numerator, denominator, out=np.full(BOOTSTRAPS, np.nan), where=denominator > 0))
                bootstrap_precision[method] = np.nanmean(np.vstack(isoform_values), axis=0)
            difference = bootstrap_precision["isoform_specific_knn"] - bootstrap_precision["pooled_chemical_knn"]
            interval = [float(np.nanpercentile(difference, 2.5)), float(np.nanpercentile(difference, 97.5))]
            expected = stored["matched_coverage_comparison"][str(coverage)]["scaffold_bootstrap_95ci"]
            max_bootstrap_difference = max(max_bootstrap_difference, *(abs(interval[i] - expected[i]) for i in (0, 1)))

    checks["label_blind_selection_invariant_to_labels"] = label_blind_masks_unchanged
    checks["all_curve_metrics_reconstructed"] = max_metric_difference <= 1e-15
    checks["all_bootstrap_intervals_reconstructed"] = max_bootstrap_difference <= 1e-15
    with TABLE.open(encoding="utf-8", newline="") as stream:
        table_rows = list(csv.DictReader(stream, delimiter="\t"))
    checks["curve_table_has_all_cells"] = len(table_rows) == len(STRATA) * 2 * len(METHODS) * len(COVERAGES)

    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "thresholds_checked": len(threshold_differences),
        "maximum_threshold_difference": max(threshold_differences),
        "maximum_curve_metric_difference": max_metric_difference,
        "maximum_bootstrap_interval_difference": max_bootstrap_difference,
        "analysis_sha256": digest(ANALYSIS),
        "table_sha256": digest(TABLE),
        "verification_script_sha256": digest(Path(__file__)),
        "interpretation": "Selection is demonstrably label-blind or development-frozen. External precision-risk results are descriptive, not a prospectively validated threshold.",
    }
    OUTPUT.mkdir()
    (OUTPUT / "validation.json").write_bytes((json.dumps(report, indent=2, ensure_ascii=False) + "\n").encode("utf-8"))
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] != "PASS":
        raise AssertionError("selective-behavior verification failed")


if __name__ == "__main__":
    main()
