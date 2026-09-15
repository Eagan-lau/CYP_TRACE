"""Descriptive fixed-human coverage-risk analysis without external-label tuning."""
from __future__ import annotations

from collections import defaultdict
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
OUTPUT = Path(os.environ.get("CYPTRACE_SELECTIVE_OUTPUT", ROOT / "human_selective_behavior_01")).resolve()
METHODS = ("pooled_chemical_knn", "isoform_specific_knn")
STRATA = ("lineage_eligible", "scaffold_and_lineage_eligible")
TARGET_COVERAGES = (1.0, 0.75, 0.5, 0.25, 0.1)
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


def stable(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def ranked_indices(rows: list[dict], method: str) -> list[int]:
    return sorted(range(len(rows)), key=lambda index: (-rows[index][method], stable(rows[index]["compound_inchikey"])))


def summarize(rows: list[dict], accepted: set[int], isoforms: list[str]) -> dict:
    per_isoform = []
    for isoform in isoforms:
        indices = [number for number, row in enumerate(rows) if row["isoform"] == isoform]
        selected = [number for number in indices if number in accepted]
        positives = sum(rows[number]["label"] for number in indices)
        selected_positives = sum(rows[number]["label"] for number in selected)
        precision = selected_positives / len(selected) if selected else None
        recall = selected_positives / positives if positives else None
        per_isoform.append({
            "isoform": isoform, "rows": len(indices), "positives": positives,
            "selected_rows": len(selected), "selected_positives": selected_positives,
            "coverage": len(selected) / len(indices) if indices else None,
            "precision": precision, "selective_risk": 1 - precision if precision is not None else None,
            "positive_recall": recall,
        })
    return {
        "rows": len(rows), "selected_rows": len(accepted),
        "micro_coverage": len(accepted) / len(rows) if rows else None,
        "macro_coverage": float(np.mean([row["coverage"] for row in per_isoform])),
        "macro_precision": float(np.mean([row["precision"] for row in per_isoform if row["precision"] is not None])),
        "macro_selective_risk": float(np.mean([row["selective_risk"] for row in per_isoform if row["selective_risk"] is not None])),
        "macro_positive_recall": float(np.mean([row["positive_recall"] for row in per_isoform if row["positive_recall"] is not None])),
        "per_isoform": per_isoform,
    }


def frozen_thresholds(development: list[dict], isoforms: list[str]) -> list[dict]:
    output = []
    for method in METHODS:
        for isoform in isoforms:
            rows = [row for row in development if row["isoform"] == isoform]
            order = ranked_indices(rows, method)
            for target in TARGET_COVERAGES:
                count = min(len(rows), max(1, math.ceil(target * len(rows))))
                threshold = float(rows[order[count - 1]][method])
                accepted = {number for number, row in enumerate(rows) if row[method] >= threshold}
                summary = summarize(rows, accepted, [isoform])["per_isoform"][0]
                output.append({
                    "method": method, "isoform": isoform, "target_development_coverage": target,
                    "threshold": threshold, "development_achieved_coverage": summary["coverage"],
                    "development_precision": summary["precision"],
                    "development_positive_recall": summary["positive_recall"],
                })
    return output


def apply_frozen(rows: list[dict], thresholds: list[dict], isoforms: list[str]) -> list[dict]:
    results = []
    threshold_map = {(row["method"], row["isoform"], row["target_development_coverage"]): row for row in thresholds}
    for method in METHODS:
        for target in TARGET_COVERAGES:
            accepted = {
                number for number, row in enumerate(rows)
                if row[method] >= threshold_map[(method, row["isoform"], target)]["threshold"]
            }
            results.append({"method": method, "target_coverage": target, **summarize(rows, accepted, isoforms)})
    return results


def matched_masks(rows: list[dict], isoforms: list[str]) -> dict[tuple[str, float], set[int]]:
    output = {}
    grouped = {isoform: [number for number, row in enumerate(rows) if row["isoform"] == isoform] for isoform in isoforms}
    for method in METHODS:
        for target in TARGET_COVERAGES:
            accepted = set()
            for isoform in isoforms:
                indices = grouped[isoform]
                order = sorted(indices, key=lambda index: (-rows[index][method], stable(rows[index]["compound_inchikey"])))
                count = min(len(order), max(1, math.ceil(target * len(order))))
                accepted.update(order[:count])
            output[(method, target)] = accepted
    return output


def bootstrap_difference(rows: list[dict], masks: dict, isoforms: list[str], stratum_number: int) -> dict:
    scaffolds = sorted({row["scaffold_group"] for row in rows})
    scaffold_index = {value: number for number, value in enumerate(scaffolds)}
    rng = np.random.default_rng(SEED + stratum_number)
    weights = rng.multinomial(len(scaffolds), np.full(len(scaffolds), 1 / len(scaffolds)), size=BOOTSTRAPS)
    output = {}
    for target in TARGET_COVERAGES:
        method_precision = {}
        for method in METHODS:
            accepted = masks[(method, target)]
            values = []
            for isoform in isoforms:
                selected_count = np.zeros(len(scaffolds), dtype=np.int64)
                selected_positive = np.zeros(len(scaffolds), dtype=np.int64)
                for number, row in enumerate(rows):
                    if row["isoform"] == isoform and number in accepted:
                        group = scaffold_index[row["scaffold_group"]]
                        selected_count[group] += 1
                        selected_positive[group] += row["label"]
                denominator = weights @ selected_count
                numerator = weights @ selected_positive
                values.append(np.divide(numerator, denominator, out=np.full(BOOTSTRAPS, np.nan), where=denominator > 0))
            method_precision[method] = np.nanmean(np.vstack(values), axis=0)
        difference = method_precision["isoform_specific_knn"] - method_precision["pooled_chemical_knn"]
        observed = (
            summarize(rows, masks[("isoform_specific_knn", target)], isoforms)["macro_precision"]
            - summarize(rows, masks[("pooled_chemical_knn", target)], isoforms)["macro_precision"]
        )
        output[str(target)] = {
            "isoform_specific_minus_pooled_macro_precision": observed,
            "scaffold_bootstrap_95ci": [float(np.nanpercentile(difference, 2.5)), float(np.nanpercentile(difference, 97.5))],
            "bootstrap_replicates": BOOTSTRAPS,
        }
    return output


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to replace {OUTPUT}")
    development = read_jsonl(DEVELOPMENT)
    external_all = read_jsonl(EXTERNAL)
    isoforms = sorted({row["isoform"] for row in external_all})
    thresholds = frozen_thresholds(development, isoforms)
    strata = {}
    table_rows = []
    for stratum_number, stratum in enumerate(STRATA):
        rows = [row for row in external_all if row["strata"][stratum]]
        masks = matched_masks(rows, isoforms)
        frozen = apply_frozen(rows, thresholds, isoforms)
        matched = []
        for method in METHODS:
            for target in TARGET_COVERAGES:
                summary = summarize(rows, masks[(method, target)], isoforms)
                record = {"method": method, "target_coverage": target, **summary}
                matched.append(record)
                table_rows.append({
                    "stratum": stratum, "curve": "label_blind_matched", "method": method,
                    "target_coverage": target, "achieved_coverage": summary["macro_coverage"],
                    "macro_precision": summary["macro_precision"], "macro_selective_risk": summary["macro_selective_risk"],
                    "macro_positive_recall": summary["macro_positive_recall"], "selected_rows": summary["selected_rows"],
                })
        for record in frozen:
            table_rows.append({
                "stratum": stratum, "curve": "development_frozen_threshold", "method": record["method"],
                "target_coverage": record["target_coverage"], "achieved_coverage": record["macro_coverage"],
                "macro_precision": record["macro_precision"], "macro_selective_risk": record["macro_selective_risk"],
                "macro_positive_recall": record["macro_positive_recall"], "selected_rows": record["selected_rows"],
            })
        strata[stratum] = {
            "rows": len(rows), "compounds": len({row["compound_inchikey"] for row in rows}),
            "scaffolds": len({row["scaffold_group"] for row in rows}),
            "development_frozen_threshold_curve": frozen,
            "label_blind_matched_coverage_curve": matched,
            "matched_coverage_comparison": bootstrap_difference(rows, masks, isoforms, stratum_number),
        }
    OUTPUT.mkdir()
    analysis = {
        "status": "COMPUTED_QC_PENDING",
        "scope": "descriptive fixed-human external selective behavior; no external retuning and no confirmatory gate",
        "methods": list(METHODS), "isoforms": isoforms, "target_coverages": list(TARGET_COVERAGES),
        "thresholds": thresholds, "strata": strata,
        "input_sha256": {"development": digest(DEVELOPMENT), "external": digest(EXTERNAL), "protocol": digest(PROTOCOL)},
        "seed": SEED, "bootstrap_replicates": BOOTSTRAPS,
    }
    (OUTPUT / "analysis.json").write_bytes((json.dumps(analysis, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8"))
    with (OUTPUT / "coverage_risk_curves.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=list(table_rows[0]))
        writer.writeheader()
        writer.writerows(table_rows)
    (OUTPUT / "input_manifest.json").write_bytes((json.dumps(analysis["input_sha256"], indent=2) + "\n").encode("utf-8"))
    print(json.dumps({
        "status": analysis["status"],
        "strata": {
            name: {
                "rows": value["rows"],
                "matched_precision_difference": value["matched_coverage_comparison"],
            }
            for name, value in strata.items()
        },
    }, indent=2))


if __name__ == "__main__":
    main()
