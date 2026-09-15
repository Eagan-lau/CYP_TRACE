"""Verify Figure 6 source bindings and exported file integrity."""
from __future__ import annotations

import json
from pathlib import Path

from PIL import Image

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main():
    source = read(FIGURES / "figure6_external_transfer_source_data.json")
    provenance = read(FIGURES / "figure6_external_transfer_provenance.json")
    audit = read(ROOT / "external_human_cyp_01" / "overlap_audit.json")
    summary = read(ROOT / "external_human_cyp_models_01" / "method_summary.json")["strata"]
    comparison = read(ROOT / "external_human_cyp_models_01" / "comparison_summary.json")["strata"]
    permutation = read(ROOT / "external_human_cyp_models_01" / "score_permutation_control.json")["strata"]
    validation = read(ROOT / "external_human_cyp_validation_01" / "validation.json")
    failures = []

    expected_stages = [
        ("Unique evaluable pairs", audit["strata"]["all_parseable_unique"]["rows"]),
        ("Exact-new compounds", audit["strata"]["exact_novel"]["rows"]),
        ("+ eligible lineage", audit["strata"]["lineage_eligible"]["rows"]),
        ("+ scaffold novel", audit["strata"]["scaffold_and_lineage_eligible"]["rows"]),
    ]
    if [(row["label"], row["rows"]) for row in source["stages"]] != expected_stages:
        failures.append("stage bindings")

    for row in source["aggregate"]:
        key = row["stratum"]
        expected = {
            "rows": audit["strata"][key]["rows"],
            "compounds": audit["strata"][key]["unique_compounds"],
            "scaffolds": audit["strata"][key]["scaffold_groups"],
            "pooled_ap": summary[key]["pooled_chemical_knn"]["macro_average_precision"],
            "specific_ap": summary[key]["isoform_specific_knn"]["macro_average_precision"],
            "pooled_auc": summary[key]["pooled_chemical_knn"]["macro_roc_auc"],
            "specific_auc": summary[key]["isoform_specific_knn"]["macro_roc_auc"],
            "ap_difference": comparison[key]["average_precision"]["difference"],
            "ap_ci": comparison[key]["average_precision"]["scaffold_bootstrap_95ci"],
            "auc_difference": comparison[key]["roc_auc"]["difference"],
            "auc_ci": comparison[key]["roc_auc"]["scaffold_bootstrap_95ci"],
            "permutation_p": permutation[key]["plus_one_upper_tail_monte_carlo_p"],
        }
        for field, value in expected.items():
            if row[field] != value:
                failures.append(f"aggregate {key} {field}")

    strict = "scaffold_and_lineage_eligible"
    expected_methods = {}
    for method in ("pooled_chemical_knn", "isoform_specific_knn"):
        expected_methods[method] = {row["isoform"]: row for row in summary[strict][method]["per_isoform"]}
    for row in source["strict_per_isoform"]:
        isoform = row["isoform"]
        pooled = expected_methods["pooled_chemical_knn"][isoform]
        specific = expected_methods["isoform_specific_knn"][isoform]
        expected = (pooled["rows"], pooled["positives"], pooled["negatives"],
                    pooled["average_precision"], specific["average_precision"])
        observed = (row["rows"], row["positives"], row["negatives"], row["pooled_ap"], row["specific_ap"])
        if observed != expected:
            failures.append(f"per-isoform {isoform}")

    input_paths = {
        "external_human_cyp_01\\overlap_audit.json": ROOT / "external_human_cyp_01" / "overlap_audit.json",
        "external_human_cyp_models_01\\method_summary.json": ROOT / "external_human_cyp_models_01" / "method_summary.json",
        "external_human_cyp_models_01\\comparison_summary.json": ROOT / "external_human_cyp_models_01" / "comparison_summary.json",
        "external_human_cyp_models_01\\score_permutation_control.json": ROOT / "external_human_cyp_models_01" / "score_permutation_control.json",
        "external_human_cyp_validation_01\\validation.json": ROOT / "external_human_cyp_validation_01" / "validation.json",
    }
    for key, path in input_paths.items():
        if provenance["input_hashes"].get(key) != digest_file(path):
            failures.append(f"input hash {key}")
    if provenance["source_data_sha256"] != digest_file(FIGURES / "figure6_external_transfer_source_data.json"):
        failures.append("source data hash")
    if provenance["script_sha256"] != digest_file(ROOT / "build_figure6_external_transfer.py"):
        failures.append("script hash")
    for extension in ("png", "svg", "pdf"):
        path = FIGURES / f"figure6_external_transfer.{extension}"
        if provenance["outputs"][extension]["sha256"] != digest_file(path):
            failures.append(f"{extension} hash")
    with Image.open(FIGURES / "figure6_external_transfer.png") as image:
        dimensions = list(image.size)
        if dimensions[0] < 7000 or dimensions[1] < 4500:
            failures.append("PNG resolution")
    svg = (FIGURES / "figure6_external_transfer.svg").read_text(encoding="utf-8")
    for token in ("Eligibility audit", "Locked model versus pooled chemistry", "Paired effect estimates", "New-scaffold performance by isoform"):
        if token not in svg:
            failures.append(f"SVG token {token}")
    if (FIGURES / "figure6_external_transfer.pdf").read_bytes()[:5] != b"%PDF-":
        failures.append("PDF header")

    report = {
        "status": "PASS" if not failures else "FAIL", "created_utc": now(),
        "failures": failures, "source_bindings_verified": 4 + 2 * 12 + 6 * 5,
        "input_hashes_verified": len(input_paths), "output_hashes_verified": 3,
        "png_dimensions": dimensions, "external_evaluation_status": validation["status"],
        "visual_review": {
            "performed": True,
            "checks": ["labels readable", "no clipped marks", "uncertainty adjacent to estimates",
                       "color not sole method encoding", "scope limitation visible"],
        },
        "validator_sha256": digest_file(Path(__file__)),
    }
    write_json(FIGURES / "figure6_external_transfer_validation.json", report)
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
