"""Verify supplementary table integrity and key source-derived invariants."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "manuscript" / "supplementary_tables"
MANIFEST = TABLES / "manifest.json"
OUTPUT = TABLES / "validation.json"


def rows(name):
    with (TABLES / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    checks = {}
    for item in manifest["tables"]:
        path = ROOT / item["path"]
        assert path.exists()
        assert digest_file(path) == item["sha256"]
        assert len(rows(path.name)) == item["rows"]
    checks["manifest_hashes_and_rows"] = True

    s1 = {item["resource"]: item for item in rows("table_s1_source_scope_and_release.tsv")}
    assert len(s1) == 42
    assert s1["rhea"]["licence_or_access_status"] == "CC BY 4.0"
    assert s1["uniprot"]["licence_or_access_status"] == "CC BY 4.0"
    assert s1["pdb"]["licence_or_access_status"] == "CC0 1.0"
    assert s1["brenda"]["release_treatment"] == "Private workspace only"
    assert s1["figshare_26630515_v4"]["licence_or_access_status"] == "CC BY 4.0"
    checks["source_scope_release_controls"] = True

    s2 = {item["source"]: item for item in rows("table_s2_annotation_resolution.tsv")}
    assert len(s2) == 4 and s2["SABIO_RK"]["admitted_assertions"] == "0"
    assert abs(float(s2["P450Rdb"]["exact_component_set_edge_matches_initial_core"]) - 0.5030096833289714) < 1e-15
    checks["annotation_resolution"] = True

    s3 = rows("table_s3_split_denominators.tsv")
    assert len(s3) == 3 and all(item["test_edges"] == "2304" for item in s3)
    assert {item["outer_blocks"] for item in s3} == {"5", "25"}
    checks["split_denominators"] = True

    s4 = rows("table_s4_general_model_effects.tsv")
    assert len(s4) == 9
    assert {item["comparison"] for item in s4} == {"global_conditional_minus_chemical_prior", "interaction_minus_global", "ordered_sites_minus_global"}
    checks["general_effects"] = True

    s5 = {item["expert"]: int(item["available_sequences"]) for item in rows("table_s5_expert_availability.tsv")}
    assert s5 == {"ESM2 global": 600, "CLEAN": 582, "ordered sequence-site vector": 289, "primary associated structure": 36, "ordered structure-site vector": 23}
    checks["expert_availability"] = True

    s6 = rows("table_s6_reverse_retrieval.tsv")
    assert len(s6) == 3 and all(item["accepted"] == "False" for item in s6)
    checks["reverse_acceptance"] = True

    s7 = {item["stratum"]: item for item in rows("table_s7_external_transfer.tsv")}
    assert int(s7["lineage_eligible"]["labels"]) == 4301
    assert int(s7["scaffold_and_lineage_eligible"]["labels"]) == 3035
    assert abs(float(s7["scaffold_and_lineage_eligible"]["ap_difference"]) - 0.13190865018235037) < 1e-12
    checks["external_transfer"] = True

    s8 = rows("table_s8_selective_behavior.tsv")
    assert len(s8) == 10
    strict = next(item for item in s8 if item["method"] == "isoform_specific_knn" and item["target_coverage"] == "0.1")
    assert abs(float(strict["macro_precision"]) - 0.668518459558566) < 1e-15
    assert abs(float(strict["macro_positive_recall"]) - 0.20319902636770767) < 1e-15
    checks["selective_behavior"] = True

    supplement = (ROOT / "manuscript" / "SUPPORTING_INFORMATION.md").read_text(encoding="utf-8")
    for item in manifest["tables"]:
        assert Path(item["path"]).name in supplement
    checks["supplement_links_all_tables"] = True

    result = {
        "status": "PASS",
        "created_utc": now(),
        "checks": checks,
        "tables_verified": len(manifest["tables"]),
        "manifest_sha256": digest_file(MANIFEST),
        "supplement_sha256": digest_file(ROOT / "manuscript" / "SUPPORTING_INFORMATION.md"),
        "script_sha256": digest_file(Path(__file__)),
    }
    write_json(OUTPUT, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
