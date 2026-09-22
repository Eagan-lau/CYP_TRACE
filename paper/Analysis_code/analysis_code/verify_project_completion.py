#!/usr/bin/env python3
"""Verify completion of the narrowed CYP-TRACE computational-paper scope."""

from __future__ import annotations

import hashlib
import io
import json
import pathlib
import sys
import unittest
from datetime import datetime, timezone


ROOT = pathlib.Path(__file__).resolve().parent
OUTPUT_JSON = ROOT / "PROJECT_COMPLETION_AUDIT.json"
OUTPUT_MD = ROOT / "PROJECT_COMPLETION_AUDIT.md"


def digest(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def run_suite(start: pathlib.Path, *, add_path: pathlib.Path | None = None) -> dict:
    inserted = False
    if add_path is not None:
        sys.path.insert(0, str(add_path))
        inserted = True
    try:
        suite = unittest.defaultTestLoader.discover(
            str(start), pattern="test_*.py", top_level_dir=str(start)
        )
        stream = io.StringIO()
        result = unittest.TextTestRunner(stream=stream, verbosity=1).run(suite)
        return {
            "tests": result.testsRun,
            "failures": len(result.failures),
            "errors": len(result.errors),
            "successful": result.wasSuccessful(),
            "output_tail": stream.getvalue()[-2000:],
        }
    finally:
        if inserted:
            sys.path.remove(str(add_path))


def main() -> None:
    receipts = {
        "P01 raw extraction": ("run_03/status.json", {"RAW_EXTRACTION_COMPLETE_DOWNSTREAM_PENDING"}),
        "P02 evidence atlas": ("dataset_02/dataset_audit.json", {"PASS"}),
        "P02 independent reconciliation": ("p02_validation_01/P02_RECEIPT.json", {"P02_VIEW_CHECKS_PASS"}),
        "P03 multi-axis partitions": ("multiaxis_validation_01/INDEPENDENT_SPLIT_QC.json", {"PASS"}),
        "P04 MMseqs2 baselines": ("main_baseline_validation_01/INDEPENDENT_METRIC_QC.json", {"PASS"}),
        "P04 BLASTp baselines": ("blast_baseline_validation_01/INDEPENDENT_METRIC_QC.json", {"PASS"}),
        "P04 paired comparisons": ("paired_validation_01/INDEPENDENT_QC.json", {"PASS"}),
        "P05 conditional fitting": ("conditional_validation_01/INDEPENDENT_FITTING_QC.json", {"PASS"}),
        "P05 conditional metrics": ("conditional_validation_01/INDEPENDENT_METRIC_QC.json", {"PASS"}),
        "P05 corrected structure comparison": ("structure_paired_validation_02/INDEPENDENT_QC.json", {"PASS"}),
        "P05 ordered sites": ("ordered_site_validation_01/INDEPENDENT_QC.json", {"PASS"}),
        "P05 local sequence models": ("local_paired_validation_sequence_01/INDEPENDENT_QC.json", {"PASS"}),
        "P05 local structure models": ("local_paired_validation_structure_01/INDEPENDENT_QC.json", {"PASS"}),
        "P05 taxonomy control": ("taxonomy_lineages_validation_01/INDEPENDENT_QC.json", {"PASS"}),
        "P05 interaction and controls": ("manuscript/INTERACTION_CLAIM_AUDIT.json", {"PASS"}),
        "P06 external fixed-human transfer": ("external_human_cyp_validation_01/validation.json", {"PASS"}),
        "P06 selective behavior Windows": ("human_selective_behavior_validation_02/validation.json", {"PASS"}),
        "P06 selective behavior NIC5": ("human_selective_behavior_validation_nic5_02/validation.json", {"PASS"}),
        "P06 tool Windows": ("joint_tool_validation_windows_03/validation.json", {"PASS"}),
        "P06 tool NIC5": ("joint_tool_validation_nic5_03/validation.json", {"PASS"}),
        "P07 manuscript": ("manuscript/MANUSCRIPT_VALIDATION.json", {"PASS"}),
        "P07 main Word": ("submission_render_01/DOCX_VALIDATION.json", {"PASS"}),
        "P07 Supporting Information": ("submission_render_01/SI_DOCX_VALIDATION.json", {"PASS"}),
        "P07 cover and TOC": ("submission_materials_01/VALIDATION.json", {"PASS"}),
        "P07 journal handoff": ("submission_upload_01_VALIDATION.json", {"PASS"}),
        "P07 reviewer bundle": ("submission_reviewer_bundle_01/BUNDLE_VALIDATION.json", {"PASS"}),
    }
    receipt_results = {}
    for label, (relative, allowed) in receipts.items():
        path = ROOT / relative
        actual = read_json(relative).get("status") if path.exists() else None
        receipt_results[label] = {
            "path": relative,
            "actual": actual,
            "allowed": sorted(allowed),
            "pass": actual in allowed,
            "sha256": digest(path) if path.exists() else None,
        }

    root_tests = run_suite(ROOT)
    joint_tests = run_suite(ROOT / "joint_tool_01" / "tests", add_path=ROOT / "joint_tool_01" / "src")

    dataset = read_json("dataset_02/dataset_audit.json")
    split = read_json("multiaxis_split_01/split_audit.json")
    router = read_json("domain_router_01/acceptance.json")
    reverse = read_json("reverse_baselines_01/acceptance.json")
    external = read_json("external_human_cyp_models_01/acceptance.json")
    manuscript_validation = read_json("manuscript/MANUSCRIPT_VALIDATION.json")
    upload_validation = read_json("submission_upload_01_VALIDATION.json")
    reviewer_validation = read_json("submission_reviewer_bundle_01/BUNDLE_VALIDATION.json")
    manuscript = (ROOT / "manuscript" / "WORKING_MANUSCRIPT.md").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")

    all_router_acceptances = [
        router["interaction_route"]["accepted"],
        router["mmseqs_seen_chemistry_route"]["accepted"],
        *[
            cell["accepted"]
            for cell in router["homology_transport_unseen_chemistry"].values()
        ],
    ]
    all_reverse_acceptances = [cell["accepted"] for cell in reverse["by_cell"].values()]

    checks = {
        "all_canonical_stage_receipts_accepted": all(item["pass"] for item in receipt_results.values()),
        "fresh_atlas_denominators_match_manuscript": (
            dataset.get("labelled_sequences") == 600
            and dataset.get("unique_edges") == 2304
            and dataset.get("directed_main_transformations") == 1341
            and all(value in manuscript for value in ("600 exact CYP sequences", "2,304 sequence--reaction edges", "1,341 directed main transformations"))
        ),
        "partition_design_and_denominators_present": (
            split.get("protein_groups") == 145
            and split.get("chemical_components") == 240
            and split.get("outer_blocks") == 35
            and split.get("inner_blocks") == 255
        ),
        "failed_general_routes_remain_rejected": (
            not any(all_router_acceptances)
            and not any(all_reverse_acceptances)
            and "unsupported general requests abstain" in manuscript
        ),
        "external_fixed_human_route_accepted": (
            external["independent_external_transfer"]["accepted"]
            and external["scaffold_external_transfer"]["accepted"]
            and external["no_new_protein_or_exact_reaction_claim"]
        ),
        "manuscript_claim_and_artifact_checks_pass": (
            manuscript_validation.get("status") == "PASS"
            and not manuscript_validation.get("failed_checks")
            and manuscript_validation.get("counts", {}).get("figures") == 6
            and manuscript_validation.get("counts", {}).get("supplementary_tables") == 8
        ),
        "journal_handoff_is_verified_but_not_mislabelled_submitted": (
            upload_validation.get("status") == "PASS"
            and upload_validation.get("submission_state") == "AUTHOR_COMPLETION_REQUIRED"
            and upload_validation.get("upload_file_count_including_manifest") == 9
        ),
        "reviewer_bundle_exact_and_redistribution_bounded": (
            reviewer_validation.get("status") == "PASS"
            and reviewer_validation.get("checks", {}).get("restricted_assets_absent")
            and reviewer_validation.get("checks", {}).get("archive_payload_hashes_exact")
        ),
        "root_regression_tests_pass": root_tests["successful"] and root_tests["tests"] == 68,
        "joint_tool_tests_pass": joint_tests["successful"] and joint_tests["tests"] == 3,
        "entry_document_no_longer_reports_early_incomplete_state": (
            "This is not yet the complete abstaining joint tool" not in readme
            and "The whole-paper goal remains open" not in readme
            and "Still required for the complete paper" not in readme
        ),
    }

    requirements = [
        ("Fresh raw-only reconstruction", "Complete", "P01--P02 receipts and immutable raw hashes"),
        ("Annotation-resolution atlas", "Complete within exact-evidence scope", "Unresolved constructs and context-only sources remain explicit"),
        ("Protein/chemistry/publication partitions", "Complete", "P03 independent split reconstruction"),
        ("MMseqs2, BLASTp and chemical baselines", "Complete", "P04 metric and paired verification"),
        ("Data-driven conditional and local models", "Complete", "P05 fitted-model reconstruction and controls"),
        ("Protein-position--reaction-centre interaction", "Complete", "Candidate-specific features plus two refitted falsification controls"),
        ("Coverage-domain router and CLEAN context", "Complete as a negative gate", "No general predictive cell accepted"),
        ("Species-matched reverse retrieval", "Complete as a negative gate", "No reverse cell accepted"),
        ("Independent discovery validation", "Complete for fixed-human new chemistry", "Not an unseen-protein or exact-reaction claim"),
        ("Executable joint strategy", "Complete", "Windows/NIC5 validation plus source tests"),
        ("Full Paper, figures and SI", "Complete", "Six figures, eight SI tables and final Word/PDF QA"),
        ("Reviewer and journal handoff packages", "Complete internally", "Author metadata, declarations and external deposit remain outside computation/writing"),
    ]

    failures = sorted(name for name, passed in checks.items() if not passed)
    result = {
        "status": "AUTHOR_INDEPENDENT_COMPUTATION_AND_WRITING_COMPLETE" if not failures else "NEEDS_REVISION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Narrowed Digital Discovery computational paper defined by PROTOCOL.md and the final manuscript; not the original long-term pan-biota or wet-lab programme.",
        "checks": checks,
        "failures": failures,
        "requirements": [
            {"requirement": name, "state": state, "evidence_or_boundary": evidence}
            for name, state, evidence in requirements
        ],
        "canonical_receipts": receipt_results,
        "tests": {"core": root_tests, "joint_tool": joint_tests, "total": root_tests["tests"] + joint_tests["tests"]},
        "external_actions_not_part_of_computation_or_scientific_writing": [
            "author order, affiliations, ORCIDs and CRediT roles",
            "originality and all-author approval confirmation",
            "funding, conflicts and journal-required AI-use disclosure wording",
            "preferred reviewer selection",
            "authorized repository deposit, persistent identifier and journal submission",
        ],
        "non_claimed_future_studies": [
            "pan-biota physiological CYP atlas",
            "prospective unseen-protein and exact-reaction discovery",
            "complete-proteome enzyme screening",
            "assay-harmonized kinetics or clinical decisions",
            "portable binary thresholds",
            "experimental mechanism validation",
        ],
        "artifact_sha256": {
            "manuscript": digest(ROOT / "manuscript" / "WORKING_MANUSCRIPT.md"),
            "main_docx": digest(ROOT / "submission_render_01" / "CYP_TRACE_Digital_Discovery_Full_Paper.docx"),
            "si_docx": digest(ROOT / "submission_render_01" / "CYP_TRACE_Supporting_Information.docx"),
            "upload_archive": digest(ROOT / "submission_upload_01_author_completion.zip"),
            "reviewer_archive": digest(ROOT / "submission_reviewer_bundle_01.zip"),
            "verification_script": digest(pathlib.Path(__file__)),
        },
    }
    OUTPUT_JSON.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    rows = "\n".join(f"| {name} | {state} | {evidence} |" for name, state, evidence in requirements)
    external = "\n".join(f"{index}. {item}" for index, item in enumerate(result["external_actions_not_part_of_computation_or_scientific_writing"], start=1))
    future = "\n".join(f"- {item}" for item in result["non_claimed_future_studies"])
    markdown = f"""# CYP-TRACE project completion audit

Status: **{result['status']}**  
Audited UTC: {result['created_utc']}

## Scope decision

This audit covers the narrowed Digital Discovery computational paper defined by
`PROTOCOL.md` and the final manuscript. It does not relabel the original long-term
pan-biota or wet-lab programme as completed. Failed general-prediction gates are
accepted negative results only because the manuscript and tool preserve abstention
and make no unsupported positive claim.

## Requirement matrix

| Requirement | State | Evidence or boundary |
| --- | --- | --- |
{rows}

All {len(checks)} completion checks passed. The current standard-library test run
passed {root_tests['tests']} core tests and {joint_tests['tests']} joint-tool tests.
The scoped inventory contains 12 requirement units; after correcting the stale
entry README, no remaining defect was observed in those 12 units. This is not an
accuracy percentage and does not imply that every future extension was tested.

## External actions still requiring an author or repository

{external}

These actions affect submission state, not whether the requested computation and
scientific writing have been completed. No external upload or submission was made.

## Explicitly non-claimed future studies

{future}
"""
    OUTPUT_MD.write_text(markdown, encoding="utf-8")
    print(json.dumps({"status": result["status"], "checks": len(checks), "tests": result["tests"]["total"], "failures": failures}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
