"""Build compact, source-backed supplementary tables for the manuscript."""
from __future__ import annotations

import csv
import json
from pathlib import Path

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "manuscript" / "supplementary_tables"


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_tsv(name, fields, rows):
    path = OUT / name
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def source_inventory():
    inventory = read(ROOT / "run_03" / "status.json")["stages"]["raw_inventory"]["resources"]
    overrides = {
        "p450rdb": ("Core sequence--reaction assertions", "yes", "Academic offline access; explicit redistribution licence not certified", "Raw records excluded; derived identifiers require final permission review"),
        "rhea": ("Directional reaction structures and identifiers", "yes", "CC BY 4.0", "Attribution required"),
        "uniprot": ("Sequence and catalytic-activity annotations", "yes", "CC BY 4.0", "Attribution required"),
        "uniprot_swissprot": ("Reviewed UniProt flat-file lineage", "yes", "CC BY 4.0", "Attribution required; not counted as independent from UniProt API"),
        "sabio_rk": ("Annotation-resolution audit; not admitted to exact core", "yes", "Redistribution terms not certified", "Raw records excluded"),
        "brenda": ("Licensed EC/context audit; not upgraded to protein--reaction truth", "yes", "Licensed user access", "Private workspace only"),
        "pdb": ("Experimentally determined structure availability", "yes", "CC0 1.0", "Attribution to entries and archive retained"),
        "rcsb_cyp_ligand_templates_v1": ("Heme/ligand geometry templates derived from PDB", "yes", "PDB archive CC0 1.0", "Derived identifiers and geometry may be released with attribution"),
        "alphafold": ("Predicted-structure availability and sensitivity context", "yes", "CC BY 4.0", "Attribution required"),
        "chebi": ("Chemical identifiers and reaction participants", "yes", "CC BY 4.0", "Attribution required"),
        "pubchem": ("Structure resolution and identity cross-checks", "yes", "Contributor-specific terms", "No blanket redistribution assumption"),
        "pfam_interpro": ("PF00067 CYP-family support", "yes", "Release terms require final record-level review", "Derived family-support flag only"),
        "ncbi_taxonomy": ("Taxonomic identifier reconciliation", "yes", "US government resource; policy review retained", "Identifiers only"),
        "crossref": ("DOI metadata reconciliation", "yes", "Metadata terms apply", "Identifiers and hashes only"),
        "europe_pmc": ("Publication metadata and full-text availability checks", "yes", "Record/article-specific rights", "No article text redistributed"),
        "pubmed": ("Publication identifier reconciliation", "yes", "NLM policies apply", "Identifiers only"),
    }
    rows = []
    default = ("Downloaded inventory; not used in primary numerical claims", "no", "Not adjudicated for release", "Excluded from manuscript redistribution bundle")
    for resource in sorted(inventory):
        role, used, licence, treatment = overrides.get(resource, default)
        values = inventory[resource]
        rows.append({
            "resource": resource,
            "raw_files": values["files"],
            "raw_bytes": values["bytes"],
            "empty_files": values["empty_files"],
            "primary_manuscript_role": role,
            "used_in_main_analysis": used,
            "licence_or_access_status": licence,
            "release_treatment": treatment,
        })
    external = read(ROOT / "external_human_cyp_01" / "acquisition_manifest.json")
    rows.extend([
        {
            "resource": "cypstrate",
            "raw_files": "not enumerated in run_03 inventory",
            "raw_bytes": "not enumerated in run_03 inventory",
            "empty_files": "not applicable",
            "primary_manuscript_role": "Nine-isoform fixed-human development labels",
            "used_in_main_analysis": "yes",
            "licence_or_access_status": "Publicly accessible source; redistribution terms not certified",
            "release_treatment": "Raw records excluded; derived model bundle carries source notice",
        },
        {
            "resource": "figshare_26630515_v4",
            "raw_files": len(external["files"]),
            "raw_bytes": sum(item["observed_bytes"] for item in external["files"]),
            "empty_files": 0,
            "primary_manuscript_role": "Independently qualified fixed-human external labels",
            "used_in_main_analysis": "yes",
            "licence_or_access_status": external["license"]["name"],
            "release_treatment": "Attribution and source DOI required",
        },
        {
            "resource": "CLEAN_v1.0.0",
            "raw_files": "external software assets",
            "raw_bytes": "not included",
            "empty_files": "not applicable",
            "primary_manuscript_role": "EC-context/pruning analysis",
            "used_in_main_analysis": "yes",
            "licence_or_access_status": "Non-exclusive research-use licence",
            "release_treatment": "Weights and licensed assets excluded",
        },
        {
            "resource": "ESM2_650M",
            "raw_files": "external model assets",
            "raw_bytes": "not included",
            "empty_files": "not applicable",
            "primary_manuscript_role": "Global protein embeddings",
            "used_in_main_analysis": "yes",
            "licence_or_access_status": "Upstream model terms apply",
            "release_treatment": "Checkpoint excluded; feature hashes and specification retained",
        },
    ])
    return rows


def evidence_resolution():
    dataset = read(ROOT / "dataset_02" / "dataset_audit.json")
    annotation = read(ROOT / "annotation_02" / "annotation_resolution_audit.json")
    metrics = ["unique_sequence_link_available", "linked_sequence_has_PF00067", "has_recorded_publication", "complete_component_set_parseable", "exact_component_set_edge_matches_initial_core"]
    rows = []
    for source in ("P450Rdb", "SABIO_RK", "SwissProt", "UniProt_API"):
        entry = {"source": source, "source_assertions": annotation["source_counts"][source]["source_assertions"], "admitted_assertions": dataset["source_counts"].get(source, {}).get("admitted_assertions", 0)}
        entry.update({metric: annotation["source_assertion_rates"][source][metric]["fraction"] for metric in metrics})
        rows.append(entry)
    return rows


def split_summary():
    split = read(ROOT / "multiaxis_split_01" / "split_audit.json")
    rows = []
    for task, values in split["task_denominators"].items():
        rows.append({"task": task, "outer_blocks": values["blocks"], "evaluable_blocks": values["evaluable_blocks"], "test_edges": values["evaluable_test_edges"], "query_panels": values["evaluable_query_panel_rows"], "unique_queries": values["unique_query_sequences"]})
    return rows


def general_effects():
    paired = read(ROOT / "paired_comparison_01" / "comparison_summary.json")["comparisons"]
    interaction = read(ROOT / "interaction_paired_sequence_01" / "comparison_summary.json")["comparisons"]
    local = read(ROOT / "local_paired_sequence_01" / "comparison_summary.json")["comparisons"]
    specifications = [
        (paired, "global_conditional_minus_chemical_prior", "esm_conditional_learned", "conditional_chemistry_prior"),
        (interaction, "interaction_minus_global", "joint_global_interaction", "global_residual"),
        (local, "ordered_sites_minus_global", "ordered_residual", "global_residual"),
    ]
    rows = []
    for source, comparison, left, right in specifications:
        for task in ("chemical_cold_all", "protein_cold_all", "double_cold_all"):
            item = next(x for x in source if x["task"] == task and x["left_method"] == left and x["right_method"] == right)
            effect = item["end_to_end"]
            rows.append({"comparison": comparison, "task": task, "left_method": left, "right_method": right, "query_panels": effect["query_panels"], "protein_groups": effect["protein_groups"], "mrr_difference": effect["difference"], "ci_low": effect["conditional_protein_bootstrap_95ci"][0], "ci_high": effect["conditional_protein_bootstrap_95ci"][1]})
    return rows


def expert_coverage():
    structure = read(ROOT / "structure_mapping_audit_02" / "audit.json")
    ordered = read(ROOT / "ordered_site_features_01" / "audit.json")
    clean = read(ROOT / "clean_core_features_01" / "audit.json")
    esm = read(ROOT / "esm_global_01" / "receipt.json")
    return [
        {"expert": "ESM2 global", "available_sequences": esm["sequences"], "denominator": 600, "availability_fraction": esm["sequences"]/600},
        {"expert": "CLEAN", "available_sequences": clean["supported_sequences"], "denominator": 600, "availability_fraction": clean["supported_sequences"]/600},
        {"expert": "ordered sequence-site vector", "available_sequences": ordered["counts"]["sequence_projected_available_sequences"], "denominator": 600, "availability_fraction": ordered["counts"]["sequence_projected_available_sequences"]/600},
        {"expert": "primary associated structure", "available_sequences": structure["primary_structure_sequences"], "denominator": 600, "availability_fraction": structure["primary_structure_sequences"]/600},
        {"expert": "ordered structure-site vector", "available_sequences": ordered["counts"]["structure_projected_available_sequences"], "denominator": 600, "availability_fraction": ordered["counts"]["structure_projected_available_sequences"]/600},
    ]


def reverse_summary():
    summaries = read(ROOT / "reverse_baselines_01" / "method_summary.json")["summaries"]
    acceptance = read(ROOT / "reverse_baselines_01" / "acceptance.json")["by_cell"]
    rows = []
    for cell, result in acceptance.items():
        method = result["primary_method"]
        source = next(x for x in summaries if x["cell"] == cell and x["method"] == method)
        rows.append({"cell": cell, "primary_method": method, "taxids": source["taxids"], "taxid_macro_mrr": source["taxid_macro_mrr"], "documented_positive_coverage": source["documented_positive_coverage"], "permutation_exact_p": result["permutation_exact_p"], "accepted": result["accepted"]})
    return rows


def external_summary():
    audit = read(ROOT / "external_human_cyp_01" / "overlap_audit.json")
    summaries = read(ROOT / "external_human_cyp_models_01" / "method_summary.json")["strata"]
    comparisons = read(ROOT / "external_human_cyp_models_01" / "comparison_summary.json")["strata"]
    rows = []
    for stratum in ("lineage_eligible", "scaffold_and_lineage_eligible"):
        rows.append({"stratum": stratum, "labels": audit["strata"][stratum]["rows"], "compounds": audit["strata"][stratum]["unique_compounds"], "scaffolds": audit["strata"][stratum]["scaffold_groups"], "pooled_macro_ap": summaries[stratum]["pooled_chemical_knn"]["macro_average_precision"], "isoform_specific_macro_ap": summaries[stratum]["isoform_specific_knn"]["macro_average_precision"], "ap_difference": comparisons[stratum]["average_precision"]["difference"], "ci_low": comparisons[stratum]["average_precision"]["scaffold_bootstrap_95ci"][0], "ci_high": comparisons[stratum]["average_precision"]["scaffold_bootstrap_95ci"][1]})
    return rows


def selective_summary():
    data = read(ROOT / "human_selective_behavior_02" / "analysis.json")["strata"]["scaffold_and_lineage_eligible"]
    comparisons = data["matched_coverage_comparison"]
    rows = []
    for item in data["label_blind_matched_coverage_curve"]:
        comparison = comparisons[str(item["target_coverage"])]
        rows.append({"method": item["method"], "target_coverage": item["target_coverage"], "actual_macro_coverage": item["macro_coverage"], "macro_precision": item["macro_precision"], "macro_selective_risk": item["macro_selective_risk"], "macro_positive_recall": item["macro_positive_recall"], "isoform_specific_minus_pooled_precision": comparison["isoform_specific_minus_pooled_macro_precision"], "comparison_ci_low": comparison["scaffold_bootstrap_95ci"][0], "comparison_ci_high": comparison["scaffold_bootstrap_95ci"][1]})
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tables = {
        "table_s1_source_scope_and_release.tsv": source_inventory(),
        "table_s2_annotation_resolution.tsv": evidence_resolution(),
        "table_s3_split_denominators.tsv": split_summary(),
        "table_s4_general_model_effects.tsv": general_effects(),
        "table_s5_expert_availability.tsv": expert_coverage(),
        "table_s6_reverse_retrieval.tsv": reverse_summary(),
        "table_s7_external_transfer.tsv": external_summary(),
        "table_s8_selective_behavior.tsv": selective_summary(),
    }
    paths = []
    for name, rows in tables.items():
        if not rows:
            raise RuntimeError(f"empty table: {name}")
        paths.append(write_tsv(name, list(rows[0]), rows))
    manifest = {
        "status": "PASS",
        "created_utc": now(),
        "script_sha256": digest_file(Path(__file__)),
        "tables": [{"path": str(path.relative_to(ROOT)), "rows": sum(1 for _ in path.open(encoding="utf-8"))-1, "sha256": digest_file(path)} for path in paths],
        "license_scope": "Release treatment is conservative and does not replace provider terms or legal review.",
    }
    write_json(OUT / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
