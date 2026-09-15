"""Independent source-binding and file-integrity checks for manuscript Figures 1--5."""
from __future__ import annotations

import json
import math
from pathlib import Path

from PIL import Image

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
PROVENANCE = FIGURES / "figures1_5_provenance.json"
OUTPUT = FIGURES / "figures1_5_validation.json"


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def close(left, right, tolerance=1e-12):
    return math.isclose(float(left), float(right), rel_tol=0, abs_tol=tolerance)


def one(rows, **match):
    hits = [item for item in rows if all(item.get(key) == value for key, value in match.items())]
    if len(hits) != 1:
        raise AssertionError(f"expected one row for {match}, found {len(hits)}")
    return hits[0]


def verify_source_bindings(checks):
    f1 = read(FIGURES / "figure1_evidence_atlas_source_data.json")
    dataset = read(ROOT / "dataset_02" / "dataset_audit.json")
    annotation = read(ROOT / "annotation_02" / "annotation_resolution_audit.json")
    source_map = {item["source"]: item for item in f1["source_assertions"]}
    for source in ("P450Rdb", "SABIO_RK", "SwissProt", "UniProt_API"):
        assert source_map[source]["source_assertions"] == annotation["source_counts"][source]["source_assertions"]
    assert source_map["P450Rdb"]["admitted_assertions"] == dataset["source_counts"]["P450Rdb"]["admitted_assertions"]
    assert source_map["SABIO_RK"]["admitted_assertions"] == 0
    for i, source in enumerate(f1["rate_sources"]):
        for j, metric in enumerate(f1["rate_metrics"]):
            assert close(f1["rates"][i][j], annotation["source_assertion_rates"][source][metric]["fraction"])
    overlap = one(dataset["cross_source_edge_overlap"], source_a="P450Rdb", source_b="UniProt_API")
    assert f1["edge_union"]["shared"] == overlap["shared_unique_edges"] == 797
    assert f1["edge_union"]["p450_unique"] + f1["edge_union"]["shared"] + f1["edge_union"]["uniprot_unique"] == dataset["unique_edges"]
    assert f1["core"] == {
        "sequences": dataset["labelled_sequences"], "edges": dataset["unique_edges"],
        "transformations": dataset["directed_main_transformations"],
        "single_pair": dataset["single_pair_transformations"], "multicomponent": dataset["multicomponent_transformations"]}
    checks["figure1_source_values"] = True

    f2 = read(FIGURES / "figure2_split_design_source_data.json")
    groups = read(ROOT / "protein_split_02" / "sequence_groups.json")["0.4"]
    split = read(ROOT / "multiaxis_split_01" / "split_audit.json")
    protein = read(ROOT / "protein_split_02" / "search_and_protein_split_audit.json")
    assert f2["protein_group_sizes"] == sorted((len(values) for values in groups.values()), reverse=True)
    assert f2["chemical_component_sizes"] == split["chemical_component_sizes_descending"]
    assert f2["protein_groups"] == len(groups) == 145 and f2["chemical_components"] == 240
    assert f2["outer_blocks"] == 35 and f2["inner_blocks"] == 255
    for observed, expected in zip(f2["folds"], protein["fold_audits"]):
        assert observed["fold"] == expected["fold"]
        assert observed["sequences"] == expected["query_sequences"]
        assert observed["groups"] == expected["query_groups"]
        assert observed["publication_purged_edges"] == expected["publication_purged_training_edges"]
        assert close(observed["max_cross_identity"], expected["maximum_detected_bilateral_cross_partition_identity"])
    task_map = {item["task"]: item for item in f2["tasks"]}
    for task, values in split["task_denominators"].items():
        assert task_map[task]["blocks"] == values["blocks"]
        assert task_map[task]["query_panels"] == values["evaluable_query_panel_rows"]
        assert task_map[task]["test_edges"] == values["evaluable_test_edges"] == 2304
    checks["figure2_source_values"] = True

    f3 = read(FIGURES / "figure3_domain_methods_source_data.json")
    router = read(ROOT / "domain_router_01" / "method_summary.json")["summaries"]
    for cell in f3["router_cells"]:
        for method, observed in cell["mrr"].items():
            if observed is None:
                assert not [item for item in router if item["cell"] == cell["cell"] and item["method"] == method]
                continue
            source_method = method
            if method == "chemistry_prior" and cell["cell"] == "represented_protein_unseen_chemistry":
                source_method = "availability_chemistry_prior"
            expected = one(router, cell=cell["cell"], method=source_method)
            assert close(observed, expected["protein_group_macro"]["rr"])
            assert cell["query_panels"] == expected["query_panels"]
    paired = read(ROOT / "paired_comparison_01" / "comparison_summary.json")["comparisons"]
    interaction = read(ROOT / "interaction_paired_sequence_01" / "comparison_summary.json")["comparisons"]
    for item in f3["global_effects"]:
        expected = one(paired, task=item["task"], left_method="esm_conditional_learned", right_method="conditional_chemistry_prior")["end_to_end"]
        assert close(item["difference"], expected["difference"])
        assert all(close(a, b) for a, b in zip(item["ci"], expected["conditional_protein_bootstrap_95ci"]))
    for item in f3["interaction_effects"]:
        expected = one(interaction, task=item["task"], left_method="joint_global_interaction", right_method="global_residual")["end_to_end"]
        assert close(item["difference"], expected["difference"])
    controls = read(ROOT / "interaction_control_comparison_01" / "sequence_summary.json")["summaries"]
    for item in f3["control_effects"]:
        expected = one(controls, control=item["control"], method="joint_global_interaction", task="chemical_cold_all")["paired_gain_comparison"]
        assert close(item["difference"], expected["difference"])
        assert all(close(a, b) for a, b in zip(item["ci"], expected["conditional_protein_bootstrap_95ci"]))
    availability = {item["label"]: item["sequences"] for item in f3["availability"]}
    assert availability == {"ESM global": 600, "CLEAN supported": 582, "Ordered sequence sites": 289, "Ordered structure sites": 23}
    checks["figure3_source_values"] = True

    f4 = read(FIGURES / "figure4_reverse_retrieval_source_data.json")
    census = read(ROOT / "reverse_candidate_census_01" / "audit.json")
    assert f4["census"]["raw_proteins"] == census["raw_proteins"] == 2447
    assert f4["census"]["admitted_candidates"] == census["admitted_candidate_sequences"] == 1803
    assert f4["census"]["admitted_core"] + f4["census"]["admitted_noncore"] == f4["census"]["admitted_candidates"]
    assert f4["census"]["query_panels"] == census["outer_query_panels"] == 4880
    reverse = read(ROOT / "reverse_baselines_01" / "method_summary.json")["summaries"]
    method_by_cell = {
        "seen_reaction__new_positive": "mmseqs_weighted_positive",
        "unseen_reaction__new_positive": "mmseqs_homology_chemical_transport",
        "unseen_reaction__represented_positive": "mmseqs_homology_chemical_transport"}
    for item in f4["performance"]:
        expected = one(reverse, cell=item["cell"], method=method_by_cell[item["cell"]])
        assert close(item["mrr"]["MMseqs2"], expected["taxid_macro_mrr"])
        assert close(item["coverage"]["MMseqs2"], expected["documented_positive_coverage"])
    accepted = read(ROOT / "reverse_baselines_01" / "acceptance.json")["by_cell"]
    for i, item in enumerate(f4["performance"]):
        assert f4["gates"][i] == [accepted[item["cell"]]["gates"][name] for name in f4["gate_names"]]
        assert f4["accepted"][i] is False
    checks["figure4_source_values"] = True

    f5 = read(FIGURES / "figure5_coverage_controls_source_data.json")
    structure = read(ROOT / "structure_mapping_audit_02" / "audit.json")
    ordered = read(ROOT / "ordered_site_features_01" / "audit.json")
    clean = read(ROOT / "clean_core_features_01" / "audit.json")
    availability = {item["label"]: item["count"] for item in f5["availability"]}
    assert availability["Core sequences"] == structure["core_sequences"] == 600
    assert availability["CLEAN supported"] == clean["supported_sequences"] == 582
    assert availability["Sequence-site vector"] == ordered["counts"]["sequence_projected_available_sequences"] == 289
    assert availability["Primary structure"] == structure["primary_structure_sequences"] == 36
    assert availability["Structure-site vector"] == ordered["counts"]["structure_projected_available_sequences"] == 23
    local = read(ROOT / "local_paired_sequence_01" / "comparison_summary.json")["comparisons"]
    for item, task in zip(f5["local_effects"][:3], ("chemical_cold_all", "protein_cold_all", "double_cold_all")):
        expected = one(local, task=task, left_method="ordered_residual", right_method="global_residual")["end_to_end"]
        assert close(item["difference"], expected["difference"])
    structure_pairs = read(ROOT / "structure_paired_02" / "comparison_summary.json")["comparisons"]
    expected = one(structure_pairs, task="protein_cold_seen", left_method="foldseek_top1", right_method="matched_mmseqs_top1")["end_to_end"]
    assert close(f5["local_effects"][3]["difference"], expected["difference"])
    pruning = read(ROOT / "clean_pruning_01" / "comparison_summary.json")["comparisons"]
    for item in f5["pruning_effects"]:
        task_lookup = {"Chemical cold": "new_protein_unseen_chemistry_with_homologue__chemical_cold_accounting", "Protein cold": "new_protein_unseen_chemistry_with_homologue__protein_cold", "Double cold": "new_protein_unseen_chemistry_with_homologue__double_cold"}
        expected = one(pruning, cell=task_lookup[item["label"]])
        for source_key, output_key in (("all", "all"), ("clean_exact_exposed", "exposed"), ("clean_exact_unexposed", "unexposed")):
            assert close(item[output_key]["difference"], expected[source_key]["difference"])
    selective = read(ROOT / "human_selective_behavior_02" / "analysis.json")["strata"]["scaffold_and_lineage_eligible"]
    assert len(f5["selective_curve"]) == len(selective["label_blind_matched_coverage_curve"]) == 10
    for observed, expected in zip(f5["selective_curve"], selective["label_blind_matched_coverage_curve"]):
        for key in ("method", "target_coverage"):
            assert observed[key] == expected[key]
        assert close(observed["coverage"], expected["macro_coverage"])
        assert close(observed["precision"], expected["macro_precision"])
        assert close(observed["positive_recall"], expected["macro_positive_recall"])
    checks["figure5_source_values"] = True


def main():
    provenance = read(PROVENANCE)
    checks = {}
    for entry in provenance["input_hashes"].values():
        assert digest_file(ROOT / entry["path"]) == entry["sha256"]
    checks["all_input_hashes"] = True
    verify_source_bindings(checks)
    dimensions = {}
    for stem, entry in provenance["figures"].items():
        source_path = ROOT / entry["source_data"]["path"]
        assert digest_file(source_path) == entry["source_data"]["sha256"]
        dimensions[stem] = {}
        for extension, output in entry["files"].items():
            path = ROOT / output["path"]
            assert path.exists() and path.stat().st_size > 1000
            assert digest_file(path) == output["sha256"]
            if extension == "png":
                with Image.open(path) as image:
                    width, height = image.size
                    assert width >= 4000 and height >= 3000
                    dimensions[stem] = {"width": width, "height": height}
            elif extension == "pdf":
                assert path.read_bytes()[:4] == b"%PDF"
            elif extension == "svg":
                assert "<svg" in path.read_text(encoding="utf-8")[:1000]
    checks["all_source_and_output_hashes"] = True
    checks["all_three_formats_present"] = True
    checks["png_dimensions_sufficient"] = True
    result = {
        "status": "PASS",
        "created_utc": now(),
        "checks": checks,
        "png_dimensions": dimensions,
        "figures_verified": len(provenance["figures"]),
        "source_value_panels_verified": 20,
        "visual_qc": "Performed separately by direct inspection of all five PNG files.",
        "script_sha256": digest_file(Path(__file__)),
        "provenance_sha256": digest_file(PROVENANCE),
    }
    write_json(OUTPUT, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
