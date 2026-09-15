"""Bind the interaction prose to independently checked numerical artefacts."""
from __future__ import annotations

import json
from pathlib import Path

from run_raw import digest_file, now, write_json


def comparison(path, task, left, right):
    payload = json.loads(path.read_text())
    return next(item for item in payload["comparisons"] if (item["task"], item["left_method"], item["right_method"]) == (task, left, right))


def close(actual, expected, tolerance=5e-9):
    return abs(float(actual) - float(expected)) <= tolerance


def run(root: Path):
    failures = []; checks = []
    required_pass = [
        "reaction_center_environment_validation_01/audit.json",
        "interaction_transform_validation_01/audit.json",
        "interaction_validation_sequence_01/INDEPENDENT_QC.json",
        "interaction_validation_structure_01/INDEPENDENT_QC.json",
        "interaction_paired_validation_sequence_01/INDEPENDENT_QC.json",
        "interaction_paired_validation_structure_01/INDEPENDENT_QC.json",
        "interaction_control_validation_01/audit.json",
        "interaction_control_comparison_01/sequence_validation.json",
        "interaction_control_comparison_01/structure_validation.json",
        "interaction_information_controls_01/sequence_validation.json",
        "interaction_information_controls_01/structure_validation.json",
        "domain_router_validation_01/INDEPENDENT_QC.json",
        "clean_ec_bridge_validation_01/INDEPENDENT_QC.json",
        "clean_core_features_validation_01/INDEPENDENT_QC.json",
        "clean_pruning_validation_01/INDEPENDENT_QC.json",
        "reverse_candidate_census_validation_01/INDEPENDENT_QC.json",
        "reverse_alignment_validation_01/INDEPENDENT_QC.json",
        "reverse_baseline_validation_01/INDEPENDENT_QC.json",
        "human_substrate_model_validation_01/INDEPENDENT_QC.json",
        "external_human_cyp_validation_01/validation.json",
    ]
    for name in required_pass:
        value = json.loads((root / name).read_text())["status"]
        checks.append({"claim": name + " status", "actual": value, "expected": "PASS"})
        if value != "PASS": failures.append(name)
    feature = json.loads((root / "reaction_center_environments_01/audit.json").read_text())
    transform = json.loads((root / "interaction_transform_01/audit.json").read_text())
    count_claims = [
        ("reaction count", feature["reaction_count"], 1341),
        ("feature cells", feature["feature_shape"][0] * feature["feature_shape"][1], 1377207),
        ("feature dimension", feature["feature_shape"][1], 1027),
        ("mapping warnings", feature["mapping_warning_count"], 21),
        ("reaction numerical rank", transform["numerical_rank"], 886),
    ]
    for label, actual, expected in count_claims:
        checks.append({"claim": label, "actual": actual, "expected": expected})
        if actual != expected: failures.append(label)

    sequence = root / "interaction_paired_sequence_01/comparison_summary.json"
    values = []
    item = comparison(sequence, "chemical_cold_all", "joint_global_interaction", "global_residual")["end_to_end"]
    values.extend([("sequence chemical-cold joint MRR", item["left_mrr"], 0.005325150342707457),
                   ("sequence chemical-cold global MRR", item["right_mrr"], 0.004880797461659627),
                   ("sequence chemical-cold increment", item["difference"], 0.00044435288104783003),
                   ("sequence chemical-cold lower", item["conditional_protein_bootstrap_95ci"][0], 0.00008264486313979849),
                   ("sequence chemical-cold upper", item["conditional_protein_bootstrap_95ci"][1], 0.000853494500037665)])
    protein = comparison(sequence, "protein_cold_all", "joint_global_interaction", "global_residual")["end_to_end"]
    unseen = comparison(sequence, "protein_cold_unseen", "joint_global_interaction", "global_residual")["end_to_end"]
    double = comparison(sequence, "double_cold_all", "joint_global_interaction", "global_residual")["end_to_end"]
    values.extend([("sequence protein-cold increment", protein["difference"], -0.007512073997929498),
                   ("sequence protein-cold-unseen increment", unseen["difference"], -0.007890386678904224),
                   ("sequence double-cold increment", double["difference"], -0.00038979645761227596)])

    for control, expected_gain in [("protein_position_permutation", -0.00010228896934159887), ("reaction_center_row_permutation", 0.0000534007891484515)]:
        workspace = "interaction_control_protein_position_permutation_01" if control.startswith("protein") else "interaction_control_reaction_center_row_permutation_01"
        control_item = comparison(root / workspace / "interaction_paired_sequence_01/comparison_summary.json", "chemical_cold_all", "joint_global_interaction", "global_residual")["end_to_end"]
        values.append((control + " chemical-cold gain", control_item["difference"], expected_gain))
    gain = json.loads((root / "interaction_control_comparison_01/sequence_summary.json").read_text())
    for control, expected, interval in [
        ("protein_position_permutation", 0.0005466418503894292, (0.0000775795556170144, 0.0011090583481439598)),
        ("reaction_center_row_permutation", 0.0003909520918993789, (-0.000033934630711868065, 0.0009022398328807044)),
    ]:
        found = next(item["paired_gain_comparison"] for item in gain["summaries"] if item["control"] == control and item["method"] == "joint_global_interaction" and item["task"] == "chemical_cold_all")
        values.extend([(control + " gain difference", found["difference"], expected),
                       (control + " gain difference lower", found["conditional_protein_bootstrap_95ci"][0], interval[0]),
                       (control + " gain difference upper", found["conditional_protein_bootstrap_95ci"][1], interval[1])])
    info = json.loads((root / "interaction_information_controls_01/sequence_summary.json").read_text())
    for right, expected in [("missing_residual", 0.001179981392224704), ("composition_residual", 0.0011143395279465036), ("lineage_residual", 0.0018581346584190502)]:
        found = next(item["end_to_end"] for item in info["comparisons"] if item["task"] == "chemical_cold_all" and item["left_method"] == "interaction_residual" and item["right_method"] == right)
        values.append(("interaction versus " + right, found["difference"], expected))
    router = json.loads((root / "domain_router_01/comparison_summary.json").read_text())
    router_index = {(item["cell"], item["left_method"], item["right_method"]): item["paired"]
                    for item in router["comparisons"]}
    routed_interaction = router_index[("represented_protein_unseen_chemistry", "joint_global_interaction", "global_residual")]
    routed_mmseqs = router_index[("new_protein_seen_chemistry_with_homologue", "mmseqs_weighted", "uniform_expectation")]
    routed_double = router_index[("new_protein_unseen_chemistry_with_homologue__double_cold", "homology_chemical_transport", "uniform_expectation")]
    values.extend([
        ("routed interaction MRR", routed_interaction["left_mrr"], 0.22958281584292115),
        ("routed global MRR", routed_interaction["right_mrr"], 0.21606525430943865),
        ("routed interaction difference", routed_interaction["difference"], 0.013517561533482526),
        ("routed interaction lower", routed_interaction["conditional_protein_bootstrap_95ci"][0], -0.018266361781993187),
        ("routed interaction upper", routed_interaction["conditional_protein_bootstrap_95ci"][1], 0.04530148484895824),
        ("routed MMseqs MRR", routed_mmseqs["left_mrr"], 0.07226223263284534),
        ("routed MMseqs random MRR", routed_mmseqs["right_mrr"], 0.01443473133722585),
        ("routed double transport minus random", routed_double["difference"], -0.003935142672678189),
    ])
    router_counts = [
        ("router metric rows", json.loads((root / "domain_router_01/method_summary.json").read_text())["metric_rows"], 6513),
        ("router decision rows", json.loads((root / "domain_router_01/routing_coverage.json").read_text())["target_subpanels"], 2217),
        ("router paired rows", router["paired_query_rows"], 4355),
        ("router control rows", router["control_query_rows"], 132),
        ("routed interaction panels", routed_interaction["query_panels"], 66),
        ("routed interaction groups", routed_interaction["protein_groups"], 4),
        ("routed MMseqs panels", routed_mmseqs["query_panels"], 93),
        ("routed MMseqs groups", routed_mmseqs["protein_groups"], 25),
    ]
    for label, actual, expected in router_counts:
        checks.append({"claim": label, "actual": actual, "expected": expected})
        if actual != expected: failures.append(label)
    acceptance = json.loads((root / "domain_router_01/acceptance.json").read_text())
    if acceptance["interaction_route"]["accepted"] or acceptance["mmseqs_seen_chemistry_route"]["accepted"] or any(
            item["accepted"] for item in acceptance["homology_transport_unseen_chemistry"].values()):
        failures.append("router_zero_acceptance")
    checks.append({"claim": "router predictive routes accepted", "actual": 0, "expected": 0})

    clean_bridge = json.loads((root / "clean_ec_bridge_01/audit.json").read_text())
    clean_features = json.loads((root / "clean_core_features_01/audit.json").read_text())
    clean_qc = json.loads((root / "clean_pruning_validation_01/INDEPENDENT_QC.json").read_text())
    clean_counts = [
        ("CLEAN mapped reactions", clean_bridge["reactions_with_direct_clean_exact_ec"], 769),
        ("CLEAN leave-query-out edges", clean_bridge["edges_with_leave_query_sequence_out_candidate_ec"], 981),
        ("CLEAN supported sequences", clean_features["supported_sequences"], 582),
        ("CLEAN unsupported overlength sequences", clean_features["unsupported_over_1022"], 18),
        ("CLEAN exact exposed sequences", clean_features["clean_exact_exposed"], 328),
        ("CLEAN exact unexposed sequences", clean_features["clean_exact_unexposed"], 272),
        ("CLEAN inner panels checked", clean_qc["inner_query_panels_checked"], 15961),
        ("CLEAN nested selections checked", clean_qc["nested_selections_checked"], 45),
        ("CLEAN outer panels checked", clean_qc["outer_query_panels_reconstructed"], 2065),
    ]
    for label, actual, expected in clean_counts:
        checks.append({"claim": label, "actual": actual, "expected": expected})
        if actual != expected: failures.append(label)
    clean_comparisons = {item["cell"]: item for item in json.loads((root / "clean_pruning_01/comparison_summary.json").read_text())["comparisons"]}
    clean_expected = {
        "new_protein_seen_chemistry_with_homologue": (-0.0025774750306021056, -0.05266251675359727, 81),
        "new_protein_unseen_chemistry_with_homologue__chemical_cold_accounting": (0.11294958728340111, -0.0311524721289497, 566),
        "new_protein_unseen_chemistry_with_homologue__protein_cold": (0.06913023331981107, -0.025496147610841863, 551),
        "new_protein_unseen_chemistry_with_homologue__double_cold": (0.12397539930878888, -0.0034009168708885504, 675),
        "represented_protein_unseen_chemistry": (-0.08823018769877974, -0.3297627397493377, 255),
    }
    for cell, (all_difference, unexposed_difference, covered_positives) in clean_expected.items():
        found = clean_comparisons[cell]
        values.extend([
            ("CLEAN all difference " + cell, found["all"]["difference"], all_difference),
            ("CLEAN unexposed difference " + cell, found["clean_exact_unexposed"]["difference"], unexposed_difference),
        ])
        checks.append({"claim": "CLEAN covered positives " + cell, "actual": found["selected_covered_positive_instances"], "expected": covered_positives})
        if found["selected_covered_positive_instances"] != covered_positives: failures.append("CLEAN covered positives " + cell)
    clean_acceptance = json.loads((root / "clean_pruning_01/acceptance.json").read_text())["by_cell"]
    accepted_clean_cells = sum(item["clean_pruning_useful"] for item in clean_acceptance.values())
    checks.append({"claim": "CLEAN pruning cells accepted", "actual": accepted_clean_cells, "expected": 0})
    if accepted_clean_cells != 0: failures.append("CLEAN zero acceptance")

    reverse_census = json.loads((root / "reverse_candidate_census_01/audit.json").read_text())
    reverse_qc = json.loads((root / "reverse_baseline_validation_01/INDEPENDENT_QC.json").read_text())
    reverse_counts = [
        ("reverse candidates", reverse_census["admitted_candidate_sequences"], 1803),
        ("reverse candidate taxids", reverse_census["candidate_taxids"], 463),
        ("reverse query panels", reverse_census["outer_query_panels"], 4880),
        ("reverse query taxids", reverse_census["outer_query_taxids"], 133),
        ("reverse query reactions", reverse_census["outer_query_reactions"], 1167),
        ("reverse positive instances", reverse_census["outer_documented_positive_instances"], 6033),
        ("reverse reconstructed metric rows", reverse_qc["metric_rows_reconstructed"], 20000),
        ("reverse reconstructed summaries", reverse_qc["method_summaries_reconstructed"], 14),
        ("reverse reconstructed comparisons", reverse_qc["comparisons_reconstructed"], 16),
    ]
    for label, actual, expected in reverse_counts:
        checks.append({"claim": label, "actual": actual, "expected": expected})
        if actual != expected: failures.append(label)
    reverse_methods = {(item["cell"], item["method"]): item for item in json.loads(
        (root / "reverse_baselines_01/method_summary.json").read_text())["summaries"]}
    reverse_comparisons = {(item["cell"], item["left_method"], item["right_method"]): item for item in json.loads(
        (root / "reverse_baselines_01/comparison_summary.json").read_text())["comparisons"]}
    reverse_controls = {(item["cell"], item["method"]): item for item in json.loads(
        (root / "reverse_baselines_01/permutation_controls.json").read_text())["controls"]}
    unseen_key = ("unseen_reaction__new_positive", "mmseqs_homology_chemical_transport")
    seen_key = ("seen_reaction__new_positive", "mmseqs_weighted_positive")
    unseen_prior = reverse_comparisons[(*unseen_key, "training_label_count_prior")]
    unseen_uniform = reverse_comparisons[(*unseen_key, "uniform_expectation")]
    values.extend([
        ("reverse unseen-new MRR", reverse_methods[unseen_key]["taxid_macro_mrr"], 0.5477818285704138),
        ("reverse unseen-new prior MRR", unseen_prior["right_mrr"], 0.48351530565643513),
        ("reverse unseen-new prior difference", unseen_prior["difference"], 0.06426652291397866),
        ("reverse unseen-new prior lower", unseen_prior["taxid_bootstrap_95ci"][0], 0.024284642654211895),
        ("reverse unseen-new prior upper", unseen_prior["taxid_bootstrap_95ci"][1], 0.10550334518654012),
        ("reverse unseen-new uniform MRR", unseen_uniform["right_mrr"], 0.5321822206398328),
        ("reverse unseen-new uniform difference", unseen_uniform["difference"], 0.015599607930580883),
        ("reverse unseen-new permutation p", reverse_controls[unseen_key]["exact_upper_tail_p"], 0.03),
        ("reverse seen-new MRR", reverse_methods[seen_key]["taxid_macro_mrr"], 0.19887325187626162),
        ("reverse seen-new coverage", reverse_methods[seen_key]["documented_positive_coverage"], 0.5660377358490566),
    ])
    accepted_reverse_cells = sum(item["accepted"] for item in json.loads(
        (root / "reverse_baselines_01/acceptance.json").read_text())["by_cell"].values())
    checks.append({"claim": "reverse cells accepted", "actual": accepted_reverse_cells, "expected": 0})
    if accepted_reverse_cells != 0: failures.append("reverse zero acceptance")

    human_audit = json.loads((root / "human_substrate_models_01/audit.json").read_text())
    human_qc = json.loads((root / "human_substrate_model_validation_01/INDEPENDENT_QC.json").read_text())
    human_counts = [
        ("human normalized rows", human_audit["normalized_rows"], 14955),
        ("human compounds", human_audit["compounds"], 1768),
        ("human scaffold components", human_audit["scaffold_components"], 1144),
        ("human outer predictions", human_audit["outer_predictions"], 14955),
        ("human fingerprint bits", human_qc["fingerprint_bits_reconstructed"], 3620864),
        ("human similarity cells", human_qc["similarity_cells_reconstructed"], 3125824),
        ("human outer scores", human_qc["outer_scores_reconstructed"], 44865),
        ("human bootstraps", human_qc["bootstrap_replicates_reconstructed"], 5000),
        ("human permutations", human_qc["score_permutations_reconstructed"], 99),
    ]
    for label, actual, expected in human_counts:
        checks.append({"claim": label, "actual": actual, "expected": expected})
        if actual != expected: failures.append(label)
    human_summary = {(item["method"], item["scope"]): item for item in json.loads(
        (root / "human_substrate_models_01/method_summary.json").read_text())["summaries"]}
    human_comparison = {item["right_method"]: item for item in json.loads(
        (root / "human_substrate_models_01/comparison_summary.json").read_text())["comparisons"]}
    human_permutation = json.loads((root / "human_substrate_models_01/score_permutation_control.json").read_text())
    values.extend([
        ("human pooled AP", human_summary[("pooled_chemical_knn", "isoform_macro")]["macro_average_precision"], 0.35417246695651744),
        ("human pooled AUC", human_summary[("pooled_chemical_knn", "isoform_macro")]["macro_roc_auc"], 0.8183502505883281),
        ("human isoform AP", human_summary[("isoform_specific_knn", "isoform_macro")]["macro_average_precision"], 0.46662820213157885),
        ("human isoform AUC", human_summary[("isoform_specific_knn", "isoform_macro")]["macro_roc_auc"], 0.843906445819373),
        ("human shrinkage AP", human_summary[("data_driven_shrinkage", "isoform_macro")]["macro_average_precision"], 0.4506021993156326),
        ("human shrinkage AUC", human_summary[("data_driven_shrinkage", "isoform_macro")]["macro_roc_auc"], 0.8438363110383285),
        ("human shrinkage vs pooled AP", human_comparison["pooled_chemical_knn"]["average_precision"]["difference"], 0.09642973235911517),
        ("human shrinkage vs pooled lower", human_comparison["pooled_chemical_knn"]["average_precision"]["scaffold_bootstrap_95ci"][0], 0.0718777712025986),
        ("human shrinkage vs pooled upper", human_comparison["pooled_chemical_knn"]["average_precision"]["scaffold_bootstrap_95ci"][1], 0.12435950707141386),
        ("human shrinkage vs isoform AP", human_comparison["isoform_specific_knn"]["average_precision"]["difference"], -0.01602600281594624),
        ("human permutation p", human_permutation["plus_one_upper_tail_monte_carlo_p"], 0.01),
    ])
    human_acceptance = json.loads((root / "human_substrate_models_01/acceptance.json").read_text())
    if not human_acceptance["isoform_condition_signal"]["accepted"] or human_acceptance["shrinkage_improves_both_components"]["accepted"]:
        failures.append("human_acceptance_boundary")
    checks.append({"claim": "human accepted claims", "actual": [True, False], "expected": [True, False]})

    external_audit = json.loads((root / "external_human_cyp_01/overlap_audit.json").read_text())
    external_qc = json.loads((root / "external_human_cyp_validation_01/validation.json").read_text())
    external_counts = [
        ("external raw rows", external_audit["counts"]["raw_rows"], 15005),
        ("external unresolved rows", external_audit["counts"]["unresolved_structure_rows"], 366),
        ("external normalized groups", external_audit["normalized_groups"], 14526),
        ("external exact pair overlaps", external_audit["counts"]["exact_pair_overlap_groups"], 9992),
        ("external exact pair agreements", external_audit["counts"]["exact_pair_agree"], 9827),
        ("external exact pair conflicts", external_audit["counts"]["exact_pair_conflict"], 158),
        ("external lineage rows", external_audit["strata"]["lineage_eligible"]["rows"], 4301),
        ("external lineage compounds", external_audit["strata"]["lineage_eligible"]["unique_compounds"], 1365),
        ("external lineage scaffolds", external_audit["strata"]["lineage_eligible"]["scaffold_groups"], 963),
        ("external strict rows", external_audit["strata"]["scaffold_and_lineage_eligible"]["rows"], 3035),
        ("external strict compounds", external_audit["strata"]["scaffold_and_lineage_eligible"]["unique_compounds"], 1001),
        ("external strict scaffolds", external_audit["strata"]["scaffold_and_lineage_eligible"]["scaffold_groups"], 844),
        ("external predictions verified", external_qc["predictions_verified"], 14519),
        ("external train fingerprint bits", external_qc["train_fingerprint_bits_verified"], 3620864),
        ("external query fingerprint bits", external_qc["external_fingerprint_bits_verified"], 6477824),
        ("external similarity cells", external_qc["similarity_cells_verified"], 5592184),
    ]
    for label, actual, expected in external_counts:
        checks.append({"claim": label, "actual": actual, "expected": expected})
        if actual != expected: failures.append(label)
    external_summary = json.loads((root / "external_human_cyp_models_01/method_summary.json").read_text())["strata"]
    external_comparison = json.loads((root / "external_human_cyp_models_01/comparison_summary.json").read_text())["strata"]
    external_permutation = json.loads((root / "external_human_cyp_models_01/score_permutation_control.json").read_text())["strata"]
    external_acceptance = json.loads((root / "external_human_cyp_models_01/acceptance.json").read_text())
    lineage = "lineage_eligible"; strict = "scaffold_and_lineage_eligible"
    values.extend([
        ("external lineage pooled AP", external_summary[lineage]["pooled_chemical_knn"]["macro_average_precision"], 0.401715096036451),
        ("external lineage specific AP", external_summary[lineage]["isoform_specific_knn"]["macro_average_precision"], 0.5392172106554336),
        ("external lineage AP difference", external_comparison[lineage]["average_precision"]["difference"], 0.1375021146189826),
        ("external lineage AP lower", external_comparison[lineage]["average_precision"]["scaffold_bootstrap_95ci"][0], 0.11463932591098845),
        ("external lineage AP upper", external_comparison[lineage]["average_precision"]["scaffold_bootstrap_95ci"][1], 0.1571825647698241),
        ("external lineage pooled AUC", external_summary[lineage]["pooled_chemical_knn"]["macro_roc_auc"], 0.5428034204927845),
        ("external lineage specific AUC", external_summary[lineage]["isoform_specific_knn"]["macro_roc_auc"], 0.6725381205666383),
        ("external lineage permutation p", external_permutation[lineage]["plus_one_upper_tail_monte_carlo_p"], 0.01),
        ("external strict pooled AP", external_summary[strict]["pooled_chemical_knn"]["macro_average_precision"], 0.41097145717927),
        ("external strict specific AP", external_summary[strict]["isoform_specific_knn"]["macro_average_precision"], 0.5428801073616204),
        ("external strict AP difference", external_comparison[strict]["average_precision"]["difference"], 0.13190865018235037),
        ("external strict AP lower", external_comparison[strict]["average_precision"]["scaffold_bootstrap_95ci"][0], 0.10806726348022931),
        ("external strict AP upper", external_comparison[strict]["average_precision"]["scaffold_bootstrap_95ci"][1], 0.15714035854241562),
        ("external strict pooled AUC", external_summary[strict]["pooled_chemical_knn"]["macro_roc_auc"], 0.550724392545219),
        ("external strict specific AUC", external_summary[strict]["isoform_specific_knn"]["macro_roc_auc"], 0.6767614650022752),
        ("external strict permutation p", external_permutation[strict]["plus_one_upper_tail_monte_carlo_p"], 0.01),
    ])
    if not external_acceptance["independent_external_transfer"]["accepted"] or not external_acceptance["scaffold_external_transfer"]["accepted"]:
        failures.append("external_human_acceptance_boundary")
    checks.append({"claim": "external human accepted claims", "actual": [True, True], "expected": [True, True]})
    for label, actual, expected in values:
        checks.append({"claim": label, "actual": actual, "expected": expected})
        if not close(actual, expected): failures.append(label)

    manuscript = (root / "manuscript/WORKING_MANUSCRIPT.md").read_text(encoding="utf-8")
    phrases = ["1,377,207", "0.00533", "0.00044", "0.00055", "0.00111", "0.00186", "−0.00789",
               "2,217", "6,513", "0.22958", "0.21606", "four protein groups", "0.07226", "no predictive route was promoted",
               "769 of", "981 of", "582 of 600", "328 had", "15,961", "2,065", "0.11295", "-0.03115", "does not qualify CLEAN pruning",
               "1,803", "4,880", "6,033", "0.54778", "0.48352", "0.53218", "0.19887", "none of the three cells passed all frozen gates",
               "14,955", "1,768", "0.35417", "0.46663", "0.45060", "0.09643", "0.01603", "0.01", "rejects cross-isoform score shrinkage",
               "15,005", "14,519", "9,992", "9,827", "158 conflicted", "4,301", "1,365", "963", "0.53922", "0.40172", "0.13750", "0.11464", "0.15718",
               "3,035", "1,001", "844", "0.54288", "0.41097", "0.13191", "0.10807", "0.15714", "externally reproduced"]
    for phrase in phrases:
        present = phrase in manuscript; checks.append({"claim": "manuscript token " + phrase, "actual": present, "expected": True})
        if not present: failures.append("manuscript:" + phrase)
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
              "checks": checks, "manuscript_sha256": digest_file(root / "manuscript/WORKING_MANUSCRIPT.md"),
              "scope": "candidate-interaction, coverage-router, nested-CLEAN, reverse-retrieval, human-substrate and external-transfer counts, selected point estimates, intervals and required QC gates"}
    write_json(root / "manuscript/INTERACTION_CLAIM_AUDIT.json", report)
    print(json.dumps({"status": report["status"], "checks": len(checks), "failures": failures}, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": run(Path(__file__).resolve().parent)
