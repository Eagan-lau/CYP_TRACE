"""Build manuscript Figures 1--5 from accepted, source-backed result files."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "figures"

INK = "#20252b"
MUTED = "#66717e"
GRID = "#d9dee5"
TEAL = "#087f8c"
BLUE = "#225ea8"
GOLD = "#c58a13"
RED = "#b44b4b"
PALE_TEAL = "#d7eef0"
PALE_GOLD = "#f4e8c9"
PALE_RED = "#f4dede"
LIGHT = "#edf0f3"

INPUT_PATHS = {
    "dataset_audit": ROOT / "dataset_02" / "dataset_audit.json",
    "annotation_audit": ROOT / "annotation_02" / "annotation_resolution_audit.json",
    "sequence_groups": ROOT / "protein_split_02" / "sequence_groups.json",
    "split_audit": ROOT / "multiaxis_split_01" / "split_audit.json",
    "protein_split_audit": ROOT / "protein_split_02" / "search_and_protein_split_audit.json",
    "router_summary": ROOT / "domain_router_01" / "method_summary.json",
    "paired_summary": ROOT / "paired_comparison_01" / "comparison_summary.json",
    "interaction_summary": ROOT / "interaction_paired_sequence_01" / "comparison_summary.json",
    "interaction_controls": ROOT / "interaction_control_comparison_01" / "sequence_summary.json",
    "esm_receipt": ROOT / "esm_global_01" / "receipt.json",
    "ordered_audit": ROOT / "ordered_site_features_01" / "audit.json",
    "clean_audit": ROOT / "clean_core_features_01" / "audit.json",
    "reverse_census": ROOT / "reverse_candidate_census_01" / "audit.json",
    "reverse_summary": ROOT / "reverse_baselines_01" / "method_summary.json",
    "reverse_acceptance": ROOT / "reverse_baselines_01" / "acceptance.json",
    "structure_audit": ROOT / "structure_mapping_audit_02" / "audit.json",
    "structure_paired": ROOT / "structure_paired_02" / "comparison_summary.json",
    "local_paired": ROOT / "local_paired_sequence_01" / "comparison_summary.json",
    "clean_pruning": ROOT / "clean_pruning_01" / "comparison_summary.json",
    "selective_analysis": ROOT / "human_selective_behavior_02" / "analysis.json",
}


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def row(rows, **match):
    hits = [entry for entry in rows if all(entry.get(key) == value for key, value in match.items())]
    if len(hits) != 1:
        raise RuntimeError(f"expected one row for {match}, found {len(hits)}")
    return hits[0]


def effect(entry, block="end_to_end"):
    value = entry[block]
    return {
        "difference": value["difference"],
        "ci": value["conditional_protein_bootstrap_95ci"],
        "left": value.get("left_mrr"),
        "right": value.get("right_mrr"),
        "query_panels": value.get("query_panels"),
        "protein_groups": value.get("protein_groups"),
    }


def collect_figure1():
    data = read(INPUT_PATHS["dataset_audit"])
    annotation = read(INPUT_PATHS["annotation_audit"])
    sources = ["P450Rdb", "SABIO_RK", "SwissProt", "UniProt_API"]
    metrics = [
        "unique_sequence_link_available",
        "linked_sequence_has_PF00067",
        "has_recorded_publication",
        "complete_component_set_parseable",
        "exact_component_set_edge_matches_initial_core",
    ]
    p450_edges = data["unique_edges_by_export"]["P450Rdb"]
    uniprot_edges = data["unique_edges_by_export"]["UniProt_API"]
    shared = row(data["cross_source_edge_overlap"], source_a="P450Rdb", source_b="UniProt_API")["shared_unique_edges"]
    return {
        "source_assertions": [
            {
                "source": source,
                "source_assertions": annotation["source_counts"][source]["source_assertions"],
                "admitted_assertions": data["source_counts"].get(source, {}).get("admitted_assertions", 0),
            }
            for source in sources
        ],
        "rate_sources": sources,
        "rate_metrics": metrics,
        "rates": [[annotation["source_assertion_rates"][source][metric]["fraction"] for metric in metrics] for source in sources],
        "edge_union": {
            "p450_unique": p450_edges - shared,
            "shared": shared,
            "uniprot_unique": uniprot_edges - shared,
            "union": data["unique_edges"],
            "swissprot_uniprot_same_lineage": True,
        },
        "core": {
            "sequences": data["labelled_sequences"],
            "edges": data["unique_edges"],
            "transformations": data["directed_main_transformations"],
            "single_pair": data["single_pair_transformations"],
            "multicomponent": data["multicomponent_transformations"],
        },
    }


def collect_figure2():
    groups = read(INPUT_PATHS["sequence_groups"])["0.4"]
    split = read(INPUT_PATHS["split_audit"])
    protein = read(INPUT_PATHS["protein_split_audit"])
    return {
        "protein_group_sizes": sorted((len(members) for members in groups.values()), reverse=True),
        "chemical_component_sizes": split["chemical_component_sizes_descending"],
        "identity_threshold": protein["primary_identity"],
        "minimum_bilateral_coverage": protein["minimum_bilateral_coverage"],
        "folds": [
            {
                "fold": item["fold"],
                "sequences": item["query_sequences"],
                "groups": item["query_groups"],
                "publication_purged_edges": item["publication_purged_training_edges"],
                "max_cross_identity": item["maximum_detected_bilateral_cross_partition_identity"],
            }
            for item in protein["fold_audits"]
        ],
        "tasks": [
            {
                "task": task,
                "blocks": values["blocks"],
                "query_panels": values["evaluable_query_panel_rows"],
                "test_edges": values["evaluable_test_edges"],
            }
            for task, values in split["task_denominators"].items()
        ],
        "outer_blocks": split["outer_blocks"],
        "inner_blocks": split["inner_blocks"],
        "protein_groups": split["protein_groups"],
        "chemical_components": split["chemical_components"],
    }


def collect_figure3():
    router = read(INPUT_PATHS["router_summary"])["summaries"]
    paired = read(INPUT_PATHS["paired_summary"])["comparisons"]
    interactions = read(INPUT_PATHS["interaction_summary"])["comparisons"]
    controls = read(INPUT_PATHS["interaction_controls"])["summaries"]
    ordered = read(INPUT_PATHS["ordered_audit"])
    clean = read(INPUT_PATHS["clean_audit"])
    esm = read(INPUT_PATHS["esm_receipt"])
    cells = [
        ("new_protein_seen_chemistry_with_homologue", "New protein\nseen chemistry"),
        ("new_protein_unseen_chemistry_with_homologue__chemical_cold_accounting", "New protein\nunseen chemistry"),
        ("new_protein_unseen_chemistry_with_homologue__double_cold", "Double cold"),
        ("represented_protein_unseen_chemistry", "Represented protein\nunseen chemistry"),
    ]
    methods = [
        ("chemistry_prior", "Chemistry prior"),
        ("homology_chemical_transport", "Homology transport"),
        ("mmseqs_top1", "MMseqs2 top-1"),
        ("global_residual", "Global residual"),
        ("joint_global_interaction", "Joint interaction"),
        ("uniform_expectation", "Uniform expectation"),
    ]
    cell_values = []
    for cell, label in cells:
        by_method = {entry["method"]: entry for entry in router if entry["cell"] == cell}
        if "chemistry_prior" not in by_method and "availability_chemistry_prior" in by_method:
            by_method["chemistry_prior"] = by_method["availability_chemistry_prior"]
        cell_values.append({
            "cell": cell,
            "label": label,
            "query_panels": next(iter(by_method.values()))["query_panels"],
            "protein_groups": next(iter(by_method.values()))["protein_groups"],
            "mrr": {key: by_method[key]["protein_group_macro"]["rr"] if key in by_method else None for key, _ in methods},
        })

    global_effects = []
    interaction_effects = []
    for task, label in (("chemical_cold_all", "Chemical cold"), ("protein_cold_all", "Protein cold"), ("double_cold_all", "Double cold")):
        global_effects.append({
            "task": task,
            "label": label,
            **effect(row(paired, task=task, left_method="esm_conditional_learned", right_method="conditional_chemistry_prior")),
        })
        interaction_effects.append({
            "task": task,
            "label": label,
            **effect(row(interactions, task=task, left_method="joint_global_interaction", right_method="global_residual")),
        })

    control_effects = []
    for control, label in (("protein_position_permutation", "Protein positions permuted"), ("reaction_center_row_permutation", "Reaction-center rows permuted")):
        entry = row(controls, control=control, method="joint_global_interaction", task="chemical_cold_all")
        comparison = entry["paired_gain_comparison"]
        control_effects.append({
            "control": control,
            "label": label,
            "primary_gain": comparison["left_mrr"],
            "control_gain": comparison["right_mrr"],
            "difference": comparison["difference"],
            "ci": comparison["conditional_protein_bootstrap_95ci"],
        })
    return {
        "router_cells": cell_values,
        "router_methods": [{"key": key, "label": label} for key, label in methods],
        "global_effects": global_effects,
        "interaction_effects": interaction_effects,
        "availability": [
            {"label": "ESM global", "sequences": esm["sequences"]},
            {"label": "CLEAN supported", "sequences": clean["supported_sequences"]},
            {"label": "Ordered sequence sites", "sequences": ordered["counts"]["sequence_projected_available_sequences"]},
            {"label": "Ordered structure sites", "sequences": ordered["counts"]["structure_projected_available_sequences"]},
        ],
        "control_effects": control_effects,
    }


def collect_figure4():
    census = read(INPUT_PATHS["reverse_census"])
    summaries = read(INPUT_PATHS["reverse_summary"])["summaries"]
    acceptance = read(INPUT_PATHS["reverse_acceptance"])["by_cell"]
    cells = [
        ("seen_reaction__new_positive", "Seen reaction\nnew positive"),
        ("unseen_reaction__new_positive", "Unseen reaction\nnew positive"),
        ("unseen_reaction__represented_positive", "Unseen reaction\nrepresented positive"),
    ]
    method_map = {
        "seen_reaction__new_positive": {
            "MMseqs2": "mmseqs_weighted_positive", "Label prior": "training_label_count_prior", "Uniform": "uniform_expectation"},
        "unseen_reaction__new_positive": {
            "MMseqs2": "mmseqs_homology_chemical_transport", "Label prior": "training_label_count_prior", "Uniform": "uniform_expectation"},
        "unseen_reaction__represented_positive": {
            "MMseqs2": "mmseqs_homology_chemical_transport", "Label prior": "training_label_count_prior", "Uniform": "uniform_expectation"},
    }
    performance = []
    for cell, label in cells:
        values = {}
        coverage = {}
        for display, method in method_map[cell].items():
            entry = row(summaries, cell=cell, method=method)
            values[display] = entry["taxid_macro_mrr"]
            coverage[display] = entry["documented_positive_coverage"]
        performance.append({
            "cell": cell, "label": label, "mrr": values, "coverage": coverage,
            "taxids": acceptance[cell]["taxids"],
        })
    gate_names = [
        "at_least_20_taxids",
        "documented_positive_coverage_at_least_0_8",
        "positive_vs_training_label_count_prior",
        "positive_vs_uniform_expectation",
        "protein_permutation_exact_p_at_most_0_05",
    ]
    return {
        "census": {
            "raw_proteins": census["raw_proteins"],
            "admitted_candidates": census["admitted_candidate_sequences"],
            "admitted_core": census["admitted_core_sequences"],
            "admitted_noncore": census["admitted_noncore_exact_pf00067_sequences"],
            "candidate_taxids": census["candidate_taxids"],
            "query_panels": census["outer_query_panels"],
            "query_taxids": census["outer_query_taxids"],
            "query_reactions": census["outer_query_reactions"],
            "positives": census["outer_documented_positive_instances"],
            "candidate_min": census["candidate_count_min"],
            "candidate_max": census["candidate_count_max"],
        },
        "performance": performance,
        "gate_names": gate_names,
        "gates": [[acceptance[cell]["gates"][gate] for gate in gate_names] for cell, _ in cells],
        "accepted": [acceptance[cell]["accepted"] for cell, _ in cells],
    }


def collect_figure5():
    structure = read(INPUT_PATHS["structure_audit"])
    ordered = read(INPUT_PATHS["ordered_audit"])
    clean = read(INPUT_PATHS["clean_audit"])
    structure_pairs = read(INPUT_PATHS["structure_paired"])["comparisons"]
    local = read(INPUT_PATHS["local_paired"])["comparisons"]
    pruning = read(INPUT_PATHS["clean_pruning"])["comparisons"]
    selective = read(INPUT_PATHS["selective_analysis"])["strata"]["scaffold_and_lineage_eligible"]
    local_effects = []
    for task, label in (("chemical_cold_all", "Chemical cold"), ("protein_cold_all", "Protein cold"), ("double_cold_all", "Double cold")):
        local_effects.append({
            "label": f"Ordered sites · {label}",
            **effect(row(local, task=task, left_method="ordered_residual", right_method="global_residual")),
        })
    local_effects.append({
        "label": "Foldseek top-1 · seen chemistry",
        **effect(row(structure_pairs, task="protein_cold_seen", left_method="foldseek_top1", right_method="matched_mmseqs_top1")),
    })
    pruning_cells = [
        ("new_protein_unseen_chemistry_with_homologue__chemical_cold_accounting", "Chemical cold"),
        ("new_protein_unseen_chemistry_with_homologue__protein_cold", "Protein cold"),
        ("new_protein_unseen_chemistry_with_homologue__double_cold", "Double cold"),
    ]
    pruning_effects = []
    for cell, label in pruning_cells:
        entry = row(pruning, cell=cell)
        pruning_effects.append({
            "label": label,
            "all": {"difference": entry["all"]["difference"], "ci": entry["all"]["conditional_protein_bootstrap_95ci"]},
            "exposed": {"difference": entry["clean_exact_exposed"]["difference"], "ci": entry["clean_exact_exposed"]["conditional_protein_bootstrap_95ci"]},
            "unexposed": {"difference": entry["clean_exact_unexposed"]["difference"], "ci": entry["clean_exact_unexposed"]["conditional_protein_bootstrap_95ci"]},
        })
    curve = selective["label_blind_matched_coverage_curve"]
    return {
        "availability": [
            {"label": "Core sequences", "count": structure["core_sequences"]},
            {"label": "CLEAN supported", "count": clean["supported_sequences"]},
            {"label": "Sequence-site vector", "count": ordered["counts"]["sequence_projected_available_sequences"]},
            {"label": "Primary structure", "count": structure["primary_structure_sequences"]},
            {"label": "Structure-site vector", "count": ordered["counts"]["structure_projected_available_sequences"]},
        ],
        "mapped_foldseek_chains": structure["mapped_foldseek_chains"],
        "ordered_positions": ordered["ordered_positions"],
        "local_effects": local_effects,
        "clean_counts": {"exposed": clean["clean_exact_exposed"], "unexposed": clean["clean_exact_unexposed"]},
        "pruning_effects": pruning_effects,
        "selective_curve": [
            {
                "method": item["method"],
                "target_coverage": item["target_coverage"],
                "coverage": item["macro_coverage"],
                "precision": item["macro_precision"],
                "positive_recall": item["macro_positive_recall"],
            }
            for item in curve
        ],
        "selective_rows": selective["rows"],
        "selective_compounds": selective["compounds"],
        "selective_scaffolds": selective["scaffolds"],
        "matched_coverage_comparison": selective["matched_coverage_comparison"],
    }


def set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.3,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "text.color": INK,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.8,
    })


def style_axis(axis, grid="y"):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(MUTED)
    axis.spines["bottom"].set_color(MUTED)
    axis.tick_params(colors=INK, labelsize=8.4)
    if grid:
        getattr(axis, f"{grid}axis").grid(True, color=GRID, linewidth=0.65)
        axis.set_axisbelow(True)


def panel_title(axis, letter, title):
    axis.set_title(f"{letter}  {title}", loc="left", fontsize=11.5, fontweight="bold", pad=9)


def save_figure(figure, stem):
    paths = {}
    for extension in ("png", "svg", "pdf"):
        path = OUTPUT / f"{stem}.{extension}"
        figure.savefig(path, dpi=450 if extension == "png" else None, bbox_inches="tight", facecolor="white")
        paths[extension] = path
    plt.close(figure)
    return paths


def draw_figure1(data):
    figure, axes = plt.subplots(2, 2, figsize=(12.8, 9.0), constrained_layout=True)
    figure.suptitle("Evidence admission exposes the resolution limits of public CYP annotations", fontsize=15.5, fontweight="bold", x=0.02, ha="left")

    axis = axes[0, 0]
    rows = data["source_assertions"]
    x = np.arange(len(rows)); width = 0.36
    axis.bar(x - width / 2, [r["source_assertions"] for r in rows], width, color=LIGHT, edgecolor=MUTED, label="Source assertions")
    bars = axis.bar(x + width / 2, [r["admitted_assertions"] for r in rows], width, color=TEAL, edgecolor=INK, linewidth=0.7, label="Admitted to core")
    axis.set_xticks(x, [r["source"].replace("_", "\n") for r in rows])
    axis.set_ylabel("Assertions")
    axis.legend(frameon=False, fontsize=8.2, ncol=2, loc="upper right")
    for bar in bars:
        axis.text(bar.get_x() + bar.get_width()/2, bar.get_height()+85, f"{int(bar.get_height()):,}", ha="center", va="bottom", fontsize=7.8, color=TEAL)
    panel_title(axis, "A", "Source records versus core admission")
    style_axis(axis)

    axis = axes[0, 1]
    matrix = np.asarray(data["rates"])
    cmap = LinearSegmentedColormap.from_list("cyp_rates", ["#f7f8fa", PALE_TEAL, TEAL])
    image = axis.imshow(matrix, vmin=0, vmax=1, aspect="auto", cmap=cmap)
    metric_labels = ["Sequence\nresolved", "PF00067", "Publication", "Chemistry\nparseable", "Exact core\nedge"]
    axis.set_xticks(np.arange(len(metric_labels)), metric_labels)
    axis.set_yticks(np.arange(len(data["rate_sources"])), [x.replace("_", " ") for x in data["rate_sources"]])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            axis.text(j, i, f"{matrix[i,j]*100:.0f}%", ha="center", va="center", fontsize=8.2,
                      color="white" if matrix[i,j] > 0.62 else INK)
    colorbar = figure.colorbar(image, ax=axis, fraction=0.045, pad=0.025)
    colorbar.set_label("Assertion rate", fontsize=8.5)
    colorbar.ax.tick_params(labelsize=8)
    panel_title(axis, "B", "Resolution is field-specific")
    axis.tick_params(labelsize=8.2)

    axis = axes[1, 0]
    union = data["edge_union"]
    values = [union["p450_unique"], union["shared"], union["uniprot_unique"]]
    labels = ["P450Rdb only", "Shared exact edge", "UniProt lineage only"]
    colors = [GOLD, TEAL, BLUE]
    left = 0
    for value, label, color in zip(values, labels, colors):
        axis.barh([0], [value], left=left, height=0.5, color=color, edgecolor="white", label=label)
        axis.text(left + value/2, 0, f"{value:,}", ha="center", va="center", color="white" if color != GOLD else INK, fontweight="bold")
        left += value
    axis.set_xlim(0, union["union"] * 1.02)
    axis.set_yticks([])
    axis.set_xlabel("Unique sequence–reaction edges")
    axis.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.34), ncol=3, fontsize=8)
    axis.text(0.5, 0.78, f"Union = {union['union']:,}", transform=axis.transAxes, ha="center", fontsize=10, fontweight="bold")
    axis.text(0.5, -0.48, "Swiss-Prot and UniProt API exports are the same UniProt lineage, not independent evidence.", transform=axis.transAxes, ha="center", color=MUTED, fontsize=8.1)
    panel_title(axis, "C", "Exact-edge union and overlap")
    style_axis(axis, grid=None)

    axis = axes[1, 1]
    axis.set_axis_off()
    panel_title(axis, "D", "Unified evidence atlas")
    core = data["core"]
    cards = [
        ("Sequences", core["sequences"], TEAL),
        ("Sequence–reaction edges", core["edges"], BLUE),
        ("Directed transformations", core["transformations"], GOLD),
        ("Multicomponent", core["multicomponent"], RED),
    ]
    positions = [(0.03, 0.55), (0.52, 0.55), (0.03, 0.12), (0.52, 0.12)]
    for (label, value, color), (x0, y0) in zip(cards, positions):
        rect = plt.Rectangle((x0, y0), 0.44, 0.30, transform=axis.transAxes, facecolor="white", edgecolor=color, linewidth=1.6)
        axis.add_patch(rect)
        axis.text(x0+0.04, y0+0.19, f"{value:,}", transform=axis.transAxes, fontsize=20, fontweight="bold", color=color)
        axis.text(x0+0.04, y0+0.07, label, transform=axis.transAxes, fontsize=9, color=INK)
    axis.text(0.03, 0.00, f"Single-pair transformations: {core['single_pair']:,}; multicomponent chemistry is retained rather than flattened.", transform=axis.transAxes, fontsize=8.1, color=MUTED)
    figure.text(0.02, -0.008, "Admission requires resolved chemistry, exact CYP-family sequence support and recorded provenance; metadata links are not treated as proof of catalytic activity.", fontsize=8.3, color=MUTED, ha="left")
    return figure


def draw_figure2(data):
    figure, axes = plt.subplots(2, 2, figsize=(12.8, 9.0), constrained_layout=True)
    figure.suptitle("Protein, chemistry and publication dependence are separated before evaluation", fontsize=15.5, fontweight="bold", x=0.02, ha="left")

    axis = axes[0, 0]
    for values, label, color, marker in (
        (data["protein_group_sizes"], "Protein groups", BLUE, "o"),
        (data["chemical_component_sizes"], "Chemical components", GOLD, "s"),
    ):
        ranks = np.arange(1, len(values)+1)
        axis.plot(ranks, values, color=color, linewidth=1.7, marker=marker, markersize=3, markevery=max(1,len(values)//15), label=label)
    axis.set_yscale("log")
    axis.set_xlabel("Component rank")
    axis.set_ylabel("Members per component (log scale)")
    axis.legend(frameon=False, fontsize=8.4)
    panel_title(axis, "A", "Connected-component concentration")
    style_axis(axis)
    axis.text(0.98, 0.95, f"Largest protein group: {data['protein_group_sizes'][0]}\nLargest chemical component: {data['chemical_component_sizes'][0]}", transform=axis.transAxes, ha="right", va="top", fontsize=8.2, color=MUTED)

    axis = axes[0, 1]
    folds = data["folds"]; x = np.arange(len(folds))
    bars = axis.bar(x, [r["sequences"] for r in folds], width=0.62, color=PALE_TEAL, edgecolor=TEAL, linewidth=1.2, label="Sequences")
    axis.set_xticks(x, [f"Fold {r['fold']+1}" for r in folds])
    axis.set_ylabel("Held-out sequences")
    axis2 = axis.twinx()
    axis2.plot(x, [r["groups"] for r in folds], color=RED, marker="D", linewidth=1.5, markersize=5, label="Protein groups")
    axis2.set_ylabel("Held-out protein groups", color=RED)
    axis2.tick_params(colors=RED, labelsize=8.4)
    axis2.spines["top"].set_visible(False); axis2.spines["right"].set_color(RED)
    for bar, r in zip(bars, folds):
        axis.text(bar.get_x()+bar.get_width()/2, bar.get_height()+2, str(r["sequences"]), ha="center", fontsize=8)
    handles = [bars, axis2.lines[0]]
    axis.legend(handles, ["Sequences", "Protein groups"], frameon=False, fontsize=8.2, loc="upper right")
    panel_title(axis, "B", "Group-preserving protein folds")
    style_axis(axis)
    axis.text(0.02, 0.94, f"Identity < {data['identity_threshold']:.1f}; bilateral coverage ≥ {data['minimum_bilateral_coverage']:.1f}", transform=axis.transAxes, fontsize=8.1, color=MUTED, va="top")

    axis = axes[1, 0]
    tasks = data["tasks"]; y = np.arange(len(tasks))[::-1]
    labels = [r["task"].replace("_", " ").title() for r in tasks]
    bars = axis.barh(y, [r["query_panels"] for r in tasks], color=[BLUE, GOLD, TEAL], edgecolor=INK, linewidth=0.6)
    axis.set_yticks(y, labels)
    axis.set_xlabel("Evaluable query panels")
    for bar, r in zip(bars, tasks):
        axis.text(bar.get_width()+14, bar.get_y()+bar.get_height()/2, f"{r['query_panels']:,} panels · {r['blocks']} blocks", va="center", fontsize=8.3)
    axis.set_xlim(0, max(r["query_panels"] for r in tasks)*1.35)
    panel_title(axis, "C", "Three prospective generalization tasks")
    style_axis(axis, grid="x")
    axis.text(0.02, -0.25, f"Each task preserves all 2,304 edges once; {data['outer_blocks']} outer and {data['inner_blocks']} nested inner blocks were materialized.", transform=axis.transAxes, fontsize=8.1, color=MUTED)

    axis = axes[1, 1]
    purge = [r["publication_purged_edges"] for r in folds]
    bars = axis.bar(x, purge, color=PALE_GOLD, edgecolor=GOLD, linewidth=1.2)
    axis.set_xticks(x, [f"Fold {r['fold']+1}" for r in folds])
    axis.set_ylabel("Training edges removed")
    for bar, value in zip(bars, purge):
        axis.text(bar.get_x()+bar.get_width()/2, value+4, str(value), ha="center", fontsize=8.2)
    panel_title(axis, "D", "Publication-overlap purge")
    style_axis(axis)
    axis.text(0.02, 0.94, f"Total removed across folds: {sum(purge):,}\nMaximum detected cross-partition identity: {max(r['max_cross_identity'] for r in folds):.3f}", transform=axis.transAxes, fontsize=8.1, color=MUTED, va="top")
    figure.text(0.02, -0.008, "Component imbalance is retained and reported; folds are not treated as exchangeable biological replicates.", fontsize=8.3, color=MUTED, ha="left")
    return figure


def draw_forest(axis, rows, title, letter, scale=1.0, xlabel="MRR difference"):
    y = np.arange(len(rows))[::-1]
    for position, item in zip(y, rows):
        estimate = item["difference"] * scale; low, high = [x*scale for x in item["ci"]]
        axis.errorbar(estimate, position, xerr=np.asarray([[estimate-low],[high-estimate]]), fmt="o", color=item.get("color", BLUE), ecolor=item.get("color", BLUE), markeredgecolor=INK, markeredgewidth=0.5, capsize=3, linewidth=1.6, markersize=5.5)
        axis.text(high + (axis.get_xlim()[1]-axis.get_xlim()[0])*0.015 if axis.get_xlim() != (0.0,1.0) else high, position, f"{estimate:+.3f}", va="center", fontsize=8)
    axis.axvline(0, color=INK, linestyle="--", linewidth=1)
    axis.set_yticks(y, [item["label"] for item in rows])
    axis.set_xlabel(xlabel)
    panel_title(axis, letter, title)
    style_axis(axis, grid="x")


def draw_figure3(data):
    figure = plt.figure(figsize=(13.8, 10.5), constrained_layout=True)
    grid = figure.add_gridspec(3, 2, height_ratios=[1.02, 1.02, 0.46], width_ratios=[1.12, 0.88])
    figure.suptitle("No single expert dominates across CYP applicability domains", fontsize=15.5, fontweight="bold", x=0.02, ha="left")

    axis = figure.add_subplot(grid[0, :])
    methods = data["router_methods"]; cells = data["router_cells"]
    matrix = np.asarray([[cell["mrr"][method["key"]] if cell["mrr"][method["key"]] is not None else np.nan
                          for method in methods] for cell in cells], dtype=float)
    masked = np.ma.masked_invalid(matrix)
    cmap = LinearSegmentedColormap.from_list("router_mrr", ["#f6f8fa", PALE_TEAL, TEAL, BLUE]).copy()
    cmap.set_bad("#d8dde3")
    image = axis.imshow(masked, vmin=0, vmax=0.24, aspect="auto", cmap=cmap)
    axis.set_xticks(np.arange(len(methods)), [m["label"].replace(" ", "\n", 1) for m in methods], fontsize=8.3)
    axis.set_yticks(np.arange(len(cells)), [c["label"].replace("\n", " ") for c in cells], fontsize=8.5)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isnan(matrix[i, j]):
                text_value, color = "undefined", MUTED
            else:
                text_value = f"{matrix[i,j]:.3f}"
                color = "white" if matrix[i,j] > 0.095 else INK
            axis.text(j, i, text_value, ha="center", va="center", fontsize=8.1, color=color)
    panel_title(axis, "A", "Router cells expose method-specific domains")
    axis.tick_params(length=0)
    axis.set_yticklabels([f"{c['label'].replace(chr(10), ' ')}\nn={c['query_panels']}; groups={c['protein_groups']}" for c in cells])
    axis.text(0.995, 1.04, "Darker cells indicate higher MRR; grey cells are outside the expert's domain.", transform=axis.transAxes, ha="right", fontsize=8.1, color=MUTED)

    axis = figure.add_subplot(grid[1, 0])
    rows = []
    for item in data["global_effects"]:
        rows.append({**item, "label": f"Global · {item['label'].replace(' cold','')}", "color": BLUE})
    for item in data["interaction_effects"]:
        rows.append({**item, "label": f"Interaction · {item['label'].replace(' cold','')}", "color": RED})
    axis.set_xlim(-0.03, 0.018)
    draw_forest(axis, rows, "Conditional gains are task-dependent", "B")

    axis = figure.add_subplot(grid[1, 1])
    avail = data["availability"]; y = np.arange(len(avail))[::-1]
    bars = axis.barh(y, [r["sequences"] for r in avail], color=[BLUE, TEAL, GOLD, RED], edgecolor=INK, linewidth=0.6)
    axis.set_yticks(y, [r["label"] for r in avail])
    axis.set_xlabel("Core sequences with representation")
    axis.set_xlim(0, 660)
    for bar, r in zip(bars, avail):
        axis.text(bar.get_width()+10, bar.get_y()+bar.get_height()/2, f"{r['sequences']}/600", va="center", fontsize=8.3)
    panel_title(axis, "C", "Expert coverage is unequal")
    style_axis(axis, grid="x")

    axis = figure.add_subplot(grid[2, :])
    short_controls = {"Protein positions permuted": "Protein-position perm.", "Reaction-center rows permuted": "Reaction-center-row perm."}
    rows = [{"label": short_controls[item["label"]], "difference": item["difference"], "ci": item["ci"], "color": TEAL if item["ci"][0] > 0 else RED} for item in data["control_effects"]]
    axis.set_xlim(-0.00015, 0.0012)
    draw_forest(axis, rows, "Interaction gain versus matched permutations", "D", scale=1.0, xlabel="Primary minus permuted gain in MRR")
    axis.text(0.02, -0.27, "The protein-position control is separated from zero; the reaction-center-row control is not. This is suggestive, not independent biological validation.", transform=axis.transAxes, fontsize=8.1, color=MUTED)
    figure.text(0.02, -0.008, "Undefined expert scores remain undefined; they are not entered as numerical zeros into a global mixture.", fontsize=8.3, color=MUTED, ha="left")
    return figure


def draw_figure4(data):
    figure, axes = plt.subplots(2, 2, figsize=(12.9, 9.2), constrained_layout=True)
    figure.suptitle("Reverse reaction-to-enzyme retrieval fails a predeclared acceptance gate", fontsize=15.5, fontweight="bold", x=0.02, ha="left")
    census = data["census"]

    axis = axes[0, 0]
    labels = ["Raw proteins", "Admitted candidates", "Non-core candidates", "Core candidates"]
    values = [census["raw_proteins"], census["admitted_candidates"], census["admitted_noncore"], census["admitted_core"]]
    bars = axis.barh(np.arange(4)[::-1], values, color=[LIGHT, TEAL, GOLD, BLUE], edgecolor=INK, linewidth=0.6)
    axis.set_yticks(np.arange(4)[::-1], labels)
    axis.set_xlabel("Sequences")
    for bar, value in zip(bars, values):
        axis.text(value+40, bar.get_y()+bar.get_height()/2, f"{value:,}", va="center", fontsize=8.5)
    axis.set_xlim(0, max(values)*1.15)
    panel_title(axis, "A", "Source-observed candidate census")
    style_axis(axis, grid="x")
    axis.text(0.02, -0.23, f"{census['candidate_taxids']:,} candidate taxids; within-taxid candidate sets range from {census['candidate_min']} to {census['candidate_max']} proteins.", transform=axis.transAxes, fontsize=8.1, color=MUTED)

    axis = axes[0, 1]
    performance = data["performance"]; x=np.arange(3); width=0.23
    styles = [("MMseqs2", TEAL, ""), ("Label prior", GOLD, "//"), ("Uniform", LIGHT, "xx")]
    for j,(method,color,hatch) in enumerate(styles):
        axis.bar(x+(j-1)*width, [r["mrr"][method] for r in performance], width, color=color, edgecolor=INK, linewidth=0.55, hatch=hatch, label=method)
    axis.set_xticks(x, [r["label"] for r in performance])
    axis.set_ylabel("Taxid-macro MRR")
    axis.legend(frameon=False, fontsize=8.2, ncol=3, loc="upper left")
    panel_title(axis, "B", "Ranking does not consistently beat simple controls")
    style_axis(axis)

    axis = axes[1, 0]
    y=np.arange(3)[::-1]
    bars=axis.barh(y,[r["coverage"]["MMseqs2"] for r in performance],color=PALE_TEAL,edgecolor=TEAL,linewidth=1.2)
    axis.axvline(0.8,color=RED,linestyle="--",linewidth=1.2,label="Coverage gate (0.8)")
    axis.set_yticks(y,[r["label"] for r in performance])
    axis.set_xlim(0,1.04); axis.set_xlabel("Documented-positive coverage")
    for bar,r in zip(bars,performance):
        axis.text(min(bar.get_width()+0.025,0.94),bar.get_y()+bar.get_height()/2,f"{bar.get_width():.3f} · {r['taxids']} taxids",va="center",fontsize=8.2)
    axis.legend(frameon=False,fontsize=8.2,loc="lower right")
    panel_title(axis,"C","Coverage and cohort size constrain interpretation")
    style_axis(axis,grid="x")

    axis=axes[1,1]
    matrix=np.asarray(data["gates"],dtype=int)
    cmap=LinearSegmentedColormap.from_list("gate",[PALE_RED,PALE_TEAL],N=2)
    axis.imshow(matrix,vmin=0,vmax=1,aspect="auto",cmap=cmap)
    gate_labels=["≥20 taxids","Coverage ≥0.8","Beats label prior","Beats uniform","Permutation p≤0.05"]
    axis.set_xticks(np.arange(5),gate_labels,rotation=30,ha="right")
    axis.set_yticks(np.arange(3),[r["label"].replace("\n"," ") for r in performance])
    for i in range(3):
        for j in range(5): axis.text(j,i,"✓" if matrix[i,j] else "×",ha="center",va="center",fontsize=13,fontweight="bold",color=TEAL if matrix[i,j] else RED)
    for i,accepted in enumerate(data["accepted"]):
        axis.text(5.0,i,"ACCEPT" if accepted else "NOT ACCEPTED",va="center",fontsize=8.2,fontweight="bold",color=TEAL if accepted else RED)
    axis.set_xlim(-0.5,6.15)
    panel_title(axis,"D","No reverse-retrieval cell clears all gates")
    axis.tick_params(labelsize=8.1)
    figure.text(0.02,-0.008,f"Evaluation covers {census['query_panels']:,} reaction–taxid panels, {census['query_reactions']:,} reactions and {census['positives']:,} documented positives; candidate panels are source-observed, not exhaustive proteomes.",fontsize=8.3,color=MUTED,ha="left")
    return figure


def draw_figure5(data):
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 9.4), constrained_layout=True)
    figure.suptitle("Coverage-aware controls define where auxiliary experts can be used", fontsize=15.5, fontweight="bold", x=0.02, ha="left")

    axis=axes[0,0]
    avail=data["availability"]; y=np.arange(len(avail))[::-1]
    bars=axis.barh(y,[r["count"] for r in avail],color=[LIGHT,TEAL,GOLD,BLUE,RED],edgecolor=INK,linewidth=0.6)
    axis.set_yticks(y,[r["label"] for r in avail]); axis.set_xlim(0,660); axis.set_xlabel("Sequences")
    for bar,r in zip(bars,avail): axis.text(r["count"]+9,bar.get_y()+bar.get_height()/2,f"{r['count']}/600",va="center",fontsize=8.2)
    panel_title(axis,"A","Representation availability is the first constraint")
    style_axis(axis,grid="x")
    axis.text(0.02,-0.24,f"{data['mapped_foldseek_chains']:,} mapped Foldseek chains; ordered vectors preserve {data['ordered_positions']} anchor positions.",transform=axis.transAxes,fontsize=8.1,color=MUTED)

    axis=axes[0,1]
    axis.set_xlim(-0.05,0.025)
    rows=[{**r,"color":BLUE if r["label"].startswith("Foldseek") else GOLD} for r in data["local_effects"]]
    draw_forest(axis,rows,"Local and structural experts do not add robust gain","B")

    axis=axes[1,0]
    effects=data["pruning_effects"]; x=np.arange(3); width=0.23
    styles=[("all","All",TEAL,""),("exposed","CLEAN-exposed",GOLD,"//"),("unexposed","CLEAN-unexposed",LIGHT,"xx")]
    for j,(key,label,color,hatch) in enumerate(styles):
        vals=[r[key]["difference"] for r in effects]
        lows=[v-r[key]["ci"][0] for v,r in zip(vals,effects)]; highs=[r[key]["ci"][1]-v for v,r in zip(vals,effects)]
        axis.bar(x+(j-1)*width,vals,width,color=color,edgecolor=INK,linewidth=0.55,hatch=hatch,label=label)
        axis.errorbar(x+(j-1)*width,vals,yerr=np.asarray([lows,highs]),fmt="none",ecolor=INK,elinewidth=0.8,capsize=2)
    axis.axhline(0,color=INK,linestyle="--",linewidth=1)
    axis.set_xticks(x,[r["label"] for r in effects]); axis.set_ylabel("Pruned minus unpruned MRR")
    axis.legend(frameon=False,fontsize=7.8,ncol=3,loc="upper center")
    panel_title(axis,"C","CLEAN pruning gain is exposure-dependent")
    style_axis(axis)
    axis.text(0.02,-0.24,f"Exact CLEAN exposure: {data['clean_counts']['exposed']} sequences; unexposed: {data['clean_counts']['unexposed']}.",transform=axis.transAxes,fontsize=8.1,color=MUTED)

    axis=axes[1,1]
    curve=data["selective_curve"]
    styles={"pooled_chemical_knn":("Pooled chemistry",MUTED,"o","--"),"isoform_specific_knn":("Isoform-specific",TEAL,"s","-")}
    for method,(label,color,marker,linestyle) in styles.items():
        rows=sorted([r for r in curve if r["method"]==method],key=lambda r:r["coverage"])
        axis.plot([r["coverage"] for r in rows],[r["precision"] for r in rows],label=label,color=color,marker=marker,linestyle=linestyle,linewidth=1.8,markersize=5,markerfacecolor="white" if method.startswith("pooled") else color,markeredgecolor=INK)
    axis.set_xlabel("Matched external coverage"); axis.set_ylabel("Macro precision")
    axis.set_xlim(0.05,1.03); axis.set_ylim(0.32,0.71)
    axis.legend(frameon=False,fontsize=8.2,loc="upper right")
    panel_title(axis,"D","Selective release enriches fixed-human predictions")
    style_axis(axis)
    strict10=data["matched_coverage_comparison"]["0.1"]
    axis.text(0.04,0.07,f"At ~10% coverage: Δ precision = {strict10['isoform_specific_minus_pooled_macro_precision']:+.3f}\n95% scaffold bootstrap {strict10['scaffold_bootstrap_95ci'][0]:+.3f} to {strict10['scaffold_bootstrap_95ci'][1]:+.3f}",transform=axis.transAxes,fontsize=8.1,color=TEAL)
    figure.text(0.02,-0.008,f"Selective analysis: {data['selective_rows']:,} fixed-human labels, {data['selective_compounds']:,} compounds and {data['selective_scaffolds']:,} scaffolds; development thresholds are not claimed as portable deployment cutoffs.",fontsize=8.3,color=MUTED,ha="left")
    return figure


def main():
    set_style()
    OUTPUT.mkdir(exist_ok=True)
    figures = {
        "figure1_evidence_atlas": (collect_figure1(), draw_figure1),
        "figure2_split_design": (collect_figure2(), draw_figure2),
        "figure3_domain_methods": (collect_figure3(), draw_figure3),
        "figure4_reverse_retrieval": (collect_figure4(), draw_figure4),
        "figure5_coverage_controls": (collect_figure5(), draw_figure5),
    }
    outputs = {}
    for stem, (data, draw) in figures.items():
        source_path = OUTPUT / f"{stem}_source_data.json"
        write_json(source_path, data)
        paths = save_figure(draw(data), stem)
        outputs[stem] = {
            "source_data": {"path": str(source_path.relative_to(ROOT)), "sha256": digest_file(source_path)},
            "files": {ext: {"path": str(path.relative_to(ROOT)), "sha256": digest_file(path)} for ext, path in paths.items()},
        }
    provenance = {
        "created_utc": now(),
        "status": "BUILT_VISUAL_QC_PENDING",
        "script_sha256": digest_file(Path(__file__)),
        "input_hashes": {key: {"path": str(path.relative_to(ROOT)), "sha256": digest_file(path)} for key, path in INPUT_PATHS.items()},
        "figures": outputs,
        "scope": "Manuscript Figures 1--5; Figure 6 is frozen and built separately.",
    }
    write_json(OUTPUT / "figures1_5_provenance.json", provenance)
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
