"""Build the source-backed external-transfer figure for the manuscript."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
AUDIT_PATH = ROOT / "external_human_cyp_01" / "overlap_audit.json"
SUMMARY_PATH = ROOT / "external_human_cyp_models_01" / "method_summary.json"
COMPARISON_PATH = ROOT / "external_human_cyp_models_01" / "comparison_summary.json"
PERMUTATION_PATH = ROOT / "external_human_cyp_models_01" / "score_permutation_control.json"
VALIDATION_PATH = ROOT / "external_human_cyp_validation_01" / "validation.json"
OUTPUT = ROOT / "figures"

INK = "#20252b"
MUTED = "#67717e"
GRID = "#d9dee5"
BLUE = "#225ea8"
GOLD = "#c58a13"
PALE_BLUE = "#d9e7f5"


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def collect_data():
    audit = read(AUDIT_PATH)
    summaries = read(SUMMARY_PATH)["strata"]
    comparisons = read(COMPARISON_PATH)["strata"]
    permutations = read(PERMUTATION_PATH)["strata"]
    validation = read(VALIDATION_PATH)
    if audit["status"] != "PASS" or validation["status"] != "PASS":
        raise RuntimeError("figure inputs have not passed independent QC")
    stages = [
        {"label": "Unique evaluable pairs", "rows": audit["strata"]["all_parseable_unique"]["rows"]},
        {"label": "Exact-new compounds", "rows": audit["strata"]["exact_novel"]["rows"]},
        {"label": "+ eligible lineage", "rows": audit["strata"]["lineage_eligible"]["rows"]},
        {"label": "+ scaffold novel", "rows": audit["strata"]["scaffold_and_lineage_eligible"]["rows"]},
    ]
    strata = [
        ("lineage_eligible", "Exact-new +\nlineage-eligible"),
        ("scaffold_and_lineage_eligible", "+ scaffold-novel"),
    ]
    aggregate = []
    for key, label in strata:
        aggregate.append({
            "stratum": key, "label": label,
            "rows": audit["strata"][key]["rows"],
            "compounds": audit["strata"][key]["unique_compounds"],
            "scaffolds": audit["strata"][key]["scaffold_groups"],
            "pooled_ap": summaries[key]["pooled_chemical_knn"]["macro_average_precision"],
            "specific_ap": summaries[key]["isoform_specific_knn"]["macro_average_precision"],
            "pooled_auc": summaries[key]["pooled_chemical_knn"]["macro_roc_auc"],
            "specific_auc": summaries[key]["isoform_specific_knn"]["macro_roc_auc"],
            "ap_difference": comparisons[key]["average_precision"]["difference"],
            "ap_ci": comparisons[key]["average_precision"]["scaffold_bootstrap_95ci"],
            "auc_difference": comparisons[key]["roc_auc"]["difference"],
            "auc_ci": comparisons[key]["roc_auc"]["scaffold_bootstrap_95ci"],
            "permutation_p": permutations[key]["plus_one_upper_tail_monte_carlo_p"],
        })
    strict = summaries["scaffold_and_lineage_eligible"]
    pooled = {row["isoform"]: row for row in strict["pooled_chemical_knn"]["per_isoform"]}
    specific = {row["isoform"]: row for row in strict["isoform_specific_knn"]["per_isoform"]}
    isoforms = []
    for isoform in sorted(pooled):
        isoforms.append({
            "isoform": isoform, "rows": pooled[isoform]["rows"],
            "positives": pooled[isoform]["positives"], "negatives": pooled[isoform]["negatives"],
            "pooled_ap": pooled[isoform]["average_precision"],
            "specific_ap": specific[isoform]["average_precision"],
        })
    return {
        "stages": stages,
        "audit_flags": {
            "raw_rows": audit["counts"]["raw_rows"],
            "unresolved_structure_rows": audit["counts"]["unresolved_structure_rows"],
            "exact_pair_overlaps": audit["counts"]["exact_pair_overlap_groups"],
            "exact_pair_label_conflicts": audit["counts"]["exact_pair_conflict"],
            "direct_exposed_source_rows": audit["counts"]["development_exposed_source_rows"],
        },
        "aggregate": aggregate,
        "strict_per_isoform": isoforms,
        "validation": {
            "status": validation["status"],
            "predictions_verified": validation["predictions_verified"],
            "similarity_cells_verified": validation["similarity_cells_verified"],
            "bootstrap_replicates_per_stratum": validation["bootstrap_replicates_per_stratum"],
            "permutations_per_stratum": validation["permutations_per_stratum"],
        },
    }


def style_axis(axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(MUTED)
    axis.spines["bottom"].set_color(MUTED)
    axis.tick_params(colors=INK, labelsize=9)


def build(data):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelcolor": INK, "axes.titlecolor": INK,
        "text.color": INK, "svg.fonttype": "none", "pdf.fonttype": 42,
    })
    figure, axes = plt.subplots(2, 2, figsize=(13.2, 9.2), constrained_layout=True)
    figure.suptitle("External evaluation of fixed-human CYP substrate transfer", fontsize=17,
                    fontweight="bold", x=0.02, ha="left")

    # A: same-grain nested eligibility attrition.
    axis = axes[0, 0]
    stage_labels = [row["label"] for row in data["stages"]][::-1]
    stage_values = [row["rows"] for row in data["stages"]][::-1]
    y = np.arange(len(stage_labels))
    colors = [BLUE, GOLD, GOLD, "#9aa5b1"]
    axis.barh(y, stage_values, color=colors, edgecolor=INK, linewidth=0.7)
    axis.set_yticks(y, stage_labels)
    axis.set_xlabel("Unique isoform–compound labels")
    axis.set_title("A  Eligibility audit", loc="left", fontweight="bold")
    axis.set_xlim(0, max(stage_values) * 1.18)
    axis.xaxis.grid(True, color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    for position, value in zip(y, stage_values):
        axis.text(value + max(stage_values) * 0.018, position, f"{value:,}", va="center", fontsize=9)
    flags = data["audit_flags"]
    axis.text(0.01, -0.27,
              f"Raw records: {flags['raw_rows']:,}; unresolved structures: {flags['unresolved_structure_rows']:,}; "
              f"exact pair overlaps: {flags['exact_pair_overlaps']:,} ({flags['exact_pair_label_conflicts']} conflicts); "
              f"direct exposed-source records: {flags['direct_exposed_source_rows']:,}.",
              transform=axis.transAxes, fontsize=8.2, color=MUTED, va="top", wrap=True)
    style_axis(axis)

    # B: paired macro AP comparison.
    axis = axes[0, 1]
    aggregate = data["aggregate"]
    y = np.arange(len(aggregate))[::-1]
    for position, row in zip(y, aggregate):
        axis.plot([row["pooled_ap"], row["specific_ap"]], [position, position], color=MUTED, linewidth=1.5, zorder=1)
        axis.scatter(row["pooled_ap"], position, s=65, facecolor="white", edgecolor=MUTED, linewidth=1.5, zorder=2)
        axis.scatter(row["specific_ap"], position, s=70, facecolor=BLUE, edgecolor=INK, linewidth=0.7, zorder=3)
        axis.text(row["pooled_ap"] - 0.008, position + 0.12, f"{row['pooled_ap']:.3f}", ha="right", fontsize=8.5, color=MUTED)
        axis.text(row["specific_ap"] + 0.008, position + 0.12, f"{row['specific_ap']:.3f}", ha="left", fontsize=8.5, color=BLUE)
    axis.set_yticks(y, [row["label"] for row in aggregate])
    axis.set_xlim(0.34, 0.60)
    axis.set_xlabel("Macro average precision")
    axis.set_title("B  Locked model versus pooled chemistry", loc="left", fontweight="bold")
    axis.xaxis.grid(True, color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    axis.scatter([], [], s=65, facecolor="white", edgecolor=MUTED, linewidth=1.5, label="Pooled chemistry")
    axis.scatter([], [], s=70, facecolor=BLUE, edgecolor=INK, linewidth=0.7, label="Isoform-specific")
    axis.legend(frameon=False, loc="center right", ncol=1, fontsize=8.5)
    style_axis(axis)

    # C: differences with paired scaffold-bootstrap intervals.
    axis = axes[1, 0]
    rows = []
    for row in aggregate:
        rows.extend([
            (row["label"].replace("\n", " ") + " · AP", row["ap_difference"], row["ap_ci"]),
            (row["label"].replace("\n", " ") + " · ROC AUC", row["auc_difference"], row["auc_ci"]),
        ])
    y = np.arange(len(rows))[::-1]
    for position, (label, estimate, interval) in zip(y, rows):
        color = BLUE if label.endswith("AP") else GOLD
        axis.errorbar(estimate, position,
                      xerr=np.asarray([[estimate - interval[0]], [interval[1] - estimate]]),
                      fmt="o", color=color, ecolor=color, markeredgecolor=INK,
                      markeredgewidth=0.6, capsize=3, markersize=6, linewidth=1.6)
        axis.text(interval[1] + 0.004, position, f"{estimate:+.3f}", va="center", fontsize=8.5)
    axis.axvline(0, color=INK, linewidth=1.0, linestyle="--")
    axis.set_yticks(y, [row[0] for row in rows])
    axis.set_xlim(-0.01, 0.18)
    axis.set_xlabel("Isoform-specific minus pooled")
    axis.set_title("C  Paired effect estimates (95% scaffold bootstrap)", loc="left", fontweight="bold")
    axis.xaxis.grid(True, color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    style_axis(axis)

    # D: strict-stratum per-isoform AP.
    axis = axes[1, 1]
    isoforms = data["strict_per_isoform"]
    y = np.arange(len(isoforms))[::-1]
    for position, row in zip(y, isoforms):
        axis.plot([row["pooled_ap"], row["specific_ap"]], [position, position], color=GRID, linewidth=2, zorder=1)
        axis.scatter(row["pooled_ap"], position, s=45, facecolor="white", edgecolor=MUTED, linewidth=1.3, zorder=2)
        axis.scatter(row["specific_ap"], position, s=50, facecolor=BLUE, edgecolor=INK, linewidth=0.6, zorder=3)
        axis.text(0.89, position, f"n={row['rows']:,}", va="center", ha="left", fontsize=8, color=MUTED)
    axis.set_yticks(y, [row["isoform"] for row in isoforms])
    axis.set_xlim(0.15, 1.00)
    axis.set_xlabel("Average precision")
    axis.set_title("D  New-scaffold performance by isoform", loc="left", fontweight="bold")
    axis.xaxis.grid(True, color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    style_axis(axis)

    figure.text(0.02, -0.012,
                "Scope: six fixed human CYP isoforms; source-reported binary labels are not assay harmonized. "
                "All model settings were locked before external files were opened. Independent reconstruction: PASS.",
                fontsize=8.5, color=MUTED, ha="left")
    return figure


def main():
    data = collect_data()
    OUTPUT.mkdir(exist_ok=True)
    write_json(OUTPUT / "figure6_external_transfer_source_data.json", data)
    figure = build(data)
    paths = {}
    for extension in ("png", "svg", "pdf"):
        path = OUTPUT / f"figure6_external_transfer.{extension}"
        figure.savefig(path, dpi=600 if extension == "png" else None,
                       bbox_inches="tight", facecolor="white")
        paths[extension] = path
    plt.close(figure)
    provenance = {
        "created_utc": now(), "title": "External evaluation of fixed-human CYP substrate transfer",
        "script_sha256": digest_file(Path(__file__)),
        "input_hashes": {str(path.relative_to(ROOT)): digest_file(path) for path in (
            AUDIT_PATH, SUMMARY_PATH, COMPARISON_PATH, PERMUTATION_PATH, VALIDATION_PATH)},
        "source_data_sha256": digest_file(OUTPUT / "figure6_external_transfer_source_data.json"),
        "outputs": {extension: {"path": str(path.relative_to(ROOT)), "sha256": digest_file(path)}
                    for extension, path in paths.items()},
        "visual_encoding": {
            "open_circle": "pooled chemical control", "filled_blue_circle": "isoform-specific model",
            "blue_interval": "average-precision difference", "gold_interval": "ROC-AUC difference",
        },
        "scope_note": "Fixed six human isoforms; curated binary labels, not assay-harmonized or unseen-protein evidence.",
    }
    write_json(OUTPUT / "figure6_external_transfer_provenance.json", provenance)
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
