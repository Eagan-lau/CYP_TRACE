#!/usr/bin/env python3
"""Redraw the CYP-TRACE narrative from frozen figure inputs at A4 / Arial 7 pt.

Run from any directory: python build_figures.py
PDF/SVG/PNG exports retain an uncropped 210 x 297 mm page. Word inserts use
an explicitly sized 170 mm crop with no change of scale or font size.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.text import Text
from matplotlib.ticker import PercentFormatter
from matplotlib.transforms import Bbox
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent / "Source_data" / "figure_inputs"
OUT = ROOT / "exports"
OUT.mkdir(exist_ok=True)
MM = 1 / 25.4
BLUE, TEAL, GOLD = "#285C8E", "#087E8B", "#B18522"
INK, GREY, LIGHT = "#202529", "#68747C", "#E8EDF0"
PURPLE = "#72578A"
FONT = font_manager.findfont(font_manager.FontProperties(family="Arial"), fallback_to_default=False)
plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7,
    "axes.labelsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "legend.title_fontsize": 7,
    "figure.titlesize": 7, "axes.titleweight": "bold",
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": GREY,
    "xtick.color": INK, "ytick.color": INK,
    "axes.linewidth": .55, "lines.linewidth": .8, "lines.markersize": 3,
    "xtick.major.width": .5, "ytick.major.width": .5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.bbox": None, "savefig.facecolor": "white",
    "axes.unicode_minus": False,
})


def js(name):
    return json.loads((INPUT / name).read_text(encoding="utf-8"))


def tsv(name):
    with (INPUT / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


D1 = js("figure1_evidence_atlas_source_data.json")
D2 = js("figure2_split_design_source_data.json")
D3 = js("figure3_domain_methods_source_data.json")
D4 = js("figure4_reverse_retrieval_source_data.json")
D5 = js("figure5_coverage_controls_source_data.json")
D6 = js("figure6_external_transfer_source_data.json")
CASES = js("chemical_case_studies.json")
REPRESENTATION = js("representation_count_provenance.json")
RECORDS = []


def text(fig, x, y, value, **kwargs):
    return fig.text(x / 210, 1 - y / 297, value, va="top", fontsize=7,
                    fontfamily="Arial", **kwargs)


def panel(fig, letter, title, x, y):
    text(fig, x, y, letter, weight="bold")
    text(fig, x + 5, y, title, weight="bold")


def axes(fig, x, y, width, height, *, off=False, grid="x"):
    ax = fig.add_axes([x / 210, 1 - (y + height) / 297, width / 210, height / 297])
    ax.spines[["top", "right"]].set_visible(False)
    if off:
        ax.set_axis_off()
    else:
        ax.set_axisbelow(True)
        if grid:
            ax.grid(axis=grid, color=LIGHT, linewidth=.5)
    return ax


def page(title):
    fig = plt.figure(figsize=(210 * MM, 297 * MM))
    # Main-figure titles belong in manuscript legends; retain supplementary titles.
    if title.startswith('Figure S'):
        text(fig, 20, 18, title, weight="bold")
    return fig


def arrow(fig, start, end, color=GREY, style="->"):
    patch = FancyArrowPatch((start[0] / 210, 1 - start[1] / 297),
                            (end[0] / 210, 1 - end[1] / 297),
                            transform=fig.transFigure, arrowstyle=style,
                            mutation_scale=7, lw=.65, color=color)
    fig.add_artist(patch)


def forest(ax, rows, colors=None, xlabel="Difference in MRR", digits=4, numbers=False):
    for i, row in enumerate(rows):
        value, (lo, hi) = row["difference"], row["ci"]
        color = colors[i] if colors else BLUE
        ax.errorbar(value, i, xerr=[[value - lo], [hi - value]], fmt="o", color=color,
                    capsize=2, markersize=3, lw=.8)
        if numbers:
            ax.text(1.02, i, f"{value:+.{digits}f}", transform=ax.get_yaxis_transform(), va="center")
    ax.axvline(0, color=GREY, lw=.65, ls="--")
    ax.set_yticks(range(len(rows)), [row["label"] for row in rows])
    ax.set_ylim(len(rows) - .45, -.55)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0, pad=4)


def molecule(fig, smiles, x, y, width, height):
    """RDKit coordinates, matplotlib bonds and Arial atom labels at 7 pt."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(smiles)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    rdDepictor.Compute2DCoords(mol)
    xy = np.array(mol.GetConformer().GetPositions())[:, :2]
    ax = axes(fig, x, y, width, height, off=True)
    for bond in mol.GetBonds():
        a, b = xy[bond.GetBeginAtomIdx()], xy[bond.GetEndAtomIdx()]
        if bond.GetBondTypeAsDouble() >= 2:
            v = b - a
            n = np.array([-v[1], v[0]]) / np.linalg.norm(v) * .06
            ax.plot([a[0] + n[0], b[0] + n[0]], [a[1] + n[1], b[1] + n[1]], color=INK, lw=.7)
            ax.plot([a[0] - n[0], b[0] - n[0]], [a[1] - n[1], b[1] - n[1]], color=INK, lw=.7)
        else:
            ax.plot([a[0], b[0]], [a[1], b[1]], color=INK, lw=.7)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C":
            continue
        label = atom.GetSymbol()
        if atom.GetTotalNumHs():
            label += "H" + (str(atom.GetTotalNumHs()) if atom.GetTotalNumHs() > 1 else "")
        ax.text(*xy[atom.GetIdx()], label, ha="center", va="center",
                color=TEAL, bbox={"facecolor": "white", "edgecolor": "none", "pad": .3})
    ax.set_aspect("equal")
    ax.set_xlim(xy[:, 0].min() - .35, xy[:, 0].max() + .35)
    ax.set_ylim(xy[:, 1].min() - .45, xy[:, 1].max() + .45)


def finish(fig, stem, panels, bottom, inputs, pdf):
    fig.canvas.draw()
    # All matplotlib text, including panel letters and annotations, is 7 pt.
    text_artists = [item for item in fig.findobj(Text) if item.get_visible() and item.get_text()]
    bad = [(item.get_text(), item.get_fontsize(), item.get_fontfamily()) for item in text_artists
           if abs(item.get_fontsize() - 7) > 1e-9 or item.get_fontfamily() != ["Arial"]]
    if bad:
        raise AssertionError(bad)
    for suffix in ("pdf", "svg"):
        fig.savefig(OUT / f"{stem}.{suffix}")
    fig.savefig(OUT / f"{stem}.png", dpi=300)
    crop = Bbox.from_extents(20 * MM, (297 - bottom) * MM, 190 * MM, (297 - 15) * MM)
    fig.savefig(OUT / f"{stem}_content.png", dpi=450, bbox_inches=crop, pad_inches=0)
    if not stem.startswith('figureS'):
        assert bottom - 15 <= 233, (stem, bottom - 15)
        # Native-size journal pages, no font scaling; A4 working pages remain above.
        for suffix in ('pdf','svg'):
            fig.savefig(OUT / f"{stem}_journal.{suffix}", bbox_inches=crop, pad_inches=0)
        fig.savefig(OUT / f"{stem}_journal.tiff", dpi=600, bbox_inches=crop,
                    pad_inches=0, pil_kwargs={'compression':'tiff_lzw'})
    pdf.savefig(fig)
    RECORDS.append({"stem": stem, "panels": panels, "page_mm": [210, 297],
                    "font": "Arial", "font_pt": 7, "visible_text_objects": len(text_artists),
                    "content_mm": [170, bottom - 15], "content_bottom_mm": bottom,
                    "inputs": inputs})
    plt.close(fig)


def figure1(pdf):
    fig = page("Figure 1   Chemical representation determines how public evidence connects")
    panel(fig, "A", "One documented main transformation, two complete reaction identities", 20, 31)
    r = CASES["reaction_identity"]
    molecule(fig, r["main_transformation"]["substrates"][0], 33, 42, 55, 17)
    molecule(fig, r["main_transformation"]["products"][0], 119, 39, 58, 23)
    arrow(fig, (94, 51), (114, 51), INK)
    text(fig, 39, 65, "Decane")
    text(fig, 133, 65, "Decan-3-ol")
    text(fig, 20, 76, "P450Rdb and UniProt: Q0D1W9; PMID 32147902")
    text(fig, 20, 83, "Oxidized flavin is encoded in different protonation/charge forms.\nParticipant-role records preserve both complete expressions.")
    panel(fig, "B", "Shared relationships depend on reaction identity", 20, 102)
    ax = axes(fig, 72, 115, 107, 22)
    shared_counts = [REPRESENTATION["complete_component_shared_relationships"], D1["edge_union"]["shared"]]
    ax.barh([0, 1], shared_counts, height=.43, color=[GREY, TEAL])
    ax.set_yticks([0, 1], ["Complete components", "Main transformation"])
    ax.invert_yaxis(); ax.set_xlim(0, 850); ax.set_xticks([0, 200, 400, 600, 800])
    ax.set_xlabel("Cross-lineage shared relationships")
    for y, n in enumerate(shared_counts):
        ax.text(n + 13, y, f"{n:,}", va="center")
    text(fig, 20, 151, "The shared count changes with representation; the source assertions are the same.")
    panel(fig, "C", "The reconstructed atlas retains source lineage", 20, 164)
    ax = axes(fig, 24, 178, 162, 13, off=True)
    union = D1["edge_union"]
    left = 0
    for n, label, color in [(union["p450_unique"], "P450Rdb only", GOLD),
                             (union["shared"], "Shared", TEAL),
                             (union["uniprot_unique"], "UniProt only", BLUE)]:
        ax.barh(0, n, left=left, color=color, height=.55, edgecolor="white", linewidth=.5)
        ax.text(left + n / 2, 0, f"{label}\n{n:,}", va="center", ha="center", color="white")
        left += n
    ax.set_xlim(0, union["union"])
    text(fig, 20, 198, "600 exact sequences   ·   2,304 sequence–reaction edges   ·   1,341 directed transformations")
    text(fig, 20, 205, "1,225 single substrate/product pairs; 116 multicomponent transformations.\nSwiss-Prot and the UniProt API are one annotation lineage.", color=GREY)
    finish(fig, "figure1_evidence_identity", 3, 218,
           ["figure1_evidence_atlas_source_data.json", "chemical_case_studies.json", "representation_count_provenance.json"], pdf)


def figure2(pdf):
    fig = page("Figure 2   Coverage and ranking answer different evaluation questions")
    panel(fig, "A", "Held-out axes define three retrospective tasks", 20, 31)
    specs = [(31, "Protein cold", 5, 600), (86, "Chemical cold", 5, 772), (141, "Double cold", 25, 772)]
    for x, label, blocks, panels in specs:
        text(fig, x, 43, label, weight="bold")
        ax = axes(fig, x + 3, 53, 25, 25, off=True)
        for row in range(4):
            for col in range(4):
                test = row == 3 if label == "Protein cold" else col == 3 if label == "Chemical cold" else row == 3 and col == 3
                excluded = label == "Double cold" and (row == 3 or col == 3) and not test
                face = TEAL if test else "white" if excluded else LIGHT
                ax.add_patch(Rectangle((col, 3 - row), 1, 1, facecolor=face,
                                       edgecolor=GREY if excluded else "white", lw=.45 if excluded else .8,
                                       hatch="//" if excluded else None))
        ax.set_xlim(0, 4); ax.set_ylim(0, 4)
        text(fig, x, 82, f"{blocks} outer blocks\n{panels} query panels")
    text(fig, 20, 95, "Rows: protein groups; columns: chemical groups. Grey: training; teal: held out; hatching: excluded.")
    text(fig, 20, 103, "145 protein groups × 240 chemical groups; recorded publication overlap is purged.\n35 outer blocks contain 255 inner blocks. Each task evaluates all 2,304 edges once.", color=GREY)
    panel(fig, "B", "Reachability and conditional ranking decompose end-to-end MRR", 20, 122)
    text(fig, 20, 132, "New protein, seen chemistry: 95 queries from 27 protein groups; 1,341 supplied candidates.")
    rows = {r["method"]: r for r in tsv("coverage_ranking_decomposition.tsv") if r["task"] == "protein_cold_seen"}
    methods = [("mmseqs_top1", "MMseqs2", BLUE), ("esm_cosine_top1", "ESM", TEAL)]
    columns = [(69, "Reachable-query weight", "weighted_query_positive_coverage_Cw", .13),
               (112, "Conditional MRR", "covered_weighted_MRR", 1),
               (155, "End-to-end MRR", "end_to_end_MRR", .1)]
    for x, title, field, limit in columns:
        text(fig, x - 2, 144, title)
        ax = axes(fig, x, 155, 30, 20)
        for i, (method, name, color) in enumerate(methods):
            value = float(rows[method][field])
            ax.hlines(i, 0, value, color=color, lw=1)
            ax.plot(value, i, "o", color=color)
            label = f"{value:.2%}" if field == "weighted_query_positive_coverage_Cw" else f"{value:.4f}"
            ax.text(value, i - .23, label, ha="center", color=color)
        ax.set_xlim(0, limit); ax.set_ylim(1.55, -.55)
        ax.set_xticks([0, limit / 2, limit])
        ax.set_yticks([0, 1], ["MMseqs2", "ESM"] if x == 69 else [])
        if field == "weighted_query_positive_coverage_Cw":
            ax.set_xticklabels(["0%", "6.5%", "13%"])
    text(fig, 104, 162, "×", ha="center")
    text(fig, 148, 162, "=", ha="center")
    text(fig, 20, 189, "Identical candidate intersection: both methods have MRR 0.7500.", weight="bold")
    text(fig, 20, 197, "Only 3 queries from 3 groups contribute 6 shared positives; 92 query panels are undefined.\nCandidate coverage is separate: MMseqs2 0.362%; ESM 0.280% of catalogue entries.")
    finish(fig, "figure2_coverage_ranking", 2, 211,
           ["figure2_split_design_source_data.json", "coverage_ranking_decomposition.tsv", "common_domain_key_comparison.tsv"], pdf)


def figure3(pdf):
    fig = page("Figure 3   Protein information is constrained by availability, controls and study exposure")
    panel(fig, "A", "Representation availability", 20, 31)
    ax = axes(fig, 55, 45, 39, 47)
    availability = D5["availability"]
    ax.barh(range(5), [r["count"] for r in availability], color=GREY, height=.45)
    ax.set_yticks(range(5), ["Global ESM", "CLEAN", "Sequence sites", "Structure", "Structure sites"])
    ax.set_xlim(0, 680); ax.invert_yaxis(); ax.set_xticks([0, 300, 600]); ax.set_xlabel("Core sequences")
    for i, r in enumerate(availability):
        ax.text(r["count"] + 10, i, str(r["count"]), va="center")
    text(fig, 20, 103, "Matched structures: 14 queries, 11 groups.\nFoldseek and MMseqs2 both MRR 0.0417.", color=GREY)
    panel(fig, "B", "Local and interaction effects", 106, 31)
    rows = []
    for kind, source in [("Sites", D5["local_effects"][:3]), ("Joint", D3["interaction_effects"])]:
        for i, r in enumerate(source):
            rows.append({**r, "label": kind + " · " + ["chemical", "protein", "double"][i]})
    ax = axes(fig, 143, 45, 43, 55)
    forest(ax, rows, [GOLD]*3 + [BLUE]*3, xlabel="MRR difference vs global")
    ax.set_xlim(-.048, .009); ax.set_xticks([-.04, -.02, 0])
    text(fig, 106, 109, "48 groups; sequence-projected cohort.\nConditional 95% protein-group intervals.", color=GREY)
    panel(fig, "C", "Interaction permutation controls", 20, 130)
    ax = axes(fig, 53, 143, 41, 28)
    control = [{**r, "label": label} for r, label in zip(D3["control_effects"], ["Position\npermutation", "Reaction-centre\npermutation"])]
    forest(ax, control, [TEAL, GREY], xlabel="Primary minus control gain")
    ax.set_xlim(-.0001, .0012); ax.set_xticks([0, .0005, .001]); ax.set_xticklabels(["0", "0.0005", "0.0010"])
    text(fig, 20, 183, "Chemical cold; conditional 95% intervals.", color=GREY)
    panel(fig, "D", "Four-group contribution sensitivity", 106, 130)
    groups = tsv("interaction_four_group_sensitivity.tsv")
    loo = tsv("interaction_leave_one_group_out.tsv")
    ax = axes(fig, 134, 143, 52, 28)
    for i, (row, omitted) in enumerate(zip(groups, loo)):
        ax.plot(float(row["paired_difference"]), i, "o", color=BLUE)
        ax.plot(float(omitted["fixed_prediction_mean_difference"]), i, "D", color=GREY, mfc="white", ms=3)
    ax.axvline(0, color=GREY, lw=.65, ls="--")
    ax.set_xlim(-.045, .07); ax.set_yticks(range(4), [f"Group {i+1}" for i in range(4)])
    ax.invert_yaxis(); ax.set_xticks([-.04, 0, .04]); ax.set_xlabel("Joint minus global MRR")
    text(fig, 106, 183, "Filled circle: group effect\nOpen diamond: mean after omitting group", color=GREY)
    panel(fig, "E", "CLEAN pruning gains depend on exact training exposure", 20, 203)
    ax = axes(fig, 55, 224, 125, 29)
    for i, row in enumerate(D5["pruning_effects"]):
        for shift, key, color, marker, fill in [(-.22, "all", GREY, "o", GREY),
                                              (0, "exposed", GOLD, "s", GOLD),
                                              (.22, "unexposed", BLUE, "o", "white")]:
            v = row[key]
            ax.errorbar(v["difference"], i + shift,
                        xerr=[[v["difference"] - v["ci"][0]], [v["ci"][1] - v["difference"]]],
                        fmt=marker, color=color, mfc=fill, ms=3, capsize=1.5, lw=.7)
    ax.axvline(0, color=GREY, ls="--", lw=.65)
    ax.set_yticks(range(3), ["Chemical cold", "Protein cold", "Double cold"])
    ax.set_ylim(2.5, -.5); ax.set_xlim(-.08, .48); ax.set_xticks([0, .2, .4])
    ax.set_xlabel("Pruned minus unpruned MRR; conditional 95% intervals")
    text(fig, 25, 213, "● All sequences    ■ Exact training overlap: 328    ○ No exact training overlap: 272", color=GREY)
    finish(fig, "figure3_protein_controls", 5, 268,
           ["figure3_domain_methods_source_data.json", "figure5_coverage_controls_source_data.json",
            "interaction_four_group_sensitivity.tsv", "interaction_leave_one_group_out.tsv"], pdf)


def figure4(pdf):
    fig = page("Figure 4   Reverse retrieval is limited by both coverage and ranking")
    text(fig, 20, 29, "1,803 candidate sequences; 4,880 reaction–taxid panels; 6,033 documented positive instances.\nCandidate sets contain 2–141 CYPs observed for the same taxid.")
    panel(fig, "A", "Ranking within the same taxid", 20, 45)
    ax = axes(fig, 64, 63, 86, 40)
    rows = D4["performance"]
    for j, (key, color, marker) in enumerate([("MMseqs2", TEAL, "o"), ("Label prior", GOLD, "s"), ("Uniform", GREY, "D")]):
        ax.plot([r["mrr"][key] for r in rows], np.arange(3) + (j-1)*.16,
                linestyle="none", marker=marker, color=color, mfc="white" if j == 2 else color, label=key)
    ax.set_yticks(range(3), ["Seen reaction\nnew positive", "Unseen reaction\nnew positive", "Unseen reaction\nrepresented positive"])
    ax.set_ylim(2.5, -.5); ax.set_xlim(0, .65); ax.set_xticks([0, .2, .4, .6]); ax.set_xlabel("Taxid-macro MRR")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=3, frameon=False, handletextpad=.4, columnspacing=1)
    text(fig, 157, 54, "Positive\ncoverage", ha="center")
    text(fig, 183, 54, "Taxids", ha="center")
    for i, row in enumerate(rows):
        y = 63 + (i + .5)/3 * 40
        text(fig, 157, y-1, f"{row['coverage']['MMseqs2']:.1%}", ha="center")
        text(fig, 183, y-1, str(row["taxids"]), ha="center")
    panel(fig, "B", "Prespecified acceptance criteria", 20, 121)
    ax = axes(fig, 66, 141, 120, 25, off=True)
    for i, flags in enumerate(D4["gates"]):
        for j, passed in enumerate(flags):
            ax.add_patch(Rectangle((j-.47, i-.38), .94, .76,
                                  facecolor="#DCECEE" if passed else LIGHT, edgecolor="white"))
            ax.text(j, i, "Yes" if passed else "No", ha="center", va="center", color=TEAL if passed else GREY)
    ax.set_xlim(-.5, 4.5); ax.set_ylim(2.5, -.5)
    labels = ["≥20 taxids", "≥80%\ncoverage", "Above\nlabel prior", "Above\nuniform", "Permutation\np ≤ 0.05"]
    for i, label in enumerate(labels):
        text(fig, 78+24*i, 131, label, ha="center")
    for i, label in enumerate(["Seen / new", "Unseen / new", "Unseen / represented"]):
        text(fig, 20, 143 + i*8.33, label)
    text(fig, 20, 176, "No reverse-retrieval cell meets all five criteria. Candidate sets are source-observed, not complete proteomes.")
    finish(fig, "figure4_reverse_retrieval", 2, 185, ["figure4_reverse_retrieval_source_data.json"], pdf)


def figure5(pdf):
    fig = page("Figure 5   Isoform-conditioned chemistry transfers to externally curated compounds")
    panel(fig, "A", "The external endpoint is qualified before testing", 20, 31)
    xs = [20, 64, 108, 152]
    names = ["Unique evaluable\nisoform–compound pairs", "Pairs on exact-new\ncompounds", "Also eligible\nsource lineage", "Also new\nscaffolds"]
    for i, (x, row, label) in enumerate(zip(xs, D6["stages"], names)):
        text(fig, x+16, 45, f"{row['rows']:,}", ha="center", color=BLUE if i == 3 else INK, weight="bold")
        text(fig, x+16, 52, label, ha="center")
        if i < 3:
            arrow(fig, (x+33, 47), (xs[i+1]-2, 47))
    text(fig, 20, 67, "Strict endpoint: 3,035 labels; 1,001 compounds; 844 scaffolds; 6 fixed human isoforms.")
    text(fig, 20, 74, "Prediction settings were locked before the external source was opened.", color=GREY)
    panel(fig, "B", "The gain is present for all six isoforms", 20, 92)
    ax = axes(fig, 43, 111, 70, 70)
    for i, row in enumerate(D6["strict_per_isoform"]):
        ax.plot([row["pooled_ap"], row["specific_ap"]], [i, i], color=GREY, lw=.7)
        ax.plot(row["pooled_ap"], i, "o", color=GREY, mfc="white", ms=3.5)
        ax.plot(row["specific_ap"], i, "o", color=TEAL, ms=3.5)
        ax.text(1.03, i, f"n={row['rows']}", transform=ax.get_yaxis_transform(), va="center")
    ax.set_yticks(range(6), [r["isoform"] for r in D6["strict_per_isoform"]])
    ax.set_xlim(0, 1); ax.set_ylim(5.5, -.5); ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xlabel("Average precision")
    text(fig, 20, 102, "○ Pooled chemistry   ● Isoform-specific", color=GREY)
    panel(fig, "C", "Paired external effect", 139, 92)
    for i, row in enumerate(D6["aggregate"]):
        top = 111 + 42*i
        title = "Exact-new + eligible lineage" if i == 0 else "Also new scaffolds"
        text(fig, 139, top, f"{title}\n{row['rows']:,} labels")
        ax = axes(fig, 142, top + 14, 42, 15)
        forest(ax, [{"difference":row["ap_difference"], "ci":row["ap_ci"], "label":""}], [BLUE if i == 0 else TEAL], xlabel="")
        ax.set_xlim(0,.18); ax.set_xticks([0,.08,.16])
        if i == 1:
            ax.set_xlabel("Difference in macro AP")
    text(fig, 20, 204, "Paired gains: +0.1375 [0.1146, 0.1572] and +0.1319 [0.1081, 0.1571].")
    text(fig, 20, 210, "Intervals: paired 95% scaffold bootstrap, 5,000 replicates. Both external gates passed.")
    text(fig, 20, 217, "This endpoint tests new compounds for represented human CYPs. It does not test new proteins or products.", color=GREY)
    finish(fig, "figure5_external_transfer", 3, 226, ["figure6_external_transfer_source_data.json"], pdf)


def figure6(pdf):
    fig = page("Figure 6   Selective release makes the supported use of predictions explicit")
    panel(fig, "A", "Higher precision retains fewer documented positives", 20, 31)
    ax = axes(fig, 40, 48, 125, 78, grid="both")
    data = tsv("selective_release_tradeoffs.tsv")
    for method, color, marker in [("pooled_chemical_knn", GREY, "o"), ("isoform_specific_knn", TEAL, "s")]:
        rows = [r for r in data if r["method"] == method]
        x = [float(r["macro_positive_recall"]) for r in rows]
        y = [float(r["macro_precision"]) for r in rows]
        ax.plot(x, y, marker=marker, color=color, mfc="white" if marker == "o" else color,
                linestyle="--" if marker == "o" else "-", label="Pooled chemistry" if marker == "o" else "Isoform-specific")
        if marker == "s":
            for i, r in enumerate(rows):
                label = f"{float(r['target_coverage']):.0%} release; {int(r['selected_rows']):,} labels"
                offsets = [(-5,-13), (5,12), (5,7), (5,7), (5,7)]
                ax.annotate(label, (x[i], y[i]), xytext=offsets[i], textcoords="offset points",
                            ha="right" if i == 0 else "left", fontsize=7)
    ax.set_xlim(0, 1.06); ax.set_ylim(.3, .73)
    ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    ax.set_xlabel("Macro recall of documented positives"); ax.set_ylabel("Macro precision among released labels")
    ax.legend(frameon=False, loc="upper right", handlelength=1.6)
    text(fig, 20, 139, "At 10% release: precision 66.9%; positive recall 20.3%; 306 released labels.")
    contrast=D5["matched_coverage_comparison"]["0.1"]
    diff=contrast["isoform_specific_minus_pooled_macro_precision"]
    lo,hi=contrast["scaffold_bootstrap_95ci"]
    text(fig, 20, 146, f"Matched 10% coverage: gain {100*diff:.1f} percentage points; paired 95% interval [{100*lo:.1f}, {100*hi:.1f}].")
    text(fig, 20, 154, "3,035 external labels, 1,001 compounds, 844 scaffolds. Coverage is descriptive; cutoffs are not portable.", color=GREY)
    panel(fig, "B", "The executable strategy exposes the evaluated scope", 20, 166)
    text(fig, 23, 183, "Human CYP +\nsupplied compound", weight="bold")
    arrow(fig, (64, 188), (79, 188))
    text(fig, 84, 180, "Isoform-specific chemical neighbours\nScores, exposure and nearest evidence\n6 isoforms externally evaluated; 3 development only")
    text(fig, 23, 216, "General CYP sequence +\nsupplied reaction", weight="bold")
    arrow(fig, (65, 221), (79, 221))
    text(fig, 84, 215, "Exact documented edge")
    arrow(fig, (135, 219), (149, 211))
    arrow(fig, (135, 221), (149, 232))
    text(fig, 153, 206, "Present\nReturn evidence", color=TEAL)
    text(fig, 153, 230, "Absent\nAbstain", color=GREY)
    finish(fig, "figure6_selective_strategy", 2, 247,
           ["selective_release_tradeoffs.tsv", "figure5_coverage_controls_source_data.json", "strategy_scope.json"], pdf)


def supplement1(pdf):
    fig = page("Figure S1   Source resolution and group concentration")
    panel(fig, "A", "Resolution fields use source assertions as denominators", 20, 31)
    ax = axes(fig, 49, 48, 130, 35, grid=None)
    arr = np.asarray(D1["rates"])
    cmap = LinearSegmentedColormap.from_list("coverage", ["#F0F3F4", TEAL])
    ax.imshow(arr, vmin=0, vmax=1, cmap=cmap, aspect="auto")
    ax.set_yticks(range(4), ["P450Rdb", "SABIO-RK", "Swiss-Prot", "UniProt API"])
    ax.set_xticks(range(5), ["Sequence", "PF00067", "Publication", "Parseable\nchemistry", "Exact initial\ncore match"])
    ax.tick_params(length=0)
    for (i,j), value in np.ndenumerate(arr):
        ax.text(j,i,f"{value:.0%}",ha="center",va="center",color="white" if value > .6 else INK)
    panel(fig, "B", "Group-size distributions", 20, 111)
    ax = axes(fig, 38, 125, 66, 51, grid="y")
    for key, color, label in [("protein_group_sizes", BLUE, "Protein groups"), ("chemical_component_sizes", GOLD, "Chemical groups")]:
        values = D2[key]
        ax.plot(np.arange(1,len(values)+1), values, color=color, label=label)
    ax.set_yscale("log"); ax.set_yticks([1,10,100], ["1","10","100"])
    ax.set_xlabel("Group rank"); ax.set_ylabel("Members per group (log scale)")
    ax.legend(frameon=False, loc="upper right")
    panel(fig, "C", "Protein folds retain group imbalance", 118, 111)
    ax = axes(fig, 132, 125, 51, 51, grid="y")
    folds = D2["folds"]
    ax.bar(range(5), [r["sequences"] for r in folds], width=.6, color=LIGHT, edgecolor=BLUE, lw=.65)
    for i, row in enumerate(folds):
        ax.text(i,row["sequences"]+3,f"{row['sequences']}\n({row['groups']} g)",ha="center")
    ax.set_xticks(range(5),range(1,6)); ax.set_ylim(0,170); ax.set_xlabel("Outer protein fold"); ax.set_ylabel("Held-out sequences")
    text(fig, 20, 192, "Protein threshold: 40% identity and ≥80% bilateral coverage. Largest groups: 138 proteins, 587 reactions.")
    text(fig, 20, 201, "Recorded publication purges by protein fold: 10, 115, 184, 115 and 182 training edges.\nThe 138-sequence fold is one homology group; folds are not independent biological replicates.", color=GREY)
    finish(fig, "figureS1_source_splits", 3, 216, ["figure1_evidence_atlas_source_data.json", "figure2_split_design_source_data.json"], pdf)


def supplement2(pdf):
    fig = page("Figure S2   Applicability domains and global conditional effects")
    panel(fig, "A", "Scores are defined only inside each expert domain", 20, 31)
    ax = axes(fig, 59, 46, 126, 58, grid=None)
    cells = D3["router_cells"]
    methods = D3["router_methods"]
    for i,row in enumerate(cells):
        for j,method in enumerate(methods):
            value=row["mrr"][method["key"]]
            color=LIGHT if value is None else "white"
            ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=color,edgecolor=LIGHT,lw=.5))
            ax.text(j,i,"NA" if value is None else f"{value:.4f}",ha="center",va="center",color=GREY if value is None else INK)
    ax.set_xlim(-.5,5.5);ax.set_ylim(3.5,-.5)
    ax.set_xticks(range(6),["Chemistry\nprior","Homology\ntransport","MMseqs2\ntop-1","Global\nresidual","Joint\ninteraction","Uniform"])
    ax.set_yticks(range(4),["New / seen\n93 panels; 25 groups", "New / unseen*\n622 panels; 138 groups", "Double cold\n725 panels; 140 groups", "Represented / unseen\n66 panels; 4 groups"])
    ax.tick_params(length=0)
    text(fig,20,117,"*Chemical-cold accounting cell. NA denotes an undefined expert; it is not a numerical score of zero.")
    panel(fig,"B","Global protein correction relative to the chemical prior",20,139)
    ax=axes(fig,68,155,104,40)
    forest(ax,D3["global_effects"],numbers=True)
    ax.set_xlim(-.005,.015);ax.set_xticks([-.005,0,.005,.010,.015])
    text(fig,20,210,"All 145 protein groups; 600 protein-cold or 772 chemical-cold / double-cold panels.\nIntervals are conditional 95% protein-group bootstraps.",color=GREY)
    finish(fig,"figureS2_domain_global",2,225,["figure3_domain_methods_source_data.json"],pdf)


def supplement3(pdf):
    fig=page("Figure S3   Chemical distance and post-evaluation sensitivity")
    panel(fig,"A","Scaffold novelty retains a range of chemical similarities",20,31)
    ax=axes(fig,38,47,68,58,grid="y")
    bins=[r for r in tsv("external_chemical_distance_bins.tsv") if r["scope"]=="unique_compounds"]
    counts=[int(r["count"]) for r in bins]
    ax.bar(np.arange(5)*.2+.1,counts,width=.196,color=TEAL)
    for x,n in zip(np.arange(5)*.2+.1,counts):
        ax.text(x,n+10,str(n),ha="center")
    ax.set_xlim(0,1);ax.set_ylim(0,680);ax.set_xticks(np.arange(0,1.01,.2))
    ax.set_xlabel("Maximum training-set Tanimoto");ax.set_ylabel("Unique external compounds")
    text(fig,117,53,"1,001 compounds\nMedian 0.318\nIQR 0.250–0.449\n95th percentile 0.636\nRange 0.143–0.882")
    panel(fig,"B","Equal-neighbour and logistic comparisons",20,129)
    ax=axes(fig,68,146,102,45)
    eq=tsv("equal_k_sensitivity.tsv")
    rows=[("Pooled k=25",float(eq[0]["pooled_macro_AP"]),GREY),
          ("Isoform k=25",float(eq[0]["isoform_specific_macro_AP"]),TEAL),
          ("Pooled k=51",float(eq[1]["pooled_macro_AP"]),GREY),
          ("Isoform k=51",float(eq[1]["isoform_specific_macro_AP"]),TEAL),
          ("Isoform logistic",float(tsv("logistic_regression_supplement.tsv")[-1]["strict_external_AP"]),PURPLE)]
    for i,(label,v,c) in enumerate(rows):
        ax.plot(v,i,"o",color=c);ax.text(v+.006,i,f"{v:.4f}",va="center",color=c)
    ax.set_yticks(range(5),[r[0] for r in rows]);ax.set_ylim(4.5,-.5);ax.set_xlim(.38,.60);ax.set_xticks([.4,.45,.5,.55,.6]);ax.set_xlabel("Macro average precision")
    text(fig,20,205,"The strict external cohort is the same for every row: 3,035 labels across six isoforms.\nThese controls were added after the external labels were opened; they do not revise the locked model.",color=GREY)
    finish(fig,"figureS3_chemical_sensitivity",2,220,
           ["external_chemical_distance_bins.tsv","external_chemical_distance_quantiles.tsv",
            "equal_k_sensitivity.tsv","logistic_regression_supplement.tsv"],pdf)


def main():
    with PdfPages(OUT / "CYP_TRACE_Main_Figures_A4.pdf") as pdf:
        for build in [figure1,figure2,figure3,figure4,figure5,figure6]:
            build(pdf)
    with PdfPages(OUT / "CYP_TRACE_Supplementary_Figures_A4.pdf") as pdf:
        for build in [supplement1,supplement2,supplement3]:
            build(pdf)
    input_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(INPUT.iterdir()) if p.is_file()}
    report={"status":"BUILT_PENDING_VISUAL_REVIEW","figure_count":len(RECORDS),
            "font_file":FONT,"font_pt":7,"page_mm":[210,297],"figures":RECORDS,
            "input_sha256":input_hashes,"source_code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT / "FIGURE_BUILD_REPORT.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"figures":len(RECORDS),"panel_counts":[r["panels"] for r in RECORDS],
                      "page_mm":[210,297],"font":"Arial 7 pt","status":report["status"]},indent=2))


if __name__=="__main__":
    main()
