"""Figure 6C: executable routes and typed outputs from recorded runtime examples."""
from __future__ import annotations

import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle
from matplotlib.transforms import Bbox
from rdkit import Chem
from rdkit.Chem import rdDepictor

import figure_primitives as p


def _query_molecule(fig, smiles, x, y, width, height):
    """Vector chemical depiction with explicit aromatic and triple bond orders."""
    molecule = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(molecule, clearAromaticFlags=True)
    rdDepictor.Compute2DCoords(molecule)
    points = np.asarray(molecule.GetConformer().GetPositions())[:, :2]
    ax = p.axes(fig, x, y, width, height, off=True)
    for bond in molecule.GetBonds():
        a, b = points[bond.GetBeginAtomIdx()], points[bond.GetEndAtomIdx()]
        vector = b - a
        normal = np.array([-vector[1], vector[0]]) / np.linalg.norm(vector) * .085
        order = int(bond.GetBondTypeAsDouble())
        offsets = {1: [0], 2: [-1, 1], 3: [-1.65, 0, 1.65]}[order]
        for offset in offsets:
            start, end = a + offset * normal, b + offset * normal
            ax.plot([start[0], end[0]], [start[1], end[1]], color=p.INK, lw=.55)
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() == "C":
            continue
        label = atom.GetSymbol()
        hydrogens = atom.GetTotalNumHs()
        if hydrogens:
            label += "H" + (str(hydrogens) if hydrogens > 1 else "")
        ax.text(*points[atom.GetIdx()], label, ha="center", va="center", color=p.INK,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": .15})
    ax.set_aspect("equal")
    ax.set_xlim(points[:, 0].min() - .5, points[:, 0].max() + .5)
    ax.set_ylim(points[:, 1].min() - .5, points[:, 1].max() + .5)


def draw(fig):
    p.panel(fig, "C", "CYP-TRACE separates scoring, evidence retrieval and abstention", 20, 114)
    example = p.js("figure6_runtime_example.json")
    scope = p.js("strategy_scope.json")
    fixed = scope["fixed_human"]
    verified = example["existing_full_software_verification"]
    chemical = example["chemical_output"]
    evidence = example["evidence_output"]["evidence"][0]
    expected = scope["runtime_validation"]
    assert verified["external_rows_recomputed"] == expected["external_scores_expected"]
    assert verified["exact_edge_keys_retrieved"] == expected["exact_edge_keys_expected"]
    assert max(verified["external_score_max_absolute_difference"].values()) == 0
    controls = example["abstention_controls"]
    assert all(row["abstained"] and row["score"] is None and not row["evidence"]
               for row in controls.values())
    assert not chemical["domain"]["exact_training_compound"]
    assert not chemical["domain"]["training_scaffold_seen"]

    # Millimetres within Panel C. One runtime container; three differently
    # coloured output objects. Arrows encode actual code paths, not chronology.
    origin_y = 122
    # Compress layout spacing, not fonts. All labels stay native Arial 7 pt.
    vertical_scale = .82
    ax = p.axes(fig, 20, origin_y, 170, 120*vertical_scale, off=True)
    ax.set_xlim(0, 170)
    ax.set_ylim(120, 0)
    navy = "#173B56"

    def label(x, y, value, color=p.INK, bold=False, ha="left", va="top"):
        return ax.text(x, y, value, ha=ha, va=va, color=color,
                       weight="bold" if bold else "normal", linespacing=1.3)

    def rect(x, y, w, h, face, edge="none", lw=.6):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=face, edgecolor=edge, lw=lw))

    def route(points, color, lw=.9):
        xs, ys = zip(*points)
        ax.plot(xs[:-1], ys[:-1], color=color, lw=lw, solid_capstyle="round")
        ax.add_patch(FancyArrowPatch(points[-2], points[-1], arrowstyle="-|>",
                                    mutation_scale=7, linewidth=lw, color=color,
                                    shrinkA=0, shrinkB=.5))

    def document(x, y, w, h, color, face):
        fold = 3
        ax.add_patch(Polygon([(x, y), (x+w-fold, y), (x+w, y+fold),
                              (x+w, y+h), (x, y+h)], closed=True,
                             facecolor=face, edgecolor=color, lw=.6))
        rect(x, y, 1, h, color)
        ax.plot([x+w-fold, x+w-fold, x+w], [y, y+fold, y+fold], color=color, lw=.5)

    label(0, 1, "ILLUSTRATIVE INPUTS", p.GREY, True)
    label(43, 1, "TWO EXECUTABLE ROUTES", p.GREY, True)
    label(120, 1, "RETURNED OUTPUTS", p.GREY, True)

    # Distributed chemical example: visual input, not an invented molecular icon.
    _query_molecule(fig, chemical["canonical_smiles"], 20, origin_y+12*vertical_scale, 33, 24*vertical_scale)
    label(16.5, 40, "SMILES + " + chemical["isoform"], bold=True, ha="center")
    label(16.5, 45, "Fixed-human query", p.GREY, ha="center")
    route([(34, 30), (43, 30)], p.TEAL)

    # One runtime enclosure gives the software visual priority over its inputs.
    rect(43, 10, 65, 91, "#F5F7F9", "#BFCAD2")
    rect(43, 10, 65, 8, navy)
    label(47, 12.5, "CYP-TRACE  /  strategy CLI", "white", True)
    rect(43, 18, 1.1, 31, p.TEAL)
    label(47, 22, "human-substrate", p.TEAL, True)
    label(47, 28, "Valid SMILES + fixed CYP")
    label(47, 34, f"Chemical kNN · k = {fixed['neighbour_count']}")
    label(47, 39, f"{fixed['training_compounds']:,} reference compounds")
    label(47, 44, f"{fixed['externally_evaluated_isoforms']} externally evaluated / "
                    f"{fixed['development_only_isoforms']} dev. CYPs", p.GREY)
    route([(108, 30), (120, 30)], p.TEAL)

    document(120, 10, 49, 39, p.TEAL, "#F0F8F7")
    label(124, 13, "RANK SCORE", p.TEAL, True)
    label(124, 20, "Specific")
    label(143, 20, f"{chemical['score']:.3f}", p.TEAL, True, ha="right")
    label(147, 20, "Pooled")
    label(166, 20, f"{chemical['pooled_chemistry_score']:.3f}", bold=True, ha="right")
    ax.plot([124, 166], [25, 25], color="#C7DFDD", lw=.55)
    label(124, 27, "Exact structure seen: no")
    label(124, 31, "Original scaffold seen: no")
    label(124, 36, f"Max Tanimoto: {chemical['domain']['maximum_training_tanimoto']:.3f} · "
                    f"{len(chemical['nearest_labelled_neighbours'])} neighbours")
    label(124, 43, "Uncalibrated; within-isoform rank.", p.GREY)

    # Exact sequence/reaction lookup: both the match and no-match exits are visible.
    label(0, 62, "FASTA + reaction", bold=True)
    label(0, 69, example["evidence_input"]["sequence_prefix"] + "…")
    label(0, 74, example["evidence_input"]["reaction_key"][:13] + "…")
    label(0, 80, f"{example['evidence_input']['sequence_length']} aa · normalized reaction", p.GREY)
    route([(34, 76), (43, 76)], p.BLUE)
    ax.plot([47, 104], [52, 52], color="#CDD7DE", lw=.65)
    label(47, 57, "reaction-screen", p.BLUE, True)
    counts = verified["evidence_build"]
    label(47, 63, f"{counts['sequences']:,} sequences · {counts['reactions']:,} reactions")
    label(47, 68, f"{counts['edges']:,} documented edges", p.GREY)
    ax.add_patch(Polygon([(75, 74), (95, 82), (75, 90), (55, 82)], closed=True,
                         facecolor="white", edgecolor=p.BLUE, lw=.85))
    label(75, 82, "Exact edge?", p.BLUE, True, ha="center", va="center")
    route([(95, 82), (113, 82), (113, 65), (120, 65)], p.BLUE)
    label(103, 77, "match", p.BLUE, ha="center")
    route([(75, 90), (75, 97), (120, 97)], p.GOLD)
    label(97, 92, "no match", p.GOLD, ha="center")

    document(120, 55, 49, 21, p.BLUE, "#F1F5FA")
    label(124, 58, "DOCUMENTED EVIDENCE", p.BLUE, True)
    label(124, 64, evidence["accessions"][0] + " · " + evidence["sources"][0])
    label(124, 70, evidence["publication_ids"][0])
    document(120, 83, 49, 18, p.GOLD, "#FBF7EE")
    label(124, 85.5, "ABSTAIN", p.GOLD, True)
    label(124, 91, "score = null; evidence = []")
    label(124, 96, "Unknown, not a negative", p.GREY)

    # These are software checks, not biological accuracies or release cutoffs.
    ax.plot([0, 170], [106, 106], color="#ADBCC7", lw=.65)
    label(0, 109, "SOFTWARE", navy, True)
    label(0, 114, "VERIFICATION", navy, True)
    label(34, 109, f"{verified['external_rows_recomputed']:,} / {expected['external_scores_expected']:,} scores", bold=True)
    label(34, 114, "Reproduced · max |Δ| = 0", p.GREY)
    label(85, 109, f"{verified['exact_edge_keys_retrieved']:,} / {expected['exact_edge_keys_expected']:,} edge keys", bold=True)
    label(85, 114, "Retrieved from local index", p.GREY)
    label(135, 109, f"{len(controls)} abstention controls", bold=True)
    label(135, 114, "New sequence / reaction", p.GREY)

    # Enlarged preview only; submission PDF/SVG exports retain their A4 page.
    fig.savefig(p.OUT / "figure6_panelC_software.png", dpi=450,
                bbox_inches=Bbox.from_extents(19/25.4, (297-221)/25.4,
                                             191/25.4, (297-112)/25.4))
