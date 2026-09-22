#!/usr/bin/env python3
"""Data-bound A4 figures. All text is native Arial 7 pt; no model is refitted."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import PercentFormatter
from rdkit import Chem
from rdkit.Chem import rdDepictor
import figure_primitives as p
import figure6_software_layout

plt.rcParams["legend.frameon"] = False

ROOT, OUT = p.ROOT, p.OUT
BLUE, TEAL, GOLD, GREY, INK, LIGHT = p.BLUE, p.TEAL, p.GOLD, p.GREY, p.INK, p.LIGHT
text, panel, axes, page, arrow, finish = p.text, p.panel, p.axes, p.page, p.arrow, p.finish


def figure1(pdf):
    fig = page("Figure 1   Reaction identity determines how evidence connects across resources")
    panel(fig, "A", "One transformation, two complete identities", 20, 31)
    r = p.CASES["reaction_identity"]
    records = {a["source"]: a for a in r["source_assertions"]}
    full_by_key = {v["full_reaction_key"]: v for v in r["full_records"]}
    uni = full_by_key[records["UniProt_API"]["normalized_full_component_reaction_keys"][0]]
    p450 = full_by_key[records["P450Rdb"]["normalized_full_component_reaction_keys"][0]]
    assert uni["full_reaction"]["products"][1] != p450["full_reaction"]["products"][1]
    assert uni["reaction_key"] == p450["reaction_key"] == r["main_reaction_key"]

    # Compact schematic: common main chemistry above the one source-specific
    # representation difference. Detailed chemical qualifications stay in the legend.
    p.molecule(fig, r["main_transformation"]["substrates"][0], 25, 42, 29, 13)
    arrow(fig, (57, 49), (66, 49), INK)
    p.molecule(fig, r["main_transformation"]["products"][0], 69, 40, 32, 17)
    text(fig, 30, 59, "Decane", ha="center")
    text(fig, 85, 59, "Decan-3-ol", ha="center")
    text(fig, 58, 65, "shared main-transformation key", ha="center", weight="bold", color=TEAL)

    # Draw the actual isoalloxazine ring fragments from the two admitted FMN
    # encodings. The ribityl-phosphate side chain is omitted because it is
    # identical and would obscure the charge-localization difference.
    def fmn_core_smiles(record):
        mol = Chem.MolFromSmiles(record["full_reaction"]["products"][1])
        Chem.Kekulize(mol, clearAromaticFlags=True)
        # Retain the first side-chain carbon so the N substituent can be shown
        # as R rather than being incorrectly capped as N-H.
        return Chem.MolFragmentToSmiles(
            mol, atomsToUse=list(range(15)) + list(range(27, 31)),
            canonical=True, kekuleSmiles=True)

    fmn_fragments = [fmn_core_smiles(uni), fmn_core_smiles(p450)]
    assert fmn_fragments[0] != fmn_fragments[1]

    def charged_fragment(smiles, x, y, width, height, color):
        mol = Chem.MolFromSmiles(smiles)
        # Parsing restores aromatic bond types. Convert the *drawn* molecule,
        # not only the serialized precursor, to explicit alternating bonds.
        Chem.Kekulize(mol, clearAromaticFlags=True)
        assert all(not bond.GetIsAromatic() for bond in mol.GetBonds())
        assert sum(b.GetBondTypeAsDouble() == 2 for b in mol.GetBonds()) == 7
        rdDepictor.Compute2DCoords(mol)
        xy = np.array(mol.GetConformer().GetPositions())[:, :2]
        ax = axes(fig, x, y, width, height, off=True)
        for bond in mol.GetBonds():
            a, b = xy[bond.GetBeginAtomIdx()], xy[bond.GetEndAtomIdx()]
            if bond.GetBondTypeAsDouble() >= 2:
                v = b - a
                n = np.array([-v[1], v[0]]) / np.linalg.norm(v) * .055
                ax.plot([a[0]+n[0], b[0]+n[0]], [a[1]+n[1], b[1]+n[1]], color=INK, lw=.65)
                ax.plot([a[0]-n[0], b[0]-n[0]], [a[1]-n[1], b[1]-n[1]], color=INK, lw=.65)
            else:
                ax.plot([a[0], b[0]], [a[1], b[1]], color=INK, lw=.65)
        charged = []
        for atom in mol.GetAtoms():
            charge = atom.GetFormalCharge()
            r_group = (atom.GetSymbol() == "C" and atom.GetDegree() == 1
                       and atom.GetNeighbors()[0].GetSymbol() == "N")
            if atom.GetSymbol() == "C" and charge == 0 and not r_group:
                continue
            label = "R" if r_group else atom.GetSymbol()
            if atom.GetTotalNumHs() and not r_group:
                label += "H" + (str(atom.GetTotalNumHs()) if atom.GetTotalNumHs() > 1 else "")
            if charge:
                label += "(-)" if charge == -1 else f"({charge:+d})"
                charged.append(atom.GetIdx())
            label_color = color if charge else INK
            if charge:
                ax.scatter([xy[atom.GetIdx(),0]], [xy[atom.GetIdx(),1]], s=80,
                           facecolor=color, alpha=.14, edgecolor=color, linewidth=.7, zorder=2)
            ax.text(*xy[atom.GetIdx()], label, ha="center", va="center", color=label_color,
                    weight="bold" if charge else "normal",
                    bbox={"facecolor":"white", "edgecolor":"none", "pad":.25}, zorder=3)
        assert len(charged) == 1
        ax.set_aspect("equal")
        ax.set_xlim(xy[:,0].min()-.35, xy[:,0].max()+.35)
        ax.set_ylim(xy[:,1].min()-.42, xy[:,1].max()+.42)

    text(fig, 56, 69, "oxidized FMN ring core", ha="center", color=GREY)
    text(fig, 91, 69, "full key", ha="center", color=GREY)
    text(fig, 21, 78, "UniProt", color=BLUE, weight="bold")
    charged_fragment(fmn_fragments[0], 33, 72, 50, 17, BLUE)
    text(fig, 91, 78, "0793...", ha="center", color=BLUE, weight="bold")
    text(fig, 21, 97, "P450Rdb", color=GOLD, weight="bold")
    charged_fragment(fmn_fragments[1], 33, 91, 50, 17, GOLD)
    text(fig, 91, 97, "ed7e...", ha="center", color=GOLD, weight="bold")

    panel(fig, "B", "Identity rule changes shared linkage", 115, 31)
    ax = axes(fig, 127, 43, 62, 34, grid="x")
    labels = ["Complete components", "Main transformation"]
    vals = [31, 797]
    ypos = np.arange(2)
    ax.barh(ypos, vals, height=.48, color=[GREY, TEAL])
    ax.set_yticks(ypos, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 850)
    ax.set_xticks([0, 400, 800])
    ax.set_xlabel("shared cross-lineage edges")
    ax.tick_params(axis="y", length=0)
    for y, value, color in zip(ypos, vals, [GREY, TEAL]):
        if value > 100:
            ax.text(value - 14, y, f"{value:,}", ha="right", va="center",
                    color="white", weight="bold")
        else:
            ax.text(value + 14, y, f"{value:,}", va="center", color=color, weight="bold")

    panel(fig, "C", "Atlas union after role projection", 115, 88)
    u = p.D1["edge_union"]
    values = [u["p450_unique"], u["shared"], u["uniprot_unique"]]
    labels = ["P450Rdb only", "Shared", "UniProt\nonly"]
    colors = [GOLD, TEAL, BLUE]
    total = sum(values)
    assert total == u["union"] == 2304
    ax = axes(fig, 115, 99, 75, 30, off=True)
    ax.set_xlim(-12, total+12); ax.set_ylim(0, 1)
    left = 0
    for value, label, color in zip(values, labels, colors):
        ax.barh(.57, value, left=left, height=.30, color=color)
        ax.text(left + value/2, .57, f"{label}\n{value:,}", ha="center", va="center",
                color="white", weight="bold")
        left += value
    ax.plot([0, total], [.92, .92], color=INK, lw=.65)
    ax.plot([0, 0], [.88, .96], color=INK, lw=.65)
    ax.plot([total, total], [.88, .96], color=INK, lw=.65)
    ax.text(total/2, .96, "Union 2,304", ha="center", va="bottom", weight="bold")
    p450_total = values[0] + values[1]
    uniprot_start = values[0]
    ax.plot([0, p450_total], [.27, .27], color=GOLD, lw=.75)
    ax.plot([0, 0], [.23, .31], color=GOLD, lw=.75)
    ax.plot([p450_total, p450_total], [.23, .31], color=GOLD, lw=.75)
    ax.text(p450_total/2, .19, "P450Rdb 1,869", ha="center", va="top", color=GOLD, weight="bold")
    ax.plot([uniprot_start, total], [.05, .05], color=BLUE, lw=.75)
    ax.plot([uniprot_start, uniprot_start], [.01, .09], color=BLUE, lw=.75)
    ax.plot([total, total], [.01, .09], color=BLUE, lw=.75)
    ax.text((uniprot_start+total)/2, .01, "UniProt 1,232", ha="center", va="top", color=BLUE, weight="bold")

    finish(fig,"figure1_evidence_identity",3,134,
           ["figure1_evidence_atlas_source_data.json","chemical_case_studies.json","representation_count_provenance.json"],pdf)


def figure2(pdf):
    fig = page("Figure 2   Top-1 label transfer separates reachability from tie-dependent retrieval")
    panel(fig,"A","Grouped holdouts define three transfer regimes",20,31)
    text(fig,190,31,"145 protein groups × 240 chemical groups",ha="right",color=GREY)
    task_map={row["task"]:row for row in p.D2["tasks"]}
    matrix_specs=[
        (28,"Protein cold*","protein_cold"),
        (85,"Chemical cold","chemical_cold"),
        (142,"Double cold","double_cold"),
    ]
    for x0,name,key in matrix_specs:
        text(fig,x0+11,40,name,ha="center",weight="bold",color=TEAL if key=="protein_cold" else INK)
        ax=axes(fig,x0,47,22,22,off=True)
        for row in range(5):
            for col in range(5):
                test=(row==4 if key=="protein_cold" else col==4 if key=="chemical_cold" else row==4 and col==4)
                excluded=key=="double_cold" and (row==4 or col==4) and not test
                ax.add_patch(Rectangle((col,4-row),1,1,
                    facecolor=TEAL if test else "white" if excluded else LIGHT,
                    edgecolor=GREY if excluded else "white",lw=.45,
                    hatch="//" if excluded else None))
        ax.set_xlim(0,5);ax.set_ylim(0,5)
        text(fig,x0+11,71,f"{task_map[key]['blocks']} outer blocks",ha="center",color=GREY)

    legend_ax=axes(fig,20,77,66,7,off=True)
    for x0,face,hatch,label in [(0,LIGHT,None,"training"),(.34,TEAL,None,"held out"),(.67,"white","//","excluded")]:
        legend_ax.add_patch(Rectangle((x0,.25),.08,.42,facecolor=face,edgecolor=GREY,lw=.45,hatch=hatch))
        legend_ax.text(x0+.10,.46,label,va="center")
    legend_ax.set_xlim(0,1);legend_ax.set_ylim(0,1)
    text(fig,190,78,"rows: protein folds · columns: chemical folds",ha="right",color=GREY)
    text(fig,190,84,"* B–C: protein cold, seen chemistry; publication overlap purged",ha="right",color=TEAL)

    panel(fig,"B","Reachability and fixed-tie retrieval decompose end-to-end MRR",20,94)
    rows={r["method"]:r for r in p.tsv("coverage_ranking_decomposition.tsv") if r["task"]=="protein_cold_seen"}
    mm=rows["mmseqs_top1"]; esm=rows["esm_cosine_top1"]
    assert int(mm["query_panels"])==int(esm["query_panels"])==95
    text(fig,190,94,"95 query panels · 27 groups · 1,341 candidates",ha="right",color=GREY)
    x_mm=float(mm["weighted_query_positive_coverage_Cw"])
    y_mm=float(mm["covered_weighted_MRR"])
    z_mm=float(mm["end_to_end_MRR"])
    x_esm=float(esm["weighted_query_positive_coverage_Cw"])
    y_esm=float(esm["covered_weighted_MRR"])
    z_esm=float(esm["end_to_end_MRR"])

    # Restore the coverage-versus-ranking view: each grey curve is a constant
    # end-to-end MRR under MRR = Cw x conditional MRR.
    ax=axes(fig,30,108,105,40,grid=False)
    x_curve=np.linspace(.03,.13,500)
    for level in (.04,.06,.08):
        y_curve=level/x_curve
        valid=(y_curve>=.20)&(y_curve<=1.0)
        ax.plot(x_curve[valid],y_curve[valid],color=LIGHT,lw=1.0,zorder=1)
        ax.text(.128,level/.128,f"{level:.2f}",ha="right",va="bottom",color=GREY)
    ax.scatter(x_mm,y_mm,s=29,color=BLUE,marker="o",zorder=3)
    ax.scatter(x_esm,y_esm,s=31,color=TEAL,marker="s",zorder=3)
    ax.annotate("MMseqs2\n7.91%, 0.8706",(x_mm,y_mm),xytext=(-7,7),
                textcoords="offset points",ha="right",va="bottom",color=BLUE,weight="bold")
    ax.annotate("ESM\n10.85%, 0.7443",(x_esm,y_esm),xytext=(-8,-6),
                textcoords="offset points",ha="right",va="top",color=TEAL,weight="bold")
    ax.set_xlim(0,.13);ax.set_ylim(.20,1.0)
    ax.set_xticks([0,.04,.08,.12]);ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0))
    ax.set_yticks([.2,.4,.6,.8,1.0])
    ax.set_xlabel("Reachable-query weight, Cw")
    ax.set_ylabel("Conditional MRR")
    ax.grid(color=LIGHT,lw=.5)
    text(fig,132,107,"grey curves: end-to-end MRR",ha="right",color=GREY)

    text(fig,150,109,"End-to-end MRR",weight="bold")
    key_ax=axes(fig,146,117,42,27,off=True)
    for yv,name,value,color,marker in [(.72,"MMseqs2",z_mm,BLUE,"o"),(.28,"ESM",z_esm,TEAL,"s")]:
        key_ax.scatter(.06,yv,s=24,color=color,marker=marker)
        key_ax.text(.15,yv,name,va="center",color=color,weight="bold")
        key_ax.text(.98,yv,f"{value:.4f}",ha="right",va="center",weight="bold")
    key_ax.plot([.15,.98],[.50,.50],color=LIGHT,lw=.6)
    key_ax.set_xlim(0,1);key_ax.set_ylim(0,1)
    text(fig,188,146,f"Δ = {z_esm-z_mm:+.4f}",ha="right",color=GREY)
    assert np.isclose(
        float(mm["weighted_query_positive_coverage_Cw"]) * float(mm["covered_weighted_MRR"]),
        float(mm["end_to_end_MRR"]),
        rtol=1e-12,
        atol=1e-15,
    )
    assert np.isclose(
        float(esm["weighted_query_positive_coverage_Cw"]) * float(esm["covered_weighted_MRR"]),
        float(esm["end_to_end_MRR"]),
        rtol=1e-12,
        atol=1e-15,
    )

    panel(fig,"C","Shared candidates have constant scores for both routes",20,160)
    common=p.tsv("common_domain_key_comparison.tsv")[0]
    nq=int(common["common_domain_query_panels"]); tq=int(common["common_domain_query_panels"])+int(common["common_undefined_query_panels"])
    np_common=int(common["common_positive_instances"]); np_total=int(common["original_target_positive_instances"])
    assert (nq,tq,np_common,np_total)==(3,95,6,277)
    esm_total=int(common["left_covered_positive_instances"])
    mm_total=int(common["right_covered_positive_instances"])
    esm_only=esm_total-np_common
    mm_only=mm_total-np_common
    neither=np_total-(esm_only+np_common+mm_only)
    assert (esm_only,np_common,mm_only,neither)==(9,6,3,259)

    text(fig,30,171,"Positive-instance support",weight="bold")
    text(fig,111,171,f"18/{np_total} covered by ≥1 · {neither} by neither",ha="right",color=GREY)
    ax=axes(fig,40,178,68,18,grid="y")
    overlap_counts=[esm_only,np_common,mm_only]
    xpos=np.arange(3)
    ax.bar(xpos,overlap_counts,width=.56,color=[TEAL,INK,BLUE])
    ax.set_xlim(-.6,2.6);ax.set_ylim(0,10.5);ax.set_yticks([0,5,10])
    ax.set_ylabel("positive instances")
    ax.set_xticks([])
    for x0,value in zip(xpos,overlap_counts):
        ax.text(x0,value+.35,str(value),ha="center",va="bottom",weight="bold")

    membership=axes(fig,40,198,68,12,off=True)
    membership.set_xlim(-.6,2.6);membership.set_ylim(-.6,1.6)
    for x0,(has_esm,has_mm) in enumerate([(1,0),(1,1),(0,1)]):
        if has_esm and has_mm:
            membership.plot([x0,x0],[0,1],color=INK,lw=1.0)
        for y0,active,color in [(1,has_esm,TEAL),(0,has_mm,BLUE)]:
            membership.scatter(x0,y0,s=24,facecolor=color if active else "white",
                               edgecolor=color if active else GREY,lw=.7,zorder=3)
    membership.text(-.72,1,"ESM",ha="right",va="center",color=TEAL,weight="bold")
    membership.text(-.72,0,"MMseqs2",ha="right",va="center",color=BLUE,weight="bold")
    membership.text(0,-.48,"only",ha="center",va="top",color=GREY)
    membership.text(1,-.48,"both",ha="center",va="top",color=GREY)
    membership.text(2,-.48,"only",ha="center",va="top",color=GREY)

    arrow(fig,(113,193),(123,193),GREY)

    text(fig,158,171,f"Only {nq}/{tq} panels · 3 groups",ha="center",weight="bold")
    audit=json.loads((ROOT.parent/'Source_data/analysis_tables/figure2_tie_audit.json').read_text())
    shared=[r for r in audit['eligible_common_panels'] if r['method']=='mmseqs_top1']
    text(fig,146,177,"Targets /\ncandidates",ha="right",color=GREY)
    text(fig,167,177,"Fixed-tie\nRR",ha="right",color=GREY)
    text(fig,187,177,"Tie-mean\nRR",ha="right",color=GREY)
    for i,row in enumerate(shared):
        y=187+i*6
        text(fig,126,y,str(i+1),ha="right",color=GREY)
        text(fig,146,y,f"{row['positives']} / {row['candidates']}",ha="right")
        text(fig,167,y,f"{row['deterministic_RR']:.3f}",ha="right")
        text(fig,187,y,f"{row['tie_averaged_RR']:.3f}",ha="right")
    text(fig,146,206,"Both MRRs",ha="right",weight="bold")
    text(fig,167,206,f"{float(common['common_domain_right_MRR']):.3f}",ha="right",weight="bold")
    mean=audit['tie_averaged']['common_domain'][0]['left_conditional_MRR']
    text(fig,187,206,f"{mean:.3f}",ha="right",weight="bold")

    finish(fig,"figure2_coverage_ranking",3,215,
           ["figure2_split_design_source_data.json","coverage_ranking_decomposition.tsv","common_domain_key_comparison.tsv"],pdf)


def figure3(pdf):
    fig=page("Figure 3   Protein-side gains contract after availability and exposure controls")

    # A is the methodological bridge: the representation counts are retained,
    # but the routes into the two downstream controls are now explicit.
    panel(fig,"A","Protein evidence is filtered before each control",20,31)
    text(fig,79,39,"available representation",ha="center",color=GREY)
    text(fig,143,39,"analysis route",ha="center",color=GREY)

    def route_node(x0,y0,w,h,title,subtitle,color):
        ax=axes(fig,x0,y0,w,h,off=True)
        ax.add_patch(Rectangle((0,0),1,1,facecolor="white",edgecolor=GREY,lw=.55))
        ax.add_patch(Rectangle((0,0),.045,1,facecolor=color,edgecolor="none"))
        ax.text(.10,.64,title,va="center",weight="bold")
        ax.text(.10,.27,subtitle,va="center",color=GREY)
        ax.set_xlim(0,1);ax.set_ylim(0,1)

    route_node(21,49,26,29,"Core CYP set","n = 600",INK)
    route_node(59,42,43,12,"Global sequence","ESM 600 · CLEAN 582",BLUE)
    route_node(59,58,43,12,"Local sequence","site vectors 289",TEAL)
    route_node(59,74,43,12,"Structure","36 models · 23 site vectors",GREY)
    route_node(116,42,56,12,"Exposure-stratified pruning","exact training overlap",GOLD)
    route_node(116,58,56,12,"Interaction robustness","66 panels · 4 protein groups",TEAL)
    route_node(116,74,72,12,"Matched retrieval: Foldseek = MMseqs2","14 queries · 11 groups · MRR 0.0417",GREY)

    for target_y in (48,64,80):
        arrow(fig,(47,63.5),(57,target_y),GREY)
    for source_y,target_y in ((48,48),(64,64),(80,80)):
        arrow(fig,(103,source_y),(114,target_y),GREY)
    text(fig,184,44,"B",ha="center",weight="bold",color=GOLD)
    arrow(fig,(173,48),(180,48),GOLD)
    text(fig,184,60,"C",ha="center",weight="bold",color=TEAL)
    arrow(fig,(173,64),(180,64),TEAL)

    panel(fig,"B","CLEAN gain is concentrated in training-exposed proteins",20,105)
    text(fig,190,105,"ΔMRR = pruned − unpruned · conditional 95% group bootstrap",ha="right",color=GREY)
    ax=axes(fig,47,121,112,42)
    for i,row in enumerate(p.D5["pruning_effects"]):
        for shift,key,color,marker,face in [(-.20,"all",GREY,"o",GREY),(0,"exposed",GOLD,"s",GOLD),(.20,"unexposed",BLUE,"o","white")]:
            v=row[key];d=v["difference"];lo,hi=v["ci"]
            ax.errorbar(d,i+shift,xerr=[[d-lo],[hi-d]],fmt=marker,color=color,mfc=face,
                        ms=3.5,capsize=1.8,lw=.75,zorder=3)
    ax.axvline(0,color=GREY,ls="--",lw=.65)
    ax.set_xlim(-.08,.48);ax.set_xticks([0,.1,.2,.3,.4]);ax.set_ylim(2.55,-.55)
    ax.set_yticks(range(3),["Chemical cold","Protein cold","Double cold"])
    ax.tick_params(axis="y",length=0,pad=5)
    ax.set_xlabel("Change in MRR")
    handles=[Line2D([],[],ls="",marker=m,color=c,mfc=face,ms=3.5,label=label) for m,c,face,label in
             [("o",GREY,GREY,"All"),("s",GOLD,GOLD,"Exact training overlap"),("o",BLUE,"white","No exact training overlap")]]
    ax.legend(handles=handles,loc="upper left",bbox_to_anchor=(0,-.23),ncol=3,
              columnspacing=1.2,handletextpad=.4)
    text(fig,163,120,"Exact\noverlap",ha="center",weight="bold",color=GOLD)
    text(fig,190,120,"No exact\noverlap",ha="right",weight="bold",color=BLUE)
    for i,row in enumerate(p.D5["pruning_effects"]):
        text(fig,162,131+i*11,f"{row['exposed']['difference']:+.3f}",color=GOLD)
        text(fig,190,131+i*11,f"{row['unexposed']['difference']:+.3f}",ha="right",color=BLUE)

    panel(fig,"C","The interaction estimate depends on four protein groups",20,184)
    text(fig,190,184,"66 query panels · 4 groups · fixed predictions",ha="right",color=GREY)
    groups=p.tsv("interaction_four_group_sensitivity.tsv");loo=p.tsv("interaction_leave_one_group_out.tsv")
    assert len(groups)==len(loo)==4
    ax=axes(fig,47,203,107,38)
    for i,(r,l) in enumerate(zip(groups,loo)):
        v=float(r["paired_difference"]);o=float(l["fixed_prediction_mean_difference"])
        ax.plot([v,o],[i,i],color=LIGHT,lw=1.3,zorder=1)
        ax.plot(v,i,"o",color=BLUE,ms=4,zorder=3)
        ax.plot(o,i,"D",mfc="white",color=GREY,ms=3.5,zorder=3)
    ax.axvline(0,color=GREY,ls="--",lw=.65)
    ax.set_xlim(-.05,.07);ax.set_xticks([-.04,-.02,0,.02,.04,.06]);ax.set_ylim(3.55,-.55)
    ax.set_yticks(range(4),[f"Group {i+1}" for i in range(4)]);ax.tick_params(axis="y",length=0,pad=5)
    ax.set_xlabel("Joint interaction minus global MRR")
    text(fig,161,197,"Group effect",weight="bold",color=BLUE)
    text(fig,190,197,"Omit group",ha="right",weight="bold",color=GREY)
    for i,(r,l) in enumerate(zip(groups,loo)):
        text(fig,161,208+i*9.25,f"{float(r['paired_difference']):+.4f}",color=BLUE)
        text(fig,190,208+i*9.25,f"{float(l['fixed_prediction_mean_difference']):+.4f}",ha="right",color=GREY)
    legend_ax=axes(fig,47,192,96,6,off=True)
    legend_ax.plot(.02,.5,"o",color=BLUE,ms=4)
    legend_ax.text(.05,.5,"group effect",va="center")
    legend_ax.plot(.35,.5,"D",mfc="white",color=GREY,ms=3.5)
    legend_ax.text(.38,.5,"mean after omitting group",va="center")
    legend_ax.set_xlim(0,1);legend_ax.set_ylim(0,1)

    # Shift the lower controls upward without changing type or mark sizes.
    for artist in fig.texts:
        x,y=artist.get_position()
        if (1-y)*297 >= 100:
            artist.set_position((x,y+10/297))
    for ax in fig.axes:
        pos=ax.get_position()
        if (1-pos.y1)*297 >= 100:
            ax.set_position([pos.x0,pos.y0+10/297,pos.width,pos.height])
    finish(fig,"figure3_protein_controls",3,242,
           ["figure5_coverage_controls_source_data.json","interaction_four_group_sensitivity.tsv","interaction_leave_one_group_out.tsv"],pdf)


def figure4(pdf):
    fig=page("Figure 4   Reverse retrieval shows signal but does not meet all study criteria")
    d=p.D4; census=d["census"]; rows=d["performance"]
    assert (census["admitted_candidates"],census["query_panels"],census["positives"])==(1803,4880,6033)

    # A defines the practical task before the paper switches from annotation
    # (protein -> reaction) to discovery (reaction + taxid -> protein ranking).
    panel(fig,"A","Reaction–taxid queries rank CYP candidates within species",20,31)
    text(fig,190,31,"1,803 candidates · 4,880 panels · 6,033 documented positives",ha="right",color=GREY)

    def task_box(x0,y0,w,h,title,subtitle,color):
        ax=axes(fig,x0,y0,w,h,off=True)
        ax.add_patch(Rectangle((0,0),1,1,facecolor="white",edgecolor=GREY,lw=.55))
        ax.add_patch(Rectangle((0,0),.045,1,facecolor=color,edgecolor="none"))
        ax.text(.10,.64,title,va="center",weight="bold")
        ax.text(.10,.26,subtitle,va="center",color=GREY)
        ax.set_xlim(0,1);ax.set_ylim(0,1)

    task_box(22,43,29,11,"Reaction","r",TEAL)
    task_box(22,59,29,11,"Species","taxid t",BLUE)
    candidate_ax=axes(fig,64,43,43,27,off=True)
    candidate_ax.add_patch(Rectangle((0,0),1,1,facecolor="white",edgecolor=GREY,lw=.55))
    candidate_ax.text(.10,.80,"Same-taxid CYP set",weight="bold",va="center")
    candidate_ax.text(.10,.59,"2–141 candidates",color=GREY,va="center")
    for i,width in enumerate((.56,.72,.48,.64,.40)):
        candidate_ax.add_patch(Rectangle((.11,.38-i*.065),width,.032,facecolor=LIGHT,edgecolor=GREY,lw=.25))
    candidate_ax.set_xlim(0,1);candidate_ax.set_ylim(0,1)

    score_ax=axes(fig,120,43,31,27,off=True)
    score_ax.add_patch(Rectangle((0,0),1,1,facecolor="white",edgecolor=GREY,lw=.55))
    score_ax.text(.10,.80,"Score candidates",weight="bold",va="center")
    for i,(name,color,marker) in enumerate((("MMseqs2",TEAL,"o"),("Label prior",GOLD,"s"),("Uniform",GREY,"D"))):
        yv=.57-i*.18
        score_ax.plot(.13,yv,marker=marker,color=color,mfc="white" if name=="Uniform" else color,ms=3)
        score_ax.text(.24,yv,name,va="center")
    score_ax.set_xlim(0,1);score_ax.set_ylim(0,1)

    rank_ax=axes(fig,164,43,25,27,off=True)
    rank_ax.add_patch(Rectangle((0,0),1,1,facecolor="white",edgecolor=GREY,lw=.55,clip_on=False))
    rank_ax.text(.10,.80,"Ranked CYPs",weight="bold",va="center")
    for i,(label,is_positive) in enumerate((("1",False),("2",True),("3",False),("…",False))):
        yv=.59-i*.14
        rank_ax.text(.12,yv,label,va="center",ha="center",color=GREY)
        rank_ax.plot([.23,.70],[yv,yv],color=TEAL if is_positive else LIGHT,lw=2.2)
        if is_positive:
            rank_ax.plot(.82,yv,"s",color=GOLD,ms=3)
    rank_ax.set_xlim(0,1);rank_ax.set_ylim(0,1)

    arrow(fig,(51,48.5),(62,54),GREY);arrow(fig,(51,64.5),(62,59),GREY)
    arrow(fig,(108,56.5),(118,56.5),GREY);arrow(fig,(152,56.5),(162,56.5),GREY)
    text(fig,64,75,"candidate universe: source-observed CYPs",color=GREY)
    text(fig,190,75,"RR = 1 / rank(first documented positive)",ha="right",weight="bold")
    text(fig,190,81,"macro-average over taxids",ha="right",color=GREY)

    panel(fig,"B","Ranking and reachability differ across evidence regimes",20,94)
    text(fig,190,94,"same candidates, positives and denominators across methods",ha="right",color=GREY)
    perf_ax=axes(fig,54,112,80,43)
    for i,row in enumerate(rows):
        values=[row["mrr"][key] for key in ("MMseqs2","Label prior","Uniform")]
        perf_ax.plot([min(values),max(values)],[i,i],color=LIGHT,lw=1.2,zorder=1)
    for shift,key,color,marker,face in [(-.16,"MMseqs2",TEAL,"o",TEAL),(0,"Label prior",GOLD,"s",GOLD),(.16,"Uniform",GREY,"D","white")]:
        perf_ax.plot([r["mrr"][key] for r in rows],np.arange(3)+shift,linestyle="none",
                     marker=marker,color=color,mfc=face,ms=3.8,label=key,zorder=3)
    perf_ax.set_yticks(range(3),[r["label"] for r in rows]);perf_ax.tick_params(axis="y",length=0,pad=5)
    perf_ax.set_ylim(2.5,-.5);perf_ax.set_xlim(0,.65);perf_ax.set_xticks([0,.2,.4,.6])
    perf_ax.set_xlabel("Taxid-macro first-positive MRR")
    perf_ax.legend(loc="upper left",bbox_to_anchor=(0,-.23),ncol=3,columnspacing=1.2,handletextpad=.4)

    text(fig,160,105,"MMseqs2 positive coverage",ha="center",weight="bold")
    cov_ax=axes(fig,147,112,28,43,grid="x")
    coverage=[r["coverage"]["MMseqs2"] for r in rows]
    cov_ax.barh(range(3),coverage,height=.40,color=TEAL)
    cov_ax.set_xlim(0,1);cov_ax.set_ylim(2.5,-.5);cov_ax.set_xticks([0,.5,1])
    cov_ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0));cov_ax.set_yticks([])
    for i,value in enumerate(coverage):
        cov_ax.text(value-.025,i,f"{value:.1%}",ha="right",va="center",color="white",weight="bold")
    text(fig,188,105,"Taxids",ha="right",weight="bold")
    for i,row in enumerate(rows):
        text(fig,188,119+i*14.35,str(row["taxids"]),ha="right")

    panel(fig,"C","No regime meets all study criteria",20,174)
    text(fig,190,174,"all five criteria required",ha="right",color=GREY)
    labels=["≥20\ntaxids","≥80%\ncoverage","95% CI > 0\nvs label prior","95% CI > 0\nvs uniform","Permutation\np ≤ 0.05"]
    matrix=axes(fig,69,195,93,27,off=True)
    matrix.set_xlim(-.6,4.6);matrix.set_ylim(2.6,-.6)
    for i,flags in enumerate(d["gates"]):
        matrix.plot([-.45,4.45],[i,i],color=LIGHT,lw=.55,zorder=0)
        for j,passed in enumerate(flags):
            matrix.plot(j,i,"o" if passed else "x",mfc=TEAL if passed else "white",
                        color=TEAL if passed else GREY,ms=4.5,mew=.9,zorder=3)
    for j,label in enumerate(labels):
        text(fig,76.8+j*17.9,184,label,ha="center")
    row_labels=["Seen / new","Unseen / new","Unseen / represented"]
    for i,label in enumerate(row_labels):
        text(fig,20,198+i*8.45,label)
    text(fig,173,184,"Met",ha="center",weight="bold")
    text(fig,190,184,"All study\ncriteria",ha="right",weight="bold")
    for i,flags in enumerate(d["gates"]):
        text(fig,173,198+i*8.45,f"{sum(flags)}/5",ha="center",weight="bold")
        text(fig,190,198+i*8.45,"Not met",ha="right",color=GREY)
    legend=axes(fig,69,225,93,7,off=True)
    legend.plot(.03,.5,"o",color=TEAL,mfc=TEAL,ms=4.5);legend.text(.07,.5,"criterion met",va="center")
    legend.plot(.36,.5,"x",color=GREY,ms=4.5,mew=.9);legend.text(.40,.5,"not met",va="center")
    legend.set_xlim(0,1);legend.set_ylim(0,1)

    finish(fig,"figure4_reverse_retrieval",3,235,["figure4_reverse_retrieval_source_data.json"],pdf)


def figure5(pdf):
    fig=page("Figure 5   Isoform information improves ranking on qualified external compounds")
    panel(fig,"A","External eligibility is determined by compound and source provenance",20,23)
    ax=axes(fig,72,36,100,24)
    for i,row in enumerate(p.D6["stages"]):
        ax.barh(i,row["rows"],height=.46,color=TEAL if i==3 else LIGHT,edgecolor=TEAL if i==3 else GREY,lw=.4)
        ax.text(row["rows"]+160,i,f"{row['rows']:,}",va="center")
    ax.set_yticks(range(4),["Evaluable pairs","After exact-identity exclusion","After source-lineage exclusion","After original scaffold exclusion"])
    ax.set_xlim(0,16200);ax.set_xticks([0,5000,10000,15000],["0","5,000","10,000","15,000"]);ax.set_ylim(3.6,-.6)
    ax.set_xlabel("Isoform–compound labels")
    panel(fig,"B","Conditioning and model choice are separate comparisons",20,82)
    isoforms=p.D6["strict_per_isoform"]
    logistic={r['isoform']:float(r['strict_external_AP']) for r in p.tsv('logistic_regression_supplement.tsv')}
    matrix=np.array([[row["positives"]/row["rows"],row["pooled_ap"],row["specific_ap"],logistic[row['isoform']]]
                     for row in isoforms])
    cmap=p.LinearSegmentedColormap.from_list("external_performance",["#F7F9FA","#CFE4E6",TEAL])
    hm=axes(fig,51,99,90,55,grid=None)
    hm.imshow(matrix,cmap=cmap,vmin=0,vmax=1,aspect="auto",interpolation="nearest")
    hm.set_xticks(range(4),["Positive\nfraction","Pooled\nkNN AP","Isoform-specific\nkNN AP","Isoform\nlogistic AP*"])
    hm.tick_params(axis="x",top=True,labeltop=True,bottom=False,labelbottom=False,length=0,pad=3)
    hm.set_yticks(range(6),[f"{row['isoform']}   (n={row['rows']})" for row in isoforms])
    hm.tick_params(axis="y",length=0,pad=4)
    hm.spines[:].set_visible(False)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value=matrix[i,j]
            hm.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor="none",edgecolor="white",lw=1.1))
            hm.text(j,i,f"{value:.3f}",ha="center",va="center",
                    color="white" if value>=.58 else INK,
                    weight="normal")
    gain=axes(fig,156,99,26,55,grid="x")
    gains=np.array([row["specific_ap"]-row["pooled_ap"] for row in isoforms])
    gain.barh(range(6),gains,height=.52,color=TEAL)
    for i,value in enumerate(gains):
        gain.text(value+.008,i,f"+{value:.3f}",va="center",ha="left",color=INK)
    gain.set_xlim(0,.35);gain.set_xticks([0,.1,.2,.3],["0",".1",".2",".3"])
    gain.set_ylim(5.5,-.5);gain.set_yticks([]);gain.set_xlabel("kNN AP gain\nSpecific − pooled")
    gain.spines["left"].set_visible(False)
    scale=axes(fig,51,160,90,3,off=True)
    gradient=np.linspace(0,1,256)[None,:]
    scale.imshow(gradient,aspect="auto",cmap=cmap,extent=(0,1,0,1))
    scale.text(0,-.75,"0",ha="center",va="top");scale.text(1,-.75,"1",ha="center",va="top")
    scale.text(.5,-.75,"cell value",ha="center",va="top")

    text(fig,51,173,"*Exploratory comparator; fitted and tuned on development data only.",color=GREY)
    panel(fig,"C","Macro-AP gain persists after original scaffold exclusion",20,183)
    aggregate=p.D6["aggregate"]
    bridge=axes(fig,68,201,81,25,grid="x")
    labels=["Identity + lineage exclusions","After original scaffold exclusion"]
    for i,row in enumerate(aggregate):
        pooled,specific=row["pooled_ap"],row["specific_ap"]
        delta=row["ap_difference"]
        bridge.barh(i,pooled,height=.48,color=GREY)
        bridge.barh(i,delta,left=pooled,height=.48,color=TEAL)
        bridge.text(pooled/2,i,f"{pooled:.3f}",ha="center",va="center",color="white")
        bridge.text(pooled+delta/2,i,f"+{delta:.3f}",ha="center",va="center",color="white",weight="bold")
        bridge.plot([specific,specific],[i-.31,i+.31],color=INK,lw=.8)
    bridge.set_xlim(0,.60);bridge.set_xticks([0,.2,.4,.6]);bridge.set_ylim(1.55,-.55)
    bridge.set_yticks(range(2),labels);bridge.tick_params(axis="y",length=0,pad=4)
    bridge.set_xlabel("Macro average precision")
    handles=[Rectangle((0,0),1,1,facecolor=GREY,edgecolor="none",label="Pooled AP"),
             Rectangle((0,0),1,1,facecolor=TEAL,edgecolor="none",label="Isoform-specific gain")]
    bridge.legend(handles=handles,loc="lower left",bbox_to_anchor=(-.02,1.06),ncol=2,
                  columnspacing=1.2,handlelength=1.2,handletextpad=.4)
    table=axes(fig,153,201,34,25,off=True)
    table.set_xlim(0,1);table.set_ylim(1.55,-.55)
    table.text(.00,-.68,"ΔAP [95% CI]",weight="bold",va="bottom")
    table.text(1.00,-.68,"Permutation P",weight="bold",ha="right",va="bottom")
    for i,row in enumerate(aggregate):
        lo,hi=row["ap_ci"]
        table.text(.00,i,f"{row['ap_difference']:.3f} [{lo:.3f}, {hi:.3f}]",va="center")
        table.text(1.00,i,f"{row['permutation_p']:.2f}",ha="right",va="center")
    text(fig,20,240,"Strict stratum: 3,035 labels; 1,001 compounds; 844 scaffolds. Paired scaffold bootstrap: 5,000 replicates.",color=GREY)
    finish(fig,"figure5_external_transfer",3,248,["figure6_external_transfer_source_data.json","logistic_regression_supplement.tsv"],pdf)


def figure6(pdf):
    fig=page("Figure 6   Selective release exposes the tradeoff between error and recovered positives")
    data=p.tsv("selective_release_tradeoffs.tsv")
    pooled=sorted([r for r in data if r["method"]=="pooled_chemical_knn"],key=lambda r:float(r["micro_coverage"]))
    specific=sorted([r for r in data if r["method"]=="isoform_specific_knn"],key=lambda r:float(r["micro_coverage"]))

    release=np.array([10,25,50,75,100],dtype=float)
    retained=[306,761,1519,2279,3035]
    risk_pooled=np.array([1-float(row["macro_precision"]) for row in pooled])*100
    risk_specific=np.array([1-float(row["macro_precision"]) for row in specific])*100
    recall_pooled=np.array([float(row["macro_positive_recall"]) for row in pooled])*100
    recall_specific=np.array([float(row["macro_positive_recall"]) for row in specific])*100
    tick_labels=[f"{int(value)}%\n{count:,}" for value,count in zip(release,retained)]

    method_handles=[
        Line2D([0],[0],color=GREY,lw=.9,ls=(0,(2.2,1.6)),marker="o",markersize=4,
               markerfacecolor="white",markeredgecolor=GREY,label="Pooled chemistry"),
        Line2D([0],[0],color=TEAL,lw=1.15,marker="s",markersize=3.6,
               markerfacecolor=TEAL,markeredgecolor=TEAL,label="Isoform-specific"),
    ]
    fig.legend(handles=method_handles,loc="center",bbox_to_anchor=(.50,1-30.5/297),
               ncol=2,columnspacing=2.1,handlelength=2.4,handletextpad=.6)

    def tradeoff_axis(x,y,width,height,pooled_values,specific_values,ylabel):
        ax=axes(fig,x,y,width,height,grid="y")
        # The narrow band marks the prespecified low-release operating point.
        ax.axvspan(7.5,12.5,color="#E7F3F3",zorder=0)
        for j,(pooled_value,specific_value) in enumerate(zip(pooled_values,specific_values)):
            ax.plot([release[j],release[j]],[pooled_value,specific_value],
                    color=GOLD if j==0 else "#D9E0E4",lw=1.0 if j==0 else .65,
                    zorder=1)
        ax.plot(release,pooled_values,color=GREY,lw=.9,ls=(0,(2.2,1.6)),
                marker="o",markersize=4,markerfacecolor="white",markeredgecolor=GREY,
                markeredgewidth=.85,zorder=2)
        ax.plot(release,specific_values,color=TEAL,lw=1.15,
                marker="s",markersize=3.6,markerfacecolor=TEAL,markeredgecolor=TEAL,
                zorder=3)
        ax.set_xlim(6,104);ax.set_ylim(0,104)
        ax.set_xticks(release,tick_labels)
        ax.set_yticks([0,25,50,75,100])
        ax.yaxis.set_major_formatter(PercentFormatter(100,decimals=0))
        ax.set_xlabel("Release fraction\n(labels retained)")
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x",pad=3)
        return ax

    panel(fig,"A","Selective risk (1 - precision)",20,23)
    risk_ax=tradeoff_axis(35,38,63,58,risk_pooled,risk_specific,"Selective risk (%)")
    # Show the full paired estimate at the highlighted 10% operating point.
    risk_ax.annotate("-26.4 pp\n95% CI -32.7 to -20.1",
                     xy=(10,(risk_pooled[0]+risk_specific[0])/2),xytext=(18,47),
                     ha="left",va="center",color=INK,
                     arrowprops={"arrowstyle":"-","color":GOLD,"lw":.7})
    risk_ax.text(10,risk_pooled[0]+4.0,f"{risk_pooled[0]:.1f}",ha="center",va="bottom",color=GREY)
    risk_ax.text(10,risk_specific[0]-4.0,f"{risk_specific[0]:.1f}",ha="center",va="top",color=TEAL,weight="bold")
    risk_ax.text(100,risk_pooled[-1]+4.0,f"{risk_pooled[-1]:.1f}",ha="right",va="bottom",color=INK)

    panel(fig,"B","Positive-label recall",111,23)
    recall_ax=tradeoff_axis(124,38,65,58,recall_pooled,recall_specific,"Positive-label recall (%)")
    recall_ax.text(10,recall_pooled[0]-4.0,f"{recall_pooled[0]:.1f}",ha="center",va="top",color=GREY)
    recall_ax.text(10,recall_specific[0]+4.0,f"{recall_specific[0]:.1f}",ha="center",va="bottom",color=TEAL,weight="bold")
    recall_ax.text(100,recall_pooled[-1]-4.0,f"{recall_pooled[-1]:.1f}",ha="right",va="top",color=INK)

    figure6_software_layout.draw(fig)
    text(fig,20,225,"A–B: 3,035 fixed-human labels · 1,001 compounds · 844 scaffolds · five evaluated release fractions.",color=GREY)
    text(fig,20,230,"Macro precision × total released is not a hit count; the 10% interval is a paired 95% scaffold bootstrap.",color=GREY)
    finish(fig,"figure6_selective_strategy",3,239,
           ["selective_release_tradeoffs.tsv","figure5_coverage_controls_source_data.json","strategy_scope.json",
            "figure6_runtime_example.json"],pdf)


def supplement4(pdf):
    original_finish=p.finish;original_page=p.page
    p.page=lambda title:original_page("Figure S4   Detailed local-site, interaction and exposure controls")
    p.finish=lambda fig,stem,panels,bottom,inputs,book:original_finish(fig,"figureS4_detailed_protein_controls",panels,bottom,inputs,book)
    try:p.figure3(pdf)
    finally:p.finish=original_finish;p.page=original_page


def main():
    with PdfPages(OUT/"CYP_TRACE_Main_Figures_A4.pdf") as pdf:
        for f in [figure1,figure2,figure3,figure4,figure5,figure6]:f(pdf)
    with PdfPages(OUT/"CYP_TRACE_Supplementary_Figures_A4.pdf") as pdf:
        for f in [p.supplement1,p.supplement2,p.supplement3,supplement4]:f(pdf)
    report={"status":"BUILT_PENDING_VISUAL_REVIEW","figure_count":len(p.RECORDS),"font_file":p.FONT,
            "font_pt":7,"page_mm":[210,297],"figures":p.RECORDS,
            "input_sha256":{x.name:hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(p.INPUT.iterdir()) if x.is_file()},
            "source_code_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "figure2_tie_audit_sha256":hashlib.sha256((ROOT.parent/'Source_data/analysis_tables/figure2_tie_audit.json').read_bytes()).hexdigest(),
            "primitive_code_sha256":hashlib.sha256((ROOT/"figure_primitives.py").read_bytes()).hexdigest(),
            "figure6_layout_code_sha256":hashlib.sha256((ROOT/"figure6_software_layout.py").read_bytes()).hexdigest(),
            "figure6_example_builder_code_sha256":hashlib.sha256((ROOT/"build_figure6_runtime_example.py").read_bytes()).hexdigest()}
    (ROOT/"FIGURE_BUILD_REPORT.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"figures":len(p.RECORDS),"panels":[r["panels"] for r in p.RECORDS]}))


if __name__=="__main__":main()
