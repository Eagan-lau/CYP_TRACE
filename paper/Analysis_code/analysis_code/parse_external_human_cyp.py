"""Normalize and audit Figshare 26630515 v4 without calculating model scores."""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold

from run_raw import digest_file, normalize_structure, now, write_json


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_additions" / "figshare_26630515_v4"
OUTPUT = ROOT / "external_human_cyp_01"
SOURCE = ROOT / "human_substrate_01" / "normalized_labels.json"
ACQUISITION = OUTPUT / "acquisition_manifest.json"
PROTOCOL = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1.md"
AMENDMENT = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1_AMENDMENT_1.md"
AMENDMENT_2 = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1_AMENDMENT_2.md"
AMENDMENT_3 = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1_AMENDMENT_3.md"
FILE_PATTERN = re.compile(r"^(CYP1A2|CYP2C9|CYP2C19|CYP2D6|CYP2E1|CYP3A4)_(training|testing)set\.csv$")
EXPOSED_SOURCES = {"Holmer": "CYPstrate", "Siyang_Tian": "Tian"}
KNOWN_EXTERNAL_SOURCES = {
    "Yamashita", "DrugBank", "CYP Knowledgebase", "Yap", "SuperCYP",
    "Tao_Zhang", "Michielan", "Mishra", "Mayo Clinic Laboratories",
}


def scaffold_for(molecule: Chem.Mol) -> str:
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    if scaffold.GetNumAtoms():
        return "ring:" + Chem.MolToSmiles(scaffold, isomericSmiles=False)
    return "acyclic_connectivity:" + Chem.MolToSmiles(molecule, isomericSmiles=False)


def standard_identity(canonical_smiles: str):
    """Return one notation and scaffold per standard InChI identity."""
    molecule = Chem.MolFromSmiles(canonical_smiles)
    inchi = Chem.MolToInchi(molecule)
    standardized = Chem.MolFromInchi(inchi)
    if standardized is None:
        return Chem.InchiToInchiKey(inchi), canonical_smiles, scaffold_for(molecule), "inchi_reconstruction_failed"
    standardized_smiles = Chem.MolToSmiles(standardized, canonical=True, isomericSmiles=True)
    return Chem.InchiToInchiKey(inchi), standardized_smiles, scaffold_for(standardized), "inchi_roundtrip"


def stratum_counts(rows, name):
    chosen = [row for row in rows if row["strata"][name]]
    counts = Counter((row["isoform"], row["label"]) for row in chosen)
    return {
        "rows": len(chosen),
        "unique_compounds": len({row["compound_inchikey"] for row in chosen}),
        "scaffold_groups": len({row["scaffold_group"] for row in chosen}),
        "by_isoform": [
            {
                "isoform": isoform,
                "positives": counts[(isoform, 1)],
                "negatives": counts[(isoform, 0)],
            }
            for isoform in sorted({row["isoform"] for row in rows})
        ],
    }


def main() -> None:
    RDLogger.DisableLog("rdApp.*")
    acquisition = json.loads(ACQUISITION.read_text(encoding="utf-8"))
    if not acquisition.get("complete"):
        raise RuntimeError("Figshare acquisition is incomplete")
    source_rows = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_compounds = {row["compound_inchikey"] for row in source_rows}
    source_representations = defaultdict(set)
    for row in source_rows:
        key, representation, _, _ = standard_identity(row["canonical_smiles"])
        source_representations[key].add(representation)
    source_scaffolds = {
        scaffold_for(Chem.MolFromSmiles(min(representations)))
        for representations in source_representations.values()
    }
    source_pairs = {(row["isoform"], row["compound_inchikey"]): row["label"] for row in source_rows}

    raw_rows = []
    groups = defaultdict(list)
    counts = Counter()
    observed_sources = set()
    files = []
    for path in sorted(RAW.glob("*set.csv")):
        match = FILE_PATTERN.fullmatch(path.name)
        if not match:
            raise ValueError(f"unexpected model-data filename: {path.name}")
        isoform, split = match.groups()
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != ["Name", "SMILES", "Label", "Source"]:
                raise ValueError(f"unexpected schema in {path.name}: {reader.fieldnames}")
            local_rows = 0
            for line, source_row in enumerate(reader, start=2):
                local_rows += 1
                counts["raw_rows"] += 1
                label_text = source_row["Label"].strip()
                label = int(label_text) if label_text in {"0", "1"} else None
                source = source_row["Source"].strip()
                observed_sources.add(source)
                structure = normalize_structure(source_row["SMILES"])
                canonical = structure["canonical_smiles"]
                inchikey = ""
                scaffold = ""
                components = 0
                if canonical:
                    molecule = Chem.MolFromSmiles(canonical)
                    inchikey, standardized_query_smiles, scaffold, standardization_status = standard_identity(canonical)
                    components = len(Chem.GetMolFrags(molecule))
                else:
                    standardized_query_smiles = ""
                    standardization_status = "unresolved"
                row = {
                    "isoform": isoform,
                    "source_split": split,
                    "source_file": path.name,
                    "source_row": line,
                    "name": source_row["Name"],
                    "source": source,
                    "raw_label": label_text,
                    "label": label,
                    **structure,
                    "compound_inchikey": inchikey,
                    "standardized_query_smiles": standardized_query_smiles,
                    "standardization_status": standardization_status,
                    "scaffold_group": scaffold,
                    "component_count": components,
                }
                raw_rows.append(row)
                if label is None:
                    counts["nonbinary_label_rows"] += 1
                if not inchikey:
                    counts["unresolved_structure_rows"] += 1
                if components > 1:
                    counts["multicomponent_rows"] += 1
                if label is not None and inchikey:
                    groups[(isoform, inchikey)].append(row)
            files.append({"file": path.name, "isoform": isoform, "source_split": split, "rows": local_rows,
                          "sha256": digest_file(path)})

    external_representations = defaultdict(set)
    for row in raw_rows:
        if row["compound_inchikey"]:
            external_representations[row["compound_inchikey"]].add(row["standardized_query_smiles"])
    chosen_representation = {key: min(values) for key, values in external_representations.items()}
    for row in raw_rows:
        key = row["compound_inchikey"]
        if not key:
            continue
        chosen = chosen_representation[key]
        row["standardized_query_smiles"] = chosen
        row["scaffold_group"] = scaffold_for(Chem.MolFromSmiles(chosen))
        if row["standardization_status"] == "inchi_reconstruction_failed":
            counts["inchi_reconstruction_failed_rows"] += 1
        if len(external_representations[key]) > 1:
            counts["key_level_representation_resolution_rows"] += 1

    unknown_sources = observed_sources - set(EXPOSED_SOURCES) - KNOWN_EXTERNAL_SOURCES
    normalized = []
    for (isoform, inchikey), members in sorted(groups.items()):
        labels = {row["label"] for row in members}
        sources = sorted({row["source"] for row in members})
        label_conflict = len(labels) != 1
        source_exposed = any(source in EXPOSED_SOURCES for source in sources)
        source_recognized = all(source in EXPOSED_SOURCES or source in KNOWN_EXTERNAL_SOURCES for source in sources)
        representative = min(members, key=lambda row: (row["source_file"], row["source_row"]))
        standardized_variants = {row["standardized_query_smiles"] for row in members}
        if len(standardized_variants) != 1:
            raise ValueError(f"standard InChI did not resolve representation for {inchikey}")
        exact_compound_overlap = inchikey in source_compounds
        exact_pair_overlap = (isoform, inchikey) in source_pairs
        label = None if label_conflict else next(iter(labels))
        source_label = source_pairs.get((isoform, inchikey))
        label_relation = (
            "not_comparable" if label is None or source_label is None else
            "agree" if label == source_label else "conflict"
        )
        exact_novel = not exact_compound_overlap
        scaffold_novel = exact_novel and representative["scaffold_group"] not in source_scaffolds
        lineage_eligible = exact_novel and source_recognized and not source_exposed
        normalized.append({
            "isoform": isoform,
            "compound_inchikey": inchikey,
            "canonical_smiles": representative["standardized_query_smiles"],
            "source_canonical_smiles_variants": sorted({row["canonical_smiles"] for row in members}),
            "scaffold_group": representative["scaffold_group"],
            "component_count": representative["component_count"],
            "label": label,
            "label_conflict": label_conflict,
            "sources": sources,
            "source_lineage_status": (
                "development_exposed" if source_exposed else
                "recognized_non_cypstrate_lineage" if source_recognized else
                "unverifiable"
            ),
            "exact_compound_overlap": exact_compound_overlap,
            "exact_pair_overlap": exact_pair_overlap,
            "cypstrate_label": source_label,
            "label_relation": label_relation,
            "source_locators": [
                {"file": row["source_file"], "split": row["source_split"], "row": row["source_row"],
                 "source": row["source"], "name": row["name"]}
                for row in members
            ],
            "strata": {
                "all_parseable_unique": not label_conflict,
                "exact_novel": exact_novel and not label_conflict,
                "scaffold_novel": scaffold_novel and not label_conflict,
                "lineage_eligible": lineage_eligible and not label_conflict,
                "scaffold_and_lineage_eligible": scaffold_novel and lineage_eligible and not label_conflict,
            },
        })
        counts["duplicate_row_excess"] += len(members) - 1
        if label_conflict:
            counts["conflicting_isoform_compound_groups"] += 1
        if exact_pair_overlap:
            counts["exact_pair_overlap_groups"] += 1
            counts["exact_pair_" + label_relation] += 1

    for row in raw_rows:
        if row["source"] in EXPOSED_SOURCES:
            counts["development_exposed_source_rows"] += 1

    with (OUTPUT / "raw_rows.jsonl").open("w", encoding="utf-8") as stream:
        for row in raw_rows:
            stream.write(json.dumps(row, ensure_ascii=False) + "\n")
    write_json(OUTPUT / "normalized_external_labels.json", normalized)
    strata = {name: stratum_counts(normalized, name) for name in (
        "all_parseable_unique", "exact_novel", "scaffold_novel",
        "lineage_eligible", "scaffold_and_lineage_eligible",
    )}
    checks = {
        "complete_acquisition": acquisition["complete"],
        "twelve_label_files": len(files) == 12,
        "all_labels_binary": counts["nonbinary_label_rows"] == 0,
        "all_sources_recognized": not unknown_sources,
        "unique_normalized_keys": len(normalized) == len({(row["isoform"], row["compound_inchikey"]) for row in normalized}),
        "no_models_or_scores_calculated": True,
    }
    audit = {
        "created_utc": now(),
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "counts": dict(counts),
        "raw_files": files,
        "observed_sources": sorted(observed_sources),
        "unknown_sources": sorted(unknown_sources),
        "normalized_groups": len(normalized),
        "strata": strata,
        "source_label_scope": "source-reported binary human CYP substrate/non-substrate; not assay harmonized",
        "source_split_used_for_independence": False,
        "models_evaluated": False,
        "rdkit_version": rdBase.rdkitVersion,
        "input_hashes": {
            "acquisition_manifest": digest_file(ACQUISITION),
            "cypstrate_normalized_labels": digest_file(SOURCE),
            "protocol": digest_file(PROTOCOL),
            "protocol_clarification_1": digest_file(AMENDMENT),
            "protocol_clarification_2": digest_file(AMENDMENT_2),
            "protocol_clarification_3": digest_file(AMENDMENT_3),
            "script": digest_file(Path(__file__)),
        },
    }
    write_json(OUTPUT / "overlap_audit.json", audit)
    print(json.dumps({
        "status": audit["status"],
        "counts": audit["counts"],
        "normalized_groups": audit["normalized_groups"],
        "observed_sources": audit["observed_sources"],
        "unknown_sources": audit["unknown_sources"],
        "strata": audit["strata"],
    }, indent=2, ensure_ascii=False))
    if audit["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
