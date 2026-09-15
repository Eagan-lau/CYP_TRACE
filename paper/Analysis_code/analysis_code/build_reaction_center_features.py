"""Build reaction-centre-derived features from RXNMapper output.

This is the explicit bridge for the next stage: translating mapped
atom-level changes into reaction-level, fixed-size descriptors that can be
joined with protein features later.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from rdkit import Chem

from run_raw import digest_file, now, write_json


def _split_reaction(smiles: str) -> Tuple[List[str], List[str]]:
    left, right = smiles.split(">>", 1)
    return ([x for x in left.split(".") if x], [x for x in right.split(".") if x])


def _mols_with_maps(smiles_list: Iterable[str]) -> List[Chem.Mol]:
    out = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid mapped component: {smiles}")
        out.append(mol)
    return out


def _collect_atom_signatures(mols: Iterable[Chem.Mol]) -> Dict[int, Dict[str, int]]:
    out = {}
    for mol in mols:
        for atom in mol.GetAtoms():
            map_id = int(atom.GetAtomMapNum())
            if map_id <= 0:
                continue
            if map_id in out:
                raise ValueError(f"Duplicate atom-map id in one side: {map_id}")
            out[map_id] = {
                "atomic_num": int(atom.GetAtomicNum()),
                "symbol": atom.GetSymbol(),
                "isotope": int(atom.GetIsotope()),
                "charge": int(atom.GetFormalCharge()),
                "hcount": int(atom.GetTotalNumHs()),
                "valence": int(atom.GetTotalValence()),
                "neighbors": int(atom.GetDegree()),
                "aromatic": bool(atom.GetIsAromatic()),
                "implicit_h": int(atom.GetNumImplicitHs()),
                "explicit_h": int(atom.GetNumExplicitHs()),
            }
    return out


def _collect_bonds(mols: Iterable[Chem.Mol]) -> Dict[Tuple[int, int], Tuple[str, int, bool]]:
    out = {}
    for mol in mols:
        for bond in mol.GetBonds():
            a = int(bond.GetBeginAtom().GetAtomMapNum())
            b = int(bond.GetEndAtom().GetAtomMapNum())
            if a <= 0 or b <= 0:
                continue
            i, j = (a, b) if a < b else (b, a)
            out[(i, j)] = (
                str(bond.GetBondType()),
                int(bond.GetBondTypeAsDouble()),
                bool(bond.GetIsAromatic()),
            )
    return out


def _symbol_counts(atoms: Iterable[int], signatures: Dict[int, Dict[str, int]]) -> Dict[str, int]:
    counts = Counter()
    for map_id in atoms:
        counts[signatures[map_id]["symbol"]] += 1
    return dict(sorted(counts.items()))


def _feature_from_row(row: Dict, reaction: Dict, rid: str) -> Dict:
    mapped_reaction = row.get("mapped_reaction", "")
    status = row.get("status", "unavailable")
    out = {
        "reaction_key": rid,
        "status": status,
        "mapped_confidence": float(row.get("mapped_confidence", 0.0)),
        "rxnmapper_version": row.get("rxnmapper_version", ""),
    }
    if status not in {"mapped", "mapped_validation_warning"} or ">>" not in mapped_reaction:
        out.update({
            "center_atom_maps": [],
            "center_atom_map_count": 0,
            "added_bonds": [],
            "removed_bonds": [],
            "changed_bonds": [],
            "atom_count_change": 0,
            "molecule_map_ok": False,
            "map_signature": None,
            "feature_vector": [0] * 24,
            "center_coverage_ratio": 0.0,
            "mapping_warning": bool(status == "mapped_validation_warning"),
            "error": row.get("error"),
        })
        return out

    raw_left, raw_right = _split_reaction(reaction["reaction_smiles"])
    maps_left, maps_right = _split_reaction(mapped_reaction)
    mapped_sub = _mols_with_maps(maps_left)
    mapped_prod = _mols_with_maps(maps_right)
    sub_sig = _collect_atom_signatures(mapped_sub)
    prod_sig = _collect_atom_signatures(mapped_prod)
    sub_bonds = _collect_bonds(mapped_sub)
    prod_bonds = _collect_bonds(mapped_prod)

    sub_maps = set(sub_sig)
    prod_maps = set(prod_sig)
    shared = sorted(sub_maps & prod_maps)
    reactant_only = sorted(sub_maps - prod_maps)
    product_only = sorted(prod_maps - sub_maps)
    removed_bonds = sorted(set(sub_bonds) - set(prod_bonds))
    added_bonds = sorted(set(prod_bonds) - set(sub_bonds))
    changed_bonds = sorted(
        (a, b) for (a, b), info in sub_bonds.items()
        if (a, b) in prod_bonds and prod_bonds[(a, b)] != info
    )
    changed_atom_maps = set()
    for map_id in shared:
        if sub_sig[map_id] != prod_sig[map_id]:
            changed_atom_maps.add(map_id)
    for a, b in removed_bonds + added_bonds + changed_bonds:
        changed_atom_maps.add(a)
        changed_atom_maps.add(b)
    changed_atom_maps = sorted(changed_atom_maps)

    total_atoms = len(sub_maps | prod_maps)
    center_coverage_ratio = (len(changed_atom_maps) / total_atoms) if total_atoms else 0.0

    changed_symbols = _symbol_counts(changed_atom_maps, {**sub_sig, **prod_sig}) if changed_atom_maps else {}
    added_symbols = _symbol_counts(product_only, prod_sig) if product_only else {}
    removed_symbols = _symbol_counts(reactant_only, sub_sig) if reactant_only else {}

    center_signature = {
        "changed_atom_count": len(changed_atom_maps),
        "changed_bond_count": len(added_bonds) + len(removed_bonds) + len(changed_bonds),
        "reactant_only_atom_count": len(reactant_only),
        "product_only_atom_count": len(product_only),
        "changed_symbols": changed_symbols,
        "added_symbols": added_symbols,
        "removed_symbols": removed_symbols,
        "raw_reactant_side_count": len(raw_left),
        "raw_product_side_count": len(raw_right),
        "mapped_atom_count": total_atoms,
    }
    # A compact dense vector for downstream linear models.
    # [center_atom_count, changed_atom_count, added_bond_count, removed_bond_count,
    #  changed_bond_count, reactant_only_count, product_only_count, coverage_ratio]
    # + total counts per main elements (C,N,O,S,P,Cl,Br,F,I,Fe,others)
    elements = ["C", "N", "O", "S", "P", "Cl", "Br", "F", "I", "B", "Se", "Fe", "Co", "Zn", "Mg", "Cu", "Ni", "others"]
    element_counts = []
    for e in elements:
        element_counts.append(changed_symbols.get(e, 0) + added_symbols.get(e, 0) - removed_symbols.get(e, 0))

    feature_vector = [
        float(len(changed_atom_maps)),
        float(len(shared)),
        float(len(added_bonds)),
        float(len(removed_bonds)),
        float(len(changed_bonds)),
        float(len(reactant_only)),
        float(len(product_only)),
        center_coverage_ratio,
    ] + [float(v) for v in element_counts]

    return {
        **out,
        "center_atom_maps": changed_atom_maps,
        "center_atom_map_count": len(changed_atom_maps),
        "center_atom_attributes": [
            {"map": int(k), **sub_sig[k]} for k in changed_atom_maps if k in sub_sig
        ],
        "added_bonds": [[int(a), int(b)] for a, b in added_bonds],
        "removed_bonds": [[int(a), int(b)] for a, b in removed_bonds],
        "changed_bonds": [[int(a), int(b)] for a, b in changed_bonds],
        "atom_count_change": len(product_only) - len(reactant_only),
        "molecule_map_ok": True,
        "map_signature": center_signature,
        "feature_vector": feature_vector,
        "center_coverage_ratio": center_coverage_ratio,
        "mapping_warning": bool(status == "mapped_validation_warning"),
    }


def run(root: Path) -> None:
    chemistry = json.loads((root / "dataset_02/core_reactions.json").read_text())
    mapping_lines = (root / "reaction_mapping_01/mapping_results.jsonl").read_text().splitlines()
    mapping = {json.loads(line)["reaction_key"]: json.loads(line) for line in mapping_lines if line.strip()}
    if len(mapping) != len(chemistry):
        raise ValueError("Mapping and chemistry counts mismatch")
    if not (root / "reaction_mapping_01/summary.json").exists():
        raise FileNotFoundError("Missing reaction_mapping_01/summary.json")
    out = root / "reaction_center_features_01"
    if out.exists():
        raise FileExistsError(out)
    out.mkdir()
    rows = []
    by_status = Counter()
    all_features = []
    for rid in sorted(chemistry):
        if rid not in mapping:
            raise KeyError(rid)
        row = mapping[rid]
        by_status[row["status"]] += 1
        feature_row = _feature_from_row(row, {"reaction_smiles": row["reaction_smiles"]}, rid)
        rows.append(feature_row)
        all_features.append(feature_row["feature_vector"])
    with (out / "reaction_center_features.jsonl").open("w", encoding="utf-8") as handle:
        for r in rows:
            handle.write(json.dumps(r, ensure_ascii=False) + "\n")
    write_json(out / "summary.json", {
        "created_utc": now(),
        "reaction_count": len(rows),
        "mapped_count": int(sum(1 for r in rows if r["status"] == "mapped")),
        "mapped_validation_warning_count": int(sum(1 for r in rows if r["status"] == "mapped_validation_warning")),
        "mapper_error_count": int(sum(1 for r in rows if r["status"] == "mapper_error")),
        "validation_error_count": int(sum(1 for r in rows if r["status"] == "mapping_validation_error")),
        "status_counts": dict(sorted(by_status.items())),
        "feature_length": len(all_features[0]) if all_features else 0,
        "feature_stats": {
            "center_atom_count_min": float(min(r["center_atom_map_count"] for r in rows)),
            "center_atom_count_max": float(max(r["center_atom_map_count"] for r in rows)),
            "center_atom_count_mean": float(sum(r["center_atom_map_count"] for r in rows) / len(rows)),
            "added_bonds_mean": float(sum(len(r["added_bonds"]) for r in rows) / len(rows)),
            "removed_bonds_mean": float(sum(len(r["removed_bonds"]) for r in rows) / len(rows)),
            "changed_bonds_mean": float(sum(len(r["changed_bonds"]) for r in rows) / len(rows)),
            "mapping_warning_count": int(by_status["mapped_validation_warning"]),
        },
        "input_manifest": {
            "dataset_02/core_reactions.json": digest_file(root / "dataset_02/core_reactions.json"),
            "reaction_mapping_01/mapping_results.jsonl": digest_file(root / "reaction_mapping_01/mapping_results.jsonl"),
            "reaction_mapping_01/summary.json": digest_file(root / "reaction_mapping_01/summary.json"),
        },
    })
    # Persist the dense matrix for easy concatenation with sequence features.
    # Keep in the same format as other modules.
    import numpy as np

    np.savez_compressed(
        out / "center_features.npz",
        reaction_keys=[rid for rid in sorted(chemistry)],
        feature_vector=np.asarray(all_features, dtype=np.float32),
    )
    print(json.dumps({
        "created": str(out),
        "rows": len(rows),
        "vector_dim": len(all_features[0]) if all_features else 0,
    }, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    run(p.parse_args().root)
