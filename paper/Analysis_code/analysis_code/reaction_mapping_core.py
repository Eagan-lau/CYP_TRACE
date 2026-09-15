"""Core utilities for RXNMapper-driven atom-level reaction-centre reconstruction."""

from __future__ import annotations

from collections import Counter
from typing import Dict, Iterable, List, Sequence, Tuple

from rdkit import Chem


REACTION_SEPARATOR = ">>"


def reaction_to_smiles(reaction: Dict) -> str:
    substrates = reaction["substrates"]
    products = reaction["products"]
    if not substrates or not products:
        raise ValueError("Reaction side cannot be empty")
    return ".".join(substrates) + REACTION_SEPARATOR + ".".join(products)


def _canonical_smiles_components(sides: Sequence[str], strip_maps: bool = False) -> List[str]:
    canonical = []
    for s in sides:
        m = Chem.MolFromSmiles(s)
        if m is None:
            raise ValueError(f"Invalid SMILES component: {s}")
        if strip_maps:
            m = Chem.Mol(m)
            for atom in m.GetAtoms():
                atom.SetAtomMapNum(0)
        canonical.append(Chem.MolToSmiles(m, canonical=True))
    return sorted(canonical)


def _map_information(mol: Chem.Mol) -> Dict[int, Dict[str, int]]:
    info = {}
    for atom in mol.GetAtoms():
        map_id = int(atom.GetAtomMapNum())
        if map_id <= 0:
            continue
        if map_id in info:
            raise ValueError(f"Duplicate atom-map id on one side: {map_id}")
        info[map_id] = {
            "atomic_num": int(atom.GetAtomicNum()),
            "isotope": int(atom.GetIsotope()),
            "charge": int(atom.GetFormalCharge()),
            "hcount": int(atom.GetTotalNumHs()),
            "symbol": atom.GetSymbol(),
        }
    return info


def _bond_signatures(mol: Chem.Mol) -> Dict[Tuple[int, int], Tuple[str, str, int]]:
    sig = {}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if int(a) <= 0 or int(b) <= 0:
            continue
        left, right = (int(a), int(b))
        if left > right:
            left, right = right, left
        sig[(left, right)] = (
            str(bond.GetBondType()),
            str(bond.GetStereo()),
            int(bond.GetBondTypeAsDouble()),
        )
    return sig


def _combine_components(smiles_list: Sequence[str]) -> Chem.Mol:
    if not smiles_list:
        return Chem.Mol()
    merged = Chem.MolFromSmiles(smiles_list[0])
    if merged is None:
        raise ValueError(f"Cannot parse component '{smiles_list[0]}'")
    for s in smiles_list[1:]:
        next_mol = Chem.MolFromSmiles(s)
        if next_mol is None:
            raise ValueError(f"Cannot parse component '{s}'")
        merged = Chem.CombineMols(merged, next_mol)
    return Chem.Mol(merged)


def _atom_center_summary(raw_sub_sig, raw_prod_sig, map_sub_sig, map_prod_sig, map_sub_bonds, map_prod_bonds):
    shared_maps = sorted(set(map_sub_sig) & set(map_prod_sig))
    sub_only = sorted(set(map_sub_sig) - set(map_prod_sig))
    prod_only = sorted(set(map_prod_sig) - set(map_sub_sig))

    added = sorted(set(map_prod_bonds) - set(map_sub_bonds))
    removed = sorted(set(map_sub_bonds) - set(map_prod_bonds))
    changed = sorted(
        key for key in (set(map_sub_bonds) & set(map_prod_bonds))
        if map_sub_bonds[key] != map_prod_bonds[key]
    )

    changing_atoms = set(sub_only + prod_only)
    for i, j in added + removed + changed:
        changing_atoms.add(i)
        changing_atoms.add(j)
    for map_id in shared_maps:
        before = map_sub_sig[map_id]
        after = map_prod_sig[map_id]
        if before != after:
            changing_atoms.add(map_id)

    property_changes = []
    for map_id in shared_maps:
        before = map_sub_sig[map_id]
        after = map_prod_sig[map_id]
        if before != after:
            property_changes.append({
                "atom_map": map_id,
                "mapped_reactant": before,
                "mapped_product": after,
            })

    return {
        "reactant_only_atom_maps": sub_only,
        "product_only_atom_maps": prod_only,
        "shared_atom_maps": shared_maps,
        "changed_center_atom_maps": sorted(changing_atoms),
        "changed_atom_attributes": property_changes,
        "added_bonds": added,
        "removed_bonds": removed,
        "changed_bonds": changed,
    }


def analyze_mapping(candidate_key: str, raw_reaction: str, mapped_reaction: str, confidence: float, rxnmapper_version: str) -> Dict:
    if not mapped_reaction:
        raise ValueError("Mapper returned empty mapped reaction")
    if REACTION_SEPARATOR not in mapped_reaction:
        raise ValueError("Mapped reaction split is invalid")

    raw_left, raw_right = raw_reaction.split(REACTION_SEPARATOR, 1)
    mapped_left, mapped_right = mapped_reaction.split(REACTION_SEPARATOR, 1)
    raw_subs = raw_left.split(".") if raw_left else []
    raw_prods = raw_right.split(".") if raw_right else []
    map_subs = mapped_left.split(".") if mapped_left else []
    map_prods = mapped_right.split(".") if mapped_right else []

    raw_sub_counter = Counter(_canonical_smiles_components(raw_subs, strip_maps=True))
    raw_prod_counter = Counter(_canonical_smiles_components(raw_prods, strip_maps=True))
    map_sub_counter = Counter(_canonical_smiles_components(map_subs, strip_maps=True))
    map_prod_counter = Counter(_canonical_smiles_components(map_prods, strip_maps=True))

    map_side_preserved = (raw_sub_counter == map_sub_counter) and (raw_prod_counter == map_prod_counter)
    raw_sub_mols = [_to_mol(s) for s in raw_subs]
    raw_prod_mols = [_to_mol(s) for s in raw_prods]
    map_sub_mols = [_to_mol(s) for s in map_subs]
    map_prod_mols = [_to_mol(s) for s in map_prods]
    map_sub_sig = _collect_map_signatures(map_sub_mols)
    map_prod_sig = _collect_map_signatures(map_prod_mols)
    map_sub_mols_combined = _combine_components(map_subs)
    map_prod_mols_combined = _combine_components(map_prods)
    map_sub_bonds = _bond_signatures(map_sub_mols_combined)
    map_prod_bonds = _bond_signatures(map_prod_mols_combined)

    raw_sub_sig = _collect_raw_signatures(raw_sub_mols)
    raw_prod_sig = _collect_raw_signatures(raw_prod_mols)
    mapped_counts, mapped_inconsistencies = _mapped_consistency_checks(map_sub_sig, map_prod_sig)
    centers = _atom_center_summary(raw_sub_sig, raw_prod_sig, map_sub_sig, map_prod_sig, map_sub_bonds, map_prod_bonds)

    shared_symbol_isotope_mismatch = sorted(mapped_inconsistencies)
    status = "mapped" if map_side_preserved and not shared_symbol_isotope_mismatch else "mapped_validation_warning"
    return {
        "reaction_key": candidate_key,
        "status": status,
        "reaction_smiles": raw_reaction,
        "mapped_reaction": mapped_reaction,
        "mapped_confidence": float(confidence),
        "rxnmapper_version": rxnmapper_version,
        "validation": {
            "sides_preserved_exact_multiset": map_side_preserved,
            "raw_unique_substrate_atom_maps": len(raw_sub_sig),
            "raw_unique_product_atom_maps": len(raw_prod_sig),
            "mapped_unique_substrate_atom_maps": len(map_sub_sig),
            "mapped_unique_product_atom_maps": len(map_prod_sig),
            "shared_map_symbol_or_isotope_mismatch": shared_symbol_isotope_mismatch,
        },
        "centers": centers,
    }


def _to_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def _collect_map_signatures(mols: Sequence[Chem.Mol]) -> Dict[int, Dict[str, int]]:
    signatures = {}
    for mol in mols:
        signatures.update(_map_information(mol))
    return signatures


def _collect_raw_signatures(mols: Sequence[Chem.Mol]) -> Dict[int, Dict[str, int]]:
    return {}


def _mapped_consistency_checks(map_sub_sig: Dict[int, Dict[str, int]], map_prod_sig: Dict[int, Dict[str, int]]) -> Tuple[int, List[int]]:
    shared = set(map_sub_sig) & set(map_prod_sig)
    mapped_count = len(shared)
    mismatches = []
    for map_id in sorted(shared):
        s = map_sub_sig[map_id]
        p = map_prod_sig[map_id]
        if s["atomic_num"] != p["atomic_num"] or s["isotope"] != p["isotope"]:
            mismatches.append(map_id)
    return mapped_count, mismatches


def map_candidate(reaction_id: str, reaction: Dict, mapper, rxnmapper_version: str) -> Dict:
    raw_reaction = reaction_to_smiles(reaction)
    try:
        mapped = mapper.get_attention_guided_atom_maps([raw_reaction])[0]
    except Exception as exc:
        return {
            "reaction_key": reaction_id,
            "status": "mapper_error",
            "reaction_smiles": raw_reaction,
            "mapped_reaction": "",
            "mapped_confidence": 0.0,
            "rxnmapper_version": rxnmapper_version,
            "error": f"{type(exc).__name__}: {exc}",
            "validation": {
                "sides_preserved_exact_multiset": False,
                "raw_unique_substrate_atom_maps": 0,
                "raw_unique_product_atom_maps": 0,
                "mapped_unique_substrate_atom_maps": 0,
                "mapped_unique_product_atom_maps": 0,
                "shared_map_symbol_or_isotope_mismatch": [],
            },
            "centers": {
                "reactant_only_atom_maps": [],
                "product_only_atom_maps": [],
                "shared_atom_maps": [],
                "changed_center_atom_maps": [],
                "changed_atom_attributes": [],
                "added_bonds": [],
                "removed_bonds": [],
                "changed_bonds": [],
            },
        }

    mapped_reaction = mapped.get("mapped_rxn", "")
    confidence = float(mapped.get("confidence", 0.0))
    try:
        return analyze_mapping(reaction_id, raw_reaction, mapped_reaction, confidence, rxnmapper_version)
    except Exception as exc:
        return {
            "reaction_key": reaction_id,
            "status": "mapping_validation_error",
            "reaction_smiles": raw_reaction,
            "mapped_reaction": mapped_reaction,
            "mapped_confidence": confidence,
            "rxnmapper_version": rxnmapper_version,
            "error": f"{type(exc).__name__}: {exc}",
            "validation": {
                "sides_preserved_exact_multiset": False,
                "raw_unique_substrate_atom_maps": 0,
                "raw_unique_product_atom_maps": 0,
                "mapped_unique_substrate_atom_maps": 0,
                "mapped_unique_product_atom_maps": 0,
                "shared_map_symbol_or_isotope_mismatch": [],
            },
            "centers": {
                "reactant_only_atom_maps": [],
                "product_only_atom_maps": [],
                "shared_atom_maps": [],
                "changed_center_atom_maps": [],
                "changed_atom_attributes": [],
                "added_bonds": [],
                "removed_bonds": [],
                "changed_bonds": [],
            },
        }
