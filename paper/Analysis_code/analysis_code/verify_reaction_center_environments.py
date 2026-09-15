"""Independent cell-level verification of reaction-centre environments."""
from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem


BINS = 512
RADII = (0, 1, 2)


def file_hash(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def token_bin(token: str) -> int:
    return int.from_bytes(hashlib.sha256(token.encode("utf-8")).digest()[:8], "big") % BINS


def component_molecules(side: str):
    output = []
    for smiles in filter(None, side.split(".")):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise AssertionError(f"unparseable component: {smiles}")
        output.append(molecule)
    return output


def map_locations(side: str):
    locations = {}
    for molecule in component_molecules(side):
        for atom in molecule.GetAtoms():
            atom_map = int(atom.GetAtomMapNum())
            if atom_map:
                if atom_map in locations:
                    raise AssertionError(f"duplicate atom map: {atom_map}")
                locations[atom_map] = (molecule, int(atom.GetIdx()))
    return locations


def neighborhood(molecule, root, radius):
    distance = {root: 0}
    pending = deque([root])
    while pending:
        atom_index = pending.popleft()
        if distance[atom_index] == radius:
            continue
        for atom in molecule.GetAtomWithIdx(atom_index).GetNeighbors():
            neighbor = int(atom.GetIdx())
            if neighbor not in distance:
                distance[neighbor] = distance[atom_index] + 1
                pending.append(neighbor)
    return sorted(distance)


def environment(molecule, root, radius):
    molecule = Chem.Mol(molecule)
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    atoms = neighborhood(molecule, root, radius)
    selected = set(atoms)
    bonds = sorted(int(bond.GetIdx()) for bond in molecule.GetBonds()
                   if int(bond.GetBeginAtomIdx()) in selected and int(bond.GetEndAtomIdx()) in selected)
    return Chem.MolFragmentToSmiles(molecule, atomsToUse=atoms, bondsToUse=bonds,
        rootedAtAtom=root, canonical=True, isomericSmiles=True,
        allBondsExplicit=True, allHsExplicit=True)


def make_tokens(side, changed, label):
    locations = map_locations(side)
    counts = Counter()
    for atom_map in sorted(set(changed)):
        if atom_map not in locations:
            continue
        molecule, root = locations[atom_map]
        for radius in RADII:
            counts[f"{label}|radius={radius}|{environment(molecule, root, radius)}"] += 1
    return counts


def reconstruct(mapping):
    vector = np.zeros(2 * BINS + 3, dtype=np.float32)
    status = mapping.get("status")
    center = mapping.get("centers") or {}
    changed = sorted(set(int(x) for x in center.get("changed_center_atom_maps", [])))
    available = status in {"mapped", "mapped_validation_warning"} and bool(changed) and ">>" in mapping.get("mapped_reaction", "")
    substrate_tokens, product_tokens = Counter(), Counter()
    if available:
        substrate, product = mapping["mapped_reaction"].split(">>", 1)
        substrate_tokens = make_tokens(substrate, changed, "substrate")
        product_tokens = make_tokens(product, changed, "product")
        for token, count in substrate_tokens.items():
            vector[token_bin(token)] += count
        for token, count in product_tokens.items():
            vector[BINS + token_bin(token)] += count
        vector[-3:] = (len(center.get("added_bonds", [])), len(center.get("removed_bonds", [])), len(center.get("changed_bonds", [])))
    return vector, available, changed, substrate_tokens, product_tokens


def run(root: Path, feature_name: str, output_name: str) -> None:
    feature_dir = root / feature_name
    output = root / output_name
    if output.exists():
        raise FileExistsError(output)
    output.mkdir()
    source_rows = [json.loads(line) for line in (root / "reaction_mapping_01/mapping_results.jsonl").read_text(encoding="utf-8").splitlines() if line.strip()]
    source = {row["reaction_key"]: row for row in source_rows}
    records = [json.loads(line) for line in (feature_dir / "environment_records.jsonl").read_text(encoding="utf-8").splitlines() if line.strip()]
    archive = np.load(feature_dir / "environment_features.npz", allow_pickle=False)
    reaction_keys = archive["reaction_keys"].tolist()
    matrix = archive["feature_vector"]
    expected_names = ([f"substrate_env_hash_{i:03d}" for i in range(BINS)] +
                      [f"product_env_hash_{i:03d}" for i in range(BINS)] +
                      ["directed_added_bond_count", "directed_removed_bond_count", "bond_order_change_count"])
    errors = []
    max_difference = 0.0
    differing_cells = 0
    if archive["feature_names"].tolist() != expected_names:
        errors.append("feature_names_mismatch")
    if reaction_keys != sorted(source) or len(records) != len(reaction_keys):
        errors.append("reaction_key_or_row_count_mismatch")
    record_by_key = {row["reaction_key"]: row for row in records}
    for index, key in enumerate(reaction_keys):
        expected, available, changed, substrate_tokens, product_tokens = reconstruct(source[key])
        difference = np.abs(expected.astype(np.float64) - matrix[index].astype(np.float64))
        max_difference = max(max_difference, float(difference.max(initial=0.0)))
        differing_cells += int(np.count_nonzero(difference))
        row = record_by_key.get(key)
        if row is None:
            errors.append(f"missing_record:{key}")
            continue
        checks = [
            row["center_available"] == available,
            row["changed_center_atom_maps"] == changed,
            row["substrate_environment_tokens"] == dict(sorted(substrate_tokens.items())),
            row["product_environment_tokens"] == dict(sorted(product_tokens.items())),
            row["feature_nonzero_count"] == int(np.count_nonzero(expected)),
            abs(row["feature_sum"] - float(expected.sum())) <= 1e-12,
        ]
        if not all(checks):
            errors.append(f"record_mismatch:{key}")
            if len(errors) > 20:
                break
    audit = json.loads((feature_dir / "audit.json").read_text(encoding="utf-8"))
    for name, digest in audit.get("output_manifest", {}).items():
        if file_hash(feature_dir / name) != digest:
            errors.append(f"output_hash_mismatch:{name}")
    report = {
        "status": "PASS" if not errors and differing_cells == 0 else "FAIL",
        "reaction_count": len(reaction_keys),
        "matrix_shape": list(matrix.shape),
        "cells_recomputed": int(matrix.size),
        "differing_cells": differing_cells,
        "maximum_absolute_difference": max_difference,
        "center_available_count_recomputed": int(sum(reconstruct(source[key])[1] for key in reaction_keys)),
        "errors": errors,
        "independent_implementation": True,
    }
    (output / "audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--feature-name", default="reaction_center_environments_01")
    parser.add_argument("--output-name", default="reaction_center_environment_validation_01")
    arguments = parser.parse_args()
    run(arguments.root.resolve(), arguments.feature_name, arguments.output_name)
