"""Build stable candidate-specific reaction-centre environment features.

The representation is label free.  Mapped atom identifiers are used only to
locate changed atoms and are cleared before any rooted environment is encoded.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
import hashlib
import json
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import platform
import sys
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from rdkit import Chem


VALID_STATUSES = {"mapped", "mapped_validation_warning"}
RADII = (0, 1, 2)
HASH_BINS_PER_SIDE = 512
HASH_NAME = "sha256_first_8_bytes_big_endian_mod_512"
FEATURE_NAMES = (
    [f"substrate_env_hash_{i:03d}" for i in range(HASH_BINS_PER_SIDE)]
    + [f"product_env_hash_{i:03d}" for i in range(HASH_BINS_PER_SIDE)]
    + ["directed_added_bond_count", "directed_removed_bond_count", "bond_order_change_count"]
)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def digest_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def parse_side(mapped_side: str) -> List[Chem.Mol]:
    mols: List[Chem.Mol] = []
    for component in (x for x in mapped_side.split(".") if x):
        mol = Chem.MolFromSmiles(component)
        if mol is None:
            raise ValueError(f"invalid mapped component: {component}")
        mols.append(mol)
    return mols


def locate_maps(mols: Iterable[Chem.Mol]) -> Dict[int, Tuple[Chem.Mol, int]]:
    found: Dict[int, Tuple[Chem.Mol, int]] = {}
    for mol in mols:
        for atom in mol.GetAtoms():
            atom_map = int(atom.GetAtomMapNum())
            if atom_map <= 0:
                continue
            if atom_map in found:
                raise ValueError(f"duplicate atom map on one side: {atom_map}")
            found[atom_map] = (mol, int(atom.GetIdx()))
    return found


def atoms_within_radius(mol: Chem.Mol, root: int, radius: int) -> List[int]:
    distance = {root: 0}
    queue = deque([root])
    while queue:
        current = queue.popleft()
        if distance[current] >= radius:
            continue
        for neighbor in mol.GetAtomWithIdx(current).GetNeighbors():
            index = int(neighbor.GetIdx())
            if index not in distance:
                distance[index] = distance[current] + 1
                queue.append(index)
    return sorted(distance)


def rooted_environment(mol: Chem.Mol, root: int, radius: int) -> str:
    clean = Chem.Mol(mol)
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    atom_indices = atoms_within_radius(clean, root, radius)
    atom_set = set(atom_indices)
    bond_indices = sorted(
        int(bond.GetIdx())
        for bond in clean.GetBonds()
        if int(bond.GetBeginAtomIdx()) in atom_set and int(bond.GetEndAtomIdx()) in atom_set
    )
    return Chem.MolFragmentToSmiles(
        clean,
        atomsToUse=atom_indices,
        bondsToUse=bond_indices,
        rootedAtAtom=root,
        canonical=True,
        isomericSmiles=True,
        allBondsExplicit=True,
        allHsExplicit=True,
    )


def stable_bin(token: str) -> int:
    digest = hashlib.sha256(token.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big", signed=False) % HASH_BINS_PER_SIDE


def environment_tokens(mapped_side: str, changed_maps: Sequence[int], side_name: str) -> Counter:
    located = locate_maps(parse_side(mapped_side))
    tokens: Counter = Counter()
    for atom_map in sorted(set(int(x) for x in changed_maps)):
        location = located.get(atom_map)
        if location is None:
            continue
        mol, root = location
        for radius in RADII:
            descriptor = rooted_environment(mol, root, radius)
            tokens[f"{side_name}|radius={radius}|{descriptor}"] += 1
    return tokens


def make_row(mapping: Dict) -> Dict:
    reaction_key = mapping["reaction_key"]
    status = mapping.get("status", "unavailable")
    centers = mapping.get("centers") or {}
    changed_maps = sorted(set(int(x) for x in centers.get("changed_center_atom_maps", [])))
    valid_mapping = status in VALID_STATUSES and ">>" in mapping.get("mapped_reaction", "")
    center_available = bool(valid_mapping and changed_maps)
    substrate_tokens: Counter = Counter()
    product_tokens: Counter = Counter()
    if center_available:
        substrate, product = mapping["mapped_reaction"].split(">>", 1)
        substrate_tokens = environment_tokens(substrate, changed_maps, "substrate")
        product_tokens = environment_tokens(product, changed_maps, "product")

    vector = np.zeros(len(FEATURE_NAMES), dtype=np.float32)
    for token, count in substrate_tokens.items():
        vector[stable_bin(token)] += float(count)
    offset = HASH_BINS_PER_SIDE
    for token, count in product_tokens.items():
        vector[offset + stable_bin(token)] += float(count)
    added = len(centers.get("added_bonds", [])) if center_available else 0
    removed = len(centers.get("removed_bonds", [])) if center_available else 0
    changed = len(centers.get("changed_bonds", [])) if center_available else 0
    vector[-3:] = (added, removed, changed)
    return {
        "reaction_key": reaction_key,
        "status": status,
        "mapping_warning": bool(status == "mapped_validation_warning"),
        "center_available": center_available,
        "changed_center_atom_maps": changed_maps,
        "changed_center_atom_count": len(changed_maps),
        "substrate_environment_tokens": dict(sorted(substrate_tokens.items())),
        "product_environment_tokens": dict(sorted(product_tokens.items())),
        "substrate_token_count": int(sum(substrate_tokens.values())),
        "product_token_count": int(sum(product_tokens.values())),
        "directed_added_bond_count": int(added),
        "directed_removed_bond_count": int(removed),
        "bond_order_change_count": int(changed),
        "feature_nonzero_count": int(np.count_nonzero(vector)),
        "feature_sum": float(vector.sum()),
        "feature_vector": vector,
    }


def collision_summary(token_sets: Dict[str, set]) -> Dict:
    result = {}
    for side, tokens in sorted(token_sets.items()):
        occupancy: Dict[int, List[str]] = defaultdict(list)
        for token in sorted(tokens):
            occupancy[stable_bin(token)].append(token)
        occupied = len(occupancy)
        colliding_bins = sum(len(values) > 1 for values in occupancy.values())
        result[side] = {
            "unique_tokens": len(tokens),
            "occupied_bins": occupied,
            "colliding_bins": colliding_bins,
            "tokens_in_colliding_bins": sum(len(values) for values in occupancy.values() if len(values) > 1),
            "maximum_unique_tokens_in_one_bin": max((len(values) for values in occupancy.values()), default=0),
        }
    return result


def run(root: Path, output_name: str) -> None:
    mapping_path = root / "reaction_mapping_01/mapping_results.jsonl"
    mapping_summary_path = root / "reaction_mapping_01/summary.json"
    core_path = root / "dataset_02/core_reactions.json"
    mappings = [json.loads(line) for line in mapping_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    core = json.loads(core_path.read_text(encoding="utf-8"))
    by_key = {row["reaction_key"]: row for row in mappings}
    if len(by_key) != len(mappings) or set(by_key) != set(core):
        raise ValueError("mapping rows and frozen candidate catalogue do not match one-to-one")
    output = root / output_name
    if output.exists():
        raise FileExistsError(output)
    output.mkdir()

    rows = []
    vectors = []
    token_sets = {"substrate": set(), "product": set()}
    for reaction_key in sorted(core):
        row = make_row(by_key[reaction_key])
        vectors.append(row.pop("feature_vector"))
        token_sets["substrate"].update(row["substrate_environment_tokens"])
        token_sets["product"].update(row["product_environment_tokens"])
        rows.append(row)

    matrix = np.asarray(vectors, dtype=np.float32)
    np.savez_compressed(
        output / "environment_features.npz",
        reaction_keys=np.asarray([row["reaction_key"] for row in rows]),
        feature_names=np.asarray(FEATURE_NAMES),
        feature_vector=matrix,
        center_available=np.asarray([row["center_available"] for row in rows], dtype=bool),
        mapping_warning=np.asarray([row["mapping_warning"] for row in rows], dtype=bool),
        changed_center_atom_count=np.asarray([row["changed_center_atom_count"] for row in rows], dtype=np.int32),
    )
    with (output / "environment_records.jsonl").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + "\n")

    center_counts = [row["changed_center_atom_count"] for row in rows]
    audit = {
        "created_utc": now(),
        "status": "FEATURE_BUILD_COMPLETE_PENDING_INDEPENDENT_VALIDATION",
        "reaction_count": len(rows),
        "feature_shape": list(matrix.shape),
        "hash_bins_per_side": HASH_BINS_PER_SIDE,
        "hash_name": HASH_NAME,
        "radii": list(RADII),
        "center_available_count": int(sum(row["center_available"] for row in rows)),
        "center_unavailable_count": int(sum(not row["center_available"] for row in rows)),
        "mapping_warning_count": int(sum(row["mapping_warning"] for row in rows)),
        "center_atom_count_min": int(min(center_counts)),
        "center_atom_count_max": int(max(center_counts)),
        "center_atom_count_mean": float(np.mean(center_counts)),
        "collision_summary": collision_summary(token_sets),
        "nonfinite_feature_cells": int(np.size(matrix) - np.isfinite(matrix).sum()),
        "input_manifest": {
            "dataset_02/core_reactions.json": digest_file(core_path),
            "reaction_mapping_01/mapping_results.jsonl": digest_file(mapping_path),
            "reaction_mapping_01/summary.json": digest_file(mapping_summary_path),
        },
        "output_manifest": {},
    }
    for name in ("environment_features.npz", "environment_records.jsonl"):
        audit["output_manifest"][name] = digest_file(output / name)
    write_json(output / "audit.json", audit)
    runtime = {
        "created_utc": now(),
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "source_code_sha256": digest_file(Path(__file__)),
        "packages": {},
    }
    for package in ("rdkit", "numpy"):
        try:
            runtime["packages"][package] = version(package)
        except PackageNotFoundError:
            runtime["packages"][package] = None
    write_json(output / "runtime.json", runtime)
    print(json.dumps({"output": str(output), **{k: audit[k] for k in ("reaction_count", "feature_shape", "center_available_count", "mapping_warning_count", "collision_summary")}}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-name", default="reaction_center_environments_01")
    args = parser.parse_args()
    run(args.root.resolve(), args.output_name)
