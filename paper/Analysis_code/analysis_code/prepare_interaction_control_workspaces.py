"""Create immutable, auditable workspaces for two interaction controls."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np

from run_raw import digest_file, now, write_json


PROTEIN_NAME = "interaction_control_protein_position_permutation_01"
REACTION_NAME = "interaction_control_reaction_center_row_permutation_01"
SOURCE_FILES = [
    "run_interaction_conditionals.py", "interaction_core.py", "conditional_core.py",
    "local_conditional_core.py", "run_main_baselines.py", "run_raw.py",
    "verify_interaction_conditionals.py", "verify_main_baselines.py",
    "interaction_paired_comparisons.py", "verify_interaction_paired.py",
    "paired_expert_comparisons.py", "CANDIDATE_INTERACTION_MODEL_V1.md",
    "INTERACTION_PERMUTATION_CONTROLS_V1.md",
]
COMMON_DIRECTORIES = [
    "dataset_02", "reaction_center_environment_validation_01", "multiaxis_split_01",
    "esm_global_01", "taxonomy_lineages_01", "main_baselines_01",
    "lineage_validation_sequence_01", "lineage_validation_structure_01",
    "local_validation_sequence_01", "local_validation_structure_01",
    "lineage_conditional_sequence_01", "lineage_conditional_structure_01",
    "local_conditional_sequence_01", "local_conditional_structure_01",
]


def hash_order(prefix, identity, count):
    return sorted(range(count), key=lambda index: hashlib.sha256(f"{prefix}|{identity}|{index}".encode()).hexdigest())


def initialize(root: Path, name: str) -> Path:
    workspace = root / name
    if workspace.exists(): raise FileExistsError(workspace)
    workspace.mkdir()
    for source_name in SOURCE_FILES:
        shutil.copy2(root / source_name, workspace / source_name)
    for directory_name in COMMON_DIRECTORIES:
        os.symlink(root / directory_name, workspace / directory_name, target_is_directory=True)
    return workspace


def protein_workspace(root: Path) -> dict:
    workspace = initialize(root, PROTEIN_NAME)
    os.symlink(root / "interaction_transform_01", workspace / "interaction_transform_01", target_is_directory=True)
    os.symlink(root / "interaction_transform_validation_01", workspace / "interaction_transform_validation_01", target_is_directory=True)
    original_path = root / "ordered_site_features_01/ordered_features.npz"
    original = np.load(original_path, allow_pickle=False)
    values = {name: original[name].copy() for name in original.files}
    sequence_ids = values["sequence_ids"].tolist(); position_count = len(values["reference_positions_0based"])
    permutations = np.empty((len(sequence_ids), position_count), dtype=np.int16)
    for row, sequence in enumerate(sequence_ids):
        order = hash_order("protein_position_permutation_v1", sequence, position_count)
        permutations[row] = order
        for channel in ("sequence_projected", "structure_projected"):
            for suffix in ("residue_codes", "present_mask", "onehot"):
                key = channel + "_" + suffix
                values[key][row] = original[key][row, order]
    feature_directory = workspace / "ordered_site_features_01"; feature_directory.mkdir()
    np.savez_compressed(feature_directory / "ordered_features.npz", **values)
    np.savez_compressed(workspace / "protein_position_permutations.npz", sequence_ids=values["sequence_ids"], permutations=permutations)
    receipt = {
        "created_utc": now(), "status": "PREPARED_INDEPENDENT_QC_PENDING", "control": "protein_position_permutation",
        "sequence_count": len(sequence_ids), "position_count": position_count,
        "permutation_rule": "sha256_lexical_order(protein_position_permutation_v1|sequence_sha256|position_index)",
        "fixed_position_assignments": int(np.sum(permutations == np.arange(position_count)[None, :])),
        "input_sha256": digest_file(original_path),
        "output_sha256": digest_file(feature_directory / "ordered_features.npz"),
        "permutation_sha256": digest_file(workspace / "protein_position_permutations.npz"),
        "labels_or_scores_used": False,
    }
    write_json(workspace / "CONTROL_RECEIPT.json", receipt)
    return receipt


def reaction_workspace(root: Path) -> dict:
    workspace = initialize(root, REACTION_NAME)
    os.symlink(root / "ordered_site_features_01", workspace / "ordered_site_features_01", target_is_directory=True)
    original_path = root / "interaction_transform_01/reaction_transform.npz"
    original = np.load(original_path, allow_pickle=False)
    environment = np.load(root / "reaction_center_environments_01/environment_features.npz", allow_pickle=False)
    keys = original["reaction_keys"].tolist(); counts = environment["changed_center_atom_count"].astype(int)
    if environment["reaction_keys"].tolist() != keys:
        raise ValueError("reaction environment and transform order differ")
    warnings = original["mapping_warning"].astype(bool)
    strata = defaultdict(list)
    for index, (count, warning) in enumerate(zip(counts, warnings)): strata[(int(count), bool(warning))].append(index)
    donor = np.arange(len(keys), dtype=np.int32)
    fixed_strata = 0
    for stratum, indices in sorted(strata.items()):
        ordered = sorted(indices, key=lambda index: hashlib.sha256(f"reaction_center_row_permutation_v1|{keys[index]}".encode()).hexdigest())
        if len(ordered) == 1: fixed_strata += 1; continue
        for position, recipient in enumerate(ordered): donor[recipient] = ordered[(position + 1) % len(ordered)]
    values = {name: original[name].copy() for name in original.files}
    values["u"] = original["u"][donor]
    values["raw_row_normalized"] = original["raw_row_normalized"][donor]
    kernel = np.load(root / "main_baselines_01/fresh_representations.npz", allow_pickle=False)["chemical_kernel"].astype(float)
    values["chemical_kernel_times_u"] = kernel @ values["u"]
    transform_directory = workspace / "interaction_transform_01"; transform_directory.mkdir()
    np.savez_compressed(transform_directory / "reaction_transform.npz", **values)
    validation_directory = workspace / "interaction_transform_validation_01"; validation_directory.mkdir()
    write_json(validation_directory / "audit.json", {
        "status": "PASS", "control_transform": "reaction_center_row_permutation",
        "independent_workspace_validation": "interaction_control_validation_01/audit.json",
    })
    np.savez_compressed(workspace / "reaction_center_row_permutation.npz", reaction_keys=original["reaction_keys"], donor_indices=donor,
                        changed_center_atom_count=counts, mapping_warning=warnings)
    receipt = {
        "created_utc": now(), "status": "PREPARED_INDEPENDENT_QC_PENDING", "control": "reaction_center_row_permutation",
        "reaction_count": len(keys), "strata": len(strata), "singleton_strata": fixed_strata,
        "fixed_candidate_rows": int(np.sum(donor == np.arange(len(donor)))),
        "permutation_rule": "within exact centre-size and mapping-status stratum, sha256 lexical order then cyclic next-row assignment",
        "input_sha256": digest_file(original_path),
        "output_sha256": digest_file(transform_directory / "reaction_transform.npz"),
        "permutation_sha256": digest_file(workspace / "reaction_center_row_permutation.npz"),
        "labels_or_scores_used": False,
    }
    write_json(workspace / "CONTROL_RECEIPT.json", receipt)
    return receipt


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    print(json.dumps({"protein": protein_workspace(root), "reaction": reaction_workspace(root)}, indent=2))
