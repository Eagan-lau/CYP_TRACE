"""Independent reconstruction of both deterministic control workspaces."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np

from run_raw import digest_file, now, write_json


def protein_check(root: Path, failures: list) -> dict:
    workspace = root / "interaction_control_protein_position_permutation_01"
    source = np.load(root / "ordered_site_features_01/ordered_features.npz", allow_pickle=False)
    control = np.load(workspace / "ordered_site_features_01/ordered_features.npz", allow_pickle=False)
    declared = np.load(workspace / "protein_position_permutations.npz", allow_pickle=False)
    cell_differences = 0; invariant_failures = 0
    for row, sequence in enumerate(source["sequence_ids"].tolist()):
        order = sorted(range(len(source["reference_positions_0based"])), key=lambda index: hashlib.sha256(f"protein_position_permutation_v1|{sequence}|{index}".encode()).hexdigest())
        if declared["permutations"][row].tolist() != order: failures.append("protein_permutation_index")
        for channel in ("sequence_projected", "structure_projected"):
            for suffix in ("residue_codes", "present_mask", "onehot"):
                key = channel + "_" + suffix; expected = source[key][row, order]
                cell_differences += int(np.count_nonzero(expected != control[key][row]))
            if int(source[channel + "_present_mask"][row].sum()) != int(control[channel + "_present_mask"][row].sum()): invariant_failures += 1
            if not np.array_equal(source[channel + "_onehot"][row].sum(axis=0), control[channel + "_onehot"][row].sum(axis=0)): invariant_failures += 1
    for key in source.files:
        if not any(key.endswith(suffix) for suffix in ("residue_codes", "present_mask", "onehot")) and not np.array_equal(source[key], control[key]): failures.append("protein_unexpected_change:" + key)
    if cell_differences: failures.append("protein_feature_cells")
    if invariant_failures: failures.append("protein_invariants")
    return {"feature_cell_differences": cell_differences, "invariant_failures": invariant_failures,
            "output_sha256": digest_file(workspace / "ordered_site_features_01/ordered_features.npz")}


def reaction_check(root: Path, failures: list) -> dict:
    workspace = root / "interaction_control_reaction_center_row_permutation_01"
    source = np.load(root / "interaction_transform_01/reaction_transform.npz", allow_pickle=False)
    control = np.load(workspace / "interaction_transform_01/reaction_transform.npz", allow_pickle=False)
    declared = np.load(workspace / "reaction_center_row_permutation.npz", allow_pickle=False)
    keys = source["reaction_keys"].tolist(); counts = declared["changed_center_atom_count"].astype(int); warnings = source["mapping_warning"].astype(bool)
    strata = defaultdict(list)
    for index, value in enumerate(zip(counts, warnings)): strata[(int(value[0]), bool(value[1]))].append(index)
    donor = np.arange(len(keys), dtype=int)
    for indices in strata.values():
        ordered = sorted(indices, key=lambda index: hashlib.sha256(f"reaction_center_row_permutation_v1|{keys[index]}".encode()).hexdigest())
        if len(ordered) > 1:
            for position, recipient in enumerate(ordered): donor[recipient] = ordered[(position + 1) % len(ordered)]
    if not np.array_equal(donor, declared["donor_indices"]): failures.append("reaction_permutation_index")
    differences = {
        "u": float(np.max(np.abs(source["u"][donor] - control["u"]))),
        "raw_row_normalized": float(np.max(np.abs(source["raw_row_normalized"][donor] - control["raw_row_normalized"]))),
    }
    kernel = np.load(root / "main_baselines_01/fresh_representations.npz", allow_pickle=False)["chemical_kernel"].astype(float)
    differences["chemical_kernel_times_u"] = float(np.max(np.abs(kernel @ control["u"].astype(float) - control["chemical_kernel_times_u"])))
    for key in source.files:
        if key not in ("u", "raw_row_normalized", "chemical_kernel_times_u") and not np.array_equal(source[key], control[key]): failures.append("reaction_unexpected_change:" + key)
    if any(value > 1e-10 for value in differences.values()): failures.append("reaction_transform_cells")
    return {"differences": differences, "strata": len(strata), "fixed_candidate_rows": int(np.sum(donor == np.arange(len(donor)))),
            "output_sha256": digest_file(workspace / "interaction_transform_01/reaction_transform.npz")}


if __name__ == "__main__":
    root = Path(__file__).resolve().parent; output = root / "interaction_control_validation_01"
    if output.exists(): raise FileExistsError(output)
    failures = []; protein = protein_check(root, failures); reaction = reaction_check(root, failures)
    for name, details in [("interaction_control_protein_position_permutation_01", protein), ("interaction_control_reaction_center_row_permutation_01", reaction)]:
        receipt = json.loads((root / name / "CONTROL_RECEIPT.json").read_text())
        if receipt["output_sha256"] != details["output_sha256"]: failures.append("receipt_hash:" + name)
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
              "protein_position_permutation": protein, "reaction_center_row_permutation": reaction,
              "independent_implementation": True, "labels_or_scores_used": False}
    output.mkdir(); write_json(output / "audit.json", report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)
