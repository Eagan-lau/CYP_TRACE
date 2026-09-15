"""Independent algebraic audit of the fixed reaction transform."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def run(root: Path, output_name: str):
    source = root / "interaction_transform_01"
    output = root / output_name
    if output.exists():
        raise FileExistsError(output)
    output.mkdir()
    declared = json.loads((source / "audit.json").read_text())
    archive = np.load(source / "reaction_transform.npz", allow_pickle=False)
    environments = np.load(root / "reaction_center_environments_01/environment_features.npz", allow_pickle=False)
    representation = np.load(root / "main_baselines_01/fresh_representations.npz", allow_pickle=False)
    raw = environments["feature_vector"].astype(np.float64)
    rownorm = raw / np.maximum(np.linalg.norm(raw, axis=1, keepdims=True), 1e-12)
    centered = rownorm - rownorm.mean(axis=0)
    scale = float(np.trace(centered @ centered.T) / len(centered))
    standardized = centered / np.sqrt(scale)
    reconstruction = (archive["u"] * archive["singular_values"][None, :]) @ archive["vt"]
    kernel_u = representation["chemical_kernel"].astype(np.float64) @ archive["u"]
    differences = {
        "normalized_rows": float(np.max(np.abs(rownorm - archive["raw_row_normalized"]))),
        "center": float(np.max(np.abs(rownorm.mean(axis=0) - archive["center"]))),
        "scale": abs(scale - float(archive["scale"])),
        "svd_reconstruction": float(np.max(np.abs(standardized - reconstruction))),
        "left_orthonormality": float(np.max(np.abs(archive["u"].T @ archive["u"] - np.eye(len(archive["singular_values"]))))),
        "right_orthonormality": float(np.max(np.abs(archive["vt"] @ archive["vt"].T - np.eye(len(archive["singular_values"]))))),
        "chemical_kernel_times_u": float(np.max(np.abs(kernel_u - archive["chemical_kernel_times_u"]))),
    }
    failures = [name for name, value in differences.items() if value > (1e-9 if name != "svd_reconstruction" else 1e-8)]
    if list(archive["reaction_keys"]) != list(environments["reaction_keys"]) or list(archive["reaction_keys"]) != list(representation["reaction_ids"]):
        failures.append("reaction_order")
    if digest(source / "reaction_transform.npz") != declared["output_sha256"]:
        failures.append("output_sha256")
    report = {
        "status": "PASS" if not failures else "FAIL",
        "reaction_count": len(archive["reaction_keys"]),
        "numerical_rank": len(archive["singular_values"]),
        "differences": differences,
        "failures": failures,
        "independent_implementation": True,
    }
    (output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-name", default="interaction_transform_validation_01")
    args = parser.parse_args()
    run(args.root.resolve(), args.output_name)
