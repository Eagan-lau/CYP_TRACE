"""Prepare the fixed catalogue SVD used by all interaction fits."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from interaction_core import normalized_rows
from run_raw import digest_file, now, write_json


def run(root: Path, output_name: str) -> None:
    output = root / output_name
    if output.exists():
        raise FileExistsError(output)
    feature_dir = root / "reaction_center_environments_01"
    feature_qc = root / "reaction_center_environment_validation_01/audit.json"
    qc = json.loads(feature_qc.read_text())
    if qc["status"] != "PASS":
        raise ValueError("reaction-centre environment QC gate")
    archive = np.load(feature_dir / "environment_features.npz", allow_pickle=False)
    representation = np.load(root / "main_baselines_01/fresh_representations.npz", allow_pickle=False)
    reaction_keys = archive["reaction_keys"]
    if list(reaction_keys) != list(representation["reaction_ids"]):
        raise ValueError("reaction catalogue order mismatch")
    raw = archive["feature_vector"].astype(np.float64)
    z = normalized_rows(raw)
    center = z.mean(axis=0)
    centered = z - center
    scale = float(np.trace(centered @ centered.T) / len(centered))
    if scale < 1e-12:
        raise ValueError("reaction feature transform has no variation")
    zs = centered / np.sqrt(scale)
    u, singular, vt = np.linalg.svd(zs, full_matrices=False)
    tolerance = float(np.finfo(np.float64).eps * max(zs.shape) * singular[0])
    keep = singular > tolerance
    u, singular, vt = u[:, keep], singular[keep], vt[keep]
    kernel = representation["chemical_kernel"].astype(np.float64)
    kernel_u = kernel @ u
    output.mkdir()
    np.savez_compressed(
        output / "reaction_transform.npz",
        reaction_keys=reaction_keys,
        feature_names=archive["feature_names"],
        raw_row_normalized=z,
        center=center,
        scale=np.asarray(scale),
        u=u,
        singular_values=singular,
        vt=vt,
        numerical_rank_tolerance=np.asarray(tolerance),
        chemical_kernel_times_u=kernel_u,
        center_available=archive["center_available"],
        mapping_warning=archive["mapping_warning"],
    )
    reconstruction = (u * singular[None, :]) @ vt
    audit = {
        "status": "TRANSFORM_COMPLETE_INDEPENDENT_QC_PENDING",
        "created_utc": now(),
        "reaction_count": len(reaction_keys),
        "raw_feature_dimension": raw.shape[1],
        "numerical_rank": len(singular),
        "numerical_rank_tolerance": tolerance,
        "row_normalization": "l2_with_1e-12_floor",
        "catalogue_centering": True,
        "catalogue_scale": scale,
        "maximum_svd_reconstruction_difference": float(np.max(np.abs(reconstruction - zs))),
        "labels_or_performance_used": False,
        "input_manifest": {
            "reaction_center_environments_01/environment_features.npz": digest_file(feature_dir / "environment_features.npz"),
            "reaction_center_environment_validation_01/audit.json": digest_file(feature_qc),
            "main_baselines_01/fresh_representations.npz": digest_file(root / "main_baselines_01/fresh_representations.npz"),
            "prepare_interaction_transform.py": digest_file(Path(__file__)),
            "interaction_core.py": digest_file(root / "interaction_core.py"),
        },
    }
    audit["output_sha256"] = digest_file(output / "reaction_transform.npz")
    write_json(output / "audit.json", audit)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-name", default="interaction_transform_01")
    args = parser.parse_args()
    run(args.root.resolve(), args.output_name)
