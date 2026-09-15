"""Independent end-to-end verification of the installable CYP-TRACE strategy."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
TOOL = ROOT / "joint_tool_01"
OUTPUT = Path(os.environ.get("CYPTRACE_TOOL_VALIDATION_OUTPUT", ROOT / "joint_tool_validation_01")).resolve()
sys.path.insert(0, str(TOOL / "src"))

from cyptrace_pipeline.evidence import build_evidence_index, load_evidence_index, screen  # noqa: E402
from cyptrace_pipeline.human import HumanSubstrateModel  # noqa: E402


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def content_digest(bundle: dict) -> str:
    copy = dict(bundle)
    copy.pop("bundle_content_sha256", None)
    return hashlib.sha256(json.dumps(copy, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def run(command: list[str], env: dict | None = None) -> subprocess.CompletedProcess:
    result = subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stdout}\n{result.stderr}")
    return result


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to replace {OUTPUT}")
    OUTPUT.mkdir()
    checks = {}

    asset_manifest = json.loads((TOOL / "ASSET_MANIFEST.json").read_text(encoding="utf-8"))
    model_path = TOOL / "src" / "cyptrace_pipeline" / "data" / "human_model_v1.json.gz"
    with gzip.open(model_path, "rt", encoding="utf-8") as stream:
        bundle = json.load(stream)
    checks["model_file_hash"] = digest(model_path) == asset_manifest["model_bundle_sha256"]
    checks["model_content_hash"] = content_digest(bundle) == bundle["bundle_content_sha256"] == asset_manifest["bundle_content_sha256"]
    checks["model_source_hash"] = digest(ROOT / "human_substrate_01" / "normalized_labels.json") == bundle["source"]["normalized_labels_sha256"]

    model = HumanSubstrateModel.load(model_path)
    external = json.loads((ROOT / "external_human_cyp_01" / "normalized_external_labels.json").read_text(encoding="utf-8"))
    structures = {row["compound_inchikey"]: row["canonical_smiles"] for row in external}
    expected_path = ROOT / "external_human_cyp_models_01" / "external_predictions.jsonl"
    tested = 0
    max_specific_difference = 0.0
    max_pooled_difference = 0.0
    route_counts = {}
    for line in expected_path.read_text(encoding="utf-8").splitlines():
        expected = json.loads(line)
        actual = model.score(structures[expected["compound_inchikey"]], expected["isoform"])
        max_specific_difference = max(max_specific_difference, abs(actual["score"] - expected["isoform_specific_knn"]))
        max_pooled_difference = max(max_pooled_difference, abs(actual["pooled_chemistry_score"] - expected["pooled_chemical_knn"]))
        route_counts[actual["route"]] = route_counts.get(actual["route"], 0) + 1
        tested += 1
    checks["all_external_rows_recomputed"] = tested == 14519
    checks["external_specific_scores_exact"] = max_specific_difference <= 1e-12
    checks["external_pooled_scores_exact"] = max_pooled_difference <= 1e-12

    evidence_path = OUTPUT / "local_evidence_index_v1.json.gz"
    evidence_build = build_evidence_index(ROOT / "dataset_02", evidence_path)
    evidence = load_evidence_index(evidence_path)
    checks["evidence_counts"] = evidence["counts"] == {"sequences": 600, "reactions": 1341, "edges": 2304}
    exact_edges = 0
    for edge_key in evidence["edges"]:
        sequence_hash, reaction_key = edge_key.split("|", 1)
        result = screen(evidence, {"q": evidence["sequences"][sequence_hash]}, [{"candidate_id": "c", "reaction_key": reaction_key}])[0]
        if result["route"] == "DOCUMENTED_EXACT_EVIDENCE" and result["evidence"]:
            exact_edges += 1
    checks["all_exact_edge_keys_retrievable"] = exact_edges == len(evidence["edges"])
    first_edge = sorted(evidence["edges"])[0]
    sequence_hash, reaction_key = first_edge.split("|", 1)
    sequence = evidence["sequences"][sequence_hash]
    replacement = "A" if sequence[0] != "A" else "G"
    novel = screen(evidence, {"novel": replacement + sequence[1:]}, [{"candidate_id": "c", "reaction_key": reaction_key}])[0]
    checks["new_sequence_abstains"] = novel["route"] == "ABSTAIN_UNVALIDATED_GENERAL_REACTION_ROUTE" and novel["score"] is None
    unknown_reaction = screen(evidence, {"known": sequence}, [{"candidate_id": "c", "reaction_key": "MAIN:not_present"}])[0]
    checks["unknown_reaction_abstains"] = unknown_reaction["route"] == "ABSTAIN_UNVALIDATED_GENERAL_REACTION_ROUTE" and unknown_reaction["score"] is None
    structure_edge = next(
        key for key in sorted(evidence["edges"])
        if evidence["reactions"][key.split("|", 1)[1]].get("single_pair")
        and len(evidence["reactions"][key.split("|", 1)[1]].get("substrates", [])) == 1
        and len(evidence["reactions"][key.split("|", 1)[1]].get("products", [])) == 1
    )
    structure_sequence_hash, structure_reaction_key = structure_edge.split("|", 1)
    reaction = evidence["reactions"][structure_reaction_key]
    structure_match = screen(evidence, {"known": evidence["sequences"][structure_sequence_hash]}, [{
        "candidate_id": "structure_pair",
        "substrate_smiles": reaction["substrates"][0],
        "product_smiles": reaction["products"][0],
    }])[0]
    checks["structure_pair_resolves_exact_reaction"] = structure_match["reaction_key"] == structure_reaction_key and structure_match["route"] == "DOCUMENTED_EXACT_EVIDENCE"

    install = OUTPUT / "install_target"
    pip = run([sys.executable, "-m", "pip", "install", "--no-deps", "--no-build-isolation", "--target", str(install), str(TOOL)])
    environment = os.environ.copy()
    inherited_pythonpath = environment.get("PYTHONPATH", "")
    environment["PYTHONPATH"] = str(install) + (os.pathsep + inherited_pythonpath if inherited_pythonpath else "")
    doctor = run([sys.executable, "-m", "cyptrace_pipeline", "doctor"], environment)
    doctor_json = json.loads(doctor.stdout)
    smoke = run([
        sys.executable, "-m", "cyptrace_pipeline", "human-substrate", "--smiles", "CCO", "--isoform", "CYP3A4"
    ], environment)
    smoke_json = json.loads(smoke.stdout)
    checks["clean_install_doctor"] = doctor_json["status"] == "PASS" and doctor_json["human_model"]["compounds"] == 1768
    checks["clean_install_prediction"] = smoke_json["results"][0]["route"] == "FIXED_HUMAN_EXTERNALLY_VALIDATED"
    (OUTPUT / "install_stdout.txt").write_text(pip.stdout + pip.stderr, encoding="utf-8")
    (OUTPUT / "doctor.json").write_text(json.dumps(doctor_json, indent=2) + "\n", encoding="utf-8")
    (OUTPUT / "human_smoke.json").write_text(json.dumps(smoke_json, indent=2) + "\n", encoding="utf-8")

    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "external_rows_recomputed": tested,
        "external_score_max_absolute_difference": {
            "isoform_specific": max_specific_difference,
            "pooled": max_pooled_difference,
        },
        "external_route_counts": route_counts,
        "evidence_build": evidence_build,
        "exact_edge_keys_retrieved": exact_edges,
        "input_sha256": {
            "human_model": digest(model_path),
            "external_predictions": digest(expected_path),
            "dataset_core_edges": digest(ROOT / "dataset_02" / "core_edges.json"),
            "dataset_core_reactions": digest(ROOT / "dataset_02" / "core_reactions.json"),
            "dataset_core_sequences": digest(ROOT / "dataset_02" / "core_sequences.fasta"),
        },
        "claim_boundary": "Fixed-human new-chemistry scores are reproduced; exact general evidence is retrievable; unseen-sequence and unseen-reaction cases abstain.",
    }
    report["validation_script_sha256"] = digest(Path(__file__))
    (OUTPUT / "validation.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if report["status"] != "PASS":
        raise AssertionError("joint tool verification failed")


if __name__ == "__main__":
    main()
