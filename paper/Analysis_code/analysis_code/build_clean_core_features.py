"""Build official CLEAN v1.0.0 features and maximum-separation calls for the fresh core."""
from __future__ import annotations

from collections import defaultdict
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import re

import numpy as np

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
OUTPUT = ROOT / "clean_core_features_01"
CACHE = OUTPUT / "esm1b_cache"
MAX_RESIDUES = 1022
DISTANCE_SCALE = 1_000_000
ASSET_HASHES = {
    "data/restricted/clean_v1_0_0/pretrained/split100.pth": "02ba4f54d44cef4e5c3ff0b94e7cafeeeaaa649cd45dc109ed701ee01b8fad99",
    "data/restricted/clean_v1_0_0/pretrained/100.pt": "a8dfe7e3401d356c771e84467d7649ee43ee80e2b688f0b4b9778536e693bcb8",
    "data/restricted/esm1b_v1/esm1b_t33_650M_UR50S.pt": "0569754efaff7dcb7e068c27367bc73f10afb4b450ea30aac30d9bc60783a8b1",
}


def fasta(path: Path) -> list[tuple[str, str]]:
    output = []
    identifier = None
    parts = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                output.append((identifier, "".join(parts)))
            identifier, parts = line[1:].split()[0], []
        else:
            if identifier is None:
                raise ValueError("sequence before FASTA identifier")
            parts.append(line)
    if identifier is not None:
        output.append((identifier, "".join(parts)))
    return output


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def parse_ecs(value: str) -> list[str]:
    output = []
    for token in re.split(r"[;|/]", str(value or "")):
        ec = token.strip()
        if ec.startswith("EC:"):
            ec = ec[3:].strip()
        if ec and ec not in output:
            output.append(ec)
    return output


def clean_model_class(path: Path):
    specification = importlib.util.spec_from_file_location("clean_official_model", path)
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot import CLEAN model")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module.LayerNormNet


def maximum_separation(distances: np.ndarray) -> int:
    """Reproduce CLEAN v1.0.0 first-gradient maximum-separation call count."""
    gamma_mean = float(np.mean(np.append(distances[1:], np.repeat(distances[-1], 10))))
    separation = np.abs(distances - gamma_mean)
    gradients = np.abs(separation[:-1] - separation[1:])
    large = np.where(gradients > np.mean(gradients))[0]
    index = int(large[0]) if len(large) else 0
    return 1 if index >= 5 else index + 1


def save_vector(path: Path, value: np.ndarray) -> None:
    temporary = path.with_suffix(".tmp.npy")
    np.save(temporary, value.astype(np.float32, copy=False), allow_pickle=False)
    temporary.replace(path)


def main() -> None:
    import torch
    import esm

    torch.set_num_threads(32)
    torch.manual_seed(0)
    OUTPUT.mkdir(exist_ok=True)
    CACHE.mkdir(exist_ok=True)
    cached_rows_at_start = len(list(CACHE.glob("*.npy")))
    final_marker = OUTPUT / "audit.json"
    if final_marker.exists():
        raise FileExistsError("completed CLEAN feature directory already exists")

    for relative, expected in ASSET_HASHES.items():
        path = PROJECT / relative
        if digest_file(path) != expected:
            raise ValueError("restricted asset hash: " + relative)

    core_fasta = ROOT / "dataset_02/core_sequences.fasta"
    records = fasta(core_fasta)
    if len(records) != 600 or len({identifier for identifier, _ in records}) != len(records):
        raise ValueError("fresh core FASTA cardinality")
    for identifier, sequence in records:
        if sequence_hash(sequence) != identifier:
            raise ValueError("fresh core sequence identifier is not its SHA256")
    records.sort()
    supported = [(identifier, sequence) for identifier, sequence in records if len(sequence) <= MAX_RESIDUES]

    old_dir = PROJECT / "data/processed/precutoff_clean_baseline_v1"
    old_manifest = list(csv.DictReader((old_dir / "query_manifest.tsv").open(newline=""), delimiter="\t"))
    old_fasta = fasta(old_dir / "unique_query_sequences.fasta")
    old_id_to_hash = {identifier: sequence_hash(sequence) for identifier, sequence in old_fasta}
    old_hash_to_index = {}
    for row in old_manifest:
        if row["clean_input_supported"].lower() != "true":
            continue
        index = int(row["clean_embedding_index"])
        sequence = row["sequence_sha256"]
        query_id = row["clean_query_id"]
        if old_id_to_hash.get(query_id) != sequence:
            raise ValueError("old CLEAN FASTA/manifest mismatch")
        if sequence in old_hash_to_index and old_hash_to_index[sequence] != index:
            raise ValueError("old CLEAN sequence index conflict")
        old_hash_to_index[sequence] = index
    old_esm = np.load(old_dir / "query_esm1b_embeddings.npy", allow_pickle=False)
    if old_esm.shape != (266, 1280) or old_esm.dtype != np.float32 or not np.isfinite(old_esm).all():
        raise ValueError("old independently verified ESM matrix")

    reused = 0
    computed = 0
    missing = []
    for identifier, _ in supported:
        target = CACHE / f"{identifier}.npy"
        if target.exists():
            value = np.load(target, allow_pickle=False)
            if value.shape != (1280,) or value.dtype != np.float32 or not np.isfinite(value).all():
                raise ValueError("invalid resumable ESM row: " + identifier)
        elif identifier in old_hash_to_index:
            save_vector(target, old_esm[old_hash_to_index[identifier]])
            reused += 1
        else:
            missing.append((identifier, dict(supported)[identifier]))
    print(json.dumps({"supported": len(supported), "cached_or_reused": len(supported) - len(missing),
                      "new_esm1b_rows": len(missing), "old_rows_newly_reused": reused}), flush=True)

    if missing:
        model_path = PROJECT / "data/restricted/esm1b_v1/esm1b_t33_650M_UR50S.pt"
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(model_path))
        model.eval().to(torch.device("cpu"))
        converter = alphabet.get_batch_converter()
        with torch.no_grad():
            for index, (identifier, sequence) in enumerate(missing, start=1):
                _, _, tokens = converter([(identifier, sequence)])
                result = model(tokens, repr_layers=[33], return_contacts=False)
                vector = result["representations"][33][0, 1:len(sequence) + 1].mean(0).cpu().numpy()
                save_vector(CACHE / f"{identifier}.npy", vector)
                computed += 1
                if index % 10 == 0 or index == len(missing):
                    print(f"ESM-1b newly computed {index}/{len(missing)}", flush=True)
        del model

    esm_matrix = np.stack([np.load(CACHE / f"{identifier}.npy", allow_pickle=False)
                           for identifier, _ in supported]).astype(np.float32, copy=False)
    np.save(OUTPUT / "esm1b_embeddings.npy", esm_matrix, allow_pickle=False)

    model_source = PROJECT / "tools/external/CLEAN_v1_0_0/src/CLEAN/model.py"
    LayerNormNet = clean_model_class(model_source)
    clean_model = LayerNormNet(512, 128, torch.device("cpu"), torch.float32)
    checkpoint = torch.load(PROJECT / "data/restricted/clean_v1_0_0/pretrained/split100.pth", map_location="cpu")
    clean_model.load_state_dict(checkpoint)
    clean_model.eval()
    with torch.no_grad():
        clean_matrix = clean_model(torch.from_numpy(esm_matrix)).cpu().numpy().astype(np.float32)
    np.save(OUTPUT / "clean_embeddings.npy", clean_matrix, allow_pickle=False)

    split_path = PROJECT / "tools/external/CLEAN_v1_0_0/data/split100.csv"
    split_rows = list(csv.DictReader(split_path.open(newline=""), delimiter="\t"))
    order = []
    seen = set()
    counts = defaultdict(int)
    clean_entries = set()
    clean_sequences = set()
    for row in split_rows:
        clean_entries.add(row["Entry"])
        clean_sequences.add(row["Sequence"])
        for ec in parse_ecs(row["EC number"]):
            if ec not in seen:
                order.append(ec); seen.add(ec)
            counts[ec] += 1
    if len(order) != 5242 or sum(counts.values()) != 241025:
        raise ValueError("official CLEAN catalog cardinality")
    train_embeddings = torch.load(PROJECT / "data/restricted/clean_v1_0_0/pretrained/100.pt", map_location="cpu")
    train_embeddings = torch.as_tensor(train_embeddings, dtype=torch.float32)
    if tuple(train_embeddings.shape) != (241025, 128):
        raise ValueError("official CLEAN training embedding shape")
    centers = np.empty((len(order), 128), dtype=np.float32)
    offset = 0
    for index, ec in enumerate(order):
        count = counts[ec]
        centers[index] = train_embeddings[offset:offset + count].mean(0).detach().numpy()
        offset += count
    np.save(OUTPUT / "ec_centers.npy", centers, allow_pickle=False)
    write_json(OUTPUT / "ec_catalog.json", [{"ec_index": index, "ec": ec,
                                               "training_assignment_count": counts[ec]}
                                              for index, ec in enumerate(order)])

    query_tensor = torch.from_numpy(clean_matrix)
    center_tensor = torch.from_numpy(centers)
    units = np.empty((len(query_tensor), len(center_tensor)), dtype=np.int32)
    with torch.no_grad():
        for start in range(0, len(query_tensor), 16):
            stop = min(start + 16, len(query_tensor))
            value = torch.cdist(query_tensor[start:stop], center_tensor).double() * DISTANCE_SCALE
            units[start:stop] = torch.round(value).clamp(0, np.iinfo(np.int32).max).numpy().astype(np.int32)
    np.save(OUTPUT / "ec_distance_units.npy", units, allow_pickle=False)

    supported_index = {identifier: index for index, (identifier, _) in enumerate(supported)}
    calls = {}
    for identifier, index in supported_index.items():
        nearest = np.argsort(units[index], kind="stable")[:10]
        count = maximum_separation(units[index, nearest].astype(np.float64))
        calls[identifier] = [{"ec": order[int(ec_index)], "distance_units": int(units[index, ec_index])}
                             for ec_index in nearest[:count]]
    write_json(OUTPUT / "maximum_separation_calls.json", calls)

    edges = json.loads((ROOT / "dataset_02/core_edges.json").read_text())
    accessions = defaultdict(set)
    for edge in edges:
        accessions[edge["sequence_sha256"]].update(edge["accessions"])
    manifest = []
    for identifier, sequence in records:
        row_accessions = sorted(accessions[identifier])
        exact_accession = any(value in clean_entries for value in row_accessions)
        exact_sequence = sequence in clean_sequences
        manifest.append({
            "sequence_sha256": identifier, "sequence_length": len(sequence),
            "clean_supported": identifier in supported_index,
            "feature_index": supported_index.get(identifier, -1),
            "accessions": row_accessions,
            "clean_exact_accession_overlap": exact_accession,
            "clean_exact_sequence_overlap": exact_sequence,
            "clean_exact_exposed": exact_accession or exact_sequence,
            "esm1b_embedding_reused_from_verified_prior_run": identifier in old_hash_to_index,
        })
    with (OUTPUT / "sequence_manifest.jsonl").open("w") as handle:
        for row in manifest:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    inputs = [core_fasta, ROOT / "dataset_02/core_edges.json", ROOT / "dataset_02/input_manifest.json",
              ROOT / "CLEAN_CANDIDATE_PRUNING_V1.md", ROOT / "clean_ec_bridge_01/audit.json",
              ROOT / "clean_ec_bridge_validation_01/INDEPENDENT_QC.json", split_path, model_source,
              PROJECT / "tools/external/CLEAN_v1_0_0/src/CLEAN/evaluate.py",
              PROJECT / "tools/external/CLEAN_v1_0_0/NON-EXCLUSIVE RESEARCH USE LICENSE FOR CLEAN SOFTWARE.pdf",
              old_dir / "query_manifest.tsv", old_dir / "unique_query_sequences.fasta",
              old_dir / "query_esm1b_embeddings.npy", PROJECT / "reports/PRECUTOFF_CLEAN_BASELINE_V1_INDEPENDENT_QC.json",
              Path(__file__)]
    for relative in ASSET_HASHES:
        inputs.append(PROJECT / relative)
    write_json(OUTPUT / "input_manifest.json", {str(path.resolve()): digest_file(path) for path in inputs})
    audit = {
        "status": "COMPLETE_INDEPENDENT_QC_PENDING", "created_utc": now(),
        "core_sequences": len(records), "supported_sequences": len(supported),
        "unsupported_over_1022": len(records) - len(supported),
        "esm1b_rows_reused_total": sum(row["esm1b_embedding_reused_from_verified_prior_run"] for row in manifest if row["clean_supported"]),
        "esm1b_rows_not_from_verified_prior_run": sum(not row["esm1b_embedding_reused_from_verified_prior_run"] for row in manifest if row["clean_supported"]),
        "esm1b_cached_rows_at_start_of_this_attempt": cached_rows_at_start,
        "esm1b_rows_newly_computed_this_run": computed,
        "clean_exact_accession_overlap": sum(row["clean_exact_accession_overlap"] for row in manifest),
        "clean_exact_sequence_overlap": sum(row["clean_exact_sequence_overlap"] for row in manifest),
        "clean_exact_exposed": sum(row["clean_exact_exposed"] for row in manifest),
        "clean_exact_unexposed": sum(not row["clean_exact_exposed"] for row in manifest),
        "official_split_rows": len(split_rows), "official_training_assignments": sum(counts.values()),
        "official_ec_labels": len(order), "maximum_separation_queries": len(calls),
        "esm1b_shape": list(esm_matrix.shape), "clean_shape": list(clean_matrix.shape),
        "center_shape": list(centers.shape), "distance_shape": list(units.shape),
        "distance_quantization_scale": DISTANCE_SCALE,
        "current_snapshot_not_historical": True, "independent_biological_validation": False,
        "clean_assets_publicly_redistributable": False,
    }
    write_json(final_marker, audit)
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
