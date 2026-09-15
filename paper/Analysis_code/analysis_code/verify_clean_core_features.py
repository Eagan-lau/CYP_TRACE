"""Independently reconstruct CLEAN projection, centres, distances, calls and exposure."""
from __future__ import annotations

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
SOURCE = ROOT / "clean_core_features_01"
OUTPUT = ROOT / "clean_core_features_validation_01"


def read_fasta(path):
    records = []
    key = None
    sequence = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if key is not None: records.append((key, sequence))
            key, sequence = line[1:].split()[0], ""
        else: sequence += line.strip()
    if key is not None: records.append((key, sequence))
    return records


def ec_tokens(value):
    return [part.strip().removeprefix("EC:").strip() for part in re.split(r"[;|/]", value) if part.strip()]


def layer_norm_net(path):
    spec = importlib.util.spec_from_file_location("independent_clean_model", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module.LayerNormNet


def independent_call_count(values):
    baseline = np.concatenate((values[1:], np.full(10, values[-1]))).mean()
    gradients = np.abs(np.diff(np.abs(values - baseline)))
    eligible = [index for index, value in enumerate(gradients) if value > gradients.mean()]
    separation_index = eligible[0] if eligible else 0
    return 1 if separation_index >= 5 else separation_index + 1


def main():
    import torch

    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    failures = []
    for name, expected in json.loads((SOURCE / "input_manifest.json").read_text()).items():
        if digest_file(Path(name)) != expected: failures.append("input_hash:" + name)
    records = sorted(read_fasta(ROOT / "dataset_02/core_sequences.fasta"))
    for identifier, sequence in records:
        if hashlib.sha256(sequence.encode()).hexdigest() != identifier: failures.append("fasta_sha256")
    manifest = [json.loads(line) for line in (SOURCE / "sequence_manifest.jsonl").read_text().splitlines()]
    if [row["sequence_sha256"] for row in manifest] != [identifier for identifier, _ in records]: failures.append("manifest_order")
    supported = [(identifier, sequence) for identifier, sequence in records if len(sequence) <= 1022]
    feature_index = {identifier: index for index, (identifier, _) in enumerate(supported)}
    for row in manifest:
        if row["clean_supported"] != (row["sequence_length"] <= 1022): failures.append("support_domain")
        if row["feature_index"] != feature_index.get(row["sequence_sha256"], -1): failures.append("feature_index")

    esm = np.load(SOURCE / "esm1b_embeddings.npy", allow_pickle=False)
    clean = np.load(SOURCE / "clean_embeddings.npy", allow_pickle=False)
    centers = np.load(SOURCE / "ec_centers.npy", allow_pickle=False)
    units = np.load(SOURCE / "ec_distance_units.npy", allow_pickle=False)
    if esm.shape != (len(supported), 1280) or clean.shape != (len(supported), 128) or centers.shape != (5242, 128) or units.shape != (len(supported), 5242): failures.append("matrix_shape")
    if not all(np.isfinite(value).all() for value in (esm, clean, centers)) or np.min(units) < 0: failures.append("matrix_values")

    model_path = PROJECT / "tools/external/CLEAN_v1_0_0/src/CLEAN/model.py"
    Model = layer_norm_net(model_path)
    model = Model(512, 128, torch.device("cpu"), torch.float32)
    model.load_state_dict(torch.load(PROJECT / "data/restricted/clean_v1_0_0/pretrained/split100.pth", map_location="cpu"))
    model.eval()
    with torch.no_grad(): recalculated_clean = model(torch.from_numpy(esm)).numpy()
    clean_max_difference = float(np.max(np.abs(recalculated_clean - clean)))
    if clean_max_difference > 1e-6: failures.append("clean_projection")

    split = list(csv.DictReader((PROJECT / "tools/external/CLEAN_v1_0_0/data/split100.csv").open(newline=""), delimiter="\t"))
    order = []
    counts = {}
    entries = set()
    sequences = set()
    for row in split:
        entries.add(row["Entry"]); sequences.add(row["Sequence"])
        for ec in ec_tokens(row["EC number"]):
            if ec not in counts: order.append(ec); counts[ec] = 0
            counts[ec] += 1
    catalog = json.loads((SOURCE / "ec_catalog.json").read_text())
    if [row["ec"] for row in catalog] != order or [row["training_assignment_count"] for row in catalog] != [counts[ec] for ec in order]: failures.append("ec_catalog")
    training = torch.as_tensor(torch.load(PROJECT / "data/restricted/clean_v1_0_0/pretrained/100.pt", map_location="cpu"), dtype=torch.float32)
    rebuilt = []
    offset = 0
    for ec in order:
        count = counts[ec]; rebuilt.append(training[offset:offset + count].mean(0).detach().numpy()); offset += count
    rebuilt_centers = np.stack(rebuilt).astype(np.float32)
    center_max_difference = float(np.max(np.abs(rebuilt_centers - centers)))
    if center_max_difference > 1e-7: failures.append("ec_centers")
    rebuilt_units = np.empty_like(units)
    with torch.no_grad():
        for start in range(0, len(clean), 16):
            stop = min(start + 16, len(clean))
            value = torch.cdist(torch.from_numpy(clean[start:stop]), torch.from_numpy(centers)).double() * 1_000_000
            rebuilt_units[start:stop] = torch.round(value).clamp(0, np.iinfo(np.int32).max).numpy().astype(np.int32)
    distance_max_difference = int(np.max(np.abs(rebuilt_units.astype(np.int64) - units.astype(np.int64))))
    if distance_max_difference != 0: failures.append("ec_distances")

    with torch.no_grad():
        whole_matrix_units = torch.round(torch.cdist(torch.from_numpy(clean), torch.from_numpy(centers)).double() * 1_000_000).clamp(0, np.iinfo(np.int32).max).numpy().astype(np.int32)
    whole_matrix_max_difference = int(np.max(np.abs(whole_matrix_units.astype(np.int64) - units.astype(np.int64))))

    calls = json.loads((SOURCE / "maximum_separation_calls.json").read_text())
    whole_matrix_call_label_disagreements = 0
    for identifier, index in feature_index.items():
        nearest = np.argsort(units[index], kind="stable")[:10]
        count = independent_call_count(units[index, nearest].astype(np.float64))
        expected = [(order[int(position)], int(units[index, position])) for position in nearest[:count]]
        observed = [(row["ec"], row["distance_units"]) for row in calls[identifier]]
        if observed != expected: failures.append("maximum_separation_call")
        alternative_nearest = np.argsort(whole_matrix_units[index], kind="stable")[:10]
        alternative_count = independent_call_count(whole_matrix_units[index, alternative_nearest].astype(np.float64))
        alternative_labels = [order[int(position)] for position in alternative_nearest[:alternative_count]]
        if alternative_labels != [row["ec"] for row in calls[identifier]]:
            whole_matrix_call_label_disagreements += 1

    edge_accessions = {}
    for edge in json.loads((ROOT / "dataset_02/core_edges.json").read_text()):
        edge_accessions.setdefault(edge["sequence_sha256"], set()).update(edge["accessions"])
    for row, (identifier, sequence) in zip(manifest, records):
        accession_overlap = any(value in entries for value in edge_accessions[identifier])
        sequence_overlap = sequence in sequences
        if row["clean_exact_accession_overlap"] != accession_overlap or row["clean_exact_sequence_overlap"] != sequence_overlap or row["clean_exact_exposed"] != (accession_overlap or sequence_overlap): failures.append("exposure")

    audit = json.loads((SOURCE / "audit.json").read_text())
    if audit["supported_sequences"] != len(supported) or audit["maximum_separation_queries"] != len(calls): failures.append("audit_counts")
    report = {
        "status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
        "core_sequences_reconstructed": len(records), "supported_sequences": len(supported),
        "ec_labels_reconstructed": len(order), "training_assignments_reconstructed": sum(counts.values()),
        "clean_projection_max_abs_difference": clean_max_difference,
        "center_max_abs_difference": center_max_difference,
        "distance_unit_max_abs_difference": distance_max_difference,
        "whole_matrix_kernel_distance_unit_max_abs_difference": whole_matrix_max_difference,
        "whole_matrix_kernel_maximum_separation_label_disagreements": whole_matrix_call_label_disagreements,
        "maximum_separation_calls_reconstructed": len(calls),
        "independent_implementation": True, "independent_biological_validation": False,
    }
    OUTPUT.mkdir(); write_json(OUTPUT / "INDEPENDENT_QC.json", report)
    print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": main()
