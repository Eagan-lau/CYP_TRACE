"""Independently reconstruct every taxonomy mapping and feature cell."""
from collections import defaultdict
import io
import json
from pathlib import Path
import tarfile
import numpy as np
from run_raw import digest_file, write_json, now


RANKS = ("domain", "phylum", "class", "order", "family", "genus")


def fields(line):
    cells = line.rstrip("\n").split("|")
    return [cell.strip() for cell in cells[:-1]]


def load_independent(path):
    parent = {}; rank = {}; scientific = {}; redirect = {}; deleted = set()
    with tarfile.open(path, "r:gz") as archive:
        members = {Path(m.name).name: m for m in archive.getmembers() if m.isfile()}
        streams = {}
        for name in ["nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"]:
            streams[name] = io.TextIOWrapper(archive.extractfile(members[name]), encoding="utf-8")
        for line in streams["nodes.dmp"]:
            row = fields(line); parent[row[0]] = row[1]; rank[row[0]] = row[2]
        for line in streams["names.dmp"]:
            row = fields(line)
            if row[3] == "scientific name": scientific[row[0]] = row[1]
        for line in streams["merged.dmp"]:
            row = fields(line); redirect[row[0]] = row[1]
        for line in streams["delnodes.dmp"]: deleted.add(fields(line)[0])
    return parent, rank, scientific, redirect, deleted


def resolved(value, parent, redirect, deleted):
    value = str(value); seen = set()
    while value not in parent and value in redirect:
        if value in seen: raise ValueError("redirect cycle")
        seen.add(value); value = redirect[value]
    if value in deleted or value not in parent: raise ValueError("unresolved " + value)
    return value


def path_to_root(value, parent):
    out = []; seen = set()
    while value not in seen:
        seen.add(value); out.append(value)
        if parent[value] == value: return out
        value = parent[value]
    raise ValueError("parent cycle")


def lca(values, parent):
    forward = [list(reversed(path_to_root(v, parent))) for v in values]
    answer = forward[0][0]
    for column in zip(*forward):
        if len(set(column)) > 1: break
        answer = column[0]
    return answer


def run(root):
    result = root / "taxonomy_lineages_01"; out = root / "taxonomy_lineages_validation_01"
    if out.exists(): raise FileExistsError(out)
    manifest = json.loads((result / "input_manifest.json").read_text()); failures = []
    for name, expected in manifest.items():
        if digest_file(root / name) != expected: failures.append("source:" + name)
    taxdump = root.parent / "data/raw/ncbi_taxonomy/new_taxdump.tar.gz"
    parent, rank, scientific, redirect, deleted = load_independent(taxdump)
    edges = json.loads((root / "dataset_02/core_edges.json").read_text()); source = defaultdict(set)
    for edge in edges: source[edge["sequence_sha256"]].update(str(t) for t in edge["taxids"] if str(t))
    records = json.loads((result / "sequence_lineages.json").read_text())
    with np.load(result / "lineage_features.npz", allow_pickle=False) as saved:
        seqs = list(saved["sequence_ids"]); rank_names = list(saved["rank_names"])
        feature_names = list(saved["feature_names"]); actual = saved["features"]
    if seqs != sorted(source) or rank_names != list(RANKS): failures.append("identifier_or_rank_order")
    if [r["sequence_sha256"] for r in records] != seqs: failures.append("record_order")
    expected_tokens = []; lca_count = 0; redirected = 0
    for sequence, record in zip(seqs, records):
        original = sorted(source[sequence], key=lambda x: (int(x) if x.isdigit() else 10**30, x))
        clean = sorted({resolved(t, parent, redirect, deleted) for t in original}, key=int)
        common = lca(clean, parent); lineage = path_to_root(common, parent); by_rank = {}
        for node in lineage:
            if rank[node] in RANKS and rank[node] not in by_rank: by_rank[rank[node]] = node
        tokens = [f"{r}={by_rank[r]}" if r in by_rank else f"{r}=MISSING" for r in RANKS]
        expected_tokens.append(tokens); lca_count += len(clean) > 1; redirected += set(original) != set(clean)
        expected_record = {"sequence_sha256": sequence, "original_taxids": original, "resolved_taxids": clean,
                           "lca_taxid": common, "lca_name": scientific.get(common, ""),
                           "rank_taxids": {r: by_rank.get(r) for r in RANKS},
                           "rank_names": {r: scientific.get(by_rank[r], "") if r in by_rank else None for r in RANKS},
                           "feature_tokens": tokens}
        if record != expected_record: failures.append("record:" + sequence)
    vocabulary = sorted({token for row in expected_tokens for token in row}); matrix = np.zeros((len(seqs), len(vocabulary)), np.uint8)
    index = {token: i for i, token in enumerate(vocabulary)}
    for i, tokens in enumerate(expected_tokens):
        for token in tokens: matrix[i, index[token]] = 1
    if feature_names != vocabulary or not np.array_equal(actual, matrix): failures.append("feature_matrix")
    audit = json.loads((result / "audit.json").read_text())
    checks = {"sequences": len(seqs), "multi_taxid_sequences": sum(len(v) > 1 for v in source.values()),
              "lca_used_sequences": lca_count, "sequences_with_redirected_taxids": redirected,
              "feature_columns": len(vocabulary), "active_features_per_sequence": len(RANKS)}
    for key, value in checks.items():
        if audit.get(key) != value: failures.append("audit:" + key)
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
              "sequence_records_reconstructed": len(records), "feature_cells_reconstructed": int(matrix.size),
              "active_cells_reconstructed": int(matrix.sum()), "taxids_resolved": len({t for v in source.values() for t in v}),
              "lca_sequences_reconstructed": lca_count, "verifier_sha256": digest_file(Path(__file__)),
              "independent_biological_validation": False}
    out.mkdir(); write_json(out / "INDEPENDENT_QC.json", report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": run(Path(__file__).resolve().parent)
