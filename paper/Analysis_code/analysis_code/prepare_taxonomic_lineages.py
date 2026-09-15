"""Build fixed six-rank NCBI-taxonomy features for the 600 core sequences."""
from collections import defaultdict
import io
import json
from pathlib import Path
import tarfile
import numpy as np
from run_raw import write_json, digest_file, now


RANKS = ("domain", "phylum", "class", "order", "family", "genus")


def dmp_fields(line):
    return [part.strip() for part in line.rstrip("\n").split("|")[:-1]]


def read_taxdump(path):
    parents = {}; ranks = {}; names = {}; merged = {}; deleted = set()
    with tarfile.open(path, "r:gz") as archive:
        members = {Path(m.name).name: m for m in archive.getmembers() if m.isfile()}
        required = {"nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"}
        if not required.issubset(members):
            raise ValueError("Taxdump lacks required members: " + ",".join(sorted(required - set(members))))
        for line in io.TextIOWrapper(archive.extractfile(members["nodes.dmp"]), encoding="utf-8"):
            fields = dmp_fields(line); parents[fields[0]] = fields[1]; ranks[fields[0]] = fields[2]
        for line in io.TextIOWrapper(archive.extractfile(members["names.dmp"]), encoding="utf-8"):
            fields = dmp_fields(line)
            if fields[3] == "scientific name": names[fields[0]] = fields[1]
        for line in io.TextIOWrapper(archive.extractfile(members["merged.dmp"]), encoding="utf-8"):
            fields = dmp_fields(line); merged[fields[0]] = fields[1]
        for line in io.TextIOWrapper(archive.extractfile(members["delnodes.dmp"]), encoding="utf-8"):
            fields = dmp_fields(line); deleted.add(fields[0])
    return parents, ranks, names, merged, deleted


def resolve_taxid(taxid, parents, merged, deleted):
    current = str(taxid); visited = set()
    while current not in parents and current in merged:
        if current in visited: raise ValueError("Taxid merge cycle: " + current)
        visited.add(current); current = merged[current]
    if current in deleted or current not in parents: raise ValueError("Unresolved taxid: " + str(taxid))
    return current


def ancestry(taxid, parents):
    result = []; current = taxid; visited = set()
    while True:
        if current in visited: raise ValueError("Taxonomy parent cycle: " + current)
        visited.add(current); result.append(current)
        parent = parents[current]
        if parent == current: break
        if parent not in parents: raise ValueError("Missing taxonomy parent: " + parent)
        current = parent
    return result


def lowest_common_ancestor(taxids, parents):
    paths = [list(reversed(ancestry(t, parents))) for t in taxids]
    common = paths[0][0]
    for level in zip(*paths):
        if len(set(level)) != 1: break
        common = level[0]
    return common


def build(root, taxdump):
    out = root / "taxonomy_lineages_01"
    if out.exists(): raise FileExistsError(out)
    edges = json.loads((root / "dataset_02/core_edges.json").read_text())
    seq_taxids = defaultdict(set)
    for edge in edges:
        seq_taxids[edge["sequence_sha256"]].update(str(t) for t in edge["taxids"] if str(t))
    seqs = sorted(seq_taxids)
    if len(seqs) != 600 or any(not seq_taxids[q] for q in seqs): raise ValueError("Expected 600 taxid-bearing sequences")
    parents, ranks, names, merged, deleted = read_taxdump(taxdump)
    records = []; token_rows = []
    for sequence in seqs:
        original = sorted(seq_taxids[sequence], key=lambda x: (int(x) if x.isdigit() else 10**30, x))
        resolved = sorted({resolve_taxid(t, parents, merged, deleted) for t in original}, key=lambda x: int(x))
        lca = lowest_common_ancestor(resolved, parents); path = ancestry(lca, parents)
        by_rank = {}
        for node in path:
            rank = ranks[node]
            if rank in RANKS and rank not in by_rank: by_rank[rank] = node
        tokens = [f"{rank}={by_rank[rank]}" if rank in by_rank else f"{rank}=MISSING" for rank in RANKS]
        token_rows.append(tokens)
        records.append({"sequence_sha256": sequence, "original_taxids": original, "resolved_taxids": resolved,
                        "lca_taxid": lca, "lca_name": names.get(lca, ""),
                        "rank_taxids": {rank: by_rank.get(rank) for rank in RANKS},
                        "rank_names": {rank: names.get(by_rank[rank], "") if rank in by_rank else None for rank in RANKS},
                        "feature_tokens": tokens})
    vocabulary = sorted({token for row in token_rows for token in row}); vi = {token: i for i, token in enumerate(vocabulary)}
    features = np.zeros((len(seqs), len(vocabulary)), dtype=np.uint8)
    for i, tokens in enumerate(token_rows):
        for token in tokens: features[i, vi[token]] = 1
    if not np.all(features.sum(1) == len(RANKS)): raise ValueError("Taxonomy feature cardinality")
    out.mkdir()
    write_json(out / "sequence_lineages.json", records)
    np.savez_compressed(out / "lineage_features.npz", sequence_ids=np.array(seqs), rank_names=np.array(RANKS),
                        feature_names=np.array(vocabulary), features=features)
    missing = {rank: sum(r["rank_taxids"][rank] is None for r in records) for rank in RANKS}
    audit = {"status": "COMPLETE_INDEPENDENT_QC_PENDING", "created_utc": now(), "sequences": len(seqs),
             "input_unique_taxids": len({t for values in seq_taxids.values() for t in values}),
             "multi_taxid_sequences": sum(len(v) > 1 for v in seq_taxids.values()),
             "sequences_with_redirected_taxids": sum(set(r["original_taxids"]) != set(r["resolved_taxids"]) for r in records),
             "feature_ranks": list(RANKS), "feature_columns": len(vocabulary), "active_features_per_sequence": len(RANKS),
             "missing_rank_counts": missing, "lca_used_sequences": sum(len(r["resolved_taxids"]) > 1 for r in records),
             "outcome_fields_used": [], "transductive_label_free_vocabulary": True,
             "independent_biological_validation": False}
    write_json(out / "audit.json", audit)
    reltax = taxdump.relative_to(root).as_posix() if taxdump.is_relative_to(root) else Path("../data/raw/ncbi_taxonomy/new_taxdump.tar.gz").as_posix()
    names_in = ["prepare_taxonomic_lineages.py", "TAXONOMIC_LINEAGE_CONTROL_V1.md", "dataset_02/core_edges.json"]
    manifest = {name: digest_file(root / name) for name in names_in}; manifest[reltax] = digest_file(taxdump)
    write_json(out / "input_manifest.json", manifest)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    build(root, root.parent / "data/raw/ncbi_taxonomy/new_taxdump.tar.gz")
