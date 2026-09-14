"""Exact evidence route and explicit abstention for unsupported generalization."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

from rdkit import Chem


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    sequence: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(sequence).upper()
            name = line[1:].split()[0]
            if not name or name in records:
                raise ValueError("FASTA identifiers must be non-empty and unique")
            sequence = []
        else:
            if name is None:
                raise ValueError("FASTA sequence encountered before an identifier")
            sequence.append(line)
    if name is not None:
        records[name] = "".join(sequence).upper()
    if not records:
        raise ValueError("FASTA contains no records")
    allowed = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*")
    for identifier, value in records.items():
        if not value or not set(value) <= allowed:
            raise ValueError(f"invalid protein sequence for {identifier}")
    return records


def build_evidence_index(dataset_dir: Path, output: Path) -> dict:
    files = {
        "sequences": dataset_dir / "core_sequences.fasta",
        "reactions": dataset_dir / "core_reactions.json",
        "edges": dataset_dir / "core_edges.json",
    }
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise ValueError("missing dataset files: " + ", ".join(missing))
    if output.exists():
        raise ValueError(f"output already exists: {output}")
    fasta = read_fasta(files["sequences"])
    sequences = {hashlib.sha256(sequence.encode("ascii")).hexdigest(): sequence for sequence in fasta.values()}
    if set(fasta) != set(sequences):
        raise ValueError("core FASTA headers do not equal sequence SHA256 values")
    reactions = json.loads(files["reactions"].read_text(encoding="utf-8"))
    edges = json.loads(files["edges"].read_text(encoding="utf-8"))
    edge_map: dict[str, list[dict]] = {}
    for edge in edges:
        key = edge["sequence_sha256"] + "|" + edge["reaction_key"]
        edge_map.setdefault(key, []).append({
            field: edge[field]
            for field in ("publication_ids", "accessions", "taxids", "sources", "source_lineages", "assertion_ids")
        })
    bundle = {
        "schema_version": "cyptrace-evidence-index-v1",
        "input_sha256": {name: _sha256(path) for name, path in files.items()},
        "counts": {"sequences": len(sequences), "reactions": len(reactions), "edges": len(edges)},
        "sequences": sequences,
        "reactions": reactions,
        "edges": edge_map,
        "claim_boundary": "Exact matches retrieve source evidence. All non-exact sequence/reaction cases abstain because no general reaction route passed the frozen development gates.",
        "public_redistribution_ready": False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(bundle, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    with output.open("wb") as raw_stream:
        with gzip.GzipFile(filename="cyptrace_evidence_v1.json", mode="wb", fileobj=raw_stream, compresslevel=9, mtime=0) as stream:
            stream.write(payload)
    return {"status": "PASS", "output": str(output.resolve()), "output_sha256": _sha256(output), **bundle["counts"]}


def load_evidence_index(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        data = json.load(stream)
    if data.get("schema_version") != "cyptrace-evidence-index-v1":
        raise ValueError("unsupported evidence index schema")
    return data


def _canonical_smiles(value: str) -> str | None:
    molecule = Chem.MolFromSmiles(value)
    return Chem.MolToSmiles(molecule, isomericSmiles=True) if molecule is not None else None


def _candidate_reaction(candidate: dict, reactions: dict) -> tuple[str | None, str]:
    supplied = candidate.get("reaction_key", "").strip()
    if supplied:
        return (supplied if supplied in reactions else None), "reaction_key"
    substrate = _canonical_smiles(candidate.get("substrate_smiles", ""))
    product = _canonical_smiles(candidate.get("product_smiles", ""))
    if substrate is None or product is None:
        return None, "invalid_smiles"
    hits = [
        key for key, row in reactions.items()
        if row.get("single_pair") and row.get("substrates") == [substrate] and row.get("products") == [product]
    ]
    return (hits[0] if len(hits) == 1 else None), "structure_pair"


def screen(index: dict, queries: dict[str, str], candidates: list[dict]) -> list[dict]:
    results = []
    for query_id, sequence in queries.items():
        sequence_sha256 = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        sequence_seen = sequence_sha256 in index["sequences"]
        for candidate in candidates:
            candidate_id = candidate["candidate_id"]
            reaction_key, resolution = _candidate_reaction(candidate, index["reactions"])
            chemistry_seen = reaction_key is not None
            evidence = index["edges"].get(sequence_sha256 + "|" + reaction_key, []) if chemistry_seen else []
            if evidence:
                route = "DOCUMENTED_EXACT_EVIDENCE"
                abstained = False
                interpretation = "Source evidence was retrieved for the exact sequence and exact normalized reaction; this is not a novel prediction."
            else:
                route = "ABSTAIN_UNVALIDATED_GENERAL_REACTION_ROUTE"
                abstained = True
                interpretation = "No validated generalization route is available. Missing evidence is UNKNOWN, not a negative label or a zero expert score."
            results.append({
                "query_id": query_id,
                "sequence_sha256": sequence_sha256,
                "candidate_id": candidate_id,
                "reaction_key": reaction_key,
                "candidate_resolution": resolution,
                "route": route,
                "abstained": abstained,
                "domain": {
                    "exact_reference_sequence": sequence_seen,
                    "exact_reaction_observed": chemistry_seen,
                    "exact_sequence_reaction_edge": bool(evidence),
                    "homology_expert_contribution": "not_scored",
                    "structure_expert_contribution": "not_scored",
                    "clean_expert_contribution": "context_only_not_scored",
                },
                "evidence": evidence,
                "score": None,
                "interpretation": interpretation,
            })
    return results
