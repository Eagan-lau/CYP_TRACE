"""Independent raw-store and panel checks for reverse-retrieval candidates."""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sqlite3

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "reverse_candidate_census_01"
OUTPUT = ROOT / "reverse_candidate_census_validation_01"


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    failures = []
    for name, expected in json.loads((SOURCE / "input_manifest.json").read_text()).items():
        if digest_file(Path(name)) != expected: failures.append("input_hash:" + name)
    saved = [json.loads(line) for line in (SOURCE / "candidate_records.jsonl").read_text().splitlines()]
    saved_index = {row["sequence_sha256"]: row for row in saved}
    fasta = {}
    key = None
    for line in (SOURCE / "candidate_sequences.fasta").read_text().splitlines():
        if line.startswith(">"):
            key = line[1:]; fasta[key] = ""
        else: fasta[key] += line.strip()
    if set(fasta) != set(saved_index): failures.append("fasta_candidate_set")
    for key, sequence in fasta.items():
        if hashlib.sha256(sequence.encode()).hexdigest() != key: failures.append("fasta_hash")

    edges = json.loads((ROOT / "dataset_02/core_edges.json").read_text())
    core = {edge["sequence_sha256"] for edge in edges}
    core_taxids = defaultdict(set)
    for edge in edges: core_taxids[edge["sequence_sha256"]].update(str(value) for value in edge["taxids"] if str(value))
    connection = sqlite3.connect((ROOT / "run_03/raw_rebuild.sqlite").resolve().as_uri() + "?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    proteins = {row["sequence_sha256"]: dict(row) for row in connection.execute("SELECT sequence_sha256,sequence,model_sequence_valid FROM proteins")}
    evidence = defaultdict(list)
    for row in connection.execute("SELECT sequence_sha256,taxid,family_evidence FROM protein_assertions"):
        evidence[row["sequence_sha256"]].append(dict(row))
    expected = set()
    for sequence, protein in proteins.items():
        exact = any(row["family_evidence"] == "PF00067" for row in evidence[sequence])
        if sequence not in core and not exact: continue
        taxids = {str(row["taxid"]) for row in evidence[sequence] if str(row["taxid"] or "")}
        if sequence in core: taxids.update(core_taxids[sequence])
        if protein["model_sequence_valid"] and len(taxids) == 1: expected.add(sequence)
    if set(saved_index) != expected: failures.append("raw_candidate_reconstruction")
    for sequence, row in saved_index.items():
        if row["sequence_length"] != len(proteins[sequence]["sequence"]): failures.append("length")

    taxid_panels = json.loads((SOURCE / "taxid_candidate_panels.json").read_text())
    reconstructed = defaultdict(list)
    for sequence, row in saved_index.items(): reconstructed[row["taxid"]].append(sequence)
    if taxid_panels != {key: sorted(value) for key, value in sorted(reconstructed.items())}: failures.append("taxid_panels")
    panels = [json.loads(line) for line in (SOURCE / "query_panels.jsonl").read_text().splitlines()]
    for panel in panels:
        candidates = taxid_panels[panel["taxid"]]
        if panel["candidate_count"] != len(candidates) or len(candidates) < 2: failures.append("panel_denominator")
        if not set(panel["positive_sequence_sha256"]).issubset(candidates): failures.append("panel_positive")
    if len({(row["block_number"], row["taxid"], row["reaction_key"]) for row in panels}) != len(panels): failures.append("duplicate_query_panel")
    blocks = json.loads((ROOT / "multiaxis_split_01/outer_inner_blocks.json").read_text())
    rebuilt_panels = {}
    for block_number, block in enumerate(blocks):
        if not block["evaluable"]: continue
        grouped = defaultdict(set)
        for edge_number in block["test_edge_indices"]:
            edge = edges[edge_number]
            candidate = saved_index.get(edge["sequence_sha256"])
            if candidate and candidate["taxid"] in {str(value) for value in edge["taxids"]}:
                grouped[(candidate["taxid"], edge["reaction_key"])].add(edge["sequence_sha256"])
        for (taxid, reaction), positives in grouped.items():
            if len(taxid_panels[taxid]) >= 2:
                rebuilt_panels[(block_number, taxid, reaction)] = sorted(positives)
    observed_panels = {(row["block_number"], row["taxid"], row["reaction_key"]): row["positive_sequence_sha256"] for row in panels}
    if observed_panels != rebuilt_panels: failures.append("query_panel_reconstruction")

    audit = json.loads((SOURCE / "audit.json").read_text())
    if audit["admitted_candidate_sequences"] != len(saved) or audit["outer_query_panels"] != len(panels): failures.append("audit_count")
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
              "raw_candidate_sequences_reconstructed": len(expected), "taxid_panels_reconstructed": len(taxid_panels),
              "outer_query_panels_checked": len(panels), "candidate_membership_reaction_blind": True,
              "independent_implementation": True, "independent_biological_validation": False}
    OUTPUT.mkdir(); write_json(OUTPUT / "INDEPENDENT_QC.json", report); connection.close()
    print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": main()
