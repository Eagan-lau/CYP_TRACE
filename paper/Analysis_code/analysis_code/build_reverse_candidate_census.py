"""Build reaction-blind exact-taxid CYP candidate panels for reverse retrieval."""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import sqlite3

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "reverse_candidate_census_01"


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    database = ROOT / "run_03/raw_rebuild.sqlite"
    dataset_manifest = json.loads((ROOT / "dataset_02/input_manifest.json").read_text())
    expected = next(value for name, value in dataset_manifest.items() if name.replace("\\", "/").endswith("run_03/raw_rebuild.sqlite"))
    if digest_file(database) != expected: raise ValueError("raw database hash")
    edges = json.loads((ROOT / "dataset_02/core_edges.json").read_text())
    blocks = json.loads((ROOT / "multiaxis_split_01/outer_inner_blocks.json").read_text())
    core = {edge["sequence_sha256"] for edge in edges}
    core_taxids = defaultdict(set)
    for edge in edges:
        core_taxids[edge["sequence_sha256"]].update(str(value) for value in edge["taxids"] if str(value))

    connection = sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    proteins = {row["sequence_sha256"]: dict(row) for row in connection.execute(
        "SELECT sequence_sha256,sequence,length,model_sequence_valid FROM proteins")}
    assertions = defaultdict(list)
    for row in connection.execute("SELECT source,locator,accession,sequence_sha256,taxid,family_evidence FROM protein_assertions"):
        assertions[row["sequence_sha256"]].append(dict(row))

    admitted = []; quarantined = []
    for sequence_sha256, protein in sorted(proteins.items()):
        records = assertions.get(sequence_sha256, [])
        exact_raw = any(row["family_evidence"] == "PF00067" for row in records)
        if sequence_sha256 not in core and not exact_raw: continue
        reasons = []
        if not protein["model_sequence_valid"]: reasons.append("invalid_sequence")
        taxids = {str(row["taxid"]) for row in records if str(row["taxid"] or "")}
        if sequence_sha256 in core: taxids.update(core_taxids[sequence_sha256])
        if len(taxids) == 0: reasons.append("taxid_missing")
        if len(taxids) > 1: reasons.append("taxid_conflict")
        record = {
            "sequence_sha256": sequence_sha256, "sequence": protein["sequence"], "sequence_length": protein["length"],
            "taxids": sorted(taxids), "taxid": next(iter(taxids)) if len(taxids) == 1 else None,
            "core_sequence": sequence_sha256 in core, "raw_exact_PF00067": exact_raw,
            "accessions": sorted({row["accession"] for row in records if row["accession"]}),
            "sources": sorted({row["source"] for row in records}),
            "family_evidence": sorted({row["family_evidence"] for row in records}),
        }
        if reasons: quarantined.append({**record, "reasons": reasons})
        else: admitted.append(record)

    panels = defaultdict(list)
    for row in admitted: panels[row["taxid"]].append(row["sequence_sha256"])
    panel_payload = {taxid: sorted(values) for taxid, values in sorted(panels.items())}
    sequence_taxid = {row["sequence_sha256"]: row["taxid"] for row in admitted}
    edge_taxid = defaultdict(set)
    for edge_number, edge in enumerate(edges):
        candidate_taxid = sequence_taxid.get(edge["sequence_sha256"])
        if candidate_taxid and candidate_taxid in {str(value) for value in edge["taxids"]}:
            edge_taxid[edge_number].add(candidate_taxid)

    query_panels = []
    excluded = defaultdict(int)
    for block_number, block in enumerate(blocks):
        if not block["evaluable"]: continue
        grouped = defaultdict(set)
        for edge_number in block["test_edge_indices"]:
            edge = edges[edge_number]
            if not edge_taxid[edge_number]:
                excluded["test_positive_without_unambiguous_candidate_taxid"] += 1
                continue
            for taxid in edge_taxid[edge_number]:
                grouped[(taxid, edge["reaction_key"])].add(edge["sequence_sha256"])
        for (taxid, reaction), positives in sorted(grouped.items()):
            candidates = panel_payload[taxid]
            if len(candidates) < 2:
                excluded["single_candidate_taxid_panel"] += 1
                continue
            if not positives.issubset(candidates): raise ValueError("positive outside candidate panel")
            query_panels.append({
                "block_number": block_number, "block_id": block["block_id"], "source_task": block["task"],
                "taxid": taxid, "reaction_key": reaction, "positive_sequence_sha256": sorted(positives),
                "candidate_count": len(candidates), "core_candidate_count": sum(value in core for value in candidates),
                "positive_count": len(positives),
            })

    OUTPUT.mkdir()
    with (OUTPUT / "candidate_records.jsonl").open("w") as handle:
        for row in admitted:
            public = {key: value for key, value in row.items() if key != "sequence"}
            handle.write(json.dumps(public, sort_keys=True) + "\n")
    with (OUTPUT / "quarantined_records.jsonl").open("w") as handle:
        for row in quarantined:
            public = {key: value for key, value in row.items() if key != "sequence"}
            handle.write(json.dumps(public, sort_keys=True) + "\n")
    with (OUTPUT / "candidate_sequences.fasta").open("w") as handle:
        for row in admitted: handle.write(f">{row['sequence_sha256']}\n{row['sequence']}\n")
    with (OUTPUT / "query_panels.jsonl").open("w") as handle:
        for row in query_panels: handle.write(json.dumps(row, sort_keys=True) + "\n")
    write_json(OUTPUT / "taxid_candidate_panels.json", panel_payload)
    panel_sizes = [len(value) for value in panel_payload.values()]
    audit = {
        "status": "COMPLETE_INDEPENDENT_QC_PENDING", "created_utc": now(),
        "raw_proteins": len(proteins), "core_sequences": len(core),
        "admitted_candidate_sequences": len(admitted), "quarantined_sequences": len(quarantined),
        "admitted_core_sequences": sum(row["core_sequence"] for row in admitted),
        "admitted_noncore_exact_pf00067_sequences": sum(not row["core_sequence"] for row in admitted),
        "candidate_taxids": len(panel_payload), "taxids_with_at_least_two_candidates": sum(value >= 2 for value in panel_sizes),
        "maximum_candidates_in_taxid": max(panel_sizes), "outer_query_panels": len(query_panels),
        "outer_query_taxids": len({row["taxid"] for row in query_panels}),
        "outer_query_reactions": len({row["reaction_key"] for row in query_panels}),
        "outer_documented_positive_instances": sum(row["positive_count"] for row in query_panels),
        "candidate_count_min": min(row["candidate_count"] for row in query_panels),
        "candidate_count_max": max(row["candidate_count"] for row in query_panels),
        "candidate_count_mean": sum(row["candidate_count"] for row in query_panels) / len(query_panels),
        "exclusions": dict(sorted(excluded.items())),
        "candidate_source_is_exhaustive_proteome": False, "performance_calculated": False,
        "independent_biological_validation": False,
    }
    write_json(OUTPUT / "audit.json", audit)
    inputs = [database, ROOT / "dataset_02/core_edges.json", ROOT / "dataset_02/input_manifest.json",
              ROOT / "multiaxis_split_01/outer_inner_blocks.json", ROOT / "multiaxis_split_01/output_checksums.json",
              ROOT / "REVERSE_RETRIEVAL_V1.md", Path(__file__)]
    write_json(OUTPUT / "input_manifest.json", {str(path.resolve()): digest_file(path) for path in inputs})
    connection.close(); print(json.dumps(audit, indent=2))


if __name__ == "__main__": main()
