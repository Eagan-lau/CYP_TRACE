"""Build a fresh, provenance-preserving reaction/sequence EC bridge for CLEAN."""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import re
import sqlite3

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
OUTPUT = ROOT / "clean_ec_bridge_01"
EXACT_EC = re.compile(r"(?<![0-9])([1-9]\d*\.\d+\.\d+\.\d+)(?![0-9])")


def ecs(text) -> set[str]:
    return set(EXACT_EC.findall(str(text or "")))


def assertion_ec(source: str, payload: dict) -> set[str]:
    if source == "UniProt_API":
        return ecs(payload.get("reaction", {}).get("ecNumber", ""))
    if source == "SwissProt":
        return ecs(payload.get("comment", ""))
    if source == "P450Rdb":
        return ecs(payload.get("EC number", ""))
    return set()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    database = ROOT / "run_03/raw_rebuild.sqlite"
    dataset_manifest = json.loads((ROOT / "dataset_02/input_manifest.json").read_text())
    expected_database = next(value for name, value in dataset_manifest.items() if name.replace("\\", "/").endswith("run_03/raw_rebuild.sqlite"))
    if digest_file(database) != expected_database:
        raise ValueError("fresh raw database hash")
    edges = json.loads((ROOT / "dataset_02/core_edges.json").read_text())
    reactions = sorted(json.loads((ROOT / "dataset_02/core_reactions.json").read_text()))
    sequences = sorted({edge["sequence_sha256"] for edge in edges})
    assertion_ids = {value for edge in edges for value in edge["assertion_ids"]}
    connection = sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    assertion_rows = {}
    for row in connection.execute("SELECT assertion_id,source,sequence_sha256,raw_json FROM assertions"):
        if row["assertion_id"] in assertion_ids:
            assertion_rows[row["assertion_id"]] = dict(row)
    if set(assertion_rows) != assertion_ids:
        raise ValueError("core assertion lookup")

    clean_split = PROJECT / "tools/external/CLEAN_v1_0_0/data/split100.csv"
    clean_ec = set()
    import csv
    with clean_split.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            clean_ec.update(ecs(row.get("EC number", "")))
    reaction_support = defaultdict(lambda: defaultdict(list))
    sequence_direct = defaultdict(set)
    edge_records = []
    source_assertions = defaultdict(int)
    for edge in edges:
        direct = set(); supports = []
        for assertion_id in edge["assertion_ids"]:
            row = assertion_rows[assertion_id]
            found = assertion_ec(row["source"], json.loads(row["raw_json"]))
            source_assertions[row["source"]] += 1
            for ec in sorted(found):
                direct.add(ec); sequence_direct[edge["sequence_sha256"]].add(ec)
                support = {"assertion_id": assertion_id, "source": row["source"],
                           "sequence_sha256": edge["sequence_sha256"]}
                reaction_support[edge["reaction_key"]][ec].append(support)
                supports.append({"ec": ec, **support})
        edge_records.append({
            "sequence_sha256": edge["sequence_sha256"], "reaction_key": edge["reaction_key"],
            "direct_exact_ecs": sorted(direct), "clean_catalog_exact_ecs": sorted(direct & clean_ec),
            "supports": supports,
        })

    accession_sequences = defaultdict(set)
    for edge in edges:
        for accession in edge["accessions"]:
            accession_sequences[accession].add(edge["sequence_sha256"])
    brenda_sequence = defaultdict(set); brenda_links = 0
    brenda_path = ROOT / "brenda_01/brenda_assertions.jsonl"
    for line in brenda_path.read_text().splitlines():
        row = json.loads(line); ec = row["ec"]
        if ec not in clean_ec:
            continue
        for accession in row["linked_cyp_accessions"]:
            linked = accession_sequences.get(accession, set())
            if len(linked) == 1:
                brenda_sequence[next(iter(linked))].add(ec); brenda_links += 1
    sequence_context = {
        sequence: {
            "direct_reaction_exact_ecs": sorted(sequence_direct.get(sequence, set())),
            "brenda_accession_context_exact_ecs": sorted(brenda_sequence.get(sequence, set())),
            "combined_clean_catalog_exact_ecs": sorted((sequence_direct.get(sequence, set()) | brenda_sequence.get(sequence, set())) & clean_ec),
        }
        for sequence in sequences
    }
    reaction_payload = {}
    for reaction in reactions:
        mapping = {}
        for ec, supports in sorted(reaction_support.get(reaction, {}).items()):
            unique = {(item["assertion_id"], item["source"], item["sequence_sha256"]): item for item in supports}
            mapping[ec] = {
                "in_clean_catalog": ec in clean_ec,
                "supporting_sequences": sorted({item["sequence_sha256"] for item in unique.values()}),
                "supporting_sources": sorted({item["source"] for item in unique.values()}),
                "supporting_assertions": sorted({item["assertion_id"] for item in unique.values()}),
            }
        reaction_payload[reaction] = mapping

    OUTPUT.mkdir()
    with (OUTPUT / "edge_ec_records.jsonl").open("w") as handle:
        for row in edge_records: handle.write(json.dumps(row, sort_keys=True) + "\n")
    write_json(OUTPUT / "reaction_ec_support.json", reaction_payload)
    write_json(OUTPUT / "sequence_ec_context.json", sequence_context)
    mapped_reactions = [reaction for reaction, values in reaction_payload.items() if any(item["in_clean_catalog"] for item in values.values())]
    leave_one_out_pairs = 0; total_pairs = 0
    for edge in edges:
        total_pairs += 1
        if any(item["in_clean_catalog"] and any(sequence != edge["sequence_sha256"] for sequence in item["supporting_sequences"])
               for item in reaction_payload[edge["reaction_key"]].values()):
            leave_one_out_pairs += 1
    audit = {
        "status": "COMPLETE_INDEPENDENT_QC_PENDING", "created_utc": now(),
        "core_sequences": len(sequences), "core_reactions": len(reactions), "core_edges": len(edges),
        "clean_exact_ec_catalog": len(clean_ec),
        "reactions_with_direct_clean_exact_ec": len(mapped_reactions),
        "reaction_mapping_fraction": len(mapped_reactions) / len(reactions),
        "edges_with_direct_clean_exact_ec": sum(bool(row["clean_catalog_exact_ecs"]) for row in edge_records),
        "edges_with_leave_query_sequence_out_candidate_ec": leave_one_out_pairs,
        "leave_query_sequence_out_edge_fraction": leave_one_out_pairs / total_pairs,
        "sequences_with_direct_clean_exact_ec": sum(bool(value["direct_reaction_exact_ecs"]) for value in sequence_context.values()),
        "sequences_with_brenda_clean_exact_ec_context": sum(bool(value["brenda_accession_context_exact_ecs"]) for value in sequence_context.values()),
        "brenda_accession_assertion_links": brenda_links,
        "source_assertions_scanned": dict(sorted(source_assertions.items())),
        "old_integrated_label_table_used": False,
        "brenda_used_as_precise_reaction_mapping": False,
        "independent_biological_validation": False,
    }
    write_json(OUTPUT / "audit.json", audit)
    inputs = [database, ROOT / "dataset_02/core_edges.json", ROOT / "dataset_02/core_reactions.json",
              ROOT / "dataset_02/input_manifest.json", brenda_path, ROOT / "brenda_01/brenda_audit.json",
              clean_split, ROOT / "CLEAN_CANDIDATE_PRUNING_V1.md", Path(__file__)]
    write_json(OUTPUT / "input_manifest.json", {str(path.resolve()): digest_file(path) for path in inputs})
    connection.close(); print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
