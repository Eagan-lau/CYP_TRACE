"""Independent checks for the raw-source reaction/EC bridge."""
from __future__ import annotations

from collections import defaultdict
import csv
import json
from pathlib import Path
import re
import sqlite3

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "clean_ec_bridge_01"
OUTPUT = ROOT / "clean_ec_bridge_validation_01"
PATTERN = re.compile(r"(?<![0-9])([1-9]\d*\.\d+\.\d+\.\d+)(?![0-9])")


def extract(source, raw):
    payload = json.loads(raw)
    if source == "UniProt_API": text = payload.get("reaction", {}).get("ecNumber", "")
    elif source == "SwissProt": text = payload.get("comment", "")
    elif source == "P450Rdb": text = payload.get("EC number", "")
    else: text = ""
    return set(PATTERN.findall(str(text)))


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    failures = []
    for name, expected in json.loads((SOURCE / "input_manifest.json").read_text()).items():
        if digest_file(Path(name)) != expected: failures.append("input_hash:" + name)
    edges = json.loads((ROOT / "dataset_02/core_edges.json").read_text())
    declared = json.loads((SOURCE / "reaction_ec_support.json").read_text())
    clean_split = ROOT.parent / "tools/external/CLEAN_v1_0_0/data/split100.csv"
    clean_catalog = set()
    with clean_split.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            clean_catalog.update(PATTERN.findall(row.get("EC number", "")))
    database = ROOT / "run_03/raw_rebuild.sqlite"
    connection = sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    needed = {value for edge in edges for value in edge["assertion_ids"]}
    rows = {row["assertion_id"]: dict(row) for row in connection.execute(
        "SELECT assertion_id,source,sequence_sha256,raw_json FROM assertions") if row["assertion_id"] in needed}
    expected = defaultdict(lambda: defaultdict(lambda: {"sequences": set(), "sources": set(), "assertions": set()}))
    for edge in edges:
        for assertion in edge["assertion_ids"]:
            row = rows[assertion]
            for ec in extract(row["source"], row["raw_json"]):
                value = expected[edge["reaction_key"]][ec]
                value["sequences"].add(edge["sequence_sha256"]); value["sources"].add(row["source"]); value["assertions"].add(assertion)
    checked = 0
    for reaction, mapping in declared.items():
        if set(mapping) != set(expected.get(reaction, {})): failures.append("reaction_ec_set:" + reaction)
        for ec, value in mapping.items():
            source = expected[reaction][ec]
            if value["supporting_sequences"] != sorted(source["sequences"]): failures.append("sequence_support")
            if value["supporting_sources"] != sorted(source["sources"]): failures.append("source_support")
            if value["supporting_assertions"] != sorted(source["assertions"]): failures.append("assertion_support")
            if value["in_clean_catalog"] != (ec in clean_catalog): failures.append("clean_catalog_membership")
            checked += 1
    edge_rows = [json.loads(line) for line in (SOURCE / "edge_ec_records.jsonl").read_text().splitlines()]
    if len(edge_rows) != len(edges): failures.append("edge_count")
    for edge, saved in zip(edges, edge_rows):
        observed = sorted({ec for assertion in edge["assertion_ids"] for ec in extract(rows[assertion]["source"], rows[assertion]["raw_json"])})
        observed_clean = sorted(set(observed) & clean_catalog)
        if saved["sequence_sha256"] != edge["sequence_sha256"] or saved["reaction_key"] != edge["reaction_key"] or saved["direct_exact_ecs"] != observed or saved["clean_catalog_exact_ecs"] != observed_clean:
            failures.append("edge_reconstruction")
    audit = json.loads((SOURCE / "audit.json").read_text())
    mapped_reactions = sum(any(value["in_clean_catalog"] for value in mapping.values()) for mapping in declared.values())
    if audit["clean_exact_ec_catalog"] != len(clean_catalog): failures.append("clean_catalog_count")
    if audit["reactions_with_direct_clean_exact_ec"] != mapped_reactions: failures.append("mapped_reaction_count")
    if audit["edges_with_direct_clean_exact_ec"] != sum(bool(row["clean_catalog_exact_ecs"]) for row in edge_rows): failures.append("mapped_edge_count")
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
              "reaction_ec_mappings_reconstructed": checked, "edges_reconstructed": len(edges),
              "raw_assertions_reparsed": len(rows), "clean_catalog_ecs_reparsed": len(clean_catalog),
              "independent_implementation": True,
              "independent_biological_validation": False}
    OUTPUT.mkdir(); write_json(OUTPUT / "INDEPENDENT_QC.json", report); connection.close()
    print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": main()
