"""Normalize MMseqs2 and BLASTp candidate-to-core alignments for reverse retrieval."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "reverse_alignment_01"


def identifiers(path): return {line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")}


def keep_best(store, row):
    key = row[:2]
    if key not in store or row[-1] > store[key][-1]: store[key] = row


def parse_mmseqs(path, queries, targets):
    output = {}
    with path.open() as handle:
        for raw in csv.reader(handle, delimiter="\t"):
            if len(raw) != 8: raise ValueError("MMseqs2 field count")
            query, target = raw[:2]
            values = list(map(float, raw[2:]))
            identity, qcov, tcov, length, evalue, bits = values
            if query not in queries or target not in targets or not all(map(math.isfinite, values)): raise ValueError("MMseqs2 row")
            if not (0 <= identity <= 1 and 0 <= qcov <= 1 and 0 <= tcov <= 1 and length > 0 and evalue >= 0): raise ValueError("MMseqs2 range")
            keep_best(output, (query, target, identity, qcov, tcov, int(length), evalue, bits))
    return output


def parse_blast(path, queries, targets):
    output = {}
    with path.open() as handle:
        for raw in csv.reader(handle, delimiter="\t"):
            if len(raw) != 12: raise ValueError("BLASTp field count")
            query, target = raw[:2]
            pident, length, qstart, qend, sstart, send, evalue, bits, qlen, slen = map(float, raw[2:])
            values = (pident, length, qstart, qend, sstart, send, evalue, bits, qlen, slen)
            if query not in queries or target not in targets or not all(map(math.isfinite, values)): raise ValueError("BLASTp row")
            identity = pident / 100
            qcov = (abs(qend - qstart) + 1) / qlen
            tcov = (abs(send - sstart) + 1) / slen
            if not (0 <= identity <= 1 and 0 <= qcov <= 1 and 0 <= tcov <= 1 and length > 0 and evalue >= 0): raise ValueError("BLASTp range")
            keep_best(output, (query, target, identity, qcov, tcov, int(length), evalue, bits))
    return output


def write_table(path, mapping):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for key in sorted(mapping): writer.writerow(mapping[key])


def main():
    if (OUTPUT / "audit.json").exists(): raise FileExistsError(OUTPUT / "audit.json")
    candidates = identifiers(ROOT / "reverse_candidate_census_01/candidate_sequences.fasta")
    core = identifiers(ROOT / "dataset_02/core_sequences.fasta")
    mmseqs = parse_mmseqs(OUTPUT / "mmseqs_raw.tsv", candidates, core)
    blast = parse_blast(OUTPUT / "blast_raw.tsv", candidates, core)
    write_table(OUTPUT / "mmseqs_normalized.tsv", mmseqs)
    write_table(OUTPUT / "blast_normalized.tsv", blast)
    core_candidates = candidates & core
    mmseqs_self = sum((value, value) in mmseqs for value in core_candidates)
    blast_self = sum((value, value) in blast for value in core_candidates)
    if mmseqs_self != len(core_candidates) or blast_self != len(core_candidates): raise ValueError("missing exact self hit")
    audit = {
        "status": "COMPLETE_INDEPENDENT_QC_PENDING", "created_utc": now(),
        "candidate_sequences": len(candidates), "core_targets": len(core), "core_candidate_overlap": len(core_candidates),
        "mmseqs_unique_pairs": len(mmseqs), "blast_unique_pairs": len(blast),
        "mmseqs_queries_with_hit": len({key[0] for key in mmseqs}), "blast_queries_with_hit": len({key[0] for key in blast}),
        "mmseqs_exact_self_hits": mmseqs_self, "blast_exact_self_hits": blast_self,
        "performance_calculated": False, "independent_biological_validation": False,
    }
    write_json(OUTPUT / "audit.json", audit)
    inputs = [ROOT / "reverse_candidate_census_01/candidate_sequences.fasta",
              ROOT / "reverse_candidate_census_validation_01/INDEPENDENT_QC.json",
              ROOT / "dataset_02/core_sequences.fasta", ROOT / "REVERSE_RETRIEVAL_V1.md",
              OUTPUT / "mmseqs_raw.tsv", OUTPUT / "blast_raw.tsv", OUTPUT / "mmseqs_version.txt",
              OUTPUT / "blast_version.txt", Path(__file__)]
    write_json(OUTPUT / "input_manifest.json", {str(path.resolve()): digest_file(path) for path in inputs})
    print(json.dumps(audit, indent=2))


if __name__ == "__main__": main()
