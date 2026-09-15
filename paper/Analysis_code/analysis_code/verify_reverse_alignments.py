"""Independently parse and verify candidate-to-core MMseqs2 and BLASTp outputs."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from run_raw import digest_file, now, write_json


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "reverse_alignment_01"
OUTPUT = ROOT / "reverse_alignment_validation_01"


def ids(path): return {line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")}


def independently_parse(path, blast):
    result = {}
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if blast:
                q, t = row[:2]; pident, length, qs, qe, ss, se, evalue, bits, qlen, slen = map(float, row[2:])
                qcov = (abs(qe - qs) + 1) / qlen
                tcov = (abs(se - ss) + 1) / slen
                value = (q, t, pident / 100, qcov, tcov, int(length), evalue, bits)
            else:
                q, t = row[:2]; identity, qcov, tcov, length, evalue, bits = map(float, row[2:])
                value = (q, t, identity, qcov, tcov, int(length), evalue, bits)
            key = (q, t)
            if key not in result or value[-1] > result[key][-1]: result[key] = value
    return result


def normalized(path):
    output = {}
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            value = (row[0], row[1], *map(float, row[2:]))
            output[value[:2]] = value
    return output


def main():
    if OUTPUT.exists(): raise FileExistsError(OUTPUT)
    failures = []
    for name, expected in json.loads((SOURCE / "input_manifest.json").read_text()).items():
        if digest_file(Path(name)) != expected: failures.append("input_hash:" + name)
    candidate_ids = ids(ROOT / "reverse_candidate_census_01/candidate_sequences.fasta")
    core_ids = ids(ROOT / "dataset_02/core_sequences.fasta")
    counts = {}
    for method, blast in (("mmseqs", False), ("blast", True)):
        rebuilt = independently_parse(SOURCE / f"{method}_raw.tsv", blast)
        saved = normalized(SOURCE / f"{method}_normalized.tsv")
        if set(rebuilt) != set(saved): failures.append(method + "_pair_set")
        for key, expected in rebuilt.items():
            observed = saved[key]
            if observed[:2] != expected[:2] or any(abs(float(a) - float(b)) > 1e-12 for a, b in zip(observed[2:], expected[2:])):
                failures.append(method + "_value")
        if any(q not in candidate_ids or t not in core_ids for q, t in rebuilt): failures.append(method + "_identifier")
        self_hits = sum((value, value) in rebuilt for value in candidate_ids & core_ids)
        if self_hits != len(candidate_ids & core_ids): failures.append(method + "_self_hits")
        counts[method] = {"unique_pairs": len(rebuilt), "queries_with_hit": len({key[0] for key in rebuilt}), "exact_self_hits": self_hits}
    audit = json.loads((SOURCE / "audit.json").read_text())
    if audit["mmseqs_unique_pairs"] != counts["mmseqs"]["unique_pairs"] or audit["blast_unique_pairs"] != counts["blast"]["unique_pairs"]: failures.append("audit_count")
    report = {"status": "PASS" if not failures else "FAIL", "created_utc": now(), "failures": failures,
              "candidate_sequences": len(candidate_ids), "core_targets": len(core_ids), "methods": counts,
              "independent_implementation": True, "performance_calculated": False,
              "independent_biological_validation": False}
    OUTPUT.mkdir(); write_json(OUTPUT / "INDEPENDENT_QC.json", report)
    print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == "__main__": main()
