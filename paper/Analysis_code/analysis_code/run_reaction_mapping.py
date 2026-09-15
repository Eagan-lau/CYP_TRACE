"""Run RXNMapper over all frozen core MAIN reactions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from collections import Counter

from rxnmapper import RXNMapper

from reaction_mapping_core import map_candidate
from run_raw import digest_file, now, write_json


def run(root: Path, batch_size: int) -> None:
    out = root / "reaction_mapping_01"
    if out.exists():
        raise FileExistsError(out)
    out.mkdir()

    chemistry = json.loads((root / "dataset_02" / "core_reactions.json").read_text())
    reaction_ids = sorted(chemistry)

    mapper = RXNMapper()
    rxnmapper_version = getattr(__import__("rxnmapper"), "__version__", "unknown")

    rows = []
    status_count = Counter()
    for i, rid in enumerate(reaction_ids, 1):
        result = map_candidate(rid, chemistry[rid], mapper, rxnmapper_version)
        status_count[result["status"]] += 1
        rows.append(result)
        if i % max(1, batch_size) == 0:
            print(f"mapped {i}/{len(reaction_ids)} status={result['status']}", flush=True)

    summary = {
        "created_utc": now(),
        "created_by": "run_reaction_mapping.py",
        "reaction_count": len(reaction_ids),
        "mapper": {
            "name": "RXNMapper",
            "version": rxnmapper_version,
        },
        "status_counts": dict(sorted(status_count.items())),
        "mapped_count": sum(1 for row in rows if row["status"] in {"mapped", "mapped_validation_warning"}),
        "mapped_validation_warning_count": sum(1 for row in rows if row["status"] == "mapped_validation_warning"),
        "mapper_error_count": sum(1 for row in rows if row["status"] == "mapper_error"),
        "validation_error_count": sum(1 for row in rows if row["status"] == "mapping_validation_error"),
    }

    with (out / "mapping_results.jsonl").open("w", encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps(row, ensure_ascii=False) + "\n")
    write_json(out / "summary.json", summary)
    write_json(out / "input_manifest.json", {
        "dataset_02/core_reactions.json": digest_file(root / "dataset_02" / "core_reactions.json"),
        "run_reaction_mapping.py": digest_file(root / "run_reaction_mapping.py"),
        "reaction_mapping_core.py": digest_file(root / "reaction_mapping_core.py"),
    })

    print("mapping complete", json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=100)
    args = parser.parse_args()
    run(Path(__file__).resolve().parent, args.batch_size)
