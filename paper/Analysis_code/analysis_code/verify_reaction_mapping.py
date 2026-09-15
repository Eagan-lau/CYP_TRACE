"""Independent verification for reaction-mapping outputs."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from run_raw import digest_file, now


def run(root: Path, out_name: str = "reaction_mapping_01") -> None:
    out = root / out_name
    mapping = out / "mapping_results.jsonl"
    manifest = out / "input_manifest.json"
    if not mapping.exists():
        raise FileNotFoundError(mapping)

    if manifest.exists():
        for name, expected in json.loads(manifest.read_text()).items():
            if name not in ["dataset_02/core_reactions.json", "run_reaction_mapping.py", "reaction_mapping_core.py"]:
                continue
            if name.endswith(".jsonl"):
                continue
            if digest_file(root / name) != expected:
                raise ValueError(f"Input hash mismatch: {name}")

    lines = [json.loads(line) for line in mapping.read_text().splitlines() if line.strip()]
    chemistry = json.loads((root / "dataset_02" / "core_reactions.json").read_text())
    expected_ids = sorted(chemistry)

    failures = []
    if len(lines) != len(expected_ids):
        failures.append(f"line_count:{len(lines)}")
    if [r.get("reaction_key") for r in lines] != expected_ids:
        failures.append("reaction_order")

    status = Counter(r.get("status") for r in lines)
    mapped_ok = sum(1 for r in lines if r.get("status") in {"mapped", "mapped_validation_warning"})
    for row in lines:
        required = {"reaction_key", "status", "reaction_smiles", "mapped_reaction", "mapped_confidence", "validation", "centers"}
        missing = [k for k in required if k not in row]
        if missing:
            failures.append(f"missing:{row.get('reaction_key')}:{','.join(sorted(missing))}")
            continue
        if row["status"] in {"mapped", "mapped_validation_warning"}:
            if not isinstance(row["mapped_reaction"], str) or not row["mapped_reaction"]:
                failures.append(f"empty_mapped:{row['reaction_key']}")

    report = {
        "status": "PASS" if not failures else "FAIL",
        "created_utc": now(),
        "scope": "order-preserved candidate-wise mapping with per-candidate status/validation/centres",
        "reaction_count": len(expected_ids),
        "mapped_rows": mapped_ok,
        "status_counts": dict(sorted(status.items())),
        "failures": failures,
        "independent_biological_validation": False,
        "atom_level_reaction_center_model": False,
    }
    out_report = out / "INDEPENDENT_QC.json"
    out_report.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--dir", default="reaction_mapping_01")
    args = p.parse_args()
    run(Path(__file__).resolve().parent, args.dir)
