"""Command-line interface for the evidence-bounded CYP-TRACE strategy."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

from rdkit import rdBase

from . import __version__
from .evidence import build_evidence_index, load_evidence_index, read_fasta, screen
from .human import HumanSubstrateModel


def _write(path: Path | None, result: dict) -> None:
    text = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    if path is None:
        sys.stdout.write(text)
        return
    if path.exists():
        raise ValueError(f"output already exists: {path}")
    if not path.parent.is_dir():
        raise ValueError(f"output parent directory does not exist: {path.parent}")
    path.write_text(text, encoding="utf-8")
    print(json.dumps({"status": "COMPLETED", "output": str(path.resolve())}))


def _read_human_input(path: Path) -> list[dict]:
    if path.suffix.lower() == ".jsonl":
        rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    else:
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    for number, row in enumerate(rows, 1):
        if not row.get("smiles"):
            raise ValueError(f"human input row {number} has no smiles")
    return rows


def _read_candidates(path: Path) -> list[dict]:
    if path.suffix.lower() == ".jsonl":
        rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    else:
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
    for number, row in enumerate(rows, 1):
        row.setdefault("candidate_id", f"candidate_{number}")
        if not row["candidate_id"]:
            row["candidate_id"] = f"candidate_{number}"
        if not row.get("reaction_key") and not (row.get("substrate_smiles") and row.get("product_smiles")):
            raise ValueError(f"candidate row {number} needs reaction_key or substrate_smiles plus product_smiles")
    return rows


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description="CYP-TRACE: domain-aware evidence and fixed-human substrate strategy")
    root.add_argument("--version", action="version", version=__version__)
    commands = root.add_subparsers(dest="command", required=True)
    doctor = commands.add_parser("doctor", help="Report model and runtime scope")
    doctor.add_argument("--model", type=Path)
    human = commands.add_parser("human-substrate", help="Score compounds for fixed human CYP isoforms")
    source = human.add_mutually_exclusive_group(required=True)
    source.add_argument("--smiles")
    source.add_argument("--input", type=Path, help="TSV or JSONL with smiles and optional id/isoform")
    human.add_argument("--isoform", action="append", help="Repeat for multiple isoforms; omitted means all nine")
    human.add_argument("--model", type=Path)
    human.add_argument("--output", type=Path)
    build = commands.add_parser("build-evidence-index", help="Build a portable local exact-evidence index")
    build.add_argument("--dataset-dir", required=True, type=Path)
    build.add_argument("--output", required=True, type=Path)
    reaction = commands.add_parser("reaction-screen", help="Route FASTA queries and supplied reaction candidates")
    reaction.add_argument("--fasta", required=True, type=Path)
    reaction.add_argument("--candidates", required=True, type=Path)
    reaction.add_argument("--evidence-index", required=True, type=Path)
    reaction.add_argument("--output", type=Path)
    return root


def execute(args: argparse.Namespace) -> dict:
    if args.command == "doctor":
        model = HumanSubstrateModel.load(args.model)
        return {
            "status": "PASS",
            "tool_version": __version__,
            "rdkit_version": rdBase.rdkitVersion,
            "human_model": {
                "compounds": len(model.compounds),
                "isoforms": list(model.isoforms),
                "externally_validated_isoforms": sorted(model.validated_isoforms),
                "isoform_k": model.isoform_k,
                "pooled_k": model.pooled_k,
            },
            "general_reaction_route": "exact_evidence_or_abstain",
            "missing_expert_policy": "domain_weight_zero_not_score_zero",
        }
    if args.command == "build-evidence-index":
        return build_evidence_index(args.dataset_dir, args.output)
    if args.command == "human-substrate":
        model = HumanSubstrateModel.load(args.model)
        if args.smiles is not None:
            rows = [{"id": "query_1", "smiles": args.smiles}]
        else:
            rows = _read_human_input(args.input)
        results = []
        for number, row in enumerate(rows, 1):
            requested = args.isoform or ([row["isoform"]] if row.get("isoform") else list(model.isoforms))
            for isoform in requested:
                scored = model.score(row["smiles"], isoform)
                scored["query_id"] = row.get("id") or f"query_{number}"
                results.append(scored)
        result = {
            "schema_version": "cyptrace-output-v1",
            "task": "fixed_human_substrate_ranking",
            "model_content_sha256": model.bundle["bundle_content_sha256"],
            "result_count": len(results),
            "results": results,
        }
        _write(args.output, result)
        return result
    index = load_evidence_index(args.evidence_index)
    results = screen(index, read_fasta(args.fasta), _read_candidates(args.candidates))
    result = {
        "schema_version": "cyptrace-output-v1",
        "task": "general_reaction_evidence_routing",
        "evidence_index_sha256": hashlib.sha256(args.evidence_index.read_bytes()).hexdigest(),
        "result_count": len(results),
        "results": results,
    }
    _write(args.output, result)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = execute(args)
        if args.command in {"doctor", "build-evidence-index"}:
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
        return 0
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(json.dumps({"status": "ERROR", "error_type": type(exc).__name__, "error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2
