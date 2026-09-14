"""Build deterministic model and local evidence assets for joint_tool_01."""
from __future__ import annotations

from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "human_substrate_01" / "normalized_labels.json"
OUTPUT = ROOT / "joint_tool_01" / "src" / "cyptrace_pipeline" / "data" / "human_model_v1.json.gz"
MANIFEST = ROOT / "joint_tool_01" / "ASSET_MANIFEST.json"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def content_digest(data: dict) -> str:
    copy = dict(data)
    copy.pop("bundle_content_sha256", None)
    return hashlib.sha256(json.dumps(copy, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def main() -> None:
    if OUTPUT.exists() or MANIFEST.exists():
        raise FileExistsError("refusing to replace an existing joint-tool asset")
    rows = json.loads(SOURCE.read_text(encoding="utf-8"))
    compounds = {}
    labels = defaultdict(dict)
    for row in rows:
        key = row["compound_inchikey"]
        descriptor = (row["canonical_smiles"], row["scaffold_group"])
        if key in compounds and compounds[key] != descriptor:
            raise ValueError(f"inconsistent compound representation: {key}")
        compounds[key] = descriptor
        if row["isoform"] in labels[key]:
            raise ValueError(f"duplicate isoform-compound label: {key}")
        labels[key][row["isoform"]] = int(row["label"])
    compound_rows = [
        {
            "compound_inchikey": key,
            "canonical_smiles": compounds[key][0],
            "scaffold_group": compounds[key][1],
            "labels": dict(sorted(labels[key].items())),
        }
        for key in sorted(compounds)
    ]
    bundle = {
        "schema_version": "cyptrace-human-model-v1",
        "source": {
            "name": "CYPstrate 2021 supplementary data",
            "normalized_labels_sha256": digest(SOURCE),
            "normalized_pairs": len(rows),
            "license": "CC BY 4.0 with attribution",
        },
        "algorithm": {
            "name": "similarity_weighted_knn",
            "morgan_radius": 2,
            "morgan_bits": 2048,
            "isoform_k": 25,
            "pooled_k": 51,
            "tie_break": "sha256_compound_inchikey",
        },
        "scope": {
            "development_isoforms": sorted({row["isoform"] for row in rows}),
            "externally_validated_isoforms": ["CYP1A2", "CYP2C19", "CYP2C9", "CYP2D6", "CYP2E1", "CYP3A4"],
            "externally_validated_endpoint": "new-chemistry fixed-human substrate ranking",
            "unsupported": ["unseen proteins", "products", "reaction centres", "catalytic rates", "clinical decisions"],
        },
        "compounds": compound_rows,
    }
    bundle["bundle_content_sha256"] = content_digest(bundle)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(bundle, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    with OUTPUT.open("wb") as raw_stream:
        with gzip.GzipFile(filename="human_model_v1.json", mode="wb", fileobj=raw_stream, compresslevel=9, mtime=0) as stream:
            stream.write(payload)
    manifest = {
        "status": "PASS",
        "source_sha256": digest(SOURCE),
        "model_bundle_sha256": digest(OUTPUT),
        "bundle_content_sha256": bundle["bundle_content_sha256"],
        "normalized_pairs": len(rows),
        "unique_compounds": len(compound_rows),
        "isoforms": bundle["scope"]["development_isoforms"],
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
