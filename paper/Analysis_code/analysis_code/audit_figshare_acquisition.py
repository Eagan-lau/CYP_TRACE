"""Verify the immutable Figshare v4 acquisition and write its manifest."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw_additions" / "figshare_26630515_v4"
OUTPUT = ROOT / "external_human_cyp_01"
PROTOCOL = ROOT / "EXTERNAL_HUMAN_CYP_DATASET_V1.md"
ARTICLE = RAW / "figshare_article_v4.json"


def digest(path: Path, algorithm: str) -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    metadata = json.loads(ARTICLE.read_text(encoding="utf-8"))
    if metadata.get("id") != 26630515 or metadata.get("version") != 4:
        raise RuntimeError("unexpected Figshare record or version")
    expected = metadata.get("files", [])
    expected_names = [row["name"] for row in expected]
    actual_names = sorted(
        path.name for path in RAW.iterdir()
        if path.is_file() and path.name != ARTICLE.name
    )
    if sorted(expected_names) != actual_names:
        raise RuntimeError("downloaded file set does not match API metadata")

    records = []
    failures = []
    for row in expected:
        path = RAW / row["name"]
        size = path.stat().st_size
        md5 = digest(path, "md5")
        supplied = (row.get("supplied_md5") or row.get("computed_md5") or "").lower()
        checks = {
            "size_matches": size == int(row["size"]),
            "md5_matches": bool(supplied) and md5 == supplied,
        }
        if not all(checks.values()):
            failures.append({"name": row["name"], **checks})
        records.append({
            "file_id": int(row["id"]),
            "name": row["name"],
            "download_url": row["download_url"],
            "expected_bytes": int(row["size"]),
            "observed_bytes": size,
            "supplied_md5": supplied,
            "observed_md5": md5,
            "sha256": digest(path, "sha256"),
            **checks,
        })

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "record_id": metadata["id"],
        "version": metadata["version"],
        "title": metadata["title"],
        "doi": metadata.get("doi"),
        "published_date": metadata.get("published_date"),
        "modified_date": metadata.get("modified_date"),
        "license": metadata.get("license"),
        "protocol_sha256": digest(PROTOCOL, "sha256"),
        "api_metadata_sha256": digest(ARTICLE, "sha256"),
        "file_count": len(records),
        "total_bytes": sum(row["observed_bytes"] for row in records),
        "complete": not failures,
        "failures": failures,
        "files": records,
    }
    OUTPUT.mkdir(exist_ok=True)
    destination = OUTPUT / "acquisition_manifest.json"
    destination.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "manifest": str(destination),
        "complete": manifest["complete"],
        "file_count": manifest["file_count"],
        "total_bytes": manifest["total_bytes"],
        "failures": failures,
    }, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
