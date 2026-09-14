"""Check release payload hashes without network, dependencies or mutation."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    files = []
    for line in (ROOT / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        expected, name = line.split('  ', 1)
        path = (ROOT / name).resolve()
        assert path.is_relative_to(ROOT), name
        assert path.is_file(), name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, name
        files.append(name)
    assert len(files) == len(set(files))
    metadata = json.loads((ROOT / 'RELEASE_METADATA.json').read_text())
    print(json.dumps({'status': 'PAYLOAD_HASHES_PASS', 'files': len(files),
                      'public_release_approved': metadata['public_release_approved'],
                      'scope': 'File integrity only; not a publication or rights approval'}, indent=2))


if __name__ == '__main__':
    main()
