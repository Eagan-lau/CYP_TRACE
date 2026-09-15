"""Read-only portable verification before launching the P02 cluster search."""
import hashlib
import json
from pathlib import Path, PurePosixPath

root = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''): h.update(chunk)
    return h.hexdigest()


def local_project_path(recorded):
    text = recorded.replace('\\', '/')
    marker = 'paper_rebuild_20260909/'
    if marker not in text: raise ValueError('Input outside the rebuilt project: ' + text)
    suffix = PurePosixPath(text.split(marker, 1)[1])
    if suffix.is_absolute() or '..' in suffix.parts: raise ValueError('Unsafe manifest path')
    return root.joinpath(*suffix.parts)


checks, failures = 0, []
for folder in ['dataset_02', 'named_support_01']:
    p = root / folder
    for name, expected in json.loads((p / 'output_checksums.json').read_text()).items():
        checks += 1
        if sha(p / name) != expected: failures.append(folder + '/' + name)
    for recorded, expected in json.loads((p / 'input_manifest.json').read_text()).items():
        target = local_project_path(recorded)
        checks += 1
        if sha(target) != expected: failures.append(str(target.relative_to(root)))
audit = json.loads((root / 'dataset_02/dataset_audit.json').read_text())
if audit['status'] != 'PASS' or not all(audit['checks'].values()): failures.append('dataset_audit_failed')
receipt = json.loads((root / 'p02_validation_01/P02_RECEIPT.json').read_text())
if receipt['status'] != 'P02_VIEW_CHECKS_PASS': failures.append('P02_receipt_failed')
print(json.dumps({'status': 'PASS' if not failures else 'FAIL', 'file_hash_checks': checks,
    'failures': failures, 'labelled_sequences': audit['labelled_sequences'],
    'unique_edges': audit['unique_edges'], 'model_evaluation_performed': False}, indent=2))
if failures: raise SystemExit(1)
