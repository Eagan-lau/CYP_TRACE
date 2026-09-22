"""Reconstruct human labels from bundled original public supplements."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from prepare_workspace import prepare

EXPECTED = {
    'development': ('human_substrate_01/normalized_labels.json', 14955,
                    'eab51d5b6099720d083ba9dbcd8a778d9d133161d40968dbfeb9557df1acc82f'),
    'external': ('external_human_cyp_01/normalized_external_labels.json', 14526,
                 '951d6f74ea96c48ffc803c0f76d27c89ed21f36a4c5a5960623d96ef1c955281'),
}

def main(output):
    prepare(output)
    analysis = output.resolve() / 'analysis'
    for args in [
        ['build_human_substrate.py', '--source', 'raw_additions/molecules-26-04678-s001.zip', '--output', 'human_substrate_01'],
        ['parse_external_human_cyp.py'],
    ]:
        subprocess.run([sys.executable, '-X', 'utf8', *args], cwd=analysis, check=True)
    results = {}
    for name, (relative, count, expected_sha) in EXPECTED.items():
        records = json.loads((analysis / relative).read_text(encoding='utf-8'))
        digest = hashlib.sha256(json.dumps(records, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        results[name] = {'records': len(records), 'semantic_sha256': digest,
                         'matches_original': len(records) == count and digest == expected_sha}
    status = 'PASS' if all(r['matches_original'] for r in results.values()) else 'FAIL'
    report = {'status': status, 'scope': 'Original human source normalization; no model selection', 'outputs': results}
    (output / 'HUMAN_SOURCE_ACCEPTANCE.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    return 0 if status == 'PASS' else 1

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    raise SystemExit(main(p.parse_args().output))
