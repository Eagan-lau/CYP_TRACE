"""Read-only bounded schema inspection of one raw EC entry."""
import argparse
import json
import re
import tarfile

p = argparse.ArgumentParser()
p.add_argument('archive')
p.add_argument('--ec', default='1.14.14.1')
p.add_argument('--proteins-only', action='store_true')
a = p.parse_args()
with tarfile.open(a.archive, 'r:gz') as archive:
    member = next(m for m in archive if m.name.endswith('.json'))
    lines = []
    for raw in archive.extractfile(member):
        line = raw.decode('utf-8')
        if line.startswith(f'    "{a.ec}": '): lines.append(line.split(': ', 1)[1]); continue
        if lines:
            if re.match(r'^    "[^\"]+": \{', line): break
            lines.append(line)
    value = json.loads(''.join(lines).rstrip().rstrip(','))
summary = {}
if a.proteins_only:
    print(json.dumps(value.get('protein', {}), ensure_ascii=True, indent=2))
    raise SystemExit(0)
for key, item in value.items():
    if isinstance(item, dict): summary[key] = {'size': len(item), 'sample': list(item.items())[:2]}
    elif isinstance(item, list): summary[key] = {'size': len(item), 'sample': item[:2]}
    else: summary[key] = item
print(json.dumps(summary, ensure_ascii=True, indent=2))
