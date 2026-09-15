"""Refuse changed raw coordinate inputs before the structural search."""
import json
from pathlib import Path
from run_raw import digest_file,write_json,now

root=Path(__file__).resolve().parent
raw=root.parent/'data/raw/rcsb_cyp_ligand_templates_v1/mmcif'
expected=json.loads((root/'structure_assets_01/source_manifest.json').read_text())
if {p.name for p in raw.glob('*.cif.gz')}!={r['file'] for r in expected}: raise ValueError('Raw structure file census changed')
for row in expected:
    if digest_file(raw/row['file'])!=row['sha256']: raise ValueError('Changed coordinate file: '+row['file'])
write_json(root/'foldseek_raw_01/input_receipt.json',{'status':'PASS','created_utc':now(),'files':len(expected),
    'source_manifest_sha256':digest_file(root/'structure_assets_01/source_manifest.json'),
    'raw_directory':str(raw.resolve()),'files_checked':expected})
print(f'Raw coordinate hashes verified: {len(expected)}',flush=True)
