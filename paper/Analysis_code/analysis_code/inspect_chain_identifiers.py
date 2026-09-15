"""Read raw author/label chain crosswalks for the two mismatched mappings."""
from collections import Counter
import gzip
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import json

root = Path(__file__).resolve().parent
for file in ['2c6h.cif.gz', '2cd8.cif.gz']:
    with gzip.open(root.parent / 'data/raw/rcsb_cyp_ligand_templates_v1/mmcif' / file, 'rt') as h: d = MMCIF2Dict(h)
    pairs = Counter((a, b) for a, b, seq in zip(d['_atom_site.auth_asym_id'], d['_atom_site.label_asym_id'], d['_atom_site.label_seq_id']) if seq not in ['.', '?'])
    print(json.dumps({'file': file, 'author_label_chain_polymer_atom_counts': [[a, b, n] for (a, b), n in pairs.items()]}))
