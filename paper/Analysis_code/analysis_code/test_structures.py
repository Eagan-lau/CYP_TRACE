"""Geometry smoke test on one original structure; never adds synthetic labels."""
import gzip
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
from Bio.PDB import MMCIFParser
from rebuild_structure_assets import build


class GeometryTest(unittest.TestCase):
    def test_original_1og5_against_brute_force(self):
        source = Path(__file__).resolve().parents[1] / 'data/raw/rcsb_cyp_ligand_templates_v1/mmcif/1og5.cif.gz'
        if not source.exists(): self.skipTest('Original 1OG5 fixture unavailable')
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / 'raw').mkdir()
            shutil.copyfile(source, folder / 'raw' / source.name)
            build(folder / 'raw', folder / 'output')
            audit = json.loads((folder / 'output/structure_audit.json').read_text())
            self.assertEqual(audit['errors'], [])
            rows = [json.loads(line) for line in (folder / 'output/heme_contact_residues.jsonl').read_text().splitlines()]
            self.assertGreater(len(rows), 0)
            with gzip.open(source, 'rt') as handle:
                model = next(MMCIFParser(QUIET=True).get_structure('check', handle).get_models())
            for row in rows:
                chain = model[row['protein_chain']]
                matches = [r for r in chain if r.id[1] == row['auth_residue_number'] and r.id[2].strip() == row['insertion_code'] and r.resname == row['residue_name']]
                self.assertEqual(len(matches), 1)
                residue = matches[0]
                heme = model[row['heme_chain']][tuple(row['heme_residue_id'])]
                # Independent pairwise arithmetic, not the KD-tree routine.
                distances = [float(np.linalg.norm(a.coord.astype(float) - b.coord.astype(float)))
                             for a in residue for b in heme if a.element not in {'H', 'D'} and b.element not in {'H', 'D'}]
                expected = min(distances)
                self.assertAlmostEqual(row['minimum_heme_heavy_atom_distance_A'], expected, places=9)
                self.assertLessEqual(expected, 5.0)
                self.assertEqual(row['label'], 'GEOMETRY_ONLY_NOT_CATALYTIC_TRUTH')


if __name__ == '__main__': unittest.main()
