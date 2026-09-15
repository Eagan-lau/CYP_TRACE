import importlib.util
import csv
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('raw_rebuild', Path(__file__).with_name('run_raw.py'))
raw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(raw)


class RawTests(unittest.TestCase):
    def test_atom_maps_removed_stereo_preserved(self):
        a = raw.normalize_structure('[CH3:1][C@@H:2](O)C(=O)O')['canonical_smiles']
        b = raw.normalize_structure('C[C@@H](O)C(=O)O')['canonical_smiles']
        opposite = raw.normalize_structure('C[C@H](O)C(=O)O')['canonical_smiles']
        self.assertEqual(a, b)
        self.assertNotEqual(a, opposite)

    def test_formula_not_protonation(self):
        self.assertEqual(raw.normalize_structure('CC(=O)[O-]', 'C2H4O2')['formula_status'], 'match')
        self.assertEqual(raw.normalize_structure('CCO', 'C3H6O')['formula_status'], 'heavy_atom_conflict')
        self.assertEqual(raw.normalize_structure('CCO', '(C2H4)n')['formula_status'], 'unparsed')
        self.assertEqual(raw.normalize_structure('CC(=O)[O-]', 'C2H3O2-')['formula_status'], 'match')
        self.assertIsNone(raw.formula_heavy('Fe3+'))

    def test_generic_structures_are_not_exact_molecules(self):
        self.assertIsNone(raw.normalize_structure('*O', 'HOR')['canonical_smiles'])
        self.assertEqual(raw.normalize_structure('*O', 'HOR')['structure_status'], 'generic_query_not_exact_molecule')
        self.assertIsNone(raw.formula_heavy('HR'))

    def test_missing_invalid_not_zero(self):
        self.assertIsNone(raw.normalize_structure('/')['canonical_smiles'])
        self.assertEqual(raw.normalize_structure('not_smiles')['structure_status'], 'invalid')

    def test_sequence_no_silent_ambiguity(self):
        self.assertTrue(raw.sequence_info('ACDEFGHIKL' * 20)[2])
        self.assertFalse(raw.sequence_info('ACDEFGHIKL' * 20 + 'X')[2])
        self.assertFalse(raw.sequence_info('ACD')[2])

    def test_accession_and_publication(self):
        self.assertEqual(raw.accession_tokens('P11712; A0A067LCX8 / P11712'), ['A0A067LCX8', 'P11712'])
        self.assertEqual(raw.publications('PMID:27272333', 'https://doi.org/10.1000/ABC'), ['DOI:10.1000/abc', 'PMID:27272333'])

    def test_no_overwrite_or_legacy_input(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            output = project / 'paper_rebuild_20260909/run_01'
            run = raw.RawBuild(project, output)
            old = project / 'data/processed/old.json'
            old.parent.mkdir(parents=True)
            old.write_text('{}')
            with self.assertRaises(ValueError): run.source(old)
            run.db.close()
            with self.assertRaises(FileExistsError): raw.RawBuild(project, output)

    def test_legacy_encoding_detected_without_character_loss(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'source.csv'
            path.write_bytes('Name,Label\nCYP,copyright \u00a9\n'.encode('cp1252'))
            self.assertEqual(raw.csv_encoding(path), 'cp1252')
            self.assertIn('\u00a9', path.read_text(encoding=raw.csv_encoding(path)))

    def test_end_to_end_synthetic_parser_fixture(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            folder = project / 'data/raw/p450rdb'
            folder.mkdir(parents=True)
            protein = {'Symbol': 'SYNTHETIC_FIXTURE_ONLY', 'Sequence': 'ACDEFGHIKL' * 20,
                       'Uniprot ID': 'P11712', 'Txid': '9606'}
            reaction = {'Symbol': 'SYNTHETIC_FIXTURE_ONLY', 'sequence': protein['Sequence'],
                        'Uniprot ID': 'P11712', 'Txid': '9606', 'PMID': '12345678',
                        'Substrate1': 'ethanol', 'sub_Smiles1': 'CCO', 'sub_formula1': 'C2H6O',
                        'Product1': 'acetaldehyde', 'pro_Smiles1': 'CC=O', 'pro_formula1': 'C2H4O'}
            for name, row in [('2. P450s.CSV', protein), ('1. Reactions.CSV', reaction)]:
                with (folder / name).open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerow(row)
            run = raw.RawBuild(project, project / 'paper_rebuild_20260909/test')
            run.p450rdb()
            self.assertEqual(run.db.execute('SELECT COUNT(*) FROM components').fetchone()[0], 2)
            self.assertEqual(run.db.execute('SELECT SUM(model_eligible) FROM assertions').fetchone()[0], 1)
            run.db.close()


if __name__ == '__main__':
    unittest.main()
