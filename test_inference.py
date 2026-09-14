"""Actual frozen inference agreement and synthetic evidence-route smoke checks."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import unittest

from cyptrace_pipeline.evidence import screen
from cyptrace_pipeline.human import HumanSubstrateModel

ROOT = Path(__file__).resolve().parent


class InferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = HumanSubstrateModel.load()

    def test_frozen_training_complete(self):
        rows = json.loads((ROOT / 'data/development_labels.json').read_text(encoding='utf-8'))
        self.assertEqual(len(rows), 14955)
        self.assertEqual(len({r['compound_inchikey'] for r in rows}), 1768)
        self.assertEqual(sum(len(r['labels']) for r in self.model.compounds), 14955)
        self.assertEqual(self.model.bundle['source']['normalized_labels_sha256'],
                         hashlib.sha256((ROOT / 'data/development_labels.json').read_bytes()).hexdigest())

    def test_every_retained_external_score(self):
        labels = json.loads((ROOT / 'data/external_labels.json').read_text(encoding='utf-8'))
        structures = {r['compound_inchikey']: r['canonical_smiles'] for r in labels}
        with (ROOT / 'evaluation/strict_external_scores.tsv').open(encoding='utf-8') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 3035)
        for r in rows:
            actual = self.model.score(structures[r['compound_inchikey']], r['isoform'])
            for field, key in [('score', 'isoform_specific_knn'), ('pooled_chemistry_score', 'pooled_chemical_knn')]:
                self.assertTrue(math.isclose(actual[field], float(r[key]), rel_tol=0, abs_tol=1e-12),
                                (r['compound_inchikey'], r['isoform'], key))

    def test_scope_and_invalid_input(self):
        self.assertEqual(self.model.score('CCO', 'CYP3A4')['route'], 'FIXED_HUMAN_EXTERNALLY_VALIDATED')
        self.assertEqual(self.model.score('CCO', 'CYP2A6')['route'], 'FIXED_HUMAN_DEVELOPMENT_ONLY')
        self.assertTrue(self.model.score('CCO', 'CYP999A1')['abstained'])
        self.assertTrue(self.model.score('not-smiles', 'CYP3A4')['abstained'])
        self.assertIn('exact_training_compound', self.model.score('CCO', 'CYP3A4')['domain'])

    def test_synthetic_evidence_and_abstention(self):
        # Synthetic test-only labels; no claim about any biological sequence.
        seq = 'ACDEFGHIK'
        digest = hashlib.sha256(seq.encode()).hexdigest()
        index = {'sequences': {digest: seq},
                 'reactions': {'test': {'single_pair': True, 'substrates': ['CCO'], 'products': ['CC=O']}},
                 'edges': {digest + '|test': [{'sources': ['synthetic_unit_test']}]} }
        candidates = [{'candidate_id': 'test', 'reaction_key': 'test'}]
        exact = screen(index, {'synthetic': seq}, candidates)[0]
        unknown = screen(index, {'synthetic': 'GCDEFGHIK'}, candidates)[0]
        self.assertEqual(exact['route'], 'DOCUMENTED_EXACT_EVIDENCE')
        self.assertTrue(unknown['abstained'])
        self.assertIsNone(unknown['score'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
