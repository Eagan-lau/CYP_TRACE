"""Source-semantic tests; no baseline score or tuned cutoff is involved."""
import json
from pathlib import Path
import sqlite3
import unittest
from rdkit import Chem

from reaction_roles_v1 import RoleRegistry, REFERENCE_ANCHORS, helper_equivalence_key
from build_reaction_core_v2 import experimental_publications, macro_ids, declared_directions


class ReactionRolesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = sqlite3.connect(':memory:')
        cls.db.row_factory = sqlite3.Row
        cls.db.execute('CREATE TABLE components(name TEXT, canonical_smiles TEXT, structure_status TEXT, formula_status TEXT)')
        for name in ['FMNH2', 'FMN']:
            cls.db.execute('INSERT INTO components VALUES(?,?,?,?)', (name,
                Chem.MolToSmiles(Chem.MolFromSmiles(REFERENCE_ANCHORS[name]['smiles'])), 'valid', 'match'))
        # A name/formula-consistent but redox-mislabeled source row must not win a vote.
        cls.wrong_fmnh2 = 'Cc1cc2nc3c(=O)nc([O-])nc-3n(C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-])c2cc1C'
        cls.db.execute('INSERT INTO components VALUES(?,?,?,?)', ('FMNH2', cls.wrong_fmnh2, 'valid', 'match'))
        cls.roles = RoleRegistry(cls.db)
        cls.red = REFERENCE_ANCHORS['FMNH2']['smiles']
        cls.ox = REFERENCE_ANCHORS['FMN']['smiles']

    def test_authoritative_redox_anchor_rejects_source_name_conflict(self):
        self.assertEqual(len(self.roles.conflicts), 1)
        self.assertEqual(self.roles.conflicts[0]['name'], 'FMNH2')
        self.assertNotEqual(helper_equivalence_key(self.red), helper_equivalence_key(self.ox))
        self.assertEqual(self.roles.classify(self.wrong_fmnh2), 'FMN')

    def test_carrier_and_inorganic_separation(self):
        p = self.roles.project(['CC', self.red, 'O=O'], ['CCO', self.ox, 'O', '[H+]'])
        self.assertEqual(p['substrates'], ['CC'])
        self.assertEqual(p['products'], ['CCO'])
        self.assertEqual(len(p['removed_participants']), 5)
        self.assertFalse(p['problems'])

    def test_unpaired_helper_blocks_admission(self):
        p = self.roles.project(['CC', self.red], ['CCO'])
        self.assertIn('unpaired_known_auxiliary_remains', p['problems'])
        self.assertIn(self.red, p['substrates'])

    def test_unequal_stoichiometry_not_erased(self):
        p = self.roles.project(['CC', self.red, self.red], ['CCO', self.ox])
        self.assertIn('unequal_paired_carrier_stoichiometry', p['problems'])
        self.assertEqual(p['substrates'].count(self.red), 2)
        self.assertEqual(p['products'].count(self.ox), 1)

    def test_carbon_coproduct_and_main_multiplicity_preserved(self):
        p = self.roles.project(['CNC', 'CNC', 'O=O'], ['CN', 'C=O', 'O=C=O', 'O'])
        self.assertEqual(p['substrates'].count('CNC'), 2)
        self.assertEqual(p['products'], ['C=O', 'CN', 'O=C=O'])
        self.assertFalse(p['single_pair'])
        self.assertFalse(p['problems'])

    def test_main_charge_is_not_neutralized(self):
        p = self.roles.project(['CC(=O)[O-]'], ['CC(=O)O'])
        self.assertEqual(p['substrates'], ['CC(=O)[O-]'])
        self.assertEqual(p['products'], ['CC(=O)O'])

    def test_main_stereo_not_erased(self):
        a = self.roles.project(['C[C@H](O)C(=O)O'], ['CC(=O)C(=O)O'])
        b = self.roles.project(['C[C@@H](O)C(=O)O'], ['CC(=O)C(=O)O'])
        self.assertNotEqual(a['reaction_key'], b['reaction_key'])

    def test_identity_and_carrier_only_invalid(self):
        self.assertIn('identity_main_transformation', self.roles.project(['CC'], ['CC'])['problems'])
        p = self.roles.project([self.red], [self.ox])
        self.assertIn('empty_substrate_after_auxiliary_separation', p['problems'])
        self.assertIn('empty_product_after_auxiliary_separation', p['problems'])

    def test_reverse_redox_pair_and_direction_preserved(self):
        a = self.roles.project(['CC', self.red], ['CCO', self.ox])
        b = self.roles.project(['CCO', self.ox], ['CC', self.red])
        self.assertFalse(b['problems'])
        self.assertNotEqual(a['reaction_key'], b['reaction_key'])

    def test_evidence_does_not_leak_from_direction_to_reaction(self):
        raw = {'comment': 'CATALYTIC ACTIVITY: Reaction=x; Evidence={ECO:0000250|PubMed:123}; PhysiologicalDirection=left-to-right; Xref=Rhea:RHEA:46309; Evidence={ECO:0000269|PubMed:456};'}
        self.assertEqual(experimental_publications('SwissProt', raw), [])
        raw['comment'] = raw['comment'].replace('ECO:0000250', 'ECO:0000269')
        self.assertEqual(experimental_publications('SwissProt', raw), ['PMID:123'])
        self.assertEqual(declared_directions('SwissProt', raw), {'46309': 'left-to-right'})

    def test_api_evidence_source_and_macro_ids(self):
        raw = {'reaction': {'evidences': [
            {'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '123'},
            {'evidenceCode': 'ECO:0000250', 'source': 'PubMed', 'id': '456'},
            {'evidenceCode': 'ECO:0000269', 'source': 'DOI', 'id': '789'}],
            'reactionCrossReferences': [{'id': 'RHEA:46308'}, {'id': 'RHEA-COMP:11964'}]},
            'physiologicalReactions': [{'evidences': [{'evidenceCode': 'ECO:0000269', 'source': 'PubMed', 'id': '999'}]}]}
        self.assertEqual(experimental_publications('UniProt_API', raw), ['PMID:123'])
        self.assertEqual(macro_ids('UniProt_API', raw), ['11964'])
        self.assertEqual(macro_ids('SwissProt', {'comment': 'Rhea:RHEA-   COMP:11965'}), ['11965'])


if __name__ == '__main__': unittest.main()
