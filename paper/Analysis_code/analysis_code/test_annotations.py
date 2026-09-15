import unittest
from audit_annotations import rhea_links, exact_reaction


class AnnotationTest(unittest.TestCase):
    def test_master_and_direction_stay_distinct(self):
        raw = {'reaction': {'reactionCrossReferences': [{'database': 'Rhea', 'id': 'RHEA:12345'}, {'database': 'Rhea', 'id': 'RHEA-COMP:999'}]},
               'physiologicalReactions': [{'reactionCrossReference': {'database': 'Rhea', 'id': 'RHEA:12346'}}]}
        self.assertEqual(rhea_links('UniProt_API', raw), (['12345'], ['12346']))

    def test_reverse_and_multiplicity_not_collapsed(self):
        self.assertNotEqual(exact_reaction('CO>>C=O')['reaction_key'], exact_reaction('C=O>>CO')['reaction_key'])
        self.assertNotEqual(exact_reaction('CO.CO>>C=O')['reaction_key'], exact_reaction('CO>>C=O')['reaction_key'])
        self.assertIsNone(exact_reaction('*CO>>CC=O')['reaction_key'])


if __name__ == '__main__': unittest.main()
