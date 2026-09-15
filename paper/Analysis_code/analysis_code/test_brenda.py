import unittest
from build_brenda import typed_row


class BrendaTest(unittest.TestCase):
    def test_ec_context_not_expanded(self):
        row = typed_row('1.2.3.4', 'reaction', 0, {'value': 'a = b'},
                        {'1': {'source': 'uniprot', 'accessions': ['P00001']}}, {}, {'P00001'})
        self.assertEqual(row['linked_cyp_accessions'], [])
        self.assertEqual(row['scope'], 'EC_CONTEXT_NO_PROTEIN_LINK')

    def test_multiple_and_non_substrate_records_stay_typed(self):
        proteins = {p: {'source': 'uniprot', 'accessions': [a]} for p, a in [('1', 'P00001'), ('2', 'P00002')]}
        row = typed_row('1.2.3.4', 'inhibitor', 0, {'proteins': ['1', '2'], 'references': ['4'], 'comment': 'mutant F87A'},
                        proteins, {'4': {'pmid': 123}}, {'P00001', 'P00002'})
        self.assertEqual(row['scope'], 'MULTIPLE_SOURCE_PROTEIN_LINKS_NOT_EXPANDED')
        self.assertEqual(row['publication_ids'], ['PMID:123'])
        self.assertFalse(row['eligible_exact_reaction_positive'])
        self.assertTrue(row['variant_text_flag'])


if __name__ == '__main__': unittest.main()
