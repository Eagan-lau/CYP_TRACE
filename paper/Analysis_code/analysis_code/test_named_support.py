import unittest
from build_named_support import brenda_equation, name_key, resolve_names


class NameResolutionTest(unittest.TestCase):
    def test_ions_and_reversibility(self):
        terms, flag, problem = brenda_equation('ethanol + NAD+ = acetaldehyde + NADH + H+ {ir}')
        self.assertEqual(terms, (['ethanol', 'NAD+'], ['acetaldehyde', 'NADH', 'H+']))
        self.assertEqual(flag, 'ir')
        self.assertIsNone(problem)

    def test_case_and_stereodescriptors_are_preserved(self):
        self.assertNotEqual(name_key('(R)-x'), name_key('(r)-x'))
        self.assertNotEqual(name_key('Test'), name_key('test'))
        self.assertEqual(name_key(' a  b '), 'a b')

    def test_ambiguous_and_missing_names_never_imputed(self):
        r = resolve_names(['x', 'y', '?'], {'x': {'smiles': 'CC', 'basis': 'test'}}, {'y': {}})
        self.assertEqual(r['structures'], ['CC'])
        self.assertEqual(len(r['unresolved_terms']), 2)
        self.assertEqual(r['unresolved_terms'][0]['reason'], 'ambiguous_source_name')

    def test_integer_coefficient_preserved(self):
        r = resolve_names(['2 x'], {'x': {'smiles': 'CC', 'basis': 'test'}}, {})
        self.assertEqual(r['structures'], ['CC', 'CC'])
        self.assertEqual(r['resolved_terms'][0]['coefficient'], 2)

    def test_incomplete_equation_not_completed(self):
        self.assertEqual(brenda_equation('x')[2], 'not_one_explicit_equation')
        self.assertEqual(brenda_equation('x = y = z')[2], 'not_one_explicit_equation')


if __name__ == '__main__': unittest.main()
