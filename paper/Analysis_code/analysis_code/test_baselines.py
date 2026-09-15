import unittest
import numpy as np
from split_baselines import Union, assign_folds, rank_metrics, chemical_arrays


class BaselineTests(unittest.TestCase):
    def test_groups_are_transitive_and_deterministic(self):
        uf = Union(['a', 'b', 'c', 'd'])
        uf.join('a', 'b'); uf.join('b', 'c')
        self.assertEqual(sorted(map(len, uf.groups().values())), [1, 3])
        self.assertEqual(assign_folds(uf.groups(), 2), assign_folds(uf.groups(), 2))

    def test_unscored_positive_retains_full_denominator(self):
        result = rank_metrics(np.array([2., 1., 0.]), np.array([True, True, False]), [2], np.arange(3))
        self.assertEqual(result['rr'], 0)
        self.assertIsNone(result['conditional_rr'])
        self.assertEqual(result['positive_count'], 1)
        self.assertEqual(result['covered_positives'], 0)

    def test_shared_ties_and_multi_positive_recall(self):
        result = rank_metrics(np.ones(4), np.ones(4, dtype=bool), [1, 3], np.arange(4))
        self.assertEqual(result['rr'], .5)
        self.assertEqual(result['recall_at_1'], 0)
        self.assertEqual(result['recall_at_5'], 1)

    def test_direction_is_not_lost_in_kernel(self):
        data = {'a': {'substrates': ['CCO'], 'products': ['CC=O']},
                'b': {'substrates': ['CC=O'], 'products': ['CCO']}}
        kernel, scaffolds = chemical_arrays(data, ['a', 'b'],
            {'fingerprint': {'bits_per_side': 2048, 'radius': 2, 'include_chirality': True}})
        self.assertEqual(kernel[0, 0], 1)
        self.assertLess(kernel[0, 1], 1)
        self.assertEqual(scaffolds['a'], scaffolds['b'])


if __name__ == '__main__': unittest.main()
