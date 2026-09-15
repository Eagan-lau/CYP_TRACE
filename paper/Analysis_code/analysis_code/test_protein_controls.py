import unittest
import numpy as np
from protein_control_core import permute_profiles


class Controls(unittest.TestCase):
    def test_no_heldout_profile_moves(self):
        y = np.array([[1, 0, 0], [0, 0, 0], [1, 1, 0], [0, 0, 1.]])
        changed, receipt = permute_profiles(y, np.random.default_rng(13))
        self.assertTrue(np.array_equal(changed[1], y[1]))
        self.assertEqual(set(receipt['donor_indices']), {0, 2, 3})
        self.assertTrue(np.array_equal(changed.sum(0), y.sum(0)))
    def test_balanced_chemical_backbone_invariant(self):
        y = np.array([[1, 1, 0], [1, 0, 0], [0, 0, 1.], [0, 0, 0]])
        changed, _ = permute_profiles(y, np.random.default_rng(15))
        balance = lambda m: (m / np.maximum(m.sum(1, keepdims=True), 1)).sum(0)
        self.assertTrue(np.array_equal(balance(y), balance(changed)))
    def test_singleton_identity_reported(self):
        changed, receipt = permute_profiles([[0, 0], [1, 0]], np.random.default_rng(1))
        self.assertEqual(receipt['fixed_points'], 1)
    def test_reproducible(self):
        y = np.eye(10)
        a = permute_profiles(y, np.random.default_rng(111))
        b = permute_profiles(y, np.random.default_rng(111))
        self.assertEqual(a[1], b[1])
    def test_invalid_incidence_rejected(self):
        with self.assertRaises(ValueError): permute_profiles([[1, -1]], np.random.default_rng(1))


if __name__ == '__main__': unittest.main()
