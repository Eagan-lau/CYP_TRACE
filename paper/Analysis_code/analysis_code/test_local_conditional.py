import unittest
import numpy as np
from local_conditional_core import local_descriptions, joint_objective, select_joint


class LocalConditionalTests(unittest.TestCase):
    def setUp(self):
        self.lp = np.log(np.array([[.6, .3, .1], [.2, .4, .4], [.4, .2, .4], [.3, .3, .4]]))
        self.d = np.array([[[1., 0.], [-1, 1], [.3, -.4]], [[-.2, .6], [.5, -.3], [0., 0.]],
                           [[.1, -.3], [.2, .7], [-.3, .1]], [[.7, -.2], [.3, .1], [-.4, .2]]])
        self.targets = [[0], [1], [2], [0]]
        self.rows = [{'group': str(i), 'query': str(i)} for i in range(4)]
    def test_gradient(self):
        f = joint_objective(self.lp, self.d, self.targets, self.rows); x = np.array([.3, .7]); _, grad = f(x)
        numeric = []
        for i in range(2):
            step = np.zeros(2); step[i] = 1e-6
            numeric.append((f(x + step)[0] - f(x - step)[0]) / 2e-6)
        self.assertTrue(np.allclose(grad, numeric, atol=1e-8))
    def test_zero_residuals_stay_zero(self):
        choice = select_joint(self.lp, np.zeros((4, 3)), np.zeros((4, 3)), self.targets, self.rows)
        self.assertEqual(choice['lambdas'], [0., 0.])
    def test_joint_not_worse_than_single_or_zero(self):
        choice = select_joint(self.lp, self.d[..., 0], self.d[..., 1], self.targets, self.rows)
        self.assertLessEqual(choice['selected_loss'], choice['best_single_or_zero_loss'] + 1e-10)
        self.assertTrue(all(v >= 0 for v in choice['lambdas']))
    def test_single_group_no_claim(self):
        rows = [{'group': 'a', 'query': str(i)} for i in range(4)]
        self.assertEqual(select_joint(self.lp, self.d[..., 0], self.d[..., 1], self.targets, rows)['lambdas'], [0., 0.])
    def test_order_matters_composition_does_not(self):
        x = np.zeros((2, 2, 21)); x[0, 0, 0] = x[0, 1, 1] = x[1, 0, 1] = x[1, 1, 0] = 1
        f = local_descriptions(x, np.ones((2, 2), bool), np.eye(2))
        self.assertFalse(np.array_equal(f['ordered'][0], f['ordered'][1]))
        self.assertTrue(np.array_equal(f['composition'][0], f['composition'][1]))
        self.assertTrue(np.all(f['missing'] == 0))
    def test_missing_not_unknown(self):
        x = np.zeros((1, 2, 21)); x[0, 0, 20] = 1
        f = local_descriptions(x, [[True, False]], [[1.]])
        self.assertTrue(np.array_equal(f['missing'][0] > 0, [True, False, False, True]))


if __name__ == '__main__': unittest.main()
