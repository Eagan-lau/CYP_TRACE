import unittest
import numpy as np
from conditional_core import fit_residual,predict_delta,normalized_features,select_lambda,group_weights

class ConditionalTests(unittest.TestCase):
    def test_zero_features_give_zero_residual(self):
        Y=np.eye(4); K=np.eye(4); P=np.zeros((4,3))
        fit=fit_residual(P,Y,K)
        self.assertTrue(fit['zero_residual']); self.assertTrue(np.all(predict_delta(fit,P)==0))
    def test_query_features_cannot_affect_fitting(self):
        rng=np.random.default_rng(2); P=normalized_features(rng.normal(size=(8,4))); Y=np.zeros((8,3)); Y[:6]=rng.integers(0,2,(6,3)); Y[:6,0]=1
        a=fit_residual(P,Y,np.eye(3)); P[6:]=999; b=fit_residual(P,Y,np.eye(3))
        np.testing.assert_array_equal(a['coefficients'],b['coefficients']); np.testing.assert_array_equal(a['prior_log'],b['prior_log'])
    def test_harmful_delta_selects_zero(self):
        prior=np.log(np.full((2,2),.5)); delta=np.array([[-1,1],[-1,1]])
        r=select_lambda(prior,delta,[[0],[0]],[{'group':'a','query':'a'},{'group':'b','query':'b'}])
        self.assertEqual(r['lambda'],0.)
    def test_helpful_delta_finite_optimum(self):
        prior=np.log(np.full((4,2),.5)); delta=np.tile([1.,-1.],(4,1))
        r=select_lambda(prior,delta,[[0],[0],[0],[1]],[{'group':str(i),'query':str(i)} for i in range(4)])
        self.assertAlmostEqual(r['lambda'],np.log(3)/2,places=6)
    def test_group_panel_weighting(self):
        rows=[{'group':'a','query':'1'},{'group':'a','query':'1'},{'group':'a','query':'2'},{'group':'b','query':'3'}]
        np.testing.assert_allclose(group_weights(rows),[.125,.125,.25,.5])
    def test_one_group_is_not_identifiable(self):
        p=np.log(np.full((2,2),.5)); d=np.tile([1.,-1.],(2,1))
        r=select_lambda(p,d,[[0],[1]],[{'group':'a','query':'1'},{'group':'a','query':'2'}])
        self.assertEqual(r['lambda'],0.)

if __name__=='__main__': unittest.main()
