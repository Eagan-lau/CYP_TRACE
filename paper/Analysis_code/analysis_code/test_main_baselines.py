import itertools
import unittest
import numpy as np
from run_main_baselines import random_expectations,rank_result,train_scores


class MainBaselineTest(unittest.TestCase):
    def test_uniform_matches_exhaustive_permutations(self):
        for n,k in [(4,1),(4,2),(4,4)]:
            truths=set(range(k)); observed=[]
            for permutation in itertools.permutations(range(n)):
                ranks=[i+1 for i,v in enumerate(permutation) if v in truths]
                observed.append(1/min(ranks))
            self.assertAlmostEqual(random_expectations(n,k)['rr'],np.mean(observed))

    def test_domain_mask_not_score_zero(self):
        m=rank_result(np.array([4.,0.]),np.array([True,False]),[1],np.arange(2))
        self.assertEqual(m['rr'],0); self.assertEqual(m['covered_positives'],0)
        self.assertIsNone(m['conditional_rr'])
        n=rank_result(np.array([4.,0.]),np.ones(2,dtype=bool),[1],np.arange(2))
        self.assertEqual(n['rr'],.5)

    def test_homology_cannot_score_unseen_labels_but_transport_can(self):
        Y=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,0.]])
        H=np.array([[10.,1.,0.],[1.,10.,0.],[2.,1.,10.]])
        kernel=np.array([[1.,.1,.8],[.1,1.,.3],[.8,.3,1.]])
        scores,masks,best=train_scores(Y,kernel,H,[2],np.array([0,1]))
        self.assertEqual(best,[0]); self.assertFalse(masks['mmseqs_weighted'][0,2])
        self.assertTrue(masks['homology_chemical_transport'][0,2])
        self.assertAlmostEqual(scores['homology_chemical_transport'][0,2],(2*.8+.3)/3)

    def test_no_homology_means_missing_expert_not_fabricated_scores(self):
        Y=np.array([[1.,0.],[0.,0.]])
        scores,masks,best=train_scores(Y,np.eye(2),np.eye(2),[1],np.array([0]))
        self.assertEqual(best,[None]); self.assertFalse(masks['homology_chemical_transport'].any())
        self.assertTrue(masks['balanced_chemistry_prior'].all())

    def test_exact_protein_context_allowed_only_from_training_labels(self):
        Y=np.array([[1.,0.],[0.,0.]])
        scores,masks,best=train_scores(Y,np.eye(2),np.eye(2)*10,[0],np.array([0]))
        self.assertEqual(best,[0]); self.assertTrue(masks['mmseqs_top1'][0,0])
        self.assertFalse(masks['mmseqs_top1'][0,1])


if __name__=='__main__': unittest.main()
