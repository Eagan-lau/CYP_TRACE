import unittest
import numpy as np
from extract_esm_features import windows

class WindowTests(unittest.TestCase):
    def test_no_truncation(self):
        for n in [1,500,1022,1023,1788,2100,5000]:
            counts=np.zeros(n,dtype=int)
            spans=list(windows(n))
            for a,b in spans:
                self.assertLessEqual(b-a,1022); counts[a:b]+=1
            self.assertTrue(np.all(counts>0)); self.assertEqual(spans[-1][1],n)
    def test_short_sequence(self):
        self.assertEqual(list(windows(500)),[(0,500)])
    def test_invalid_overlap(self):
        with self.assertRaises(ValueError): list(windows(10,10,10))

if __name__=='__main__': unittest.main()
