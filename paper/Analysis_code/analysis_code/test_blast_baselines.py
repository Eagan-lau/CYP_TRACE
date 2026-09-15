import unittest
from run_blast_baselines import blast_row

class BlastTests(unittest.TestCase):
    def test_coordinate_coverage(self):
        r='q t 80 5 3 6 1 5 1e-10 20 10 20 AC-DE ACXDE'.split()
        q,t,identity,qcov,tcov,e,bits=blast_row(r)
        self.assertEqual(identity,.8); self.assertEqual(qcov,.4); self.assertEqual(tcov,.25)
    def test_bad_gap_span(self):
        with self.assertRaises(ValueError): blast_row('q t 80 5 3 7 1 5 1e-10 20 10 20 AC-DE ACXDE'.split())
    def test_bad_percent(self):
        with self.assertRaises(ValueError): blast_row('q t 120 5 1 5 1 5 1e-10 20 10 20 ACDEF ACDEF'.split())

if __name__=='__main__': unittest.main()
