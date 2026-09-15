"""Small regression checks for the portable exploratory logistic analysis."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from sklearn.metrics import average_precision_score

from reproduce_logistic import average_precision, ap_template, weighted_ap, choose, development_arrays


class LogisticTests(unittest.TestCase):
    def test_tied_ap(self):
        y = np.array([0, 1, 1, 0, 1])
        s = np.array([.2, .2, .7, .1, .1])
        self.assertAlmostEqual(average_precision(y, s), average_precision_score(y, s), places=14)
        w = np.array([0., 2., 3., 1., 0.])
        self.assertAlmostEqual(weighted_ap(ap_template(y, s), w),
                               average_precision_score(y, s, sample_weight=w), places=14)

    def test_no_positives_is_undefined(self):
        self.assertTrue(np.isnan(average_precision([0, 0], [.1, .3])))

    def test_selection_ties(self):
        rows = [dict(mean_fold_AP=.5, C=c, class_weight=w)
                for c in (1., .1, .01) for w in ("balanced", None)]
        best = choose(rows)
        self.assertEqual(best["C"], .01)
        self.assertIsNone(best["class_weight"])

    def invalid_data(self, rows, pattern):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "labels.json"
            path.write_text(json.dumps(rows), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, pattern):
                development_arrays(path)

    def test_duplicate_pair_rejected(self):
        row = dict(isoform="CYP3A4", compound_inchikey="a", label=1,
                   canonical_smiles="CCO", scaffold_group="a", development_fold=0)
        self.invalid_data([row, row], "Duplicate")

    def test_scaffold_leak_rejected(self):
        rows = [dict(isoform="CYP3A4", compound_inchikey=str(i), label=i,
                     canonical_smiles="CCO", scaffold_group="shared", development_fold=i)
                for i in range(2)]
        self.invalid_data(rows, "Scaffold crosses")


if __name__ == "__main__":
    unittest.main()
