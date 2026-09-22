import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'paper/Analysis_code/analysis_code'))
import build_clean_core_features as clean
from rebuild_core import verify_inputs

class PortabilityTests(unittest.TestCase):
    def test_clean_cache_is_optional(self):
        with tempfile.TemporaryDirectory() as folder:
            mapping,array,inputs=clean.optional_prior_cache(Path(folder))
            self.assertEqual(mapping,{})
            self.assertEqual(array.shape,(0,1280))
            self.assertEqual(inputs,[])

    def test_partial_cache_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder)/'query_manifest.tsv').touch()
            with self.assertRaises(ValueError): clean.optional_prior_cache(Path(folder))

    def test_missing_raw_inputs_are_explicit(self):
        with tempfile.TemporaryDirectory() as folder:
            rows=verify_inputs(Path(folder))
            self.assertEqual(len(rows),27)
            self.assertTrue(all(r['status']=='missing' for r in rows))

if __name__=='__main__': unittest.main()
