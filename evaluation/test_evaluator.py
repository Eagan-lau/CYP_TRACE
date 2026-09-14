"""Independent small cases for domain accounting and input contracts."""
import copy
import itertools
import math
from pathlib import Path
import unittest
from evaluate_candidates import evaluate,read_table,expected_tied_rr,rr


class ContractTests(unittest.TestCase):
    def setUp(self):
        self.rows=read_table(Path(__file__).with_name('example_candidates.tsv'))

    def test_decomposition_and_common_domain(self):
        result=evaluate(self.rows)
        hom=result['methods']['homology'];esm=result['methods']['embedding']
        self.assertEqual(hom['weighted_query_coverage'],.5)
        self.assertEqual(hom['end_to_end_MRR'],.5)
        self.assertEqual(hom['conditional_MRR'],1)
        self.assertEqual(esm['end_to_end_MRR'],.75)
        self.assertEqual(result['common_domain'][0]['eligible_query_panels'],1)
        self.assertEqual(result['common_domain'][0]['difference'],-.5)

    def test_missing_domain_is_undefined(self):
        rows=copy.deepcopy(self.rows)
        for r in rows:
            if r['method']=='homology':
                r['applicability']=0;r['score']=None
        result=evaluate(rows)
        self.assertEqual(result['methods']['homology']['end_to_end_MRR'],0)
        self.assertIsNone(result['methods']['homology']['conditional_MRR'])
        self.assertIsNone(result['common_domain'][0]['difference'])

    def test_mismatched_denominators_rejected(self):
        with self.assertRaises(ValueError):evaluate(self.rows[:-1])

    def test_conflicting_labels_rejected(self):
        rows=copy.deepcopy(self.rows);rows[0]['positive']=0
        with self.assertRaises(ValueError):evaluate(rows)

    def test_repeated_panels_do_not_overweight_queries(self):
        extra=[{**r,'panel':'replicate_panel'} for r in self.rows if r['query']=='query1']
        result=evaluate(self.rows+extra)
        self.assertEqual(result['methods']['embedding']['end_to_end_MRR'],.75)

    def test_tied_expectation_by_exhaustive_positive_positions(self):
        for size in range(1,9):
            for positives in range(1,size+1):
                for higher in (0,2,7):
                    exact=math.fsum(1/(higher+min(p)) for p in itertools.combinations(range(1,size+1),positives))/math.comb(size,positives)
                    self.assertAlmostEqual(expected_tied_rr(higher,size,positives),exact,places=14)

    def test_tied_expectation_multiple_score_blocks(self):
        rows=[dict(candidate=str(i),score=s,tie_key=str(i),positive=p)
              for i,(s,p) in enumerate([(3,0),(3,0),(2,0),(2,1),(2,1),(1,1)])]
        self.assertAlmostEqual(rr(rows,{r['candidate'] for r in rows},'average'),(2/3)/3+(1/3)/4)

    def test_average_ignores_tie_identifiers(self):
        rows=copy.deepcopy(self.rows)
        before=evaluate(rows,ties='average')
        for r in rows:r['tie_key']=''.join(chr(255-ord(c)) for c in r['tie_key'])
        self.assertEqual(before,evaluate(rows,ties='average'))

    def test_average_preserves_coverage(self):
        det=evaluate(self.rows);avg=evaluate(self.rows,ties='average')
        for m in det['methods']:
            for name in ('weighted_query_coverage','weighted_candidate_coverage','weighted_positive_instance_coverage'):
                self.assertEqual(det['methods'][m][name],avg['methods'][m][name])
            self.assertLess(avg['methods'][m]['decomposition_error'],1e-12)


if __name__=='__main__':unittest.main()
