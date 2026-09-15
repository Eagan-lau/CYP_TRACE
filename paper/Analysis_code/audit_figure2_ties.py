"""Fixed-score, exact tie-permutation sensitivity for Figure 2; no refitting."""
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'reviewer_reproduction'))
from evaluate_candidates import evaluate,read_table,rr


def main():
    path=ROOT/'reviewer_reproduction/figure2_candidates.tsv.gz'
    rows=read_table(path)
    deterministic=evaluate(rows)
    expected=json.loads((ROOT/'reviewer_reproduction/figure2_expected.json').read_text())
    expected.pop('input_sha256')
    assert deterministic==expected, 'Historical deterministic results changed'
    average=evaluate(rows,ties='average')
    assert {r['score'] for r in rows if r['applicability']}=={1.0}
    grouped=defaultdict(dict)
    for r in rows:grouped[(r['query'],r['panel'])].setdefault(r['method'],[]).append(r)
    per_panel=[]
    for key,methods in grouped.items():
        masks={m:{r['candidate'] for r in records if r['applicability']} for m,records in methods.items()}
        common=set.intersection(*masks.values())
        for domain in ('own','common'):
            for method,records in methods.items():
                allowed=masks[method] if domain=='own' else common
                positive=sum(r['positive'] for r in records if r['candidate'] in allowed)
                per_panel.append(dict(query=key[0],panel=key[1],group=records[0]['group'],method=method,
                    domain=domain,candidates=len(allowed),positives=positive,
                    deterministic_RR=rr(records,allowed),tie_averaged_RR=rr(records,allowed,'average'),
                    conditional_defined=bool(positive)))
    common=average['common_domain'][0]
    assert common['eligible_query_panels']==3 and common['protein_groups']==3
    assert common['difference']==0
    for method in average['methods']:
        for name in ('weighted_query_coverage','weighted_candidate_coverage','weighted_positive_instance_coverage'):
            assert average['methods'][method][name]==deterministic['methods'][method][name]
    summary=dict(status='PASS',scope='Post-evaluation exact averaging of within-score-block permutations; no model, eligibility, mask or primary result changed',
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        evaluator_sha256=hashlib.sha256((ROOT/'reviewer_reproduction/evaluate_candidates.py').read_bytes()).hexdigest(),
        all_applicable_scores_one=True,shared_tie_keys_verified=True,
        deterministic=deterministic,tie_averaged=average,
        eligible_common_panels=[r for r in per_panel if r['domain']=='common' and r['conditional_defined']])
    (ROOT/'tables/figure2_tie_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (ROOT/'tables/figure2_tie_per_panel.tsv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(per_panel[0]),delimiter='\t');writer.writeheader();writer.writerows(per_panel)
    (ROOT/'reviewer_reproduction/figure2_tie_average_expected.json').write_text(json.dumps({**average,'input_sha256':summary['input_sha256']},indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ('status','all_applicable_scores_one','shared_tie_keys_verified','eligible_common_panels')},indent=2))
    print(json.dumps({m:{k:v for k,v in r.items() if k in ('end_to_end_MRR','conditional_MRR','weighted_query_coverage')} for m,r in average['methods'].items()},indent=2))


if __name__=='__main__':main()
