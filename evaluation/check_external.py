"""Recompute fixed-score macro AP and early release using only Python stdlib."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
METHODS=['isoform_specific_knn','pooled_chemical_knn','isoform_logistic']

def average_precision(rows,method):
    # Tie-aware threshold AP, not an arbitrary strict score tie ordering.
    buckets={}
    for r in rows:
        score=float(r[method]);buckets.setdefault(score,[]).append(int(r['label']))
    positive=sum(int(r['label']) for r in rows)
    if not positive:raise ValueError('Undefined AP without positives')
    tp=n=0;result=0.
    for score in sorted(buckets,reverse=True):
        y=buckets[score];new=sum(y);tp+=new;n+=len(y)
        result+=new/positive*tp/n
    return result

def main():
    with (ROOT/'strict_external_scores.tsv').open(encoding='utf-8') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==3035
    isoforms=sorted({r['isoform'] for r in rows})
    expected=dict(zip(METHODS,[.5428801073616204,.41097145717927014,.5341034089854914]))
    output={}
    for method in METHODS:
        values=[average_precision([r for r in rows if r['isoform']==iso],method) for iso in isoforms]
        macro=sum(values)/len(values)
        assert abs(macro-expected[method])<1e-12,(method,macro)
        release=[]
        for fraction in (.1,.25,.5,.75,1.0):
            precisions=[];recalls=[];released=hits=0
            for iso in isoforms:
                local=[r for r in rows if r['isoform']==iso]
                ranked=sorted(local,key=lambda r:(-float(r[method]),hashlib.sha256(r['compound_inchikey'].encode()).hexdigest()))
                selected=ranked[:math.ceil(fraction*len(local))]
                positive=sum(int(r['label']) for r in local);hit=sum(int(r['label']) for r in selected)
                precisions.append(hit/len(selected));recalls.append(hit/positive)
                released+=len(selected);hits+=hit
            release.append(dict(fraction=fraction,released=released,reported_positive_hits=hits,
                isoform_macro_precision=sum(precisions)/6,isoform_macro_recall=sum(recalls)/6))
        output[method]=dict(macro_AP=macro,release=release)
    with (ROOT.parent/'paper/Tables/Table_S13a_selective_release_tradeoffs.tsv').open(encoding='utf-8') as handle:
        expected_release=list(csv.DictReader(handle,delimiter='\t'))
    checked=0
    for expected_row in expected_release:
        method=expected_row['method']
        actual=next(r for r in output[method]['release'] if r['fraction']==float(expected_row['target_coverage']))
        assert actual['released']==int(expected_row['selected_rows'])
        for a,b in [('isoform_macro_precision','macro_precision'),('isoform_macro_recall','macro_positive_recall')]:
            assert math.isclose(actual[a],float(expected_row[b]),rel_tol=0,abs_tol=1e-12),(method,a)
        checked+=1
    print(json.dumps(dict(status='PASS',labels=3035,isoforms=6,release_rows_checked=checked,methods=output),indent=2))

if __name__=='__main__':main()
