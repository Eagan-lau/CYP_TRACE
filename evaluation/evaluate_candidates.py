"""Evaluate supplied candidate tables with explicit applicability and hierarchy.

Standard-library-only CLI. Scores outside applicability must be empty. Candidate
and target sets must agree across methods; unknown alternatives are not negatives.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path

REQUIRED = {'method','query','candidate','score','applicability','positive','group'}


def read_table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    rows = []
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle,delimiter='\t')
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError('Missing fields: '+str(sorted(REQUIRED-set(reader.fieldnames or []))))
        for number,r in enumerate(reader,2):
            for name in ('method','query','candidate','group'):
                if not r[name]:
                    raise ValueError(f'Row {number}: empty {name}')
            for name in ('applicability','positive'):
                if r[name] not in ('0','1'):
                    raise ValueError(f'Row {number}: {name} must be 0 or 1')
                r[name] = int(r[name])
            if r['applicability']:
                r['score'] = float(r['score'])
                if not math.isfinite(r['score']):
                    raise ValueError(f'Row {number}: non-finite applicable score')
            else:
                if r['score'] not in ('','NA','null'):
                    raise ValueError(f'Row {number}: domain-external score must be empty, not zero')
                r['score'] = None
            r.setdefault('panel','default')
            r.setdefault('exposure','unspecified')
            r.setdefault('tie_key',hashlib.sha256(r['candidate'].encode()).hexdigest())
            if not r['panel'] or not r['tie_key'] or not r['exposure']:
                raise ValueError(f'Row {number}: empty panel, exposure or tie key')
            rows.append(r)
    return rows


def hierarchy(panel_records):
    groups = defaultdict(lambda:defaultdict(list))
    for key,r in panel_records.items():
        groups[r['group']][r['query']].append(key)
    weights = {}
    for queries in groups.values():
        for panels in queries.values():
            for panel in panels:
                weights[panel] = 1/len(groups)/len(queries)/len(panels)
    return weights


def expected_tied_rr(higher, size, positives):
    """Exact E[1 / first-positive rank] over uniform within-block permutations.

    The first positive occupies position j with probability
    choose(size-j, positives-1) / choose(size, positives). This recurrence
    avoids large factorials. Higher-score blocks contain no positives.
    """
    if positives == 0:
        return 0.0
    if not (higher >= 0 and 1 <= positives <= size):
        raise ValueError('Invalid tie-block counts')
    probability = positives / size
    terms = []
    for position in range(1, size-positives+2):
        terms.append(probability / (higher+position))
        if position < size-positives+1:
            probability *= (size-positives-position+1) / (size-position)
    return math.fsum(terms)


def rr(records, allowed, ties='deterministic'):
    ranked = sorted([r for r in records if r['candidate'] in allowed],
                    key=lambda r:(-r['score'],r['tie_key']))
    if ties == 'deterministic':
        return next((1/i for i,r in enumerate(ranked,1) if r['positive']),0.0)
    higher = 0
    for _, block_iter in itertools.groupby(ranked, key=lambda r:r['score']):
        block = list(block_iter)
        positives = sum(r['positive'] for r in block)
        if positives:
            return expected_tied_rr(higher, len(block), positives)
        higher += len(block)
    return 0.0


def summarize(records):
    if not records:
        return {'query_panels':0,'protein_groups':0,'end_to_end_MRR':None,
                'weighted_query_coverage':None,'conditional_MRR':None}
    weights=hierarchy(records)
    cw=sum(weights[k] for k,r in records.items() if r['covered_positives'])
    end=sum(weights[k]*r['rr'] for k,r in records.items())
    conditional=end/cw if cw else None
    return dict(query_panels=len(records),unique_queries=len({r['query'] for r in records.values()}),
        protein_groups=len({r['group'] for r in records.values()}),
        weighted_candidate_coverage=sum(weights[k]*r['covered_candidates']/r['candidates'] for k,r in records.items()),
        weighted_positive_instance_coverage=sum(weights[k]*r['covered_positives']/r['positives'] for k,r in records.items()),
        micro_positive_instance_coverage=sum(r['covered_positives'] for r in records.values())/sum(r['positives'] for r in records.values()),
        weighted_query_coverage=cw,end_to_end_MRR=end,conditional_MRR=conditional,
        decomposition_error=abs(end-cw*conditional) if conditional is not None else None,
        positive_instances=sum(r['positives'] for r in records.values()),
        covered_positive_instances=sum(r['covered_positives'] for r in records.values()))


def evaluate(rows, ties='deterministic'):
    if ties not in ('deterministic','average'):
        raise ValueError('Unknown tie policy')
    if not rows:
        raise ValueError('Empty table')
    grouped=defaultdict(lambda:defaultdict(dict))
    query_meta={}
    for r in rows:
        key=(r['query'],r['panel'])
        if r['candidate'] in grouped[r['method']][key]:
            raise ValueError('Duplicate method/query/panel/candidate record')
        grouped[r['method']][key][r['candidate']]=r
        meta=(r['group'],r['exposure'])
        if r['query'] in query_meta and query_meta[r['query']]!=meta:
            raise ValueError('Group or exposure changes within query')
        query_meta[r['query']]=meta
    methods=sorted(grouped)
    baseline=grouped[methods[0]]
    for method in methods:
        if set(grouped[method])!=set(baseline):
            raise ValueError('Methods do not contain identical query panels')
        for key, candidates in grouped[method].items():
            if set(candidates)!=set(baseline[key]):
                raise ValueError('Candidate denominators differ across methods')
            if not sum(r['positive'] for r in candidates.values()):
                raise ValueError('Every panel needs at least one documented target')
            if len({r['tie_key'] for r in candidates.values()}) != len(candidates):
                raise ValueError('Tie keys must be unique within a panel')
            for c,r in candidates.items():
                if (r['positive'],r['tie_key'])!=(baseline[key][c]['positive'],baseline[key][c]['tie_key']):
                    raise ValueError('Target labels or tie rules differ across methods')
    output={}; panels={}
    for method in methods:
        records={}
        for key, mapping in grouped[method].items():
            local=list(mapping.values()); first=local[0]
            allowed={r['candidate'] for r in local if r['applicability']}
            records[key]=dict(query=key[0],panel=key[1],group=first['group'],exposure=first['exposure'],
                candidates=len(local),positives=sum(r['positive'] for r in local),
                covered_candidates=len(allowed),covered_positives=sum(r['positive'] for r in local if r['applicability']),
                rr=rr(local,allowed,ties))
        panels[method]=records
        output[method]=summarize(records)
        output[method]['exposure_strata']={label:summarize({k:r for k,r in records.items() if r['exposure']==label})
                                           for label in sorted({r['exposure'] for r in records.values()})}
    common=[]
    for left,right in itertools.combinations(methods,2):
        records={left:{},right:{}}
        shared_targets=0
        for key in baseline:
            a=grouped[left][key];b=grouped[right][key]
            allowed={c for c,r in a.items() if r['applicability'] and b[c]['applicability']}
            positive=sum(a[c]['positive'] for c in allowed)
            if not positive:
                continue
            shared_targets+=positive
            for method,mapping in ((left,a),(right,b)):
                prior=panels[method][key]
                records[method][key]={**prior,'candidates':len(allowed),'covered_candidates':len(allowed),
                    'positives':positive,'covered_positives':positive,'rr':rr(list(mapping.values()),allowed,ties)}
        # Common-domain estimand rebuilds the same panel/query/group hierarchy
        # on eligible units, matching the article's three-group comparison.
        a=summarize(records[left]);b=summarize(records[right])
        common.append(dict(left=left,right=right,eligible_query_panels=len(records[left]),
            undefined_query_panels=len(baseline)-len(records[left]),
            protein_groups=a['protein_groups'],common_positive_instances=shared_targets,
            left_conditional_MRR=a['conditional_MRR'],right_conditional_MRR=b['conditional_MRR'],
            difference=(a['conditional_MRR']-b['conditional_MRR']) if records[left] else None))
    return dict(schema_version=1,rows=len(rows),methods=output,common_domain=common,
        weighting='equal groups, equal queries within group, equal panels within query; conditional MRR renormalizes original weights; common domain rebuilds hierarchy on eligible units',
        ties=('ascending supplied tie_key; defaults to SHA256(candidate)' if ties == 'deterministic'
              else 'exact expectation over uniform permutations within each equal-score block'),
        positive_definition='documented target; positive=0 means other supplied alternative, not certified negative',
        score_calibration='not performed',uncertainty='not estimated by this evaluator')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--ties',choices=['deterministic','average'],default='deterministic')
    args=parser.parse_args()
    result=evaluate(read_table(args.input),ties=args.ties)
    result['input_sha256']=hashlib.sha256(args.input.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
