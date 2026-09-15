"""Reconstruct common-domain rankings and paired group bootstrap via counts."""
from collections import defaultdict
import json
import math
from pathlib import Path
import numpy as np
from run_raw import stable,digest_file,write_json,now
from verify_main_baselines import uniform

root=Path(__file__).resolve().parent; result=root/'paired_comparison_01'; output=root/'paired_validation_01/INDEPENDENT_QC.json'
if output.exists(): raise FileExistsError(output)
read=lambda p:json.loads(p.read_text())
declared=read(result/'comparison_summary.json')
rows=[json.loads(line) for line in (result/'paired_query_records.jsonl').read_text().splitlines()]
rids=sorted(read(root/'dataset_02/core_reactions.json')); tie=[stable(r) for r in rids]
methods={}; source_rows={}
for folder in ['main_baselines_01','conditional_model_01']:
    for line in (root/folder/'per_query_metrics.jsonl').read_text().splitlines():
        r=json.loads(line); methods[r['method']]=folder
        source_rows[(r['method'],r['task'],r['block_number'],r['sequence_sha256'])]=r
failures=[]; checked=0; maxdiff=0.
for f,h in read(result/'input_manifest.json').items():
    if digest_file(root/f)!=h: failures.append('input:'+f)
byblock=defaultdict(list)
for r in rows: byblock[r['block_number']].append(r)
for bn,members in byblock.items():
    arrays={f:np.load(root/f/f'block_{bn:03d}_scores.npz',allow_pickle=False) for f in set(methods.values())}
    cache={}
    for row in members:
        sources=[]; masks=[]; scores=[]
        for side in ['left','right']:
            method=row[side+'_method']; r=source_rows[(method,row['task'],bn,row['sequence_sha256'])]
            if abs(r['rr']-row[side+'_rr'])>1e-12 or r['covered_positives']!=row[side+'_covered_positives']: failures.append('end_to_end')
            sources.append(r)
            if method=='uniform_expectation': scores.append(None); masks.append(np.ones(len(rids),dtype=bool)); continue
            data=arrays[methods[method]]; qi=r['query_index_in_block']
            score=data[method]; mask=data['domain_'+method]
            scores.append(score[qi] if score.ndim==2 else score); masks.append(mask[qi] if mask.ndim==2 else mask)
        if sources[0]['target_reaction_indices']!=sources[1]['target_reaction_indices']: failures.append('unmatched_positives')
        common=[i for i in range(len(rids)) if masks[0][i] and masks[1][i]]; common_set=set(common)
        truth=[i for i in sources[0]['target_reaction_indices'] if i in common_set]
        if len(common)!=row['common_candidate_count'] or len(truth)!=row['common_positive_count']: failures.append('common_count')
        if bool(truth)!=row['common_reranking_defined']: failures.append('undefined_reranking')
        for side,score in zip(['left','right'],scores):
            if not truth:
                if row['common_'+side+'_rr'] is not None: failures.append('undefined_not_null')
                continue
            if score is None: actual=uniform(len(common),len(truth))['rr']
            else:
                key=(row[side+'_method'],row['sequence_sha256'],tuple(common))
                if key not in cache:
                    order=sorted(common,key=lambda i:(-float(score[i]),tie[i]))
                    cache[key]={i:j+1 for j,i in enumerate(order)}
                actual=1/min(cache[key][i] for i in truth)
            difference=abs(actual-row['common_'+side+'_rr']); maxdiff=max(maxdiff,difference)
            if difference>1e-12: failures.append('common_rank')
        checked+=1
    for data in arrays.values(): data.close()
    print(f'Paired QC block {bn+1}/35',flush=True)
rng=np.random.Generator(np.random.PCG64(declared['bootstrap_seed'])); summaries=0; intervals=0
for comparison in declared['comparisons']:
    selected=[r for r in rows if (r['task'],r['left_method'],r['right_method'])==(comparison['task'],comparison['left_method'],comparison['right_method'])]
    for scope,left,right,subset in [('end_to_end','left_rr','right_rr',selected),
        ('common_domain_reranking','common_left_rr','common_right_rr',[r for r in selected if r['common_reranking_defined']])]:
        expected=comparison[scope]; groups=defaultdict(lambda:defaultdict(list))
        for r in subset: groups[r['protein_group']][r['sequence_sha256']].append(r)
        values=[]
        for g,queries in sorted(groups.items()):
            a=math.fsum(math.fsum(r[left] for r in panels)/len(panels) for panels in queries.values())/len(queries)
            b=math.fsum(math.fsum(r[right] for r in panels)/len(panels) for panels in queries.values())/len(queries)
            values.append([a,b])
        n=len(values)
        if n!=expected['protein_groups'] or len(subset)!=expected['query_panels']: failures.append('aggregation_denominator')
        if n:
            v=np.array(values); difference=v[:,0]-v[:,1]
            for field,actual in [('left_mrr',math.fsum(v[:,0])/n),('right_mrr',math.fsum(v[:,1])/n),('difference',math.fsum(difference)/n)]:
                if abs(actual-expected[field])>1e-12: failures.append('aggregate:'+field)
        if n>=2:
            draws=rng.integers(n,size=(declared['bootstrap_replicates'],n))
            counts=np.zeros((len(draws),n),dtype=int)
            np.add.at(counts,(np.repeat(np.arange(len(draws)),n),draws.ravel()),1)
            boot=np.sort(counts@difference/n)
            ci=[]
            for probability in [.025,.975]:
                position=probability*(len(boot)-1); lower=int(math.floor(position)); upper=int(math.ceil(position))
                ci.append(float(boot[lower]+(position-lower)*(boot[upper]-boot[lower])))
            if not np.allclose(ci,expected['conditional_protein_bootstrap_95ci'],rtol=0,atol=1e-12): failures.append('bootstrap_interval')
            intervals+=1
        elif expected['conditional_protein_bootstrap_95ci'] is not None: failures.append('undefined_interval')
        summaries+=1
receipt={'status':'PASS' if not failures else 'FAIL','created_utc':now(),'failures':failures,'paired_rows_recomputed':checked,
    'summary_scopes_recomputed':summaries,'conditional_bootstrap_intervals_recomputed':intervals,'maximum_rank_difference':maxdiff,
    'resampling_verifier':'Multiplicity-count weighted means versus indexed group draws; independently sorted candidate ranks',
    'conditional_interval_scope_unchanged':True,'chemical_publication_or_refit_uncertainty_certified':False,
    'verifier_sha256':digest_file(Path(__file__))}
write_json(output,receipt); print(json.dumps(receipt,indent=2))
if failures: raise SystemExit(1)
