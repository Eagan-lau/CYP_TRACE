"""Recalculate all metrics/aggregates without importing the baseline scorer."""
import argparse
from collections import defaultdict
from functools import lru_cache
import json
import math
from pathlib import Path
import numpy as np
from run_raw import write_json,digest_file,now,stable

FIELDS=['rr','mean_positive_rr','recall_at_1','recall_at_5','recall_at_10','hit_at_1','hit_at_5','hit_at_10','ndcg_at_10']


@lru_cache(maxsize=5000)
def uniform(n,k):
    denominator=math.comb(n,k)
    rr=math.fsum(math.comb(n-r,k-1)/denominator/r for r in range(1,n-k+2))
    values={'rr':rr,'mean_positive_rr':math.fsum(1/r for r in range(1,n+1))/n}
    for top in [1,5,10]:
        t=min(top,n)
        values['recall_at_'+str(top)]=t/n
        values['hit_at_'+str(top)]=1-(math.comb(n-t,k)/denominator if n-t>=k else 0)
    values['ndcg_at_10']=(k/n)*math.fsum(1/math.log2(r+1) for r in range(1,min(n,10)+1))/math.fsum(1/math.log2(r+1) for r in range(1,min(k,10)+1))
    return values


def verify(root,output):
    if output.exists(): raise FileExistsError(output)
    result=root/'main_baselines_01'
    read=lambda p:json.loads(p.read_text())
    ev=read(result/'evaluation.json'); blocks=read(root/'multiaxis_split_01/outer_inner_blocks.json')
    edges=read(root/'dataset_02/core_edges.json'); chem=read(root/'dataset_02/core_reactions.json')
    seqs=sorted({e['sequence_sha256'] for e in edges}); rids=sorted(chem)
    si={s:i for i,s in enumerate(seqs)}; ri={r:i for i,r in enumerate(rids)}
    ids_hash=[stable(r) for r in rids]
    checks={}; failures=[]; maxdiff=0.
    for recorded,expected in read(result/'input_manifest.json').items():
        suffix=recorded.replace('\\','/').split('paper_rebuild_20260909/',1)[-1]
        checks['input:'+suffix]=digest_file(root/suffix)==expected
    rows=[json.loads(line) for line in (result/'per_query_metrics.jsonl').read_text().splitlines()]
    byblock=defaultdict(list)
    for row in rows: byblock[row['block_number']].append(row)
    reps=np.load(result/'fresh_representations.npz',allow_pickle=False)
    H=reps['homology_bits']; kernel=reps['chemical_kernel']; recomputed=[]; fitted_score_checks=0
    for number,selected in sorted(byblock.items()):
        block=blocks[number]; source=np.load(result/f'block_{number:03d}_scores.npz',allow_pickle=False)
        if list(source['reaction_ids'])!=rids: failures.append('candidate_order:'+str(number))
        Y=np.zeros((len(seqs),len(rids)),dtype=float)
        for i in block['train_edge_indices']:
            e=edges[i]; Y[si[e['sequence_sha256']],ri[e['reaction_key']]]=1
        frequency=Y.sum(0); balanced=Y/np.maximum(Y.sum(1,keepdims=True),1)
        training=np.where(Y.sum(1)>0)[0]; known=set(np.where(frequency>0)[0]); truths=defaultdict(set)
        for i in block['test_edge_indices']:
            e=edges[i]; truths[e['sequence_sha256']].add(ri[e['reaction_key']])
        expected={'training_frequency':frequency,'chemical_nearest':kernel[list(sorted(known))].max(0),
            'chemistry_prior':frequency@kernel,'balanced_chemistry_prior':balanced.sum(0)@kernel}
        for method,values in expected.items():
            if not np.allclose(values,source[method],rtol=2e-5,atol=1e-5): failures.append(f'fitted:{number}:{method}')
            fitted_score_checks+=1
        qs=list(source['query_ids']); sample=sorted(set([0,len(qs)-1,len(qs)//2]))
        for qi in sample:
            h=H[si[qs[qi]]].astype(float); weights=h[training]
            supported=[s for s in training if h[s]>0]
            best=min(supported,key=lambda s:(-h[s],s)) if supported else None
            nearest=Y[best] if best is not None else np.zeros(len(rids))
            transport=(weights@balanced[training]/max(float(weights.sum()),1))@kernel
            for method,values in [('mmseqs_top1',nearest),('mmseqs_weighted',h@Y),('homology_chemical_transport',transport)]:
                if not np.allclose(values,source[method][qi],rtol=2e-5,atol=1e-5): failures.append(f'fitted:{number}:{qi}:{method}')
                fitted_score_checks+=1
        ranks_cache={}
        for row in selected:
            q=row['sequence_sha256']; method=row['method']; qi=row['query_index_in_block']; n=len(rids)
            positives=truths[q]
            if row['task']=='protein_cold_seen': positives=positives&known
            if row['task']=='protein_cold_unseen': positives=positives-known
            if sorted(positives)!=row['target_reaction_indices']: failures.append('truth:'+str(number)+':'+q)
            if row['candidate_count']!=n: failures.append('denominator:'+str(number))
            if method=='uniform_expectation':
                values=uniform(n,len(positives)); covered=n; covered_positive=len(positives); conditional=values['rr']
            else:
                cache_key=(method,qi)
                if cache_key not in ranks_cache:
                    scores=source[method]; scores=scores[qi] if scores.ndim==2 else scores
                    domain=source['domain_'+method]; domain=domain[qi] if domain.ndim==2 else domain
                    ranked=sorted([i for i,x in enumerate(domain) if x],key=lambda i:(-float(scores[i]),ids_hash[i]))
                    ranks_cache[cache_key]=({i:rank+1 for rank,i in enumerate(ranked)},len(ranked))
                ranks,covered=ranks_cache[cache_key]
                positions=[ranks.get(i,math.inf) for i in positives]; best=min(positions)
                values={'rr':0. if math.isinf(best) else 1/best,'mean_positive_rr':math.fsum(1/r for r in positions)/len(positions)}
                for top in [1,5,10]:
                    hits=sum(r<=top for r in positions)
                    values['recall_at_'+str(top)]=hits/len(positions); values['hit_at_'+str(top)]=float(hits>0)
                ideal=math.fsum(1/math.log2(i+1) for i in range(1,min(10,len(positions))+1))
                values['ndcg_at_10']=math.fsum(1/math.log2(r+1) for r in positions if r<=10)/ideal
                covered_positive=sum(math.isfinite(r) for r in positions); conditional=values['rr'] if covered_positive else None
            for field,value in values.items():
                difference=abs(value-row[field]); maxdiff=max(maxdiff,difference)
                if difference>1e-10: failures.append(f'metric:{number}:{method}:{field}')
            if covered!=row['covered_candidates'] or covered_positive!=row['covered_positives'] or len(positives)!=row['positive_count']:
                failures.append(f'coverage:{number}:{method}')
            if (conditional is None)!=(row['conditional_rr'] is None): failures.append(f'conditional_null:{number}:{method}')
            recomputed.append({**row,**values})
        source.close()
    # Independent aggregation: query-panel means, query means within protein group, then groups.
    table=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for row in recomputed: table[(row['task'],row['method'])][row['protein_group']][row['sequence_sha256']].append(row)
    aggregates=0
    for (task,method),groups in table.items():
        target=ev['aggregate'][task][method]
        for field in FIELDS:
            groupmeans=[]
            for queries in groups.values():
                querymeans=[math.fsum(r[field] for r in members)/len(members) for members in queries.values()]
                groupmeans.append(math.fsum(querymeans)/len(querymeans))
            actual=math.fsum(groupmeans)/len(groupmeans)
            if abs(actual-target['protein_group_macro_'+field])>1e-10: failures.append(f'aggregate:{task}:{method}:{field}')
        aggregates+=1
    checks['all_metric_rows_recomputed']=len(recomputed)==len(rows)==ev['metric_rows']
    checks['all_expected_outer_blocks_present']=len(byblock)==sum(b['evaluable'] for b in blocks)
    checks['all_method_and_task_keys_unique']=len(rows)==len({(r['task'],r['block_id'],r['sequence_sha256'],r['method']) for r in rows})
    checks['no_arithmetic_or_score_reconciliation_failures']=not failures
    summary={'status':'PASS' if all(checks.values()) else 'FAIL','created_utc':now(),'checks':checks,
        'failures':failures,'metric_rows_recomputed':len(recomputed),'method_task_aggregates_recomputed':aggregates,
        'fitted_score_vector_checks':fitted_score_checks,'maximum_metric_difference':maxdiff,
        'verification_scope':'All ranking metrics/aggregates and selected per-query train-only score vectors; not independent biology',
        'verifier_sha256':digest_file(Path(__file__)),'independent_biological_validation':False,'whole_paper_complete':False}
    write_json(output,summary); print(json.dumps(summary,indent=2))
    if not all(checks.values()): raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();verify(a.root,a.output)
