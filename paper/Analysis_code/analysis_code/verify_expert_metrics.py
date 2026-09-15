"""Independent dictionary/list implementation for saved expert scores and metrics."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import numpy as np
from run_raw import write_json,digest_file,now,stable

FIELDS=['rr','mean_positive_rr','recall_at_1','recall_at_5','recall_at_10','hit_at_1','hit_at_5','hit_at_10','ndcg_at_10']


def verify(root,result,output):
    if output.exists(): raise FileExistsError(output)
    read=lambda p:json.loads(p.read_text())
    ev=read(result/'evaluation.json'); edges=read(root/'dataset_02/core_edges.json')
    blocks=read(root/'multiaxis_split_01/outer_inner_blocks.json'); rids=sorted(read(root/'dataset_02/core_reactions.json'))
    ri={r:i for i,r in enumerate(rids)}; seqs=sorted({e['sequence_sha256'] for e in edges}); si={s:i for i,s in enumerate(seqs)}
    rows=[json.loads(line) for line in (result/'per_query_metrics.jsonl').read_text().splitlines()]
    checks={}; problems=[]; byblock=defaultdict(list); maxdiff=0.; fitchecks=0
    for name,h in read(result/'input_manifest.json').items():
        relative=name.replace('\\','/').split('paper_rebuild_20260909/')[-1]
        checks['input:'+relative]=digest_file(root/relative)==h
    for row in rows: byblock[row['block_number']].append(row)
    H=None
    if (result/'blast_homology.npz').exists():
        H=np.load(result/'blast_homology.npz',allow_pickle=False)['homology_bits'].astype(float)
        K=np.load(root/'main_baselines_01/fresh_representations.npz',allow_pickle=False)['chemical_kernel'].astype(float)
    table=defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
    for bn,members in byblock.items():
        b=blocks[bn]; source=np.load(result/f'block_{bn:03d}_scores.npz',allow_pickle=False)
        if list(source['reaction_ids'])!=rids: problems.append('candidate_order:'+str(bn))
        Y=np.zeros((len(seqs),len(rids))); truth=defaultdict(set)
        for idx in b['train_edge_indices']:
            e=edges[idx]; Y[si[e['sequence_sha256']],ri[e['reaction_key']]]=1
        for idx in b['test_edge_indices']:
            e=edges[idx]; truth[e['sequence_sha256']].add(ri[e['reaction_key']])
        seen=set(np.where(Y.sum(0)>0)[0]); training=np.where(Y.sum(1)>0)[0]
        qs=list(source['query_ids']); rankcache={}
        if H is not None:
            h=H[[si[q] for q in qs]]; mass=h@Y
            balanced=Y/np.maximum(Y.sum(1,keepdims=True),1)
            weighted=(h[:,training]@balanced[training]/np.maximum(h[:,training].sum(1,keepdims=True),1))@K
            nearest=np.zeros_like(mass)
            for qi,q in enumerate(qs):
                candidates=[s for s in training if h[qi,s]>0]
                if candidates: nearest[qi]=Y[min(candidates,key=lambda s:(-h[qi,s],seqs[s]))]
            expected={'blast_top1':nearest,'blast_weighted':mass,'blast_chemical_transport':weighted}
            for method,value in expected.items():
                mask=value>0 if method!='blast_chemical_transport' else np.broadcast_to(h[:,training].sum(1,keepdims=True)>0,value.shape)
                if not np.allclose(source[method],value,rtol=3e-5,atol=1e-5): problems.append(f'fit:{bn}:{method}')
                if not np.array_equal(source['domain_'+method],mask): problems.append(f'domain:{bn}:{method}')
                fitchecks+=len(qs)
        for row in members:
            method=row['method']; qi=row['query_index_in_block']; q=row['sequence_sha256']; key=(method,qi)
            expectedtruth=truth[q]
            if row['task']=='protein_cold_seen': expectedtruth=expectedtruth&seen
            if row['task']=='protein_cold_unseen': expectedtruth=expectedtruth-seen
            if expectedtruth!=set(row['target_reaction_indices']) or qs[qi]!=q: problems.append('truth_or_query:'+str(bn))
            if key not in rankcache:
                score=source[method]; score=score[qi] if score.ndim==2 else score
                mask=source['domain_'+method]; mask=mask[qi] if mask.ndim==2 else mask
                order=sorted([i for i,x in enumerate(mask) if x],key=lambda i:(-float(score[i]),stable(rids[i])))
                rankcache[key]={i:j+1 for j,i in enumerate(order)}
            ranked=rankcache[key]; positions=[ranked.get(i,math.inf) for i in sorted(expectedtruth)]
            count=len(positions); best=min(positions); rr=0. if math.isinf(best) else 1/best
            values={'rr':rr,'mean_positive_rr':math.fsum(1/r for r in positions)/count}
            for k in [1,5,10]:
                hits=sum(r<=k for r in positions); values['recall_at_'+str(k)]=hits/count; values['hit_at_'+str(k)]=float(hits>0)
            values['ndcg_at_10']=math.fsum(1/math.log2(r+1) for r in positions if r<=10)/math.fsum(1/math.log2(r+1) for r in range(1,min(count,10)+1))
            for field,value in values.items():
                diff=abs(value-row[field]); maxdiff=max(maxdiff,diff)
                if diff>1e-10: problems.append(f'metric:{bn}:{method}:{field}')
            if row['candidate_count']!=len(rids) or row['covered_candidates']!=len(ranked) or row['positive_count']!=count or row['covered_positives']!=sum(math.isfinite(r) for r in positions):
                problems.append('coverage:'+str(bn))
            if (row['conditional_rr'] is None)!=math.isinf(best): problems.append('conditional_undefined:'+str(bn))
            if row['conditional_rr'] is not None and abs(row['conditional_rr']-rr)>1e-10: problems.append('conditional_value:'+str(bn))
            table[(row['task'],method)][row['protein_group']][q].append(values)
    for (task,method),groups in table.items():
        for field in FIELDS:
            gmeans=[math.fsum(math.fsum(r[field] for r in panels)/len(panels) for panels in queries.values())/len(queries) for queries in groups.values()]
            actual=math.fsum(gmeans)/len(gmeans)
            if abs(actual-ev['aggregate'][task][method]['protein_group_macro_'+field])>1e-10: problems.append('aggregate:'+task+':'+method)
    checks.update({'all_blocks_present':len(byblock)==sum(b['evaluable'] for b in blocks),'all_rows_counted':len(rows)==ev['metric_rows'],
        'unique_row_keys':len(rows)==len({(r['task'],r['block_id'],r['sequence_sha256'],r['method']) for r in rows}),
        'all_metric_and_fit_checks_pass':not problems})
    receipt={'status':'PASS' if all(checks.values()) else 'FAIL','created_utc':now(),'checks':checks,'problems':problems,
        'metric_rows_recomputed':len(rows),'aggregates_recomputed':len(table),'fitted_score_vectors_recomputed':fitchecks,
        'maximum_metric_difference':maxdiff,'verifier_sha256':digest_file(Path(__file__)),
        'scope':'All saved-score metric arithmetic; BLAST also all training score vectors; no biological independence certification'}
    write_json(output,receipt); print(json.dumps(receipt,indent=2))
    if receipt['status']!='PASS': raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--result',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); verify(Path(__file__).resolve().parent,a.result,a.output)
