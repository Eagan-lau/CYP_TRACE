"""Independently reconstruct all nested residual predictions by linear solves."""
import json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from scipy.special import logsumexp
from run_raw import write_json,digest_file,now

root=Path(__file__).resolve().parent; result=root/'conditional_model_01'
output=root/'conditional_validation_01/INDEPENDENT_FITTING_QC.json'
if output.exists(): raise FileExistsError(output)
read=lambda name:json.loads((root/name).read_text())
edges=read('dataset_02/core_edges.json'); rids=sorted(read('dataset_02/core_reactions.json'))
seqs=sorted({e['sequence_sha256'] for e in edges}); si={s:i for i,s in enumerate(seqs)}; ri={r:i for i,r in enumerate(rids)}
P=np.load(root/'esm_global_01/global_features.npz',allow_pickle=False)['features'].astype(float)
P/=np.maximum(np.sqrt(np.sum(P*P,axis=1,keepdims=True)),1e-12)
K=np.load(root/'main_baselines_01/fresh_representations.npz',allow_pickle=False)['chemical_kernel'].astype(float)
blocks=read('multiaxis_split_01/outer_inner_blocks.json'); receipts=read('conditional_model_01/block_receipts.json')
axes=read('multiaxis_split_01/axis_assignments.json')
failures=[]; max_pred_error=0.; nested_count=0; outer_count=0; lambda_checks=0


def reconstruct(block,declared):
    Y=np.zeros((len(seqs),len(rids)))
    truth=defaultdict(set)
    for idx in block['train_edge_indices']:
        e=edges[idx]; Y[si[e['sequence_sha256']],ri[e['reaction_key']]]=1
    for idx in block['test_edge_indices']:
        e=edges[idx]; truth[e['sequence_sha256']].add(ri[e['reaction_key']])
    train=np.where(Y.sum(1)>0)[0]; target=(Y[train]/Y[train].sum(1,keepdims=True))@K
    mean=target.mean(0); prob=np.maximum(mean/mean.sum(),1e-12); prob/=prob.sum()
    residual=target-mean; rms=np.sqrt(np.mean(residual**2)); center=P[train].mean(0); x=P[train]-center
    gram=x@x.T; scale=np.trace(gram)/len(train); qs=sorted(truth)
    if declared['zero_residual']:
        delta=np.zeros((len(qs),len(rids)))
    else:
        coefficients=solve(gram/scale+declared['alpha']*np.eye(len(train)),residual,assume_a='pos',check_finite=False)
        delta=((P[[si[q] for q in qs]]-center)@x.T/scale)@coefficients/rms
    return qs,truth,np.log(prob),delta


for bn,block in enumerate(blocks):
    if not block['evaluable']: continue
    info=read(f'conditional_model_01/block_{bn:03d}_inner_manifest.json')
    saved=np.load(result/f'block_{bn:03d}_inner_predictions.npz',allow_pickle=False)
    fits={r['inner_number']:r for r in info['fits']}; rows=info['rows']; truths=info['targets']
    grouped=defaultdict(list)
    for i,row in enumerate(rows): grouped[row['inner_number']].append(i)
    for j,indices in grouped.items():
        inner=block['inner_blocks'][j]
        if not set(inner['train_edge_indices']).issubset(block['train_edge_indices']) or not set(inner['test_edge_indices']).issubset(block['train_edge_indices']): failures.append('nested_pool')
        qs,target,lp,delta=reconstruct(inner,fits[j]); nested_count+=1
        if qs!=[rows[i]['query'] for i in indices]: failures.append('inner_query_order')
        for qi,i in enumerate(indices):
            if sorted(target[qs[qi]])!=truths[i] or axes['protein_groups'][qs[qi]]!=rows[i]['group']: failures.append('inner_target_or_group')
        error=float(np.max(np.abs(delta-saved['delta'][indices]))); max_pred_error=max(max_pred_error,error)
        if error>1e-5 or not np.allclose(lp,saved['log_prior'][indices],atol=1e-10,rtol=1e-10): failures.append(f'inner_prediction:{bn}:{j}:{error}')
    choice=info['selection']; value=choice['lambda']; lp=saved['log_prior']; delta=saved['delta']
    panel_counts=Counter((r['group'],r['query']) for r in rows); group_queries=Counter(g for g,q in panel_counts)
    weights=np.array([1/len(group_queries)/group_queries[r['group']]/panel_counts[(r['group'],r['query'])] for r in rows])
    true_lp=np.array([np.mean(lp[i,t]) for i,t in enumerate(truths)]); true_delta=np.array([np.mean(delta[i,t]) for i,t in enumerate(truths)])
    score=lp+value*delta; probs=np.exp(score-logsumexp(score,axis=1,keepdims=True))
    derivative=float(weights@(np.sum(probs*delta,axis=1)-true_delta))
    loss=float(weights@(logsumexp(score,axis=1)-true_lp-value*true_delta)); base=float(weights@(logsumexp(lp,axis=1)-true_lp))
    if choice['reason']=='convex_inner_optimum' and abs(derivative)>1e-6: failures.append('lambda_stationarity:'+str(bn))
    if value<0 or loss>base+1e-8: failures.append('lambda_worse_than_zero:'+str(bn))
    if abs(loss-choice.get('selected_loss',loss))>1e-8: failures.append('lambda_loss:'+str(bn))
    lambda_checks+=1
    declared=next(r['fit'] for r in receipts if r['block_number']==bn)
    qs,target,lp,delta=reconstruct(block,declared); outer_count+=1
    scorefile=np.load(result/f'block_{bn:03d}_scores.npz',allow_pickle=False)
    expected=lp[None,:]+value*delta
    error=float(np.max(np.abs(expected-scorefile['esm_conditional_learned']))); max_pred_error=max(max_pred_error,error)
    if error>1e-5: failures.append('outer_score:'+str(bn))
    print(f'Fitting QC {bn+1}/{len(blocks)}',flush=True)
receipt={'status':'PASS' if not failures else 'FAIL','created_utc':now(),'failures':failures,
    'inner_models_reconstructed':nested_count,'outer_models_reconstructed':outer_count,'lambda_optimality_checks':lambda_checks,
    'maximum_prediction_difference':max_pred_error,'implementation':'Cholesky/positive-definite linear solve versus eigen-decomposition fitting',
    'scope':'Every nested prediction, target/group identity and selected coefficient first-order optimum; GCV scientific suitability and biological validation remain separate',
    'verifier_sha256':digest_file(Path(__file__))}
write_json(output,receipt); print(json.dumps(receipt,indent=2))
if failures: raise SystemExit(1)
