"""Nested MAIN evaluation of train-only chemical backbone and ESM residual."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import time
import numpy as np
from run_raw import write_json,digest_file,stable,now
from run_main_baselines import rank_result,summarize
from conditional_core import normalized_features,fit_residual,predict_delta,select_lambda,model_receipt


def run(root,output):
    if output.exists(): raise FileExistsError(output)
    read=lambda name:json.loads((root/name).read_text())
    if read('esm_global_01/receipt.json')['status']!='PASS': raise ValueError('ESM feature gate failed')
    if read('main_baseline_validation_01/INDEPENDENT_METRIC_QC.json')['status']!='PASS': raise ValueError('MAIN metric gate failed')
    for file,h in read('multiaxis_split_01/output_checksums.json').items():
        if digest_file(root/'multiaxis_split_01'/file)!=h: raise ValueError('Split changed')
    edges=read('dataset_02/core_edges.json'); rids=sorted(read('dataset_02/core_reactions.json'))
    seqs=sorted({e['sequence_sha256'] for e in edges}); si={s:i for i,s in enumerate(seqs)}; ri={r:i for i,r in enumerate(rids)}
    axes=read('multiaxis_split_01/axis_assignments.json'); blocks=read('multiaxis_split_01/outer_inner_blocks.json')
    protein=np.load(root/'esm_global_01/global_features.npz',allow_pickle=False)
    rep=np.load(root/'main_baselines_01/fresh_representations.npz',allow_pickle=False)
    if list(protein['sequence_ids'])!=seqs or list(rep['sequence_ids'])!=seqs or list(rep['reaction_ids'])!=rids: raise ValueError('Representation order mismatch')
    if digest_file(root/'esm_global_01/global_features.npz')!=read('esm_global_01/receipt.json')['features_sha256']: raise ValueError('ESM checksum mismatch')
    P=normalized_features(protein['features']); K=rep['chemical_kernel']
    tie=np.empty(len(rids),dtype=int)
    for rank,i in enumerate(sorted(range(len(rids)),key=lambda x:stable(rids[x]))): tie[i]=rank
    def construct(block):
        Y=np.zeros((len(seqs),len(rids)),dtype=np.float64); truth=defaultdict(set)
        for idx in block['train_edge_indices']:
            e=edges[idx]; Y[si[e['sequence_sha256']],ri[e['reaction_key']]]=1
        for idx in block['test_edge_indices']:
            e=edges[idx]; truth[e['sequence_sha256']].add(ri[e['reaction_key']])
        return Y,truth
    output.mkdir(parents=True); allrows=[]; receipts=[]
    for bn,b in enumerate(blocks):
        if not b['evaluable']: continue
        started=time.time(); innerlog=[]; innerdelta=[]; innertargets=[]; innerrows=[]; fits=[]
        outertrain=set(b['train_edge_indices'])
        for j,inner in enumerate(b['inner_blocks']):
            if not inner['evaluable']: continue
            if not set(inner['train_edge_indices']).issubset(outertrain) or not set(inner['test_edge_indices']).issubset(outertrain): raise ValueError('Inner leakage')
            if set(inner['train_edge_indices'])&set(inner['test_edge_indices']): raise ValueError('Inner target leakage')
            Y,truth=construct(inner); model=fit_residual(P,Y,K); qs=sorted(truth)
            delta=predict_delta(model,P[[si[q] for q in qs]])
            innerlog.extend([model['prior_log']]*len(qs)); innerdelta.extend(delta)
            innertargets.extend([sorted(truth[q]) for q in qs])
            innerrows.extend([{'group':axes['protein_groups'][q],'query':q,'inner_number':j} for q in qs])
            fits.append({'inner_number':j,'training_edges':len(inner['train_edge_indices']),'test_edges':len(inner['test_edge_indices']),**model_receipt(model)})
        logmat=np.array(innerlog); deltamat=np.array(innerdelta)
        choice=select_lambda(logmat,deltamat,innertargets,innerrows) if innerrows else {'lambda':0.,'reason':'no_inner_predictions'}
        np.savez_compressed(output/f'block_{bn:03d}_inner_predictions.npz',log_prior=logmat,delta=deltamat)
        write_json(output/f'block_{bn:03d}_inner_manifest.json',{'rows':innerrows,'targets':innertargets,'fits':fits,'selection':choice})
        Y,truth=construct(b); model=fit_residual(P,Y,K); qs=sorted(truth); qindices=[si[q] for q in qs]
        delta=predict_delta(model,P[qindices]); selected_lambda=choice['lambda']; lp=model['prior_log']
        training=np.flatnonzero(Y.sum(1)); cosine=P[qindices]@P[training].T
        best=training[np.argmax(cosine,axis=1)]
        nearest=Y[best]; transported=(nearest/nearest.sum(1,keepdims=True))@K
        scores={'conditional_chemistry_prior':lp,'esm_conditional_learned':lp[None,:]+selected_lambda*delta,
                'esm_conditional_lambda1_diagnostic':lp[None,:]+delta,'esm_cosine_top1':nearest,'esm_cosine_chemical_transport':transported}
        masks={name:np.ones_like(score,dtype=bool) for name,score in scores.items()}; masks['esm_cosine_top1']=nearest>0
        np.savez_compressed(output/f'block_{bn:03d}_scores.npz',query_ids=np.array(qs),reaction_ids=np.array(rids),
            **scores,**{'domain_'+name:mask for name,mask in masks.items()})
        np.savez_compressed(output/f'block_{bn:03d}_model.npz',training_sequence_ids=np.array([seqs[i] for i in model['training']]),
            center=model['center'],training_features_centered=model['training_features_centered'],kernel_scale=model['kernel_scale'],
            prior_log=lp,coefficients=model['coefficients'],residual_rms=model['residual_rms'],selected_lambda=selected_lambda)
        seen=set(np.flatnonzero(Y.sum(0)))
        for qi,q in enumerate(qs):
            tasks={b['task']+'_all':truth[q]}
            if b['task']=='protein_cold': tasks.update({'protein_cold_seen':truth[q]&seen,'protein_cold_unseen':truth[q]-seen})
            for task,positives in tasks.items():
                if not positives: continue
                for method,values in scores.items():
                    score=values[qi] if values.ndim==2 else values
                    mask=masks[method][qi] if masks[method].ndim==2 else masks[method]
                    allrows.append({'task':task,'block_id':b['block_id'],'block_number':bn,'sequence_sha256':q,
                        'protein_group':axes['protein_groups'][q],'method':method,'target_reaction_indices':sorted(positives),
                        'query_index_in_block':qi,**rank_result(score,mask,positives,tie)})
        receipts.append({'block_id':b['block_id'],'block_number':bn,'selection':choice,'fit':model_receipt(model),
            'inner_fits':len(fits),'inner_rows':len(innerrows),'elapsed_seconds':time.time()-started})
        write_json(output/'progress.json',{'status':'RUNNING','completed_blocks':len(receipts),'total_blocks':len(blocks),'updated_utc':now()})
        write_json(output/'block_receipts.json',receipts)
        print(f"Block {bn+1}/{len(blocks)} lambda={selected_lambda:.6g} inner={len(innerrows)} seconds={time.time()-started:.1f}",flush=True)
    with (output/'per_query_metrics.jsonl').open('w') as stream:
        for row in allrows: stream.write(json.dumps(row,allow_nan=False)+'\n')
    write_json(output/'evaluation.json',{'created_utc':now(),'aggregate':summarize(allrows),'metric_rows':len(allrows),
        'outer_blocks':len(blocks),'lambda_zero_blocks':sum(r['selection']['lambda']==0 for r in receipts),
        'conditional_model_trained':True,'deployment_router_complete':False,'calibrated_biochemical_probability':False,
        'independent_biological_validation':False,'permutation_and_lineage_controls_complete':False,
        'statistical_inference_complete':False,'whole_paper_complete':False,'specification':'CONDITIONAL_MODEL_V1.md'})
    names=['run_conditional_model.py','conditional_core.py','CONDITIONAL_MODEL_V1.md','esm_feature_spec_v1.json',
        'dataset_02/core_edges.json','dataset_02/core_reactions.json','esm_global_01/global_features.npz','esm_global_01/receipt.json',
        'main_baselines_01/fresh_representations.npz','multiaxis_split_01/outer_inner_blocks.json','multiaxis_split_01/axis_assignments.json']
    write_json(output/'input_manifest.json',{name:digest_file(root/name) for name in names})
    write_json(output/'progress.json',{'status':'COMPLETE_ARITHMETIC_QC_PENDING','completed_blocks':len(receipts),'updated_utc':now()})
    print(json.dumps({'status':'COMPLETE_ARITHMETIC_QC_PENDING','metric_rows':len(allrows),'zero_lambda_blocks':sum(r['selection']['lambda']==0 for r in receipts)},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True); a=p.parse_args()
    run(Path(__file__).resolve().parent,a.output)
