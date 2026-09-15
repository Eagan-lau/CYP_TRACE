"""Same-catalogue, domain-masked baselines on frozen MAIN outer blocks."""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
try:
    from rdkit import Chem, DataStructs, rdBase
    from rdkit.Chem import rdFingerprintGenerator
    _rdkit_import_error = None
except Exception as exc:  # pragma: no cover
    Chem = DataStructs = rdBase = rdFingerprintGenerator = None
    _rdkit_import_error = exc
from run_raw import stable,write_json,digest_file,now

METRICS=['rr','mean_positive_rr','recall_at_1','recall_at_5','recall_at_10',
         'hit_at_1','hit_at_5','hit_at_10','ndcg_at_10']


@lru_cache(maxsize=10000)
def random_expectations(n,k):
    if not (0<k<=n): raise ValueError('Invalid candidate/positive counts')
    probability=k/n; rr=0.
    for rank in range(1,n-k+2):
        rr+=probability/rank
        if rank<n-k+1: probability*= (n-rank-k+1)/(n-rank)
    result={'rr':rr,'mean_positive_rr':sum(1/i for i in range(1,n+1))/n,
        'positive_count':k,'covered_positives':k,'candidate_count':n,'covered_candidates':n,
        'conditional_rr':rr}
    for cutoff in [1,5,10]:
        top=min(cutoff,n); none=1.
        for j in range(top): none*=max(0,n-k-j)/(n-j)
        result['recall_at_'+str(cutoff)]=top/n
        result['hit_at_'+str(cutoff)]=1-none
    discounts=[1/math.log2(i+1) for i in range(1,min(n,10)+1)]
    result['ndcg_at_10']=(k/n)*sum(discounts)/sum(discounts[:min(k,10)])
    return result


def rank_result(scores,mask,truth,tie):
    truth=np.asarray(sorted(set(truth)),dtype=int)
    if not len(truth): raise ValueError('No positives for this retrieval query')
    if not np.all(np.isfinite(scores[mask])): raise ValueError('Non-finite defined expert score')
    selected=np.flatnonzero(mask)
    order=selected[np.lexsort((tie[selected],-scores[selected]))]
    ranks=np.full(len(scores),np.inf); ranks[order]=np.arange(1,len(order)+1)
    tr=ranks[truth]; best=float(np.min(tr)); rr=0. if not math.isfinite(best) else 1/best
    result={'rr':rr,'mean_positive_rr':float(np.mean(1/tr)),'positive_count':len(tr),
        'covered_positives':int(np.sum(mask[truth])),'candidate_count':len(scores),
        'covered_candidates':len(order),'conditional_rr':rr if math.isfinite(best) else None}
    for k in [1,5,10]:
        result['recall_at_'+str(k)]=float(np.mean(tr<=k)); result['hit_at_'+str(k)]=float(np.any(tr<=k))
    observed=tr[tr<=10]
    dcg=float(np.sum(1/np.log2(observed+1)))
    ideal=sum(1/math.log2(r+1) for r in range(1,min(len(tr),10)+1))
    result['ndcg_at_10']=dcg/ideal
    return result


def chemical_kernel(chemistry,rids,spec):
    if Chem is None:
        raise ModuleNotFoundError('rdkit is required for chemical kernel construction') from _rdkit_import_error
    width=spec['fingerprint']['bits_per_side']
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=spec['fingerprint']['radius'],
        fpSize=width,includeChirality=spec['fingerprint']['include_chirality'])
    matrix=np.zeros((len(rids),width*2),dtype=np.float32)
    for i,r in enumerate(rids):
        for j,side in enumerate(['substrates','products']):
            mol=Chem.MolFromSmiles('.'.join(chemistry[r][side]))
            if mol is None: raise ValueError('Invalid MAIN chemistry')
            bits=gen.GetFingerprint(mol); a=np.zeros(width,dtype=np.uint8)
            DataStructs.ConvertToNumpyArray(bits,a); matrix[i,j*width:(j+1)*width]=a
    intersection=matrix@matrix.T
    mass=matrix.sum(axis=1); union=mass[:,None]+mass[None,:]-intersection
    return np.divide(intersection,union,out=np.zeros_like(intersection),where=union>0),matrix


def homology_matrix(file,seqs,spec):
    index={s:i for i,s in enumerate(seqs)}; bits=np.zeros((len(seqs),len(seqs)),dtype=np.float32)
    n=0
    with file.open() as stream:
        for r in csv.reader(stream,delimiter='\t'):
            if len(r)!=8 or not all(r): raise ValueError('Invalid MMseqs row')
            q,t=r[:2]; identity,qcov,tcov,length,evalue,score=map(float,r[2:])
            if q not in index or t not in index: raise ValueError('Unknown MMseqs sequence')
            if not all(math.isfinite(v) for v in [identity,qcov,tcov,length,evalue,score]): raise ValueError('Invalid MMseqs number')
            if not (0<=identity<=1 and 0<=qcov<=1 and 0<=tcov<=1): raise ValueError('MMseqs fractions invalid')
            if evalue<=spec['homology_evalue_maximum'] and min(qcov,tcov)>=spec['homology_minimum_bilateral_coverage']:
                bits[index[q],index[t]]=max(bits[index[q],index[t]],score)
            n+=1
    return bits,n


def train_scores(Y,kernel,H,queries,training_sequence_indices):
    n=Y.shape[1]; frequency=Y.sum(axis=0)
    balanced=Y/np.maximum(Y.sum(axis=1,keepdims=True),1)
    Hq=H[np.asarray(queries)]
    mass=Hq@Y
    nearest=np.zeros_like(mass); best_proteins=[]
    for j in range(len(queries)):
        eligible=[int(s) for s in training_sequence_indices if Hq[j,s]>0]
        best=min(eligible,key=lambda s:(-float(Hq[j,s]),s)) if eligible else None
        if best is not None: nearest[j]=Y[best]
        best_proteins.append(best)
    relevant=Hq[:,training_sequence_indices]
    norm=np.maximum(relevant.sum(axis=1,keepdims=True),1.)
    transport=(relevant@balanced[training_sequence_indices]/norm)@kernel
    trainlabels=np.flatnonzero(frequency)
    methods={
        'hash_tie_control':np.ones(n,dtype=np.float32),
        'training_frequency':frequency,
        'chemical_nearest':kernel[trainlabels].max(axis=0),
        'chemistry_prior':frequency@kernel,
        'balanced_chemistry_prior':balanced.sum(axis=0)@kernel,
        'mmseqs_top1':nearest,'mmseqs_weighted':mass,'homology_chemical_transport':transport}
    masks={name:np.ones_like(score,dtype=bool) for name,score in methods.items()}
    masks['mmseqs_top1']=nearest>0; masks['mmseqs_weighted']=mass>0
    masks['homology_chemical_transport']=np.broadcast_to(relevant.sum(axis=1,keepdims=True)>0,transport.shape).copy()
    return methods,masks,best_proteins


def summarize(rows):
    result={}
    for task in sorted({r['task'] for r in rows}):
        result[task]={}
        for method in sorted({r['method'] for r in rows}):
            selected=[r for r in rows if r['task']==task and r['method']==method]
            if not selected: continue
            # A protein repeated in chemical panels contributes once after its panel mean.
            queries=defaultdict(list)
            for r in selected: queries[r['sequence_sha256']].append(r)
            group_values=defaultdict(list)
            for q,members in queries.items():
                group_values[members[0]['protein_group']].append({m:float(np.mean([r[m] for r in members])) for m in METRICS})
            groupmean={g:{m:float(np.mean([r[m] for r in members])) for m in METRICS} for g,members in group_values.items()}
            result[task][method]={'query_panel_rows':len(selected),'unique_queries':len(queries),'protein_groups':len(group_values),
                **{'protein_group_macro_'+m:float(np.mean([r[m] for r in groupmean.values()])) for m in METRICS},
                'micro_positive_coverage':sum(r['covered_positives'] for r in selected)/sum(r['positive_count'] for r in selected),
                'query_panel_mean_candidate_coverage':float(np.mean([r['covered_candidates']/r['candidate_count'] for r in selected])),
                'uncertainty_not_estimated':True}
    return result


def run(dataset,similarity,splits,output):
    if output.exists(): raise FileExistsError(output)
    root=Path(__file__).resolve().parent
    specpath=root/'main_baseline_spec_v1.json'; spec=json.loads(specpath.read_text())
    audit=json.loads((splits/'split_audit.json').read_text())
    if audit['status']!='PASS' or not all(audit['checks'].values()): raise ValueError('Split gate not passed')
    if json.loads((root/'multiaxis_validation_01/INDEPENDENT_SPLIT_QC.json').read_text())['status']!='PASS': raise ValueError('Independent split QC not passed')
    for file,expected in json.loads((splits/'output_checksums.json').read_text()).items():
        if digest_file(splits/file)!=expected: raise ValueError('Frozen split changed: '+file)
    chemistry=json.loads((dataset/'core_reactions.json').read_text()); rids=sorted(chemistry)
    edges=json.loads((dataset/'core_edges.json').read_text()); seqs=sorted({e['sequence_sha256'] for e in edges})
    si={s:i for i,s in enumerate(seqs)}; ri={r:i for i,r in enumerate(rids)}
    blocks=json.loads((splits/'outer_inner_blocks.json').read_text())
    axes=json.loads((splits/'axis_assignments.json').read_text())
    output.mkdir(parents=True)
    write_json(output/'status.json',{'status':'RUNNING','started_utc':now(),'model_selection_performed':False})
    print('Building fresh MAIN chemical and homology representations',flush=True)
    kernel,fingerprints=chemical_kernel(chemistry,rids,spec)
    H,nalign=homology_matrix(similarity/'all_vs_all.tsv',seqs,spec)
    np.savez_compressed(output/'fresh_representations.npz',reaction_ids=np.array(rids),sequence_ids=np.array(seqs),
                        chemical_kernel=kernel,chemical_fingerprints=fingerprints,homology_bits=H)
    tie=np.empty(len(rids),dtype=int)
    for rank,i in enumerate(sorted(range(len(rids)),key=lambda x:stable(rids[x]))): tie[i]=rank
    rows=[]; receipts=[]
    for bno,b in enumerate(blocks):
        if not b['evaluable']:
            receipts.append({'block_id':b['block_id'],'status':'NOT_EVALUABLE','reason':b['unevaluable_reason']}); continue
        Y=np.zeros((len(seqs),len(rids)),dtype=np.float32)
        for i in b['train_edge_indices']:
            e=edges[i]; Y[si[e['sequence_sha256']],ri[e['reaction_key']]]=1
        queries=sorted({edges[i]['sequence_sha256'] for i in b['test_edge_indices']})
        trainseq=np.flatnonzero(Y.sum(axis=1)); trainreactions=set(np.flatnonzero(Y.sum(axis=0)))
        scores,masks,best=train_scores(Y,kernel,H,[si[q] for q in queries],trainseq)
        target=defaultdict(set)
        for i in b['test_edge_indices']: target[edges[i]['sequence_sha256']].add(ri[edges[i]['reaction_key']])
        np.savez_compressed(output/f'block_{bno:03d}_scores.npz',query_ids=np.array(queries),reaction_ids=np.array(rids),
            **scores,**{'domain_'+m:mask for m,mask in masks.items()})
        for qi,q in enumerate(queries):
            truth=target[q]
            tasks={b['task']+'_all':truth}
            if b['task']=='protein_cold':
                tasks.update({'protein_cold_seen':truth&trainreactions,'protein_cold_unseen':truth-trainreactions})
            for task,positives in tasks.items():
                if not positives: continue
                for method in spec['methods']:
                    if method=='uniform_expectation': metric=dict(random_expectations(len(rids),len(positives)))
                    else:
                        score=scores[method][qi] if scores[method].ndim==2 else scores[method]
                        mask=masks[method][qi] if masks[method].ndim==2 else masks[method]
                        metric=rank_result(score,mask,positives,tie)
                    rows.append({'task':task,'block_id':b['block_id'],'block_number':bno,'sequence_sha256':q,
                        'protein_group':axes['protein_groups'][q],'method':method,
                        'target_reaction_indices':sorted(positives),'query_index_in_block':qi,**metric})
        receipts.append({'block_id':b['block_id'],'block_number':bno,'status':'COMPLETE','train_edges':len(b['train_edge_indices']),
            'test_edges':len(b['test_edge_indices']),'query_sequences':len(queries),
            'best_homology_training_sequence_ids':[seqs[i] if i is not None else None for i in best]})
        write_json(output/'status.json',{'status':'RUNNING','completed_outer_blocks':len(receipts),'total_outer_blocks':len(blocks),'updated_utc':now()})
        print(f"Block {bno+1}/{len(blocks)} {b['block_id']}: {len(queries)} queries, {len(b['train_edge_indices'])} training edges",flush=True)
    with (output/'per_query_metrics.jsonl').open('w') as stream:
        for row in rows: stream.write(json.dumps(row,allow_nan=False)+'\n')
    write_json(output/'block_receipts.json',receipts)
    evaluation={'created_utc':now(),'aggregate':summarize(rows),'metric_rows':len(rows),'alignment_rows':nalign,
        'sequences':len(seqs),'reaction_labels':len(rids),'core_edges':len(edges),'outer_blocks':len(blocks),
        'rdkit_version':'unavailable' if rdBase is None else rdBase.rdkitVersion,
        'numpy_version':np.__version__,'specification':spec,
        'conditional_model_trained':False,'independent_validation':False,'whole_paper_complete':False,
        'statistical_scope':'Descriptive development baselines; dependence-aware effect inference and independent evidence pending'}
    write_json(output/'evaluation.json',evaluation)
    inputs=[dataset/'core_edges.json',dataset/'core_reactions.json',similarity/'all_vs_all.tsv',
        splits/'outer_inner_blocks.json',splits/'axis_assignments.json',splits/'output_checksums.json',specpath,Path(__file__),root/'run_raw.py']
    write_json(output/'input_manifest.json',{str(p.resolve()):digest_file(p) for p in inputs})
    write_json(output/'status.json',{'status':'MAIN_BASELINES_COMPLETE_ARITHMETIC_QC_PENDING','completed_utc':now(),
        'outer_blocks':len(receipts),'conditional_model_trained':False,'independent_validation':False,'whole_paper_complete':False})
    print(json.dumps({'status':'MAIN_BASELINES_COMPLETE_ARITHMETIC_QC_PENDING','metric_rows':len(rows)},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--dataset',required=True,type=Path)
    p.add_argument('--similarity',required=True,type=Path); p.add_argument('--splits',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path); a=p.parse_args(); run(a.dataset,a.similarity,a.splits,a.output)
