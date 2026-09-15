"""Matched-domain reranking and conditional paired protein-group uncertainty."""
from collections import defaultdict
import json
from pathlib import Path
import numpy as np
from run_raw import write_json,digest_file,stable,now
from run_main_baselines import random_expectations,rank_result

PAIRS=[('esm_conditional_learned','conditional_chemistry_prior'),
       ('esm_cosine_chemical_transport','conditional_chemistry_prior'),
       ('esm_cosine_chemical_transport','homology_chemical_transport'),
       ('esm_cosine_top1','mmseqs_top1'),('esm_cosine_top1','mmseqs_weighted'),
       ('esm_conditional_learned','mmseqs_weighted'),('esm_cosine_chemical_transport','uniform_expectation')]


def group_summary(rows,left,right,rng):
    grouped=defaultdict(lambda:defaultdict(list))
    for r in rows: grouped[r['protein_group']][r['sequence_sha256']].append(r)
    means=[]; group_values=[]
    for group,queries in sorted(grouped.items()):
        a=np.mean([np.mean([r[left] for r in panels]) for panels in queries.values()])
        b=np.mean([np.mean([r[right] for r in panels]) for panels in queries.values()])
        group_values.append({'protein_group':group,'left_mean':float(a),'right_mean':float(b),'difference':float(a-b)})
        means.append([a,b])
    values=np.array(means); n=len(means)
    if not n: return {'query_panels':0,'unique_queries':0,'protein_groups':0,'left_mrr':None,'right_mrr':None,'difference':None,'conditional_protein_bootstrap_95ci':None,'groups':[]}
    difference=values[:,0]-values[:,1]
    interval=None
    if n>=2:
        boot=difference[rng.integers(n,size=(5000,n))].mean(1)
        interval=list(map(float,np.quantile(boot,[.025,.975])))
    return {'query_panels':len(rows),'unique_queries':len({r['sequence_sha256'] for r in rows}),'protein_groups':n,
        'left_mrr':float(values[:,0].mean()),'right_mrr':float(values[:,1].mean()),'difference':float(difference.mean()),
        'conditional_protein_bootstrap_95ci':interval,'groups':group_values}


def run(root):
    output=root/'paired_comparison_01'
    if output.exists(): raise FileExistsError(output)
    for file in ['main_baseline_validation_01/INDEPENDENT_METRIC_QC.json','conditional_validation_01/INDEPENDENT_METRIC_QC.json',
                 'conditional_validation_01/INDEPENDENT_FITTING_QC.json']:
        if json.loads((root/file).read_text())['status']!='PASS': raise ValueError('Unverified model input')
    folders=['main_baselines_01','conditional_model_01']; records={}; method_folder={}
    for folder in folders:
        for line in (root/folder/'per_query_metrics.jsonl').read_text().splitlines():
            row=json.loads(line); method=row['method']; method_folder[method]=folder
            records[(method,row['task'],row['block_number'],row['sequence_sha256'])]=row
    rids=sorted(json.loads((root/'dataset_02/core_reactions.json').read_text())); n=len(rids)
    tie=np.empty(n,dtype=int)
    for j,i in enumerate(sorted(range(n),key=lambda i:stable(rids[i]))): tie[i]=j
    output.mkdir(); rows=[]
    for bn in range(35):
        files={f:np.load(root/f/f'block_{bn:03d}_scores.npz',allow_pickle=False) for f in folders}
        def values(method,row):
            if method=='uniform_expectation': return None,np.ones(n,dtype=bool)
            data=files[method_folder[method]]; qi=row['query_index_in_block']
            if str(data['query_ids'][qi])!=row['sequence_sha256'] or list(data['reaction_ids'])!=rids: raise ValueError('Saved score order mismatch')
            score=data[method]; score=score[qi] if score.ndim==2 else score
            mask=data['domain_'+method]; mask=mask[qi] if mask.ndim==2 else mask
            return score,mask
        for left,right in PAIRS:
            leftrows=[r for (m,t,b,q),r in records.items() if m==left and b==bn]
            for a in leftrows:
                b=records[(right,a['task'],bn,a['sequence_sha256'])]
                if a['target_reaction_indices']!=b['target_reaction_indices'] or a['candidate_count']!=b['candidate_count'] or a['protein_group']!=b['protein_group']:
                    raise ValueError('Unmatched pair denominator')
                sa,ma=values(left,a); sb,mb=values(right,b); common=ma&mb
                target=[i for i in a['target_reaction_indices'] if common[i]]
                ca=cb=None
                if target:
                    ca=rank_result(sa,common,target,tie)['rr']
                    cb=random_expectations(int(common.sum()),len(target))['rr'] if right=='uniform_expectation' else rank_result(sb,common,target,tie)['rr']
                rows.append({'task':a['task'],'block_number':bn,'sequence_sha256':a['sequence_sha256'],'protein_group':a['protein_group'],
                    'left_method':left,'right_method':right,'left_rr':a['rr'],'right_rr':b['rr'],'full_candidate_count':n,
                    'target_positives':a['positive_count'],'left_covered_positives':a['covered_positives'],'right_covered_positives':b['covered_positives'],
                    'left_candidate_count':int(ma.sum()),'right_candidate_count':int(mb.sum()),'common_candidate_count':int(common.sum()),
                    'common_positive_count':len(target),'common_left_rr':ca,'common_right_rr':cb,'common_reranking_defined':bool(target)})
        for data in files.values(): data.close()
        print(f'Paired comparison block {bn+1}/35',flush=True)
    rng=np.random.default_rng(20260909); summary=[]
    for task in sorted({r['task'] for r in rows}):
        for left,right in PAIRS:
            selected=[r for r in rows if r['task']==task and r['left_method']==left and r['right_method']==right]
            common=[r for r in selected if r['common_reranking_defined']]
            summary.append({'task':task,'left_method':left,'right_method':right,
                'end_to_end':group_summary(selected,'left_rr','right_rr',rng),
                'common_domain_reranking':group_summary(common,'common_left_rr','common_right_rr',rng),
                'coverage':{'original_target_positive_instances':sum(r['target_positives'] for r in selected),
                    'left_covered_positive_instances':sum(r['left_covered_positives'] for r in selected),
                    'right_covered_positive_instances':sum(r['right_covered_positives'] for r in selected),
                    'common_positive_instances':sum(r['common_positive_count'] for r in selected),
                    'common_undefined_query_panels':len(selected)-len(common),
                    'mean_common_candidates_all_query_panels':float(np.mean([r['common_candidate_count'] for r in selected]))}})
    with (output/'paired_query_records.jsonl').open('w') as stream:
        for r in rows: stream.write(json.dumps(r,allow_nan=False)+'\n')
    write_json(output/'comparison_summary.json',{'created_utc':now(),'status':'COMPUTED_VERIFICATION_PENDING','comparisons':summary,
        'comparisons_count':len(summary),'paired_query_rows':len(rows),'bootstrap_seed':20260909,'bootstrap_replicates':5000,
        'interval_scope':'Conditional on fitted models, fixed catalogue, observed evidence and partitions; protein-group resampling only',
        'chemical_and_publication_dependence_addressed':False,'confirmatory_hypothesis_test':False,'independent_validation':False})
    files=['PAIRED_COMPARISON_PROTOCOL_V1.md','paired_expert_comparisons.py','main_baselines_01/per_query_metrics.jsonl',
        'conditional_model_01/per_query_metrics.jsonl','conditional_validation_01/INDEPENDENT_FITTING_QC.json']
    write_json(output/'input_manifest.json',{f:digest_file(root/f) for f in files})


if __name__=='__main__': run(Path(__file__).resolve().parent)
