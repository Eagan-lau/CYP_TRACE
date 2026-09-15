"""Fresh BLAST expert, retaining the frozen catalogue, partitions and metrics."""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import numpy as np
from run_raw import write_json,digest_file,now,stable
from run_main_baselines import train_scores,rank_result,summarize

METHODS={'mmseqs_top1':'blast_top1','mmseqs_weighted':'blast_weighted',
         'homology_chemical_transport':'blast_chemical_transport'}


def blast_row(row):
    if len(row)!=14: raise ValueError('BLAST row needs exactly 14 fields')
    q,t=row[:2]
    ident,length,qs,qe,ts,te,evalue,bits,qlen,tlen=map(float,row[2:12])
    if not all(math.isfinite(v) for v in [ident,length,qs,qe,ts,te,evalue,bits,qlen,tlen]): raise ValueError('Nonfinite BLAST number')
    if not (0<=ident<=100 and 1<=qs<=qe<=qlen and 1<=ts<=te<=tlen and evalue>=0 and bits>=0): raise ValueError('Invalid BLAST value')
    if any(not v.is_integer() for v in [length,qs,qe,ts,te,qlen,tlen]): raise ValueError('Noninteger coordinate')
    if len(row[12])!=length or len(row[13])!=length: raise ValueError('Aligned sequence length mismatch')
    if len(row[12].replace('-',''))!=qe-qs+1 or len(row[13].replace('-',''))!=te-ts+1: raise ValueError('Gapped coordinate mismatch')
    return q,t,ident/100,(qe-qs+1)/qlen,(te-ts+1)/tlen,evalue,bits


def read_blast(file,seqs,spec,raw_sequences=None):
    si={s:i for i,s in enumerate(seqs)}; H=np.zeros((len(seqs),len(seqs)),dtype=np.float32)
    pairs=set(); selfhits=set(); partial_self=[]; rows=0; qualified=0
    with file.open() as stream:
        for row in csv.reader(stream,delimiter='\t'):
            q,t,ident,qcov,tcov,e,bits=blast_row(row)
            if q not in si or t not in si: raise ValueError('Unknown BLAST sequence')
            if (q,t) in pairs: raise ValueError('Duplicate query-target BLAST pair')
            pairs.add((q,t)); rows+=1
            if raw_sequences is not None:
                for sid,start,end,aligned,total in [(q,int(row[4]),int(row[5]),row[12],int(row[10])),
                                                  (t,int(row[6]),int(row[7]),row[13],int(row[11]))]:
                    if len(raw_sequences[sid])!=total or raw_sequences[sid][start-1:end].upper()!=aligned.replace('-','').upper():
                        raise ValueError('BLAST alignment does not reconstruct source sequence')
            if q==t:
                if ident!=1: raise ValueError('Non-identical BLAST self alignment')
                if qcov!=1 or tcov!=1: partial_self.append({'sequence_sha256':q,'qcov':qcov,'tcov':tcov})
                selfhits.add(q)
            if e<=spec['homology_evalue_maximum'] and min(qcov,tcov)>=spec['homology_minimum_bilateral_coverage']:
                H[si[q],si[t]]=bits; qualified+=1
    if selfhits!=set(seqs): raise ValueError('Missing exact self hits')
    return H,{'rows':rows,'qualified_pairs':qualified,'present_self_hits':len(selfhits),
              'complete_self_hits':len(selfhits)-len(partial_self),'partial_self_hits':partial_self,
              'partial_self_policy':'Local HSPs can omit terminal residues with SEG enabled; retain measured coverage, do not impute full coverage',
              'all_aligned_residues_reconstructed':raw_sequences is not None,
              'duplicate_pairs':0,'coverage_definition':'ungapped coordinate span / full sequence length'}


def run(root,output):
    if output.exists(): raise FileExistsError(output)
    read=lambda name:json.loads((root/name).read_text())
    spec=read('main_baseline_spec_v1.json')
    if read('main_baseline_validation_01/INDEPENDENT_METRIC_QC.json')['status']!='PASS': raise ValueError('MAIN gate failed')
    for file,h in read('multiaxis_split_01/output_checksums.json').items():
        if digest_file(root/'multiaxis_split_01'/file)!=h: raise ValueError('Split hash mismatch')
    aln=root/'blast_main_01/all_vs_all.tsv'
    if digest_file(aln)!=(root/'blast_main_01/output_alignment.sha256').read_text().split()[0]: raise ValueError('BLAST transfer hash mismatch')
    if digest_file(root/'dataset_02/core_sequences.fasta')!=(root/'blast_main_01/input_fasta.sha256').read_text().split()[0]: raise ValueError('BLAST input hash mismatch')
    edges=read('dataset_02/core_edges.json'); rids=sorted(read('dataset_02/core_reactions.json'))
    seqs=sorted({e['sequence_sha256'] for e in edges}); si={s:i for i,s in enumerate(seqs)}; ri={r:i for i,r in enumerate(rids)}
    blocks=read('multiaxis_split_01/outer_inner_blocks.json'); axes=read('multiaxis_split_01/axis_assignments.json')
    reps=np.load(root/'main_baselines_01/fresh_representations.npz',allow_pickle=False)
    if list(reps['reaction_ids'])!=rids or list(reps['sequence_ids'])!=seqs: raise ValueError('Feature order mismatch')
    from Bio import SeqIO
    raw_sequences={r.id:str(r.seq) for r in SeqIO.parse(root/'dataset_02/core_sequences.fasta','fasta')}
    kernel=reps['chemical_kernel']; H,audit=read_blast(aln,seqs,spec,raw_sequences)
    output.mkdir(parents=True)
    np.savez_compressed(output/'blast_homology.npz',sequence_ids=np.array(seqs),homology_bits=H)
    tie=np.empty(len(rids),dtype=int)
    for rank,i in enumerate(sorted(range(len(rids)),key=lambda x:stable(rids[x]))): tie[i]=rank
    rows=[]
    for number,b in enumerate(blocks):
        if not b['evaluable']: continue
        Y=np.zeros((len(seqs),len(rids)),dtype=np.float32); truths=defaultdict(set)
        for i in b['train_edge_indices']:
            e=edges[i]; Y[si[e['sequence_sha256']],ri[e['reaction_key']]]=1
        for i in b['test_edge_indices']:
            e=edges[i]; truths[e['sequence_sha256']].add(ri[e['reaction_key']])
        qs=sorted(truths); seen=set(np.flatnonzero(Y.sum(0))); training=np.flatnonzero(Y.sum(1))
        oldscores,oldmasks,best=train_scores(Y,kernel,H,[si[q] for q in qs],training)
        scores={new:oldscores[old] for old,new in METHODS.items()}; masks={new:oldmasks[old] for old,new in METHODS.items()}
        np.savez_compressed(output/f'block_{number:03d}_scores.npz',query_ids=np.array(qs),reaction_ids=np.array(rids),
            **scores,**{'domain_'+name:mask for name,mask in masks.items()})
        for qi,q in enumerate(qs):
            tasks={b['task']+'_all':truths[q]}
            if b['task']=='protein_cold': tasks.update({'protein_cold_seen':truths[q]&seen,'protein_cold_unseen':truths[q]-seen})
            for task,truth in tasks.items():
                if not truth: continue
                for method in scores:
                    rows.append({'task':task,'block_id':b['block_id'],'block_number':number,'sequence_sha256':q,
                        'protein_group':axes['protein_groups'][q],'method':method,'target_reaction_indices':sorted(truth),
                        'query_index_in_block':qi,**rank_result(scores[method][qi],masks[method][qi],truth,tie)})
        print(f'BLAST block {number+1}/{len(blocks)}',flush=True)
    with (output/'per_query_metrics.jsonl').open('w') as stream:
        for row in rows: stream.write(json.dumps(row,allow_nan=False)+'\n')
    write_json(output/'evaluation.json',{'created_utc':now(),'aggregate':summarize(rows),'metric_rows':len(rows),
        'search_audit':audit,'methods':list(METHODS.values()),'specification':'main_baseline_spec_v1.json (identical domain thresholds, separate BLAST bits)',
        'same_candidate_count':len(rids),'same_outer_blocks':len(blocks),'independent_biological_validation':False,
        'arithmetic_qc_pending':True,'whole_paper_complete':False})
    inputs=['run_blast_baselines.py','run_main_baselines.py','main_baseline_spec_v1.json','dataset_02/core_edges.json',
        'dataset_02/core_reactions.json','blast_main_01/all_vs_all.tsv','multiaxis_split_01/outer_inner_blocks.json',
        'multiaxis_split_01/axis_assignments.json','main_baselines_01/fresh_representations.npz']
    write_json(output/'input_manifest.json',{name:digest_file(root/name) for name in inputs})
    print(json.dumps({'status':'COMPLETE_METRIC_QC_PENDING','metric_rows':len(rows),'search_audit':audit},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); run(Path(__file__).resolve().parent,a.output)
