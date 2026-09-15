"""Recalculate Figures 1, 3 and 4 estimates from assertion/panel-level records.

These are fixed-prediction checks, not model refitting or raw-source acquisition.
"""
from pathlib import Path
from collections import defaultdict,Counter
import argparse,gzip,hashlib,json,math
import numpy as np
ROOT=Path(__file__).resolve().parent/'Source_data/general_metrics'
def read(path): return json.loads(path.read_text())
def rows(path):
 with gzip.open(path,'rt',encoding='utf8') as h:return [json.loads(s) for s in h]
def close(a,b):
 if not math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-12): raise AssertionError((a,b))
def group_pairs(records,left,right):
 groups=defaultdict(lambda:defaultdict(list))
 for r in records:groups[r['protein_group']][r['sequence_sha256']].append([r[left],r[right]])
 return np.array([np.mean([np.mean(v,axis=0) for v in groups[g].values()],axis=0) for g in sorted(groups)])
def main():
 report={}
 records=rows(ROOT/'identity_assertion_audit.jsonl.gz')
 admitted=[r for r in records if r['eligible']]
 lineage=defaultdict(set);full=defaultdict(set)
 for r in admitted:
  lineage[r['source_lineage']].add((r['sequence_sha256'],r['main_reaction_keys'][0]))
  for key in r['full_reaction_keys']:full[r['source_lineage']].add((r['sequence_sha256'],key))
 assert len(lineage['P450Rdb'])==1869 and len(lineage['UniProt'])==1232
 assert len(lineage['P450Rdb']|lineage['UniProt'])==2304
 assert len(lineage['P450Rdb']&lineage['UniProt'])==797
 assert len(full['P450Rdb']&full['UniProt'])==31
 report['figure1']={'qualified_edges':2304,'shared_MAIN_edges':797,'shared_full_component_edges':31,'assertions_audited':len(records)}
 data=rows(ROOT/'clean_pruning_01/outer_paired_records.jsonl.gz')
 expected=read(ROOT/'clean_pruning_01/comparison_summary.json')
 count=0
 for item in expected['comparisons']:
  cell=item['cell'];base=[r for r in data if r['cell']==cell]
  for key,suffix,predicate in [('all','all',lambda r:True),('clean_exact_exposed','exposed',lambda r:r['clean_exact_exposed']),('clean_exact_unexposed','unexposed',lambda r:not r['clean_exact_exposed'])]:
   chosen=[r for r in base if predicate(r)]; exp=item[key]
   assert len(chosen)==exp['query_panels']
   if not chosen: continue
   array=group_pairs(chosen,'selected_rr','none_rr');delta=array[:,0]-array[:,1]
   close(array[:,0].mean(),exp['selected_mrr']);close(array[:,1].mean(),exp['unpruned_mrr']);close(delta.mean(),exp['difference'])
   assert len(array)==exp['protein_groups']
   if len(array)>=2:
    seed=(20260914+int.from_bytes(hashlib.sha256((cell+'|'+suffix).encode()).digest()[:8],'little'))%(2**63-1)
    rng=np.random.default_rng(seed);draws=rng.integers(len(array),size=(5000,len(array)))
    ci=np.quantile(delta[draws].mean(axis=1),[.025,.975])
    for a,b in zip(ci,exp['conditional_protein_bootstrap_95ci']):close(a,b)
   count+=1
 report['figure3_pruning']={'strata_recomputed':count,'panel_rows':len(data),'bootstrap_replicates':5000}
 for folder,label in [('interaction_paired_sequence_01','figure3_interaction'),('structure_paired_02','figure3_matched_structure')]:
  data=rows(ROOT/folder/'paired_query_records.jsonl.gz');expected=read(ROOT/folder/'comparison_summary.json')
  count=0
  for item in expected['comparisons']:
   chosen=[r for r in data if all(r[k]==item[k] for k in ('task','left_method','right_method'))]
   exp=item['end_to_end']; assert len(chosen)==exp['query_panels']
   if not chosen:continue
   array=group_pairs(chosen,'left_rr','right_rr')
   close(array[:,0].mean(),exp['left_mrr']);close(array[:,1].mean(),exp['right_mrr']);close((array[:,0]-array[:,1]).mean(),exp['difference'])
   assert len(array)==exp['protein_groups'];count+=1
  report[label]={'comparisons_recomputed':count,'panel_rows':len(data),'intervals':'stored intervals not resampled in this subcheck'}
 data=rows(ROOT/'reverse_baselines_01/per_query_metrics.jsonl.gz');expected=read(ROOT/'reverse_baselines_01/method_summary.json')
 for exp in expected['summaries']:
  chosen=[r for r in data if r['cell']==exp['cell'] and r['method']==exp['method']]
  by_taxid=defaultdict(list)
  for r in chosen:by_taxid[r['taxid']].append(r['rr'])
  close(np.mean([np.mean(v) for v in by_taxid.values()]),exp['taxid_macro_mrr'])
  close(sum(r['covered_positives'] for r in chosen)/sum(r['positive_count'] for r in chosen),exp['documented_positive_coverage'])
  close(np.mean([r['covered_candidates']/r['candidate_count'] for r in chosen]),exp['mean_candidate_coverage'])
  assert len(chosen)==exp['query_panels'] and len(by_taxid)==exp['taxids']
 report['figure4']={'method_regime_estimates':len(expected['summaries']),'panel_metric_rows':len(data)}
 report['status']='PASS';return report
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args()
 report=main();text=json.dumps(report,indent=2)
 if args.output:args.output.write_text(text+'\n',encoding='utf8')
 print(text)
