"""Fresh two-sided chemical components and publication-aware nested split plans.

No model fitting, scores, historical labels or old feature matrices are inputs.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import json
from pathlib import Path
import sqlite3

from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import MurckoScaffold
from run_raw import stable, digest_file, now, write_json


class Components:
    def __init__(self, values): self.parent = {x:x for x in values}
    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x
    def join(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b: self.parent[max(a,b)] = min(a,b)
    def groups(self):
        buckets = defaultdict(list)
        for value in sorted(self.parent): buckets[self.find(value)].append(value)
        return {stable('|'.join(values)): values for values in buckets.values()}


@lru_cache(maxsize=20000)
def fragment_keys(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError('Invalid core chemistry')
    keys = []
    for fragment in Chem.GetMolFrags(mol, asMols=True):
        if not any(a.GetAtomicNum() == 6 for a in fragment.GetAtoms()): continue
        scaffold = MurckoScaffold.GetScaffoldForMol(fragment)
        key = 'ring:' + Chem.MolToSmiles(scaffold,isomericSmiles=False) if scaffold.GetNumAtoms() else (
            'acyclic_connectivity:' + Chem.MolToSmiles(fragment,isomericSmiles=False))
        keys.append(key)
    return tuple(sorted(set(keys)))


def chemical_components(chemistry):
    keys, key_reactions = {}, defaultdict(set)
    uf = Components(chemistry)
    for rid, row in chemistry.items():
        keys[rid] = sorted({k for side in ['substrates','products'] for s in row[side] for k in fragment_keys(s)})
        if not keys[rid]: raise ValueError('Main transformation has no carbon-containing key')
        for key in keys[rid]: key_reactions[key].add(rid)
    for values in key_reactions.values():
        values = sorted(values)
        for value in values[1:]: uf.join(values[0], value)
    return keys, uf.groups(), key_reactions


def fold_map(groups, allowed, maximum):
    active = {g: sorted(set(members) & allowed) for g,members in groups.items()}
    active = {g: members for g,members in active.items() if members}
    k = min(maximum,len(active))
    if k == 0: return {}, 0
    loads, mapping = [0]*k, {}
    for group,members in sorted(active.items(),key=lambda item:(-len(item[1]),item[0])):
        fold = min(range(k),key=lambda x:(loads[x],x))
        loads[fold] += len(members)
        mapping.update({entity:fold for entity in members})
    return mapping,k


def publication_units(db, edges):
    ids = {p for e in edges for p in e['publication_ids']}
    links = []
    for r in db.execute("SELECT assertion_id,publication_ids,locator,source_sha256 FROM assertions WHERE source='P450Rdb'"):
        pubs = sorted(set(json.loads(r['publication_ids'])))
        pmids = [p for p in pubs if p.startswith('PMID:')]
        dois = [p for p in pubs if p.startswith('DOI:')]
        if len(pmids) == len(dois) == 1:
            ids.update(pmids+dois)
            links.append({'assertion_id':r['assertion_id'],'locator':r['locator'],
                'source_sha256':r['source_sha256'],'PMID':pmids[0],'DOI':dois[0]})
    uf = Components(ids)
    for link in links: uf.join(link['PMID'],link['DOI'])
    groups = uf.groups()
    mapping = {p:g for g,members in groups.items() for p in members}
    conflicts = {g:members for g,members in groups.items() if sum(p.startswith('PMID:') for p in members)>1}
    return mapping,groups,links,conflicts


def make_block(task, pool, edges, unit_sets, protein_folds, chemical_folds, pf, cf,
               scaffold_keys, protein_groups, chemical_groups, candidate_count):
    def is_query(e):
        return ((pf is None or protein_folds[e['sequence_sha256']]==pf)
                and (cf is None or chemical_folds[e['reaction_key']]==cf))
    test = [i for i in pool if is_query(edges[i])]
    # Protein-held-out publication scope includes all edges of those query proteins.
    pub_scope = [i for i in pool if protein_folds[edges[i]['sequence_sha256']]==pf] if pf is not None else test
    test_units = set().union(*(unit_sets[i] for i in pub_scope)) if pub_scope else set()
    before = [i for i in pool if (pf is None or protein_folds[edges[i]['sequence_sha256']]!=pf)
              and (cf is None or chemical_folds[edges[i]['reaction_key']]!=cf)]
    train = [i for i in before if not unit_sets[i] & test_units]
    queries = {edges[i]['sequence_sha256'] for i in test}
    train_sequences = {edges[i]['sequence_sha256'] for i in train}
    train_reactions = {edges[i]['reaction_key'] for i in train}
    test_reactions = {edges[i]['reaction_key'] for i in test}
    train_keys = {k for r in train_reactions for k in scaffold_keys[r]}
    test_keys = {k for r in test_reactions for k in scaffold_keys[r]}
    overlap_units = test_units & {u for i in train for u in unit_sets[i]}
    protein_overlap = ({protein_groups[s] for s in queries} & {protein_groups[s] for s in train_sequences}) if pf is not None else set()
    chemical_overlap = ({chemical_groups[r] for r in test_reactions} & {chemical_groups[r] for r in train_reactions}) if cf is not None else set()
    checks = {'disjoint_edge_indices':not set(test)&set(train), 'no_publication_unit_overlap':not overlap_units,
              'required_protein_groups_disjoint':not protein_overlap,
              'required_chemical_groups_disjoint':not chemical_overlap,
              'required_two_sided_scaffold_keys_disjoint':cf is None or not test_keys & train_keys}
    if not all(checks.values()): raise ValueError('Split exclusion invariant failed')
    reason = 'empty_test_rectangle' if not test else ('publication_or_axis_purge_left_no_training_edges' if not train else None)
    seen = {edges[i]['sequence_sha256'] for i in test if edges[i]['reaction_key'] in train_reactions}
    unseen = {edges[i]['sequence_sha256'] for i in test if edges[i]['reaction_key'] not in train_reactions}
    return {'task':task,'protein_fold':pf,'chemical_fold':cf,'train_edge_indices':train,'test_edge_indices':test,
        'publication_scope_edge_indices':pub_scope,'train_edges_before_publication_purge':len(before),
        'publication_purged_train_edges':len(before)-len(train),'train_edges':len(train),'test_edges':len(test),
        'query_sequences':len(queries),'query_protein_groups':len({protein_groups[s] for s in queries}),
        'training_sequences':len(train_sequences),'training_reactions':len(train_reactions),
        'training_protein_groups':len({protein_groups[s] for s in train_sequences}),
        'test_chemical_groups':len({chemical_groups[r] for r in test_reactions}),
        'training_chemical_groups':len({chemical_groups[r] for r in train_reactions}),
        'queries_with_seen_positive':len(seen),'queries_with_unseen_positive':len(unseen),
        'candidate_count':candidate_count,'checks':checks,'evaluable':reason is None,'unevaluable_reason':reason}


def build(raw_run, dataset, protein_split, output):
    if output.exists(): raise FileExistsError(f'Refusing to overwrite {output}')
    root = Path(__file__).resolve().parent
    audit = json.loads((dataset/'dataset_audit.json').read_text())
    protein_audit = json.loads((protein_split/'search_and_protein_split_audit.json').read_text())
    if audit['status']!='PASS' or protein_audit['status']!='PASS': raise ValueError('Prerequisite audit not passed')
    dbpath = raw_run/'raw_rebuild.sqlite'
    if digest_file(dbpath)!=audit['fresh_raw_sha256']: raise ValueError('Raw hash mismatch')
    chemistry = json.loads((dataset/'core_reactions.json').read_text())
    edges = json.loads((dataset/'core_edges.json').read_text())
    old_plan = json.loads((protein_split/'protein_fold_plan.json').read_text())
    pgroups = json.loads((protein_split/'sequence_groups.json').read_text())['0.4']
    pfolds,pindex = old_plan['protein_folds'],old_plan['protein_groups']
    skeys,cgroups,key_reactions = chemical_components(chemistry)
    cfolds,nc = fold_map(cgroups,set(chemistry),5)
    cindex = {r:g for g,members in cgroups.items() for r in members}
    db = sqlite3.connect(dbpath.resolve().as_uri()+'?mode=ro',uri=True); db.row_factory=sqlite3.Row
    units,ugroups,links,conflicts = publication_units(db,edges); db.close()
    unit_sets = [set(units[p] for p in e['publication_ids']) for e in edges]
    pool=list(range(len(edges))); blocks=[]; npfold=old_plan['folds']
    for task in ['protein_cold','chemical_cold','double_cold']:
        pvalues=range(npfold) if task!='chemical_cold' else [None]
        cvalues=range(nc) if task!='protein_cold' else [None]
        for pf in pvalues:
            for cf in cvalues:
                block=make_block(task,pool,edges,unit_sets,pfolds,cfolds,pf,cf,skeys,pindex,cindex,len(chemistry))
                block['block_id']=f'{task}:p{pf}:c{cf}'
                outertrain=block['train_edge_indices']
                ips,ni=fold_map(pgroups,{edges[i]['sequence_sha256'] for i in outertrain},3)
                ics,nj=fold_map(cgroups,{edges[i]['reaction_key'] for i in outertrain},3)
                inner=[]
                if block['evaluable'] and (task=='chemical_cold' or ni>=2) and (task=='protein_cold' or nj>=2):
                    for ip in range(ni) if task!='chemical_cold' else [None]:
                        for ic in range(nj) if task!='protein_cold' else [None]:
                            ib=make_block(task,outertrain,edges,unit_sets,ips,ics,ip,ic,skeys,pindex,cindex,len(chemistry))
                            ib['block_id']=block['block_id']+f':inner:p{ip}:c{ic}'
                            inner.append(ib)
                block['inner_blocks']=inner
                block['inner_evaluable_blocks']=sum(b['evaluable'] for b in inner)
                block['inner_split_unavailable_reason']=None if inner else 'outer_block_or_group_counts_insufficient'
                blocks.append(block)
    denominators={}
    for task in ['protein_cold','chemical_cold','double_cold']:
        selected=[b for b in blocks if b['task']==task]
        counted=Counter(i for b in selected for i in b['test_edge_indices'])
        denominators[task]={'blocks':len(selected),'evaluable_blocks':sum(b['evaluable'] for b in selected),
            'test_edges':sum(b['test_edges'] for b in selected),'distinct_test_edges':len(counted),
            'each_core_edge_assigned_once':counted==Counter(pool),
            'evaluable_test_edges':sum(b['test_edges'] for b in selected if b['evaluable']),
            'evaluable_query_panel_rows':sum(b['query_sequences'] for b in selected if b['evaluable']),
            'unique_query_sequences':len({edges[i]['sequence_sha256'] for b in selected for i in b['test_edge_indices']}),
            'query_panel_rows_not_independent_replicates':True}
    largest=max(cgroups,key=lambda g:len(cgroups[g]))
    largest_set=set(cgroups[largest])
    bridges=[{'key':key,'reaction_count':len(values),'in_largest_component':bool(values&largest_set)}
             for key,values in sorted(key_reactions.items(),key=lambda item:(-len(item[1]),item[0]))]
    nested_safe=all(set(ib['train_edge_indices']+ib['test_edge_indices'])<=set(b['train_edge_indices'])
                    and not set(ib['train_edge_indices']+ib['test_edge_indices'])&set(b['test_edge_indices'])
                    for b in blocks for ib in b['inner_blocks'])
    checks={'all_core_edges_preserved_each_outer_task':all(v['each_core_edge_assigned_once'] for v in denominators.values()),
        'all_outer_invariants_pass':all(all(b['checks'].values()) for b in blocks),
        'all_inner_invariants_pass':all(all(i['checks'].values()) for b in blocks for i in b['inner_blocks']),
        'inner_records_only_from_outer_training':nested_safe,
        'each_reaction_has_exactly_one_chemical_component':sum(map(len,cgroups.values()))==len(chemistry),
        'protein_folds_unchanged':pfolds==old_plan['protein_folds'],
        'raw_store_hash_unchanged':digest_file(dbpath)==audit['fresh_raw_sha256']}
    output.mkdir(parents=True)
    write_json(output/'chemical_components.json',cgroups)
    write_json(output/'reaction_scaffold_keys.json',skeys)
    write_json(output/'chemical_bridge_profile.json',bridges)
    write_json(output/'publication_units.json',{'mapping':units,'groups':ugroups,'explicit_raw_alias_links':links,
        'multi_PMID_units_for_conservative_purging_not_asserted_same_paper':conflicts,'alias_completeness_certified':False})
    write_json(output/'outer_inner_blocks.json',blocks)
    write_json(output/'axis_assignments.json',{'protein_folds':pfolds,'chemical_folds':cfolds,
        'protein_groups':pindex,'chemical_groups':cindex,'protein_fold_count':npfold,'chemical_fold_count':nc})
    write_json(output/'edge_index.json',[{'index':i,**e,'publication_unit_ids':sorted(unit_sets[i])} for i,e in enumerate(edges)])
    summary={'created_utc':now(),'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
        'sequences':len(pfolds),'core_edges':len(edges),'reaction_labels':len(chemistry),
        'protein_groups':len(pgroups),'chemical_components':len(cgroups),'chemical_folds':nc,
        'largest_chemical_component_reactions':len(largest_set),'largest_chemical_component_id':largest,
        'chemical_component_sizes_descending':sorted(map(len,cgroups.values()),reverse=True),
        'chemical_key_count':len(key_reactions),'raw_single_DOI_single_PMID_alias_links':len(links),
        'unique_explicit_alias_pairs':len({(r['PMID'],r['DOI']) for r in links}),
        'multi_PMID_provenance_units':len(conflicts),'task_denominators':denominators,
        'outer_blocks':len(blocks),'inner_blocks':sum(len(b['inner_blocks']) for b in blocks),
        'unevaluable_outer_blocks':[{'id':b['block_id'],'reason':b['unevaluable_reason']} for b in blocks if not b['evaluable']],
        'unevaluable_inner_blocks':sum(not i['evaluable'] for b in blocks for i in b['inner_blocks']),
        'rdkit_version':rdBase.rdkitVersion,'historical_exposure':'DEVELOPMENT_EXPOSED_OR_UNVERIFIED',
        'model_performance_evaluated':False,'independent_validation':False,
        'scope':'Nested multi-axis development split feasibility; not independent biological validation'}
    write_json(output/'split_audit.json',summary)
    inputs=[dbpath,dataset/'core_edges.json',dataset/'core_reactions.json',dataset/'dataset_audit.json',
        protein_split/'protein_fold_plan.json',protein_split/'sequence_groups.json',
        root/'MULTIAXIS_PROTOCOL_V1.md',root/'run_raw.py',Path(__file__)]
    write_json(output/'input_manifest.json',{str(p.resolve()):digest_file(p) for p in inputs})
    write_json(output/'output_checksums.json',{p.name:digest_file(p) for p in sorted(output.iterdir()) if p.is_file()})
    print(json.dumps({k:v for k,v in summary.items() if k!='chemical_component_sizes_descending'},indent=2))
    if not all(checks.values()): raise RuntimeError('Multi-axis split audit failed')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--raw-run',type=Path,required=True); p.add_argument('--dataset',type=Path,required=True)
    p.add_argument('--protein-split',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); build(a.raw_run,a.dataset,a.protein_split,a.output)
