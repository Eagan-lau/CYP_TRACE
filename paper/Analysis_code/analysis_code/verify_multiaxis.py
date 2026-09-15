"""Independent graph and row-set reconciliation, not a call to split construction."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import numpy as np
from scipy.sparse.csgraph import connected_components
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from run_raw import write_json,digest_file,now


def verify(root,output):
    if output.exists(): raise FileExistsError(output)
    def read(name): return json.loads((root/name).read_text())
    chemistry=read('dataset_02/core_reactions.json'); ids=sorted(chemistry); index={r:i for i,r in enumerate(ids)}
    declared=read('multiaxis_split_01/chemical_components.json')
    keys=read('multiaxis_split_01/reaction_scaffold_keys.json'); axes=read('multiaxis_split_01/axis_assignments.json')
    edges=read('multiaxis_split_01/edge_index.json'); blocks=read('multiaxis_split_01/outer_inner_blocks.json')
    recomputed={}
    for rid,c in chemistry.items():
        values=set()
        for side in ['substrates','products']:
            for source in c[side]:
                for mol in Chem.GetMolFrags(Chem.MolFromSmiles(source),asMols=True):
                    if not any(a.GetAtomicNum()==6 for a in mol.GetAtoms()): continue
                    scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=mol,includeChirality=False)
                    values.add('ring:'+scaffold if scaffold else 'acyclic_connectivity:'+Chem.MolToSmiles(mol,isomericSmiles=False))
        recomputed[rid]=sorted(values)
    graph=np.eye(len(ids),dtype=np.uint8); bykey=defaultdict(list)
    for rid,values in recomputed.items():
        for key in values: bykey[key].append(index[rid])
    for values in bykey.values(): graph[np.ix_(values,values)]=1
    n,labels=connected_components(graph,directed=False)
    actual=sorted(tuple(ids[i] for i in np.where(labels==g)[0]) for g in range(n))
    checks={'scaffold_keys_match_alternate_API':recomputed==keys,
        'components_match_SciPy_graph':actual==sorted(tuple(sorted(x)) for x in declared.values()),
        'input_and_output_hashes_valid':all(digest_file(root/'multiaxis_split_01'/p)==h for p,h in read('multiaxis_split_01/output_checksums.json').items())}
    problems=[]
    def check_block(b,pool):
        train=set(b['train_edge_indices']); test=set(b['test_edge_indices']); pubscope=set(b['publication_scope_edge_indices'])
        if not (train|test|pubscope)<=pool or train&test: problems.append(b['block_id']+':record_scope')
        trainpub={p for i in train for p in edges[i]['publication_unit_ids']}
        testpub={p for i in pubscope for p in edges[i]['publication_unit_ids']}
        if trainpub&testpub: problems.append(b['block_id']+':publication_overlap')
        if b['protein_fold'] is not None:
            qg={axes['protein_groups'][edges[i]['sequence_sha256']] for i in test}
            tg={axes['protein_groups'][edges[i]['sequence_sha256']] for i in train}
            if qg&tg: problems.append(b['block_id']+':protein_group_overlap')
        if b['chemical_fold'] is not None:
            qk={k for i in test for k in keys[edges[i]['reaction_key']]}
            tk={k for i in train for k in keys[edges[i]['reaction_key']]}
            if qk&tk: problems.append(b['block_id']+':two_sided_key_overlap')
        if bool(train and test)!=b['evaluable']: problems.append(b['block_id']+':evaluable_flag')
        if len(train)!=b['train_edges'] or len(test)!=b['test_edges']: problems.append(b['block_id']+':counts')
    pool=set(range(len(edges)))
    for b in blocks:
        check_block(b,pool)
        for inner in b['inner_blocks']: check_block(inner,set(b['train_edge_indices']))
    for task in ['protein_cold','chemical_cold','double_cold']:
        checks[task+'_edge_denominator_recomputed']=Counter(i for b in blocks if b['task']==task for i in b['test_edge_indices'])==Counter(pool)
    checks['all_outer_and_inner_record_checks']=not problems
    summary={'status':'PASS' if all(checks.values()) else 'FAIL','created_utc':now(),'checks':checks,'problems':problems,
        'chemical_components':int(n),'outer_blocks':len(blocks),'inner_blocks':sum(len(b['inner_blocks']) for b in blocks),
        'verifier_sha256':digest_file(Path(__file__)),'model_performance_evaluated':False,'independent_biological_validation':False}
    write_json(output,summary); print(json.dumps(summary,indent=2))
    if not all(checks.values()): raise SystemExit(1)


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--output',type=Path,required=True); a=p.parse_args(); verify(a.root,a.output)
