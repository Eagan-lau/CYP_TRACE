"""Qualify sequence-compatible structural availability without activity labels."""
from collections import defaultdict,Counter
import csv
import gzip
import json
from pathlib import Path
import math
from Bio import SeqIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from run_raw import write_json,digest_file,stable,now
from run_blast_baselines import blast_row

root=Path(__file__).resolve().parent; out=root/'structure_mapping_audit_02'
if out.exists(): raise FileExistsError(out)
read=lambda name:json.loads((root/name).read_text())
raw=root.parent/'data/raw/rcsb_cyp_ligand_templates_v1/mmcif'
metadata={}
for r in read('structure_assets_01/source_manifest.json'):
    path=raw/r['file']
    if digest_file(path)!=r['sha256']: raise ValueError('Changed raw structure')
    with gzip.open(path,'rt') as handle: d=MMCIF2Dict(handle)
    values=[]
    for key in ['_refine.ls_d_res_high','_em_3d_reconstruction.resolution']:
        for value in d.get(key,[]):
            try:
                number=float(value)
                if math.isfinite(number) and number>0: values.append(number)
            except ValueError: pass
    metadata[r['file']]={'experimental_methods':d.get('_exptl.method',[]),'resolution_A':min(values) if values else None,
        'initial_deposition_date':d.get('_pdbx_database_status.recvd_initial_deposition_date',[]),
        'source_sha256':r['sha256']}
aliases=read('structure_mapping_01/sequence_to_chain_assets.json')
chains=[a for members in aliases.values() for a in members]
byfile=defaultdict(list)
for a in chains: byfile[a['file'].removesuffix('.cif.gz')].append(a)
fs_alias={}; unmatched=[]
for r in SeqIO.parse(root/'foldseek_raw_01/structure_sequences.fasta','fasta'):
    sid=stable(str(r.seq)); stems=[s for s in byfile if r.id==s or r.id.startswith(s+'_')]
    candidates=[a for stem in stems for a in byfile[stem] if a['sequence_sha256']==sid]
    if '_' in r.id:
        candidates=[a for a in candidates if r.id==a['file'].removesuffix('.cif.gz')+'_'+a['chain']]
    if len(candidates)==1:
        a=candidates[0]; key=(a['file'],a['chain'])
        if key in fs_alias: raise ValueError('Multiple Foldseek identifiers for one parsed chain')
        fs_alias[key]=r.id
    else: unmatched.append({'foldseek_identifier':r.id,'sequence_sha256':sid,'candidate_chains':len(candidates),'sequence_length':len(r.seq)})
query={r.id:str(r.seq) for r in SeqIO.parse(root/'dataset_02/core_sequences.fasta','fasta')}
target={r.id:str(r.seq) for r in SeqIO.parse(root/'structure_mapping_01/resolved_chain_sequences.fasta','fasta')}
aln=root/'structure_blast_01/core_to_resolved_chains.tsv'
expected=(root/'structure_blast_01/output.sha256').read_text().split()[0]
if digest_file(aln)!=expected: raise ValueError('Structure BLAST checksum mismatch')
matches=defaultdict(list); counts=Counter()
with aln.open() as stream:
    for row in csv.reader(stream,delimiter='\t'):
        q,t,identity,qcov,tcov,evalue,bits=blast_row(row); counts['alignment_rows']+=1
        if q not in query or t not in target: raise ValueError('Unexpected alignment identifier')
        for sid,seqs,start,end,aligned,total in [(q,query,int(row[4]),int(row[5]),row[12],int(row[10])),
                                             (t,target,int(row[6]),int(row[7]),row[13],int(row[11]))]:
            if len(seqs[sid])!=total or seqs[sid][start-1:end].upper()!=aligned.replace('-','').upper(): raise ValueError('Alignment cannot reconstruct source')
        if evalue>0.001 or qcov<0.8 or tcov<0.95 or identity<0.95: continue
        qpos=int(row[4])-1; tpos=int(row[6])-1; mapping=[]
        for qa,ta in zip(row[12],row[13]):
            if qa!='-' and ta!='-': mapping.append([qpos,tpos,qa==ta])
            qpos+=qa!='-'; tpos+=ta!='-'
        for a in aliases[t]:
            key=(a['file'],a['chain'])
            if key not in fs_alias: counts['qualified_aliases_unmapped_to_foldseek']+=1; continue
            matches[q].append({'target_sequence_sha256':t,'file':a['file'],'chain':a['chain'],'foldseek_identifier':fs_alias[key],
                'identity':identity,'query_coverage':qcov,'resolved_chain_coverage':tcov,'evalue':evalue,'bits':bits,
                'primary_exact_aligned_identity':identity==1,'complete_sequence_identity':query[q]==target[t],
                'aligned_position_map_0based':mapping,'experimental_metadata':metadata[a['file']],
                'assayed_construct_certified':False})


def representative(candidates):
    if not candidates: return None
    coverage=max((r['query_coverage'],r['resolved_chain_coverage']) for r in candidates)
    tied=[r for r in candidates if (r['query_coverage'],r['resolved_chain_coverage'])==coverage]
    methodsets={tuple(r['experimental_metadata']['experimental_methods']) for r in tied}
    resolutions=[r['experimental_metadata']['resolution_A'] for r in tied]
    if len(methodsets)==1 and all(r is not None for r in resolutions):
        best=min(resolutions); tied=[r for r in tied if r['experimental_metadata']['resolution_A']==best]
    return min(tied,key=lambda r:(r['file'],r['chain']))


selected={}; census=[]
for q in sorted(query):
    candidates=matches[q]; primary=[r for r in candidates if r['primary_exact_aligned_identity']]
    exact=[r for r in candidates if r['complete_sequence_identity']]
    best=representative(primary)
    if best is not None: selected[q]=best
    census.append({'sequence_sha256':q,'primary_available':bool(primary),'primary_chain_count':len(primary),
        'identity95_sensitivity_available':bool(candidates),'complete_sequence_geometry_available':bool(exact),
        'primary_selected_structure':None if best is None else best['foldseek_identifier']})
out.mkdir()
write_json(out/'raw_experimental_metadata.json',metadata)
write_json(out/'foldseek_chain_reconciliation.json',{'mapped_chains':len(fs_alias),'unmatched':unmatched,
    'policy':'Unmatched Foldseek entries are excluded from query/label transfer until adjudicated; no silent identity repair'})
write_json(out/'qualified_chain_associations.json',dict(matches)); write_json(out/'primary_representatives.json',selected)
write_json(out/'availability_census.json',census)
summary={'status':'ASSOCIATIONS_COMPUTED_NOT_MODEL_EVALUATION','created_utc':now(),'counts':dict(counts),'core_sequences':len(query),
    'primary_structure_sequences':len(selected),'identity95_sensitivity_sequences':sum(r['identity95_sensitivity_available'] for r in census),
    'complete_sequence_geometry_sequences':sum(r['complete_sequence_geometry_available'] for r in census),'mapped_foldseek_chains':len(fs_alias),
    'unmatched_foldseek_entries':len(unmatched),'all_blast_aligned_residues_reconstructed':True,
    'selection_uses_reaction_labels_or_scores':False,'independent_biological_validation':False,'structural_functional_performance_evaluated':False}
write_json(out/'audit.json',summary)
names=['audit_structure_mapping_v2.py','STRUCTURE_EXPERT_PROTOCOL_V1.md','structure_assets_01/source_manifest.json',
    'structure_mapping_01/sequence_to_chain_assets.json','structure_mapping_01/resolved_chain_sequences.fasta',
    'dataset_02/core_sequences.fasta','structure_blast_01/core_to_resolved_chains.tsv','foldseek_raw_01/structure_sequences.fasta']
write_json(out/'input_manifest.json',{name:digest_file(root/name) for name in names})
print(json.dumps(summary,indent=2))
