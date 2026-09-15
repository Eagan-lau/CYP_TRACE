"""Reconstruct the UniProt/Rhea subset using frozen qualification decisions.

Requires Python >=3.10 and RDKit. Does not recreate admission screening against
the omitted P450Rdb assertions or claim an independent evaluation dataset.
"""
from pathlib import Path
from collections import defaultdict
import argparse,gzip,hashlib,json,sys
ROOT=Path(__file__).resolve().parent

def read(path): return json.loads(path.read_text(encoding='utf8'))
def records(path):
 with gzip.open(path,'rt',encoding='utf8') as h: return [json.loads(s) for s in h]
def stable(s): return hashlib.sha256(s.encode()).hexdigest()

def rebuild(output):
 from rdkit import Chem, rdBase
 if output.exists(): raise FileExistsError('Choose a new output directory: '+str(output))
 source=ROOT/'Source_data/open_biological_core'
 members=records(source/'qualification.jsonl.gz')
 admitted=[r for r in members if r['eligible']]
 assertions={r['assertion_id']:r for r in records(source/'uniprot_assertions.jsonl.gz')}
 proteins={r['sequence_sha256']:r for r in records(source/'proteins.jsonl.gz')}
 rhea={r['rhea_id']:r for r in records(source/'rhea_reactions.jsonl.gz')}
 projections=read(source/'rhea_main_projections.json')
 sys.path.insert(0,str(ROOT/'Analysis_code/analysis_code'))
 from audit_annotations import exact_reaction
 from build_reaction_core_v2 import experimental_publications, rhea_links
 edges=defaultdict(list);reactions={};checked_rhea={}
 for row in admitted:
  assert row['source'] in ('SwissProt','UniProt_API')
  raw=json.loads(assertions[row['assertion_id']]['raw_json'])
  pubs=experimental_publications(row['source'],raw)
  assert pubs and pubs==row['admissible_publication_ids']
  primary,physio=rhea_links(row['source'],raw)
  assert set(row['directional_rhea_ids']) <= set(primary+physio)
  assert stable(proteins[row['sequence_sha256']]['sequence'])==row['sequence_sha256']
  full_keys=set()
  for rid in row['directional_rhea_ids']:
   if rid not in checked_rhea: checked_rhea[rid]=exact_reaction(rhea[rid]['reaction_smiles'])
   parsed=checked_rhea[rid];key=parsed['reaction_key'];assert key
   full_keys.add(key)
   proj=projections[key]
   # Rejoin removed participants: a MAIN projection must conserve every
   # participant of its source Rhea reaction, including multiplicity.
   for side,plural in [('substrate','substrates'),('product','products')]:
    restored=proj[plural]+[p['smiles'] for p in proj['removed_participants'] if p['side']==side]
    assert sorted(restored)==sorted(parsed['sides'][side]),(rid,side)
    assert all(Chem.MolFromSmiles(s) is not None for s in proj[plural])
   identity={'substrate':sorted(proj['substrates']),'product':sorted(proj['products'])}
   assert proj['reaction_key']=='MAIN:'+stable(json.dumps(identity,sort_keys=True))
   assert proj['reaction_key']==row['main_reaction_keys'][0]
   reactions[proj['reaction_key']]={k:proj[k] for k in ('substrates','products','single_pair')}
  assert full_keys==set(row['full_reaction_keys'])
  edges[(row['sequence_sha256'],row['main_reaction_keys'][0])].append(row)
 expected=read(source/'expected_core.json')
 assert stable(json.dumps(sorted(edges),separators=(',',':')))==expected['edge_identity_sha256']
 assembled=[]
 for (seq,key),rows in sorted(edges.items()):
  union=lambda key:sorted({x for r in rows for x in r[key]})
  assembled.append({'sequence_sha256':seq,'reaction_key':key,
    'sources':sorted({r['source'] for r in rows}),'source_lineages':['UniProt'],
    'publication_ids':union('admissible_publication_ids'),'assertion_ids':sorted(r['assertion_id'] for r in rows),
    'accessions':sorted({r['accession'] for r in rows}),'taxids':sorted({r['taxid'] for r in rows}),
    'full_reaction_keys':union('full_reaction_keys'),
    'source_scope':'open_UniProt_subset_of_development_core','historical_exposure':'DEVELOPMENT_EXPOSED_OR_UNVERIFIED'})
 assert (len(assembled),len(proteins),len(reactions))==(expected['edges'],expected['sequences'],expected['reactions'])
 output.mkdir(parents=True)
 for name,obj in [('core_edges.json',assembled),('core_reactions.json',reactions)]:
  (output/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf8')
 (output/'core_sequences.fasta').write_text(''.join('>'+seq+'\n'+proteins[seq]['sequence']+'\n' for seq in sorted(proteins)),encoding='utf8')
 report={'status':'PASS','edges':len(edges),'sequences':len(proteins),'reactions':len(reactions),'reparsed_Rhea_reactions':len(checked_rhea),
 'verified_admitted_assertions':len(admitted),'rdkit':rdBase.rdkitVersion,'qualification':'frozen; not re-screened against mixed-source inputs'}
 (output/'REBUILD_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
 return report

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
 print(json.dumps(rebuild(args.output.resolve()),indent=2))
