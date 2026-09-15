"""Test the real open biological index; this is a retrieval fidelity test."""
from pathlib import Path
import argparse,json,sys,tempfile
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'inference/src'))
from cyptrace_pipeline.evidence import build_evidence_index,load_evidence_index,screen

def main(core):
 with tempfile.TemporaryDirectory(prefix='cyptrace-open-evidence-') as d:
  target=Path(d)/'index.json.gz';build_evidence_index(core,target);index=load_evidence_index(target)
  count=0
  for edge_key,evidence in index['edges'].items():
   seq,key=edge_key.split('|',1)
   result=screen(index,{'query':index['sequences'][seq]},[{'candidate_id':'target','reaction_key':key}])[0]
   assert not result['abstained'] and result['evidence']==evidence
   assert result['score'] is None
   assert all(set(x['sources']) <= {'SwissProt','UniProt_API'} for x in evidence)
   assert all(x['publication_ids'] for x in evidence)
   count+=1
  seq=next(iter(index['sequences'].values()))
  missing=screen(index,{'query':seq},[{'candidate_id':'missing','reaction_key':'MAIN:not_present'}])[0]
  assert missing['abstained'] and missing['score'] is None and missing['evidence']==[]
  novel=screen(index,{'novel':'M'*100},[{'candidate_id':'known','reaction_key':next(iter(index['reactions']))}])[0]
  assert novel['abstained'] and novel['score'] is None and novel['evidence']==[]
  assert count==1232
  return {'status':'PASS','documented_biological_edges_retrieved':count,'unknown_reaction_and_sequence_checks':2,'scope':'Exact evidence fidelity; not prospective prediction validation'}
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--core',type=Path,required=True);a=p.parse_args();print(json.dumps(main(a.core),indent=2))
