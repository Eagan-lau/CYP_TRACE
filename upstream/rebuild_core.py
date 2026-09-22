"""Reconstruct the complete CYP core from the 27 checksum-identified source files."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
UP = ROOT/'upstream'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def canonical_sha(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def verify_inputs(raw_root):
    result=[]
    for row in json.loads((UP/'core_input_manifest.json').read_text()):
        path=raw_root/row['path']
        # Windows-origin .csv names are accepted on case-sensitive systems.
        if not path.is_file() and path.parent.is_dir():
            matches=[p for p in path.parent.iterdir() if p.name.lower()==path.name.lower()]
            if len(matches)==1:
                path=matches[0]
        state='missing'
        if path.is_file():
            state='match' if path.stat().st_size==row['bytes'] and digest(path)==row['sha256'] else 'hash_mismatch'
        result.append({'input':row['path'],'status':state})
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-root',type=Path,required=True,help='Directory containing p450rdb/, rhea/ and uniprot/.')
    parser.add_argument('--output',type=Path,required=True,help='New directory; existing paths are not overwritten.')
    parser.add_argument('--check-inputs-only',action='store_true')
    args=parser.parse_args()
    checks=verify_inputs(args.raw_root.resolve())
    if not all(r['status']=='match' for r in checks):
        print(json.dumps({'status':'INPUTS_REQUIRED','inputs':checks},indent=2))
        raise SystemExit(2)
    if args.check_inputs_only:
        print(json.dumps({'status':'PASS_INPUTS','files':len(checks)},indent=2)); return
    output=args.output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    code=output/'analysis'
    code.mkdir()
    for name in ('run_raw.py','build_dataset.py','build_reaction_core_v2.py','audit_annotations.py','reaction_roles_v1.py'):
        shutil.copy2(ROOT/'paper/Analysis_code/analysis_code'/name,code/name)
    for name in ('PROTOCOL.md','REACTION_RULES_V1.md'):
        shutil.copy2(UP/'configuration'/name,code/name)
    sys.path.insert(0,str(code))
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    from run_raw import RawBuild
    from build_dataset import build as initial
    from build_reaction_core_v2 import build as core
    run=RawBuild(output,code/'run_03')
    run.raw=args.raw_root.resolve()
    for stage in ('p450rdb','rhea','uniprot_api','swissprot'):
        getattr(run,stage)()
    run.finish()
    initial(code/'run_03',code/'dataset_01')
    core(code/'run_03',code/'dataset_01',code/'dataset_02')
    result=code/'dataset_02'
    edges=json.loads((result/'core_edges.json').read_text())
    chemistry=json.loads((result/'core_reactions.json').read_text())
    semantic=sorted([{k:v for k,v in e.items() if k!='assertion_ids'} for e in edges],key=lambda e:(e['sequence_sha256'],e['reaction_key']))
    observed={'sequences':len({e['sequence_sha256'] for e in edges}),'reactions':len(chemistry),'edges':len(edges),
              'semantic_edges_sha256':canonical_sha(semantic),'reaction_structures_sha256':canonical_sha(chemistry),
              'fasta_lf_sha256':hashlib.sha256((result/'core_sequences.fasta').read_bytes().replace(b'\r\n',b'\n')).hexdigest()}
    expected=json.loads((UP/'expected_core.json').read_text())
    report={'status':'PASS' if observed==expected else 'FAIL','observed':observed,
            'matches':{k:observed[k]==v for k,v in expected.items()},'raw_inputs':checks,
            'environment':{'python':platform.python_version(),'platform':platform.system(),
                           'rdkit':importlib.metadata.version('rdkit'),'biopython':importlib.metadata.version('biopython')}}
    (output/'CORE_ACCEPTANCE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if report['status']!='PASS': raise SystemExit(1)

if __name__=='__main__': main()
