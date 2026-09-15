"""Check the extracted source-data archive using Python's standard library.

With --software-root, also execute the separately supplied saved-score checks.
This checks the stated reproduction subset, not general-CYP model refitting.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile

def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def within(root, name):
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError('Path escapes archive root: ' + name)
    return path

def require(condition, message):
    if not condition:
        raise ValueError(message)

def check(root, software=None):
    manifest = rows(root/'FILE_MANIFEST.tsv')
    require(len({r['path'] for r in manifest}) == len(manifest), 'Duplicate manifest path')
    for record in manifest:
        path=within(root,record['path'])
        require(path.is_file(), 'Missing file: '+record['path'])
        data=path.read_bytes()
        require(len(data)==int(record['bytes']), 'Byte count: '+record['path'])
        require(hashlib.sha256(data).hexdigest()==record['sha256'], 'Checksum: '+record['path'])
    index=rows(root/'Tables/Table_index.tsv')
    for item in index:
        actual=rows(within(root,item['file']))
        require(len(actual)==int(item['rows']), 'Table row count: '+item['file'])
    figure_map=rows(root/'Tables/Figure_to_data_map.tsv')
    for item in figure_map:
        if item['saved_input'].startswith('Source_data/'):
            require(within(root,item['saved_input']).is_file(), 'Figure source: '+item['saved_input'])
        for name in item['tables'].split(';'):
            require(within(root,'Tables/'+name.strip()).is_file(), 'Figure table: '+name)
    external=rows(root/'Tables/Table_S7_table_s7_external_transfer.tsv')
    for item in external:
        delta=float(item['isoform_specific_macro_ap'])-float(item['pooled_macro_ap'])
        require(math.isclose(delta,float(item['ap_difference']),abs_tol=1e-12), 'External AP difference')
    strict=next(r for r in external if r['stratum']=='scaffold_and_lineage_eligible')
    require((int(strict['labels']),int(strict['compounds']),int(strict['scaffolds']))==(3035,1001,844), 'Strict endpoint denominators')
    release=rows(root/'Tables/Table_S13a_selective_release_tradeoffs.tsv')
    for method in ('pooled_chemical_knn','isoform_specific_knn'):
        selected={float(r['target_coverage']):r for r in release if r['method']==method}
        require(set(selected)=={0.1,0.25,0.5,0.75,1.0}, 'Release fractions')
        require([int(selected[f]['selected_rows']) for f in (0.1,0.25,0.5,0.75,1.0)]==[306,761,1519,2279,3035], 'Release denominators')
    logistic=rows(root/'Tables/Table_S14b_logistic_regression_supplement.tsv')
    macro=next(r for r in logistic if r['isoform']=='MACRO')
    individual=[float(r['strict_external_AP']) for r in logistic if r['isoform']!='MACRO']
    require(len(individual)==6 and math.isclose(sum(individual)/6,float(macro['strict_external_AP']),abs_tol=1e-12), 'Logistic macro AP')
    report={'source_checks':'PASS','manifest_files':len(manifest),'supplementary_tables':len(index),
            'figure_map_rows':len(figure_map),'strict_external_macro_AP':float(strict['isoform_specific_macro_ap']),
            'logistic_macro_AP':float(macro['strict_external_AP']),
            'upstream_general_model_refitting':'NOT_RUN_REQUIRES_HISTORICAL_MIXED_SOURCE_INPUTS',
            'additional_available_checks':['reproduce_general_metrics.py','rebuild_open_core.py','test_open_evidence.py'],
            'software_checks':'NOT_REQUESTED'}
    if software:
        software=software.resolve()
        require(software.is_dir(), 'Software root is not a directory')
        commands=[
            ['verify_package.py'],
            ['-m','unittest','discover','-s','evaluation','-p','test_evaluator.py','-v'],
            ['evaluation/evaluate_candidates.py','evaluation/figure2_candidates.tsv.gz'],
            ['evaluation/evaluate_candidates.py','evaluation/figure2_candidates.tsv.gz','--ties','average'],
            ['evaluation/check_external.py'],
        ]
        with tempfile.TemporaryDirectory(prefix='cyptrace-reviewer-') as scratch:
            for i,cmd in enumerate(commands):
                if cmd[0]=='evaluation/evaluate_candidates.py':
                    cmd += ['--output',str(Path(scratch)/f'figure2_{i}.json')]
                proc=subprocess.run([sys.executable,'-X','utf8',*cmd],cwd=software,capture_output=True,text=True,encoding='utf8',timeout=240)
                require(proc.returncode==0, 'Software check failed: '+' '.join(cmd)+'\n'+proc.stdout+'\n'+proc.stderr)
                print('PASS '+' '.join(cmd[:4]),file=sys.stderr)
        report['software_checks']='PASS_5_COMMANDS'
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--software-root',type=Path)
    args=parser.parse_args()
    try:
        result=check(args.root.resolve(),args.software_root)
    except (OSError,ValueError,KeyError,StopIteration,subprocess.TimeoutExpired) as exc:
        print('FAIL: '+str(exc),file=sys.stderr)
        return 1
    print(json.dumps(result,indent=2))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
