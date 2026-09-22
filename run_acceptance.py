"""Offline acceptance of the declared public reproduction scope.

Install requirements-reproduction.txt and ./inference first. This runner uses
the current Python interpreter, never a historical workspace or cluster.
Generated outputs must go to a new directory; package inputs are unchanged.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main(output):
    output=output.resolve()
    output.mkdir(parents=True,exist_ok=False)
    # Prove that tuning can run with a data root that has no external labels,
    # saved predictions or selected-parameter table.
    isolated=output/'development_only'
    (isolated/'data').mkdir(parents=True)
    shutil.copy2(ROOT/'data/development_labels.json',isolated/'data/development_labels.json')
    commands=[
        ('payload_integrity',['verify_package.py']),
        ('source_tables',['paper/reviewer_checks.py']),
        ('portable_upstream_units',['-m','unittest','discover','-s','upstream','-p','test_portability.py','-v']),
        ('human_raw_rebuild',['upstream/rebuild_human_sources.py','--output',str(output/'human_raw')]),
        ('evaluator_units',['-m','unittest','discover','-s','evaluation','-p','test_evaluator.py','-v']),
        ('figure2_fixed',['evaluation/evaluate_candidates.py','evaluation/figure2_candidates.tsv.gz','--output',str(output/'figure2_fixed.json')]),
        ('figure2_ties',['evaluation/evaluate_candidates.py','evaluation/figure2_candidates.tsv.gz','--ties','average','--output',str(output/'figure2_average.json')]),
        ('external_AP_and_release',['evaluation/check_external.py']),
        ('general_metrics',['paper/reproduce_general_metrics.py','--output',str(output/'general_metrics.json')]),
        ('external_intervals',['paper/reproduce_external_intervals.py','--output',str(output/'external_intervals.json')]),
        ('open_core_rebuild',['paper/rebuild_open_core.py','--output',str(output/'open_core')]),
        ('open_evidence_lookup',['paper/test_open_evidence.py','--core',str(output/'open_core')]),
        ('inference_doctor',['-m','cyptrace_pipeline','doctor']),
        ('all_external_knn_scores',['test_inference.py']),
        ('human_asset_rebuild',['rebuild_human_bundle.py','--output',str(output/'rebuilt_human_model.json.gz')]),
        ('logistic_units',['-m','unittest','discover','-s','paper','-p','test_logistic_reproduction.py','-v']),
        ('logistic_development_selection',['paper/reproduce_logistic.py','tune','--root',str(isolated),'--output',str(output/'selection')]),
        ('logistic_external_evaluation',['paper/reproduce_logistic.py','evaluate','--selection',str(output/'selection/logistic_selection.json'),'--output',str(output/'logistic')]),
        ('payload_integrity_after',['verify_package.py']),
    ]
    env=dict(os.environ)
    env.pop('PYTHONPATH',None)
    env['PYTHONNOUSERSITE']='1'
    env['PYTHONDONTWRITEBYTECODE']='1'
    env['PYTHONIOENCODING']='utf-8'
    report=dict(status='RUNNING',python=platform.python_version(),platform=platform.platform(),
                versions={p:importlib.metadata.version(p) for p in ('numpy','rdkit','scikit-learn','scipy','joblib','threadpoolctl')},
                payload_manifest_sha256=sha(ROOT/'SHA256SUMS.txt'),
                runner_sha256=sha(Path(__file__)),commands=[],
                scope='Public reconstruction and fixed-prediction checks; not full mixed-source upstream refitting',
                not_run=['Complete mixed-source core (separately tested; see upstream/core_acceptance.json)',
                         'General-CYP protein searches/features/model refitting',
                         'Original nested kNN selection',
                         'Final author-edited figure layout recreation'],
                logistic_tuning_data_root_contains_external_data=False)
    for name,args in commands:
        start=time.monotonic()
        try:
            result=subprocess.run([sys.executable,'-X','utf8',*args],cwd=ROOT,env=env,
                                  capture_output=True,text=True,encoding='utf8',timeout=900)
            code=result.returncode
            content=result.stdout+'\n'+result.stderr
        except subprocess.TimeoutExpired as exc:
            code=124
            content='Timed out after 900 seconds\n'+str(exc)
        # Keep downloadable receipts portable and omit local username/path text.
        content=content.replace(str(output),'<OUTPUT>').replace(str(ROOT),'<PACKAGE>')
        log=output/(name+'.log')
        log.write_text(content,encoding='utf8')
        clean_args=[a.replace(str(output),'<OUTPUT>').replace(str(ROOT),'<PACKAGE>') for a in args]
        report['commands'].append(dict(name=name,command=['python',*clean_args],exit_code=code,
                                       seconds=round(time.monotonic()-start,3),log=log.name,log_sha256=sha(log)))
        print(('PASS ' if code==0 else 'FAIL ')+name,flush=True)
        if code:
            print(content[-3000:],flush=True)
            report['status']='FAIL'
            break
    if len(report['commands'])==len(commands) and all(r['exit_code']==0 for r in report['commands']):
        report['status']='PASS_DECLARED_SCOPE'
    report['outputs']={p.relative_to(output).as_posix():sha(p) for p in sorted(output.rglob('*'))
                       if p.is_file() and p.name!='ACCEPTANCE.json' and 'development_only' not in p.parts}
    (output/'ACCEPTANCE.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(report['status'],flush=True)
    return 0 if report['status']=='PASS_DECLARED_SCOPE' else 1

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    raise SystemExit(main(parser.parse_args().output))
