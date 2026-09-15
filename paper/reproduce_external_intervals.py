"""Run the original fixed-score scaffold bootstrap for both S7 cohorts."""
from pathlib import Path
import argparse
import csv
import json
import math
import sys

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'Analysis_code/analysis_code'))
from evaluate_external_human_cyp import bootstrap

def main(output):
    all_rows=[json.loads(line) for line in (ROOT.parent/'data/external_predictions.jsonl').read_text(encoding='utf8').splitlines() if line]
    with (ROOT/'Tables/Table_S7_table_s7_external_transfer.tsv').open(encoding='utf8') as handle:
        expected=list(csv.DictReader(handle,delimiter='\t'))
    results={}
    for record in expected:
        stratum=record['stratum']
        rows=[r for r in all_rows if r['strata'][stratum]]
        assert len(rows)==int(record['labels'])
        result=bootstrap(rows,sorted({r['isoform'] for r in rows}),stratum)
        ci=result['average_precision']['scaffold_bootstrap_95ci']
        assert math.isclose(ci[0],float(record['ci_low']),rel_tol=0,abs_tol=1e-12)
        assert math.isclose(ci[1],float(record['ci_high']),rel_tol=0,abs_tol=1e-12)
        results[stratum]=result
    report=dict(status='PASS',cohorts=results,scope='Original fixed-score bootstrap; no model refitting')
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    main(args.output)
