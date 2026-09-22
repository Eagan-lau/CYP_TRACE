"""Prepare a new portable analysis workspace without historical local directories."""
import argparse
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[1]

def prepare(destination):
    destination=destination.resolve()
    destination.mkdir(parents=True,exist_ok=False)
    analysis=destination/'analysis'
    analysis.mkdir()
    for folder in (ROOT/'paper/Analysis_code/analysis_code',ROOT/'paper/Analysis_code/protocols',ROOT/'upstream/configuration'):
        for source in folder.iterdir():
            if source.is_file() and source.suffix in {'.py','.md','.json','.sh'}:
                shutil.copy2(source,analysis/source.name)
    # Original raw human inputs are already part of the public release.
    raw=analysis/'raw_additions'
    raw.mkdir()
    shutil.copy2(ROOT/'paper/Source_data/raw_development/molecules-26-04678-s001.zip',raw/'molecules-26-04678-s001.zip')
    external=raw/'figshare_26630515_v4'
    external.mkdir()
    for source in (ROOT/'paper/Source_data/raw_external').glob('*.csv'):
        shutil.copy2(source,external/source.name)
    (analysis/'external_human_cyp_01').mkdir()
    shutil.copy2(ROOT/'upstream/configuration/external_acquisition_manifest.json',analysis/'external_human_cyp_01/acquisition_manifest.json')
    (destination/'data/raw').mkdir(parents=True)
    (destination/'tools').mkdir()
    print(json.dumps({'status':'WORKSPACE_PREPARED','analysis_directory':str(analysis),
                      'external_label_files':len(list(external.glob('*.csv')))},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    prepare(p.parse_args().output)
