"""Download the version-4 original label files and reject hash/version drift.

Unused supplier fingerprints are available only with --include-unused.
Existing files are verified, never silently overwritten.
"""
from pathlib import Path
import argparse,csv,hashlib,urllib.request
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parent
def verified(data,row):
 return (len(data)==int(row['bytes']) and hashlib.md5(data).hexdigest()==row['official_md5'] and hashlib.sha256(data).hexdigest()==row['local_sha256'])
def download(row,destination):
 name=row['filename'];assert Path(name).name==name
 p=destination/name
 if p.exists():
  if not verified(p.read_bytes(),row):raise ValueError('Existing input differs from frozen source: '+name)
  return name+' VERIFIED_EXISTING'
 request=urllib.request.Request(row['download_url'],headers={'User-Agent':'CYP-TRACE-reproduction/1.0'})
 with urllib.request.urlopen(request,timeout=90) as response:data=response.read()
 if not verified(data,row):raise ValueError('Downloaded input differs from frozen source: '+name)
 with p.open('xb') as h:h.write(data)
 return name+' DOWNLOADED_AND_VERIFIED'
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--include-unused',action='store_true');args=parser.parse_args()
 with (ROOT/'Tables/Table_S1b_External_file_downloads.tsv').open(encoding='utf8') as h:records=list(csv.DictReader(h,delimiter='\t'))
 records=[r for r in records if args.include_unused or r['analytical_use']=='Original labels']
 args.output.mkdir(parents=True,exist_ok=True)
 with ThreadPoolExecutor(max_workers=4) as pool:
  for result in pool.map(lambda r:download(r,args.output),records):print(result,flush=True)
