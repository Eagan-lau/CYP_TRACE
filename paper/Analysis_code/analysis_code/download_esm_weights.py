"""Download official public model assets to an isolated, resumable directory."""
import json
from pathlib import Path
import subprocess
from run_raw import digest_file, write_json, now

root = Path(__file__).resolve().parent
spec = json.loads((root/'esm_feature_spec_v1.json').read_text())
out = root/'raw_additions/pretrained_esm2_650m'
out.mkdir(parents=True, exist_ok=True)
receipt = {}
for kind in ['model', 'regression']:
    url = spec[kind+'_url']; dest = out/url.rsplit('/',1)[-1]
    expected = spec[kind+'_expected_bytes']
    if not dest.exists():
        partial = dest.with_suffix(dest.suffix+'.partial')
        subprocess.run(['curl','--fail','--location','--retry','5','--retry-delay','5',
            '--connect-timeout','30','--continue-at','-','--output',str(partial),url],check=True)
        if partial.stat().st_size != expected: raise ValueError('Downloaded size mismatch')
        partial.rename(dest)
    if dest.stat().st_size != expected: raise ValueError('Existing asset size mismatch')
    receipt[kind]={'url':url,'bytes':expected,'sha256':digest_file(dest)}
write_json(out/'download_receipt.json',{'status':'COMPLETE','created_utc':now(),'assets':receipt,
    'authenticity_scope':'TLS download from the official published model URL; computed checksum is not an independently publisher-signed checksum'})
print(json.dumps(receipt,indent=2),flush=True)
