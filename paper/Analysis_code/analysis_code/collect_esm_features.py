"""Validate all sequence shards and collect fresh global means without labels."""
import json
from pathlib import Path
import hashlib
import numpy as np
from Bio import SeqIO
from run_raw import write_json,digest_file,now
from extract_esm_features import windows

root=Path(__file__).resolve().parent; source=root/'esm_features_01'; output=root/'esm_global_01'
if output.exists(): raise FileExistsError(output)
spec=json.loads((root/'esm_feature_spec_v1.json').read_text()); spec_hash=digest_file(root/'esm_feature_spec_v1.json')
fasta=root/'dataset_02/core_sequences.fasta'; seqs={r.id:str(r.seq) for r in SeqIO.parse(fasta,'fasta')}
declared={}; receipts=[]
for shard in range(4):
    r=json.loads((source/f'shard_{shard:02d}_receipt.json').read_text())
    if r['status']!='COMPLETE' or r['spec_sha256']!=spec_hash or r['fasta_sha256']!=digest_file(fasta): raise ValueError('Incomplete or inconsistent shard')
    receipts.append(r)
    for entry in r['sequences']:
        sid=entry['sequence_sha256']
        if sid in declared: raise ValueError('Duplicate shard sequence')
        declared[sid]=entry
if set(declared)!=set(seqs): raise ValueError('Shard coverage mismatch')
vectors=[]; lengths=[]; windowed=[]; maxdiff=0.
for sid in sorted(seqs):
    file=source/'sequences'/(sid+'.npz'); data=np.load(file,allow_pickle=False)
    if digest_file(file)!=declared[sid]['file_sha256']: raise ValueError('Feature hash mismatch')
    if hashlib.sha256(seqs[sid].encode()).hexdigest()!=sid: raise ValueError('Sequence hash mismatch')
    if str(data['sequence_sha256'])!=sid or str(data['spec_sha256'])!=spec_hash: raise ValueError('Feature identity mismatch')
    span=np.array(list(windows(len(seqs[sid]),spec['maximum_window_residues'],spec['window_overlap_residues'])))
    if not np.array_equal(span,data['window_spans']) or not np.all(data['window_counts']>0): raise ValueError('Window coverage failure')
    residues=data['residue_features']; mean=data['mean_features']
    if residues.shape!=(len(seqs[sid]),spec['dimension']) or mean.shape!=(spec['dimension'],): raise ValueError('Feature shape mismatch')
    if not np.all(np.isfinite(residues)) or not np.all(np.isfinite(mean)): raise ValueError('Nonfinite feature')
    difference=float(np.max(np.abs(residues.astype(np.float32).mean(0)-mean))); maxdiff=max(maxdiff,difference)
    if difference>0.002: raise ValueError('Stored mean differs beyond float16 rounding tolerance')
    vectors.append(mean); lengths.append(len(seqs[sid]))
    if len(span)>1: windowed.append(sid)
output.mkdir()
np.savez_compressed(output/'global_features.npz',sequence_ids=np.array(sorted(seqs)),features=np.array(vectors),lengths=np.array(lengths))
write_json(output/'receipt.json',{'status':'PASS','created_utc':now(),'sequences':len(vectors),'dimensions':spec['dimension'],
    'windowed_sequences':windowed,'full_residue_coverage':True,'maximum_float16_mean_reconciliation_difference':maxdiff,
    'spec_sha256':spec_hash,'fasta_sha256':digest_file(fasta),'features_sha256':digest_file(output/'global_features.npz'),
    'checkpoint_sha256':receipts[0]['download_receipt']['assets']['model']['sha256'],'pretraining_exposure_certified':False,
    'supervised_model_trained':False,'per_sequence_file_hashes':{s:e['file_sha256'] for s,e in declared.items()}})
print(json.dumps({'status':'PASS','sequences':len(vectors),'windowed':len(windowed),'max_mean_difference':maxdiff},indent=2))
