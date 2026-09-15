"""Fresh, full-residue ESM2 features; resumable independent sequence shards."""
import argparse
import hashlib
import json
from pathlib import Path
import os
import time
import numpy as np
from run_raw import write_json, digest_file, now


def windows(length, maximum=1022, overlap=256):
    if length <= 0 or not 0 <= overlap < maximum: raise ValueError('Invalid window parameters')
    start = 0
    while True:
        end = min(start+maximum, length)
        yield start,end
        if end == length: return
        start = end-overlap


def main(shard, shards):
    import torch
    import esm
    if not 0 <= shard < shards: raise ValueError('Invalid shard')
    root = Path(__file__).resolve().parent
    specfile = root/'esm_feature_spec_v1.json'; spec = json.loads(specfile.read_text())
    fasta = root/'dataset_02/core_sequences.fasta'
    from Bio import SeqIO
    seqs = {r.id:str(r.seq) for r in SeqIO.parse(fasta,'fasta')}
    for sid,seq in seqs.items():
        if hashlib.sha256(seq.encode()).hexdigest()!=sid: raise ValueError('Sequence hash mismatch')
    weightdir=root/'raw_additions/pretrained_esm2_650m'
    download=json.loads((weightdir/'download_receipt.json').read_text())
    checkpoint=weightdir/(spec['model']+'.pt')
    if digest_file(checkpoint)!=download['assets']['model']['sha256']: raise ValueError('Model hash mismatch')
    out=root/'esm_features_01'; (out/'sequences').mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(int(os.environ.get('SLURM_CPUS_PER_TASK','8')))
    torch.set_num_interop_threads(1); torch.manual_seed(0)
    model,alphabet=esm.pretrained.load_model_and_alphabet_local(checkpoint)
    model.eval(); converter=alphabet.get_batch_converter()
    receipts=[]
    selected=sorted(seqs)[shard::shards]
    for n,sid in enumerate(selected):
        dest=out/'sequences'/(sid+'.npz'); seq=seqs[sid]
        if dest.exists():
            old=np.load(dest,allow_pickle=False)
            if str(old['sequence_sha256'])!=sid or str(old['spec_sha256'])!=digest_file(specfile):
                raise ValueError('Existing feature provenance mismatch')
            if old['residue_features'].shape != (len(seq),spec['dimension']): raise ValueError('Existing feature shape mismatch')
            receipts.append({'sequence_sha256':sid,'length':len(seq),'file_sha256':digest_file(dest),'resumed':True})
            continue
        started=time.time(); sums=np.zeros((len(seq),spec['dimension']),dtype=np.float32)
        counts=np.zeros(len(seq),dtype=np.int32); spans=list(windows(len(seq),spec['maximum_window_residues'],spec['window_overlap_residues']))
        with torch.inference_mode():
            for start,end in spans:
                _,_,tokens=converter([(sid,seq[start:end])])
                result=model(tokens,repr_layers=[spec['representation_layer']],return_contacts=False)
                emb=result['representations'][spec['representation_layer']][0,1:end-start+1].cpu().numpy()
                if emb.shape != (end-start,spec['dimension']): raise ValueError('Token count mismatch')
                sums[start:end]+=emb; counts[start:end]+=1
        if not np.all(counts>0): raise ValueError('Truncated sequence')
        residues=sums/counts[:,None]
        if not np.all(np.isfinite(residues)): raise ValueError('Invalid ESM feature')
        with dest.with_suffix('.npz.partial').open('wb') as stream:
            np.savez_compressed(stream,sequence_sha256=np.array(sid),spec_sha256=np.array(digest_file(specfile)),
                checkpoint_sha256=np.array(download['assets']['model']['sha256']),length=np.array(len(seq)),
                residue_features=residues.astype(np.float16),mean_features=residues.mean(0).astype(np.float32),
                window_spans=np.array(spans),window_counts=counts)
        dest.with_suffix('.npz.partial').rename(dest)
        receipts.append({'sequence_sha256':sid,'length':len(seq),'windows':len(spans),
            'elapsed_seconds':time.time()-started,'file_sha256':digest_file(dest),'resumed':False})
        print(f'Shard {shard}: {n+1}/{len(selected)} residues={len(seq)} windows={len(spans)} seconds={time.time()-started:.1f}',flush=True)
        write_json(out/f'shard_{shard:02d}_progress.json',{'status':'RUNNING','completed':len(receipts),'expected':len(selected),'updated_utc':now()})
    write_json(out/f'shard_{shard:02d}_receipt.json',{'status':'COMPLETE','created_utc':now(),'shard':shard,'shards':shards,
        'torch_version':torch.__version__,'esm_version':esm.__version__,'fasta_sha256':digest_file(fasta),
        'spec_sha256':digest_file(specfile),'download_receipt':download,'sequences':receipts,'supervised_model_trained':False})


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--shard',type=int,required=True); p.add_argument('--shards',type=int,default=4)
    a=p.parse_args(); main(a.shard,a.shards)
