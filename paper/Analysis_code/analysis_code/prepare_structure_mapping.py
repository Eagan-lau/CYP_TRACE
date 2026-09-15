"""Build a deduplicated resolved-chain search reference from fresh raw geometry."""
from collections import defaultdict
import json
from pathlib import Path
from run_raw import stable,write_json,digest_file,now

root=Path(__file__).resolve().parent; output=root/'structure_mapping_01'
if output.exists(): raise FileExistsError(output)
source=root/'structure_assets_01/chain_assets.jsonl'; aliases=defaultdict(list); sequences={}
for line in source.read_text().splitlines():
    row=json.loads(line); sid=row['sequence_sha256']; seq=row['sequence']
    if stable(seq)!=sid: raise ValueError('Chain sequence hash mismatch')
    if sid in sequences and sequences[sid]!=seq: raise ValueError('Conflicting sequence identity')
    sequences[sid]=seq; aliases[sid].append({k:v for k,v in row.items() if k!='sequence'})
output.mkdir()
with (output/'resolved_chain_sequences.fasta').open('w') as stream:
    for sid,seq in sorted(sequences.items()): stream.write('>'+sid+'\n'+seq+'\n')
write_json(output/'sequence_to_chain_assets.json',dict(aliases))
write_json(output/'input_receipt.json',{'status':'REFERENCE_READY_NOT_QUERY_MAPPING','created_utc':now(),
    'chain_records':sum(map(len,aliases.values())),'unique_resolved_sequences':len(sequences),
    'source_sha256':digest_file(source),'fasta_sha256':digest_file(output/'resolved_chain_sequences.fasta'),
    'construct_identity_certified':False,'biological_labels_added':0})
print(json.dumps({'chain_records':sum(map(len,aliases.values())),'unique_sequences':len(sequences)},indent=2))
