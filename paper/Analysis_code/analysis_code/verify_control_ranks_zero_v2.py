"""Replicate 000: independent transport/shuffle scores and all rank metrics."""
from collections import defaultdict
import json
import math
from pathlib import Path
import numpy as np
from run_raw import stable, digest_file, write_json, now

root = Path(__file__).resolve().parent; result = root / 'protein_controls_01/replicate_000'
output = root / 'protein_control_validation_01/REPLICATE_000_RANK_QC_V2.json'
if output.exists(): raise FileExistsError(output)
read = lambda p: json.loads(p.read_text())
if read(root / 'protein_control_validation_01/REPLICATE_000_FITTING_QC.json')['status'] != 'PASS': raise ValueError('Residual-fit gate')
edges = read(root / 'dataset_02/core_edges.json'); rids = sorted(read(root / 'dataset_02/core_reactions.json'))
seqs = sorted({e['sequence_sha256'] for e in edges}); si = {s: i for i, s in enumerate(seqs)}; ri = {r: i for i, r in enumerate(rids)}
P = np.load(root / 'esm_global_01/global_features.npz', allow_pickle=False)['features'].astype(float); P /= np.linalg.norm(P, axis=1, keepdims=True)
K = np.load(root / 'main_baselines_01/fresh_representations.npz', allow_pickle=False)['chemical_kernel'].astype(float)
blocks = read(root / 'multiaxis_split_01/outer_inner_blocks.json'); receipts = read(result / 'block_receipts.json')
metrics = [json.loads(s) for s in (result / 'per_query_metrics.jsonl').read_text().splitlines()]
byblock = defaultdict(list)
for r in metrics: byblock[r['block_number']].append(r)
failures = []; maxdiff = 0.; scorevectors = 0; checked = 0
for bn, b in enumerate(blocks):
    archive = np.load(result / f'block_{bn:03d}_scores.npz', allow_pickle=False)
    d = {name: archive[name] for name in archive.files}
    archive.close()
    original = np.load(root / 'conditional_model_01' / f'block_{bn:03d}_scores.npz', allow_pickle=False)
    truth = defaultdict(set); Y = np.zeros((len(seqs), len(rids)))
    for i in b['train_edge_indices']:
        e = edges[i]; Y[si[e['sequence_sha256']], ri[e['reaction_key']]] = 1
    for i in b['test_edge_indices']:
        e = edges[i]; truth[e['sequence_sha256']].add(ri[e['reaction_key']])
    qs = sorted(truth); train = sorted({si[edges[i]['sequence_sha256']] for i in b['train_edge_indices']}); seen = set(np.where(Y.sum(0) > 0)[0])
    donor = receipts[bn]['outer_assignment']['donor_indices']; qperm = receipts[bn]['query_assignment']
    nearest = []
    for q in qs:
        best = min(range(len(train)), key=lambda j: (-float(P[si[q]] @ P[train[j]]), train[j]))
        labels = np.flatnonzero(Y[donor[best]]); nearest.append(K[labels].mean(0))
    expected = {'train_profile_permuted_esm_transport': np.array(nearest),
                'query_shuffled_residual': original['esm_conditional_learned'][qperm],
                'query_shuffled_esm_transport': original['esm_cosine_chemical_transport'][qperm]}
    for name, value in expected.items():
        diff = float(np.max(np.abs(value - d[name]))); maxdiff = max(maxdiff, diff)
        if diff > 1e-10: failures.append('nonresidual_score:' + name)
        scorevectors += len(qs)
    if list(d['query_ids']) != qs or list(d['reaction_ids']) != rids: failures.append('saved_order')
    cache = {}
    for row in byblock[bn]:
        q = row['sequence_sha256']; qi = row['query_index_in_block']; name = row['method']; positive = truth[q]
        if row['task'] == 'protein_cold_seen': positive = positive & seen
        elif row['task'] == 'protein_cold_unseen': positive = positive - seen
        if sorted(positive) != row['target_reaction_indices']: failures.append('target')
        if (name, qi) not in cache:
            order = sorted(range(len(rids)), key=lambda j: (-float(d[name][qi, j]), stable(rids[j])))
            cache[(name, qi)] = {j: k + 1 for k, j in enumerate(order)}
        positions = [cache[(name, qi)][i] for i in sorted(positive)]
        values = {'rr': 1 / min(positions), 'mean_positive_rr': math.fsum(1 / i for i in positions) / len(positions)}
        for k in [1, 5, 10]:
            hits = sum(i <= k for i in positions); values['hit_at_' + str(k)] = float(hits > 0); values['recall_at_' + str(k)] = hits / len(positions)
        values['ndcg_at_10'] = math.fsum(1 / math.log2(i + 1) for i in positions if i <= 10) / math.fsum(1 / math.log2(i + 1) for i in range(1, min(len(positions), 10) + 1))
        for field, value in values.items():
            if abs(row[field] - value) > 1e-12: failures.append('rank:' + field)
        checked += 1
    original.close()
    print(f'Independent rank QC block {bn + 1}/35', flush=True)
report = {'status': 'PASS' if not failures else 'FAIL', 'created_utc': now(), 'failures': failures,
          'replicate': 0, 'nonresidual_score_vectors_reconstructed': scorevectors, 'rank_metric_rows_recomputed': checked,
          'maximum_nonresidual_score_difference': maxdiff, 'other_98_individual_scores_recomputed_here': False,
          'verifier_sha256': digest_file(Path(__file__))}
write_json(output, report); print(json.dumps(report, indent=2))
if failures: raise SystemExit(1)
