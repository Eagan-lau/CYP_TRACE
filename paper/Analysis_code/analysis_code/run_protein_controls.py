"""Nested train-profile permutations and fixed-model query-identity ablations."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import time
import numpy as np
from conditional_core import normalized_features, fit_residual, predict_delta, select_lambda, model_receipt
from protein_control_core import permute_profiles
from run_main_baselines import rank_result, summarize
from run_raw import write_json, digest_file, stable, now


def run(root, replicate):
    if not 0 <= replicate < 99: raise ValueError('Frozen replicate range is 0..98')
    out = root / 'protein_controls_01' / f'replicate_{replicate:03d}'
    if out.exists(): raise FileExistsError(out)
    read = lambda n: json.loads((root / n).read_text())
    if read('conditional_validation_01/INDEPENDENT_FITTING_QC.json')['status'] != 'PASS': raise ValueError('Fit QC gate')
    for n, h in read('multiaxis_split_01/output_checksums.json').items():
        if digest_file(root / 'multiaxis_split_01' / n) != h: raise ValueError('Changed frozen partition')
    edges = read('dataset_02/core_edges.json'); rids = sorted(read('dataset_02/core_reactions.json'))
    seqs = sorted({e['sequence_sha256'] for e in edges}); si = {s: i for i, s in enumerate(seqs)}; ri = {r: i for i, r in enumerate(rids)}
    blocks = read('multiaxis_split_01/outer_inner_blocks.json'); axes = read('multiaxis_split_01/axis_assignments.json')
    protein = np.load(root / 'esm_global_01/global_features.npz', allow_pickle=False)
    rep = np.load(root / 'main_baselines_01/fresh_representations.npz', allow_pickle=False)
    if list(protein['sequence_ids']) != seqs or list(rep['sequence_ids']) != seqs or list(rep['reaction_ids']) != rids: raise ValueError('ID order')
    P = normalized_features(protein['features']); K = rep['chemical_kernel']; n = len(rids)
    tie = np.empty(n, dtype=int)
    for j, i in enumerate(sorted(range(n), key=lambda i: stable(rids[i]))): tie[i] = j
    def construct(block):
        Y = np.zeros((len(seqs), n)); truth = defaultdict(set)
        for i in block['train_edge_indices']:
            e = edges[i]; Y[si[e['sequence_sha256']], ri[e['reaction_key']]] = 1
        for i in block['test_edge_indices']:
            e = edges[i]; truth[e['sequence_sha256']].add(ri[e['reaction_key']])
        return Y, truth
    def fit_with_control(block, bn, stage):
        Y, truth = construct(block)
        changed, permutation = permute_profiles(Y, np.random.default_rng(np.random.SeedSequence([20260909, replicate, bn, stage])))
        model = fit_residual(P, changed, K)
        active = Y[Y.sum(1) > 0]; prior = (active / active.sum(1, keepdims=True)) @ K.astype(float)
        prior = prior.mean(0); prior = np.maximum(prior / prior.sum(), 1e-12); prior /= prior.sum()
        difference = float(np.max(np.abs(np.log(prior) - model['prior_log'])))
        if difference > 1e-10: raise ValueError('Chemical backbone changed')
        return Y, changed, truth, model, {**permutation, 'maximum_chemical_prior_log_difference': difference}
    out.mkdir(parents=True); rows = []; receipts = []; started = time.time()
    for bn, b in enumerate(blocks):
        innerlog = []; innerdelta = []; targets = []; records = []; fits = []
        for j, inner in enumerate(b['inner_blocks']):
            if not inner['evaluable']: continue
            Y, changed, truth, model, assignment = fit_with_control(inner, bn, j + 1)
            qs = sorted(truth); delta = predict_delta(model, P[[si[q] for q in qs]])
            innerlog.extend([model['prior_log']] * len(qs)); innerdelta.extend(delta); targets.extend([sorted(truth[q]) for q in qs])
            records.extend([{'group': axes['protein_groups'][q], 'query': q, 'inner_number': j} for q in qs])
            fits.append({'inner_number': j, 'assignment': assignment, **model_receipt(model)})
        logmat = np.array(innerlog); deltamat = np.array(innerdelta)
        choice = select_lambda(logmat, deltamat, targets, records) if records else {'lambda': 0., 'reason': 'no_inner_predictions'}
        if replicate == 0:
            np.savez_compressed(out / f'block_{bn:03d}_inner_predictions.npz', log_prior=logmat, delta=deltamat)
            write_json(out / f'block_{bn:03d}_inner_manifest.json', {'rows': records, 'targets': targets, 'selection': choice})
        Y, changed, truth, model, assignment = fit_with_control(b, bn, 10000)
        qs = sorted(truth); qidx = [si[q] for q in qs]; delta = predict_delta(model, P[qidx])
        training = np.flatnonzero(changed.sum(1)); best = training[np.argmax(P[qidx] @ P[training].T, axis=1)]
        near = changed[best]
        actual = np.load(root / 'conditional_model_01' / f'block_{bn:03d}_scores.npz', allow_pickle=False)
        if list(actual['query_ids']) != qs: raise ValueError('Original query order changed')
        query_assignment = np.random.default_rng(np.random.SeedSequence([20260909, replicate, bn, 20000])).permutation(len(qs))
        scores = {'train_profile_permuted_residual': model['prior_log'][None, :] + choice['lambda'] * delta,
                  'train_profile_permuted_esm_transport': (near / near.sum(1, keepdims=True)) @ K,
                  'query_shuffled_residual': actual['esm_conditional_learned'][query_assignment],
                  'query_shuffled_esm_transport': actual['esm_cosine_chemical_transport'][query_assignment]}
        if replicate == 0:
            np.savez_compressed(out / f'block_{bn:03d}_scores.npz', query_ids=np.array(qs), reaction_ids=np.array(rids), **scores)
        seen = set(np.flatnonzero(Y.sum(0))); mask = np.ones(n, dtype=bool)
        for qi, q in enumerate(qs):
            tasks = {b['task'] + '_all': truth[q]}
            if b['task'] == 'protein_cold': tasks.update({'protein_cold_seen': truth[q] & seen, 'protein_cold_unseen': truth[q] - seen})
            for task, positives in tasks.items():
                if not positives: continue
                for name, values in scores.items():
                    rows.append({'replicate': replicate, 'task': task, 'block_number': bn, 'block_id': b['block_id'], 'sequence_sha256': q,
                                 'protein_group': axes['protein_groups'][q], 'method': name, 'target_reaction_indices': sorted(positives),
                                 'query_index_in_block': qi, **rank_result(values[qi], mask, positives, tie)})
        receipts.append({'block_number': bn, 'inner_fits': fits, 'selection': choice, 'outer_fit': model_receipt(model),
                         'outer_assignment': assignment, 'query_assignment': query_assignment.tolist(),
                         'query_fixed_points': int(np.sum(query_assignment == np.arange(len(qs)))), 'query_sequences': len(qs)})
        actual.close()
        write_json(out / 'progress.json', {'status': 'RUNNING', 'replicate': replicate, 'completed_blocks': len(receipts), 'updated_utc': now()})
        print(f'replicate={replicate} block={bn + 1}/35 lambda={choice["lambda"]:.5g}', flush=True)
    with (out / 'per_query_metrics.jsonl').open('w') as h:
        for r in rows: h.write(json.dumps(r, allow_nan=False) + '\n')
    write_json(out / 'block_receipts.json', receipts)
    write_json(out / 'evaluation.json', {'status': 'COMPUTED_QC_PENDING', 'created_utc': now(), 'replicate': replicate, 'metric_rows': len(rows),
               'aggregate': summarize(rows), 'inner_refits': sum(len(r['inner_fits']) for r in receipts), 'outer_refits': len(receipts),
               'lambda_zero_blocks': sum(r['selection']['lambda'] == 0 for r in receipts), 'elapsed_seconds': time.time() - started,
               'formal_exchangeability_or_pvalue_claim': False, 'independent_biological_validation': False})
    names = ['run_protein_controls.py', 'protein_control_core.py', 'conditional_core.py', 'PROTEIN_INFORMATION_CONTROLS_V1.md',
             'dataset_02/core_edges.json', 'dataset_02/core_reactions.json', 'multiaxis_split_01/outer_inner_blocks.json',
             'multiaxis_split_01/axis_assignments.json', 'esm_global_01/global_features.npz', 'main_baselines_01/fresh_representations.npz',
             'conditional_model_01/input_manifest.json', 'conditional_validation_01/INDEPENDENT_FITTING_QC.json']
    write_json(out / 'input_manifest.json', {f: digest_file(root / f) for f in names})
    write_json(out / 'progress.json', {'status': 'COMPLETE_QC_PENDING', 'replicate': replicate, 'completed_blocks': len(receipts), 'updated_utc': now()})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--replicate', type=int, required=True); args = p.parse_args()
    run(Path(__file__).resolve().parent, args.replicate)
