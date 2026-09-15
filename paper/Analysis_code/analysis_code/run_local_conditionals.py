"""Nested local/global protein descriptions on identical available cohorts."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import time
import numpy as np
from conditional_core import fit_residual, predict_delta, select_lambda, model_receipt
from local_conditional_core import local_descriptions, select_joint
from run_main_baselines import rank_result, random_expectations, summarize, train_scores
from run_raw import write_json, digest_file, stable, now


def run(root, channel):
    out = root / ('local_conditional_' + channel + '_01'); prefix = channel + '_projected'
    if out.exists(): raise FileExistsError(out)
    read = lambda n: json.loads((root / n).read_text())
    if read('ordered_site_validation_01/INDEPENDENT_QC.json')['status'] != 'PASS': raise ValueError('Local feature QA gate')
    if read('conditional_validation_01/INDEPENDENT_FITTING_QC.json')['status'] != 'PASS': raise ValueError('Global fit QA gate')
    for file, h in read('multiaxis_split_01/output_checksums.json').items():
        if digest_file(root / 'multiaxis_split_01' / file) != h: raise ValueError('Frozen split changed')
    edges = read('dataset_02/core_edges.json'); rids = sorted(read('dataset_02/core_reactions.json'))
    seqs = sorted({e['sequence_sha256'] for e in edges}); si = {q: i for i, q in enumerate(seqs)}; ri = {r: i for i, r in enumerate(rids)}
    blocks = read('multiaxis_split_01/outer_inner_blocks.json'); axes = read('multiaxis_split_01/axis_assignments.json')
    site = np.load(root / 'ordered_site_features_01/ordered_features.npz', allow_pickle=False)
    esm = np.load(root / 'esm_global_01/global_features.npz', allow_pickle=False)
    rep = np.load(root / 'main_baselines_01/fresh_representations.npz', allow_pickle=False)
    if any(list(d['sequence_ids']) != seqs for d in [site, esm, rep]) or list(rep['reaction_ids']) != rids: raise ValueError('Feature order')
    available = site[prefix + '_available']; availset = {q for q, a in zip(seqs, available) if a}
    features = local_descriptions(site[prefix + '_onehot'], site[prefix + '_present_mask'], esm['features'])
    K = rep['chemical_kernel'].astype(float); H = rep['homology_bits'].astype(float)
    tie = np.empty(len(rids), dtype=int)
    for j, i in enumerate(sorted(range(len(rids)), key=lambda i: stable(rids[i]))): tie[i] = j
    def construct(block):
        Y = np.zeros((len(seqs), len(rids))); truth = defaultdict(set); fullseen = set()
        for i in block['train_edge_indices']:
            e = edges[i]; fullseen.add(ri[e['reaction_key']])
            if e['sequence_sha256'] in availset: Y[si[e['sequence_sha256']], ri[e['reaction_key']]] = 1
        for i in block['test_edge_indices']:
            e = edges[i]
            if e['sequence_sha256'] in availset: truth[e['sequence_sha256']].add(ri[e['reaction_key']])
        return Y, truth, fullseen
    out.mkdir(); allrows = []; receipts = []
    for bn, block in enumerate(blocks):
        started = time.time(); Y, truth, seen = construct(block); qs = sorted(truth); training = np.flatnonzero(Y.sum(1))
        receipt = {'block_number': bn, 'block_id': block['block_id'], 'query_ids': qs, 'training_sequence_ids': [seqs[i] for i in training],
                   'retained_train_edges': int(Y.sum()), 'retained_test_edges': sum(map(len, truth.values())),
                   'full_training_seen_labels': len(seen), 'available_training_seen_labels': int(np.sum(Y.sum(0) > 0))}
        if not len(training) or not qs:
            receipts.append({**receipt, 'status': 'NOT_EVALUABLE', 'reason': 'empty_available_training_or_query'}); continue
        logs = []; deltas = {f: [] for f in features}; innerrecords = []; targets = []; innerfits = []
        for j, inner in enumerate(block['inner_blocks']):
            if not set(inner['train_edge_indices']).issubset(block['train_edge_indices']) or not set(inner['test_edge_indices']).issubset(block['train_edge_indices']): raise ValueError('Inner leakage')
            iy, it, _ = construct(inner); iq = sorted(it)
            if not np.any(iy) or not iq:
                innerfits.append({'inner_number': j, 'status': 'NOT_EVALUABLE', 'reason': 'empty_available_training_or_query'}); continue
            models = {f: fit_residual(P, iy, K) for f, P in features.items()}; lp = models['global']['prior_log']
            for f, model in models.items():
                if not np.allclose(lp, model['prior_log'], rtol=0, atol=1e-12): raise ValueError('Feature-specific chemical backbone')
                deltas[f].extend(predict_delta(model, features[f][[si[q] for q in iq]]))
            logs.extend([lp] * len(iq)); targets.extend([sorted(it[q]) for q in iq])
            innerrecords.extend([{'group': axes['protein_groups'][q], 'query': q, 'inner_number': j} for q in iq])
            innerfits.append({'inner_number': j, 'status': 'COMPLETE', 'training_sequence_ids': [seqs[i] for i in np.flatnonzero(iy.sum(1))],
                              'query_ids': iq, 'fits': {f: model_receipt(m) for f, m in models.items()}})
        logmat = np.array(logs).reshape(-1, len(rids)); delta_mats = {f: np.array(d).reshape(-1, len(rids)) for f, d in deltas.items()}
        choices = {f: select_lambda(logmat, d, targets, innerrecords) if innerrecords else {'lambda': 0., 'reason': 'no_inner_predictions'} for f, d in delta_mats.items()}
        joint = select_joint(logmat, delta_mats['global'], delta_mats['ordered'], targets, innerrecords) if innerrecords else {'lambdas': [0., 0.], 'reason': 'no_inner_predictions', 'joint_optimizer_accepted': False}
        np.savez_compressed(out / f'block_{bn:03d}_inner_predictions.npz', log_prior=logmat, **{'delta_' + f: d for f, d in delta_mats.items()})
        write_json(out / f'block_{bn:03d}_inner_manifest.json', {'rows': innerrecords, 'targets': targets, 'fits': innerfits, 'selections': choices, 'joint_selection': joint})
        models = {f: fit_residual(P, Y, K) for f, P in features.items()}; lp = models['global']['prior_log']; qidx = [si[q] for q in qs]
        ds = {f: predict_delta(m, features[f][qidx]) for f, m in models.items()}
        scores = {'availability_chemistry_prior': lp}
        for f in features: scores[f + '_residual'] = lp[None, :] + choices[f]['lambda'] * ds[f]
        scores['joint_global_ordered'] = lp[None, :] + joint['lambdas'][0] * ds['global'] + joint['lambdas'][1] * ds['ordered']
        baseline, domains, _ = train_scores(Y, K, H, qidx, training)
        masks = {m: np.ones_like(v, dtype=bool) for m, v in scores.items()}
        for m in ['mmseqs_top1', 'mmseqs_weighted', 'homology_chemical_transport']:
            scores['matched_' + m] = baseline[m]; masks['matched_' + m] = domains[m]
        for f, name in [('global', 'matched_esm'), ('ordered', 'ordered_nearest')]:
            similarity = features[f][qidx] @ features[f][training].T
            best = training[np.argmax(similarity, axis=1)]; near = Y[best]
            scores[name + '_chemical_transport'] = (near / near.sum(1, keepdims=True)) @ K
            masks[name + '_chemical_transport'] = np.ones_like(near, dtype=bool)
            if f == 'global': scores[name + '_top1'] = near; masks[name + '_top1'] = near > 0
        np.savez_compressed(out / f'block_{bn:03d}_scores.npz', query_ids=np.array(qs), reaction_ids=np.array(rids),
                            **scores, **{'domain_' + m: mask for m, mask in masks.items()})
        saved = {'training_sequence_ids': np.array([seqs[i] for i in training]), 'prior_log': lp, 'joint_lambdas': np.array(joint['lambdas'])}
        for f, model in models.items():
            for k in ['center', 'training_features_centered', 'kernel_scale', 'residual_rms', 'coefficients']:
                saved[f + '__' + k] = model[k]
            saved[f + '__lambda'] = choices[f]['lambda']
        np.savez_compressed(out / f'block_{bn:03d}_models.npz', **saved)
        for qi, q in enumerate(qs):
            tasks = {block['task'] + '_all': truth[q]}
            if block['task'] == 'protein_cold': tasks.update({'protein_cold_seen': truth[q] & seen, 'protein_cold_unseen': truth[q] - seen})
            for task, positive in tasks.items():
                if not positive: continue
                for method in list(scores) + ['uniform_expectation']:
                    if method == 'uniform_expectation': metric = dict(random_expectations(len(rids), len(positive)))
                    else:
                        value = scores[method]; mask = masks[method]
                        metric = rank_result(value[qi] if value.ndim == 2 else value, mask[qi] if mask.ndim == 2 else mask, positive, tie)
                    allrows.append({'cohort': channel, 'task': task, 'block_number': bn, 'block_id': block['block_id'], 'sequence_sha256': q,
                                    'protein_group': axes['protein_groups'][q], 'method': method, 'query_index_in_block': qi,
                                    'target_reaction_indices': sorted(positive), **metric})
        receipts.append({**receipt, 'status': 'COMPLETE', 'fits': {f: model_receipt(m) for f, m in models.items()},
                         'selections': choices, 'joint_selection': joint, 'inner_completed_blocks': sum(r['status'] == 'COMPLETE' for r in innerfits),
                         'inner_unavailable_blocks': sum(r['status'] != 'COMPLETE' for r in innerfits), 'inner_rows': len(innerrecords),
                         'elapsed_seconds': time.time() - started})
        write_json(out / 'block_receipts.json', receipts)
        write_json(out / 'progress.json', {'status': 'RUNNING', 'completed_outer_entries': len(receipts), 'updated_utc': now()})
        print(f'{channel} block={bn + 1}/35 queries={len(qs)} joint={joint["lambdas"]} seconds={time.time() - started:.1f}', flush=True)
    with (out / 'per_query_metrics.jsonl').open('w') as h:
        for r in allrows: h.write(json.dumps(r, allow_nan=False) + '\n')
    write_json(out / 'block_receipts.json', receipts)
    ev = {'status': 'COMPLETE_INDEPENDENT_QC_PENDING', 'created_utc': now(), 'cohort': channel, 'available_sequences': len(availset),
          'available_protein_groups': len({axes['protein_groups'][q] for q in availset}), 'core_sequences': len(seqs), 'candidate_reactions': len(rids),
          'metric_rows': len(allrows), 'eligible_outer_blocks': sum(r['status'] == 'COMPLETE' for r in receipts), 'aggregate': summarize(allrows),
          'feature_dimensions': {f: P.shape[1] for f, P in features.items()},
          'joint_optimizer_fallback_blocks': sum(not r['joint_selection']['joint_optimizer_accepted'] for r in receipts if r['status'] == 'COMPLETE'),
          'independent_biological_validation': False, 'atom_level_reaction_center_model': False, 'final_router_or_tool_complete': False}
    write_json(out / 'evaluation.json', ev)
    names = ['run_local_conditionals.py', 'local_conditional_core.py', 'conditional_core.py', 'run_main_baselines.py', 'LOCAL_CONDITIONAL_MODEL_V1.md',
             'dataset_02/core_edges.json', 'dataset_02/core_reactions.json', 'ordered_site_features_01/ordered_features.npz',
             'ordered_site_validation_01/INDEPENDENT_QC.json', 'esm_global_01/global_features.npz', 'main_baselines_01/fresh_representations.npz',
             'multiaxis_split_01/outer_inner_blocks.json', 'multiaxis_split_01/axis_assignments.json']
    write_json(out / 'input_manifest.json', {n: digest_file(root / n) for n in names})
    write_json(out / 'progress.json', {'status': ev['status'], 'outer_entries': len(receipts), 'updated_utc': now()})
    print(json.dumps({k: v for k, v in ev.items() if k != 'aggregate'}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--channel', choices=['sequence', 'structure'], required=True); args = p.parse_args()
    run(Path(__file__).resolve().parent, args.channel)
