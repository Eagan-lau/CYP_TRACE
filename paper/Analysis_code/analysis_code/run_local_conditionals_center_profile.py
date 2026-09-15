"""Nested local/global protein models with an additional reaction-center profile channel."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import time

import numpy as np

from conditional_core import fit_residual, predict_delta, select_lambda, model_receipt
from local_conditional_core import local_descriptions, select_joint
from run_main_baselines import rank_result, random_expectations, summarize, train_scores
from run_raw import digest_file, now, stable, write_json


def _center_profile_lookup(root, reaction_ids):
    data = np.load(root / 'reaction_center_features_01/center_features.npz', allow_pickle=False)
    keys = [str(x) for x in data['reaction_keys']]
    vectors = np.asarray(data['feature_vector'], dtype=float)
    if vectors.ndim != 2:
        raise ValueError('Reaction center feature matrix must be 2D')
    if len(keys) != len(reaction_ids):
        raise ValueError('Reaction-center vector count does not match reaction count')
    key_to_index = {k: i for i, k in enumerate(keys)}
    order = []
    for rid in reaction_ids:
        if rid not in key_to_index:
            raise KeyError(f'Missing reaction center features: {rid}')
        order.append(key_to_index[rid])
    return vectors[np.asarray(order, dtype=int)]


def _sequence_center_profiles(Y, center_features):
    seq_profile = np.zeros((Y.shape[0], center_features.shape[1]), dtype=float)
    counts = Y.sum(1, keepdims=True)
    observed = counts[:, 0] > 0
    if np.any(observed):
        seq_profile[observed] = Y[observed] @ center_features / counts[observed]
    return seq_profile


def _smoothed_center_features(query_indices, train_indices, sequence_center_profiles, homology_matrix):
    if len(query_indices) == 0 or len(train_indices) == 0:
        return np.zeros((len(query_indices), sequence_center_profiles.shape[1]), dtype=float)
    weights = homology_matrix[np.ix_(query_indices, train_indices)]
    weights = np.maximum(weights, 0.0)
    train_has_profile = np.linalg.norm(sequence_center_profiles[train_indices], axis=1) > 0
    if not np.any(train_has_profile):
        return np.zeros((len(query_indices), sequence_center_profiles.shape[1]), dtype=float)
    weights[:, ~train_has_profile] = 0.0
    denom = np.sum(weights, axis=1, keepdims=True)
    valid = denom[:, 0] > 0
    if not np.any(valid):
        return np.zeros((len(query_indices), sequence_center_profiles.shape[1]), dtype=float)
    smoothed = np.zeros_like(weights @ sequence_center_profiles[train_indices])
    smoothed[valid] = (weights[valid] @ sequence_center_profiles[train_indices]) / denom[valid]
    return smoothed


def run(root, channel):
    out = root / ('local_conditional_' + channel + '_center_profile_01')
    if out.exists():
        raise FileExistsError(out)
    read = lambda n: json.loads((root / n).read_text())
    if read('ordered_site_validation_01/INDEPENDENT_QC.json')['status'] != 'PASS':
        raise ValueError('Local feature QA gate')
    if read('conditional_validation_01/INDEPENDENT_FITTING_QC.json')['status'] != 'PASS':
        raise ValueError('Global fit QA gate')
    for file, hash_value in read('multiaxis_split_01/output_checksums.json').items():
        if digest_file(root / 'multiaxis_split_01' / file) != hash_value:
            raise ValueError('Frozen split changed')

    edges = read('dataset_02/core_edges.json')
    reaction_ids = sorted(read('dataset_02/core_reactions.json'))
    sequences = sorted({e['sequence_sha256'] for e in edges})
    sequence_index = {q: i for i, q in enumerate(sequences)}
    reaction_index = {r: i for i, r in enumerate(reaction_ids)}
    blocks = read('multiaxis_split_01/outer_inner_blocks.json')
    axes = read('multiaxis_split_01/axis_assignments.json')

    site = np.load(root / 'ordered_site_features_01/ordered_features.npz', allow_pickle=False)
    esm = np.load(root / 'esm_global_01/global_features.npz', allow_pickle=False)
    rep = np.load(root / 'main_baselines_01/fresh_representations.npz', allow_pickle=False)
    if any(list(d['sequence_ids']) != sequences for d in [site, esm, rep]):
        raise ValueError('Feature order')
    if list(rep['reaction_ids']) != reaction_ids:
        raise ValueError('Reaction ID order')

    available = site[(channel + '_projected_available')]
    available_set = {q for q, flag in zip(sequences, available) if flag}
    features = local_descriptions(site[channel + '_projected_onehot'], site[channel + '_projected_present_mask'], esm['features'])
    chemical_kernel = rep['chemical_kernel'].astype(float)
    homology_bits = rep['homology_bits'].astype(float)
    center_features = _center_profile_lookup(root, reaction_ids)
    tie = np.empty(len(reaction_ids), dtype=int)
    for j, i in enumerate(sorted(range(len(reaction_ids)), key=lambda i: stable(reaction_ids[i]))):
        tie[i] = j

    def construct(block):
        y = np.zeros((len(sequences), len(reaction_ids)))
        truth = defaultdict(set)
        full_seen = set()
        for edge_index in block['train_edge_indices']:
            edge = edges[edge_index]
            full_seen.add(reaction_index[edge['reaction_key']])
            if edge['sequence_sha256'] in available_set:
                y[sequence_index[edge['sequence_sha256']], reaction_index[edge['reaction_key']]] = 1
        for edge_index in block['test_edge_indices']:
            edge = edges[edge_index]
            if edge['sequence_sha256'] in available_set:
                truth[edge['sequence_sha256']].add(reaction_index[edge['reaction_key']])
        return y, truth, full_seen

    out.mkdir()
    all_rows = []
    receipts = []
    for block_number, block in enumerate(blocks):
        started = time.time()
        y, truth, seen = construct(block)
        query_ids = sorted(truth)
        training = np.flatnonzero(y.sum(1))
        receipt = {
            'block_number': block_number,
            'block_id': block['block_id'],
            'query_ids': query_ids,
            'training_sequence_ids': [sequences[i] for i in training],
            'retained_train_edges': int(y.sum()),
            'retained_test_edges': sum(len(v) for v in truth.values()),
            'full_training_seen_labels': len(seen),
            'available_training_seen_labels': int(np.sum(y.sum(0) > 0)),
        }
        if not len(training) or not query_ids:
            receipts.append({**receipt, 'status': 'NOT_EVALUABLE', 'reason': 'empty_available_training_or_query'})
            continue

        logs = []
        deltas = {f: [] for f in features}
        deltas['center_profile'] = []
        inner_rows = []
        targets = []
        inner_fits = []
        sequence_center_profiles = _sequence_center_profiles(y, center_features)

        for inner_number, inner in enumerate(block['inner_blocks']):
            inner_y, inner_truth, _ = construct(inner)
            inner_query_ids = sorted(inner_truth)
            if not np.any(inner_y) or not inner_query_ids:
                inner_fits.append({'inner_number': inner_number, 'status': 'NOT_EVALUABLE', 'reason': 'empty_available_training_or_query'})
                continue

            inner_training = np.flatnonzero(inner_y.sum(1))
            inner_seq_center_profiles = _sequence_center_profiles(inner_y, center_features)
            center_train_inner = _smoothed_center_features(
                inner_training,
                inner_training,
                inner_seq_center_profiles,
                homology_bits,
            )
            center_train_full_inner = np.zeros((len(sequences), center_features.shape[1]), dtype=float)
            center_train_full_inner[inner_training] = center_train_inner
            inner_feature_views = {
                'global': features['global'],
                'ordered': features['ordered'],
                'missing': features['missing'],
                'composition': features['composition'],
                'center_profile': center_train_full_inner,
            }
            models = {f: fit_residual(P, inner_y, chemical_kernel) for f, P in inner_feature_views.items()}
            lp = models['global']['prior_log']
            for f, model in models.items():
                if not np.allclose(lp, model['prior_log'], rtol=0, atol=1e-12):
                    raise ValueError('Feature-specific chemical backbone')

            center_query_inner = _smoothed_center_features(
                [sequence_index[q] for q in inner_query_ids],
                inner_training,
                inner_seq_center_profiles,
                homology_bits,
            )
            inner_query_features = {
                'global': features['global'][[sequence_index[q] for q in inner_query_ids]],
                'ordered': features['ordered'][[sequence_index[q] for q in inner_query_ids]],
                'missing': features['missing'][[sequence_index[q] for q in inner_query_ids]],
                'composition': features['composition'][[sequence_index[q] for q in inner_query_ids]],
                'center_profile': center_query_inner,
            }
            for feature_name, model in models.items():
                deltas[feature_name].extend(predict_delta(model, inner_query_features[feature_name]))

            logs.extend([lp] * len(inner_query_ids))
            targets.extend([sorted(inner_truth[q]) for q in inner_query_ids])
            inner_rows.extend([{'group': axes['protein_groups'][q], 'query': q, 'inner_number': inner_number} for q in inner_query_ids])
            inner_fits.append({
                'inner_number': inner_number,
                'status': 'COMPLETE',
                'training_sequence_ids': [sequences[i] for i in inner_training],
                'query_ids': inner_query_ids,
                'fits': {f: model_receipt(model) for f, model in models.items()},
            })

        log_matrix = np.array(logs).reshape(-1, len(reaction_ids))
        delta_matrices = {f: np.array(v).reshape(-1, len(reaction_ids)) for f, v in deltas.items()}
        selections = {f: select_lambda(log_matrix, d, targets, inner_rows) if inner_rows else {'lambda': 0., 'reason': 'no_inner_predictions'}
                      for f, d in delta_matrices.items()}
        joint = select_joint(
            log_matrix,
            delta_matrices['global'],
            delta_matrices['ordered'],
            targets,
            inner_rows,
        ) if inner_rows else {'lambdas': [0., 0.], 'reason': 'no_inner_predictions', 'joint_optimizer_accepted': False}
        np.savez_compressed(
            out / f'block_{block_number:03d}_inner_predictions.npz',
            log_prior=log_matrix,
            **{'delta_' + f: d for f, d in delta_matrices.items()},
        )
        write_json(out / f'block_{block_number:03d}_inner_manifest.json', {
            'rows': inner_rows,
            'targets': targets,
            'fits': inner_fits,
            'selections': selections,
            'joint_selection': joint,
        })

        center_train_outer = _smoothed_center_features(
            training,
            training,
            sequence_center_profiles,
            homology_bits,
        )
        center_train_full_outer = np.zeros((len(sequences), center_features.shape[1]), dtype=float)
        center_train_full_outer[training] = center_train_outer
        train_features = {
            'global': features['global'],
            'ordered': features['ordered'],
            'missing': features['missing'],
            'composition': features['composition'],
            'center_profile': center_train_full_outer,
        }
        models = {f: fit_residual(P, y, chemical_kernel) for f, P in train_features.items()}
        lp = models['global']['prior_log']
        qidx = [sequence_index[q] for q in query_ids]
        center_query_outer = _smoothed_center_features(
            qidx,
            training,
            sequence_center_profiles,
            homology_bits,
        )
        query_features = {
            'global': features['global'][qidx],
            'ordered': features['ordered'][qidx],
            'missing': features['missing'][qidx],
            'composition': features['composition'][qidx],
            'center_profile': center_query_outer,
        }
        delta_outer = {f: predict_delta(model, query_features[f]) for f, model in models.items()}

        scores = {'availability_chemistry_prior': lp}
        for f in features:
            scores[f + '_residual'] = lp[None, :] + selections[f]['lambda'] * delta_outer[f]
        scores['center_profile_residual'] = lp[None, :] + selections['center_profile']['lambda'] * delta_outer['center_profile']
        scores['joint_global_ordered'] = lp[None, :] + joint['lambdas'][0] * delta_outer['global'] + joint['lambdas'][1] * delta_outer['ordered']

        baseline, domains, _ = train_scores(y, chemical_kernel, homology_bits, qidx, training)
        masks = {m: np.ones_like(v, dtype=bool) for m, v in scores.items()}
        for m in ['mmseqs_top1', 'mmseqs_weighted', 'homology_chemical_transport']:
            scores['matched_' + m] = baseline[m]
        for m in ['mmseqs_top1', 'mmseqs_weighted', 'homology_chemical_transport']:
            masks['matched_' + m] = domains[m]

        for feature_name, method_name in [('global', 'matched_esm'), ('ordered', 'ordered_nearest')]:
            similarity = features[feature_name][qidx] @ features[feature_name][training].T
            nearest = training[np.argmax(similarity, axis=1)]
            near = y[nearest]
            scores[method_name + '_chemical_transport'] = (near / near.sum(1, keepdims=True)) @ chemical_kernel
            masks[method_name + '_chemical_transport'] = np.ones_like(near, dtype=bool)
            if feature_name == 'global':
                scores[method_name + '_top1'] = near
                masks[method_name + '_top1'] = near > 0

        np.savez_compressed(
            out / f'block_{block_number:03d}_scores.npz',
            query_ids=np.array(query_ids),
            reaction_ids=np.array(reaction_ids),
            **scores,
            **{'domain_' + m: mask for m, mask in masks.items()},
        )

        saved = {'training_sequence_ids': np.array([sequences[i] for i in training]), 'prior_log': lp, 'joint_lambdas': np.array(joint['lambdas'])}
        for feature_name, model in models.items():
            for key in ['center', 'training_features_centered', 'kernel_scale', 'residual_rms', 'coefficients']:
                saved[feature_name + '__' + key] = model[key]
            saved[feature_name + '__lambda'] = selections[feature_name]['lambda']

        np.savez_compressed(out / f'block_{block_number:03d}_models.npz', **saved)
        for query_index, sequence in enumerate(query_ids):
            tasks = {block['task'] + '_all': truth[sequence]}
            if block['task'] == 'protein_cold':
                tasks.update({
                    'protein_cold_seen': truth[sequence] & seen,
                    'protein_cold_unseen': truth[sequence] - seen,
                })
            for task_name, positive in tasks.items():
                if not positive:
                    continue
                for method in list(scores) + ['uniform_expectation']:
                    if method == 'uniform_expectation':
                        metric = dict(random_expectations(len(reaction_ids), len(positive)))
                    else:
                        value = scores[method]
                        mask = masks[method]
                        metric = rank_result(
                            value[query_index] if value.ndim == 2 else value,
                            mask[query_index] if mask.ndim == 2 else mask,
                            positive,
                            tie,
                        )
                    all_rows.append({
                        'cohort': channel,
                        'task': task_name,
                        'block_number': block_number,
                        'block_id': block['block_id'],
                        'sequence_sha256': sequence,
                        'protein_group': axes['protein_groups'][sequence],
                        'method': method,
                        'query_index_in_block': query_index,
                        'target_reaction_indices': sorted(positive),
                        **metric,
                    })

        receipts.append({
            **receipt,
            'status': 'COMPLETE',
            'fits': {f: model_receipt(m) for f, m in models.items()},
            'selections': selections,
            'joint_selection': joint,
            'inner_completed_blocks': sum(r['status'] == 'COMPLETE' for r in inner_fits),
            'inner_unavailable_blocks': sum(r['status'] != 'COMPLETE' for r in inner_fits),
            'inner_rows': len(inner_rows),
            'elapsed_seconds': time.time() - started,
        })
        write_json(out / 'block_receipts.json', receipts)
        write_json(out / 'progress.json', {'status': 'RUNNING', 'completed_outer_entries': len(receipts), 'updated_utc': now()})
        print(f'{channel} block={block_number + 1}/35 queries={len(query_ids)} joint={joint["lambdas"]} seconds={time.time() - started:.1f}', flush=True)

    with (out / 'per_query_metrics.jsonl').open('w', encoding='utf-8') as handle:
        for row in all_rows:
            handle.write(json.dumps(row, allow_nan=False) + '\n')

    write_json(out / 'block_receipts.json', receipts)
    evaluation = {
        'status': 'COMPLETE_INDEPENDENT_QC_PENDING',
        'created_utc': now(),
        'cohort': channel,
        'available_sequences': len(available_set),
        'available_protein_groups': len({axes['protein_groups'][q] for q in available_set}),
        'core_sequences': len(sequences),
        'candidate_reactions': len(reaction_ids),
        'metric_rows': len(all_rows),
        'eligible_outer_blocks': sum(r['status'] == 'COMPLETE' for r in receipts),
        'aggregate': summarize(all_rows),
        'feature_dimensions': {**{f: P.shape[1] for f, P in features.items()}, 'center_profile': center_features.shape[1]},
        'joint_optimizer_fallback_blocks': sum(not r['joint_selection']['joint_optimizer_accepted'] for r in receipts if r['status'] == 'COMPLETE'),
        'independent_biological_validation': False,
        'atom_level_reaction_center_model': True,
        'final_router_or_tool_complete': False,
    }
    write_json(out / 'evaluation.json', evaluation)

    inputs = [
        'run_local_conditionals.py',
        'run_local_conditionals_center_profile.py',
        'local_conditional_core.py',
        'conditional_core.py',
        'run_main_baselines.py',
        'LOCAL_CONDITIONAL_MODEL_V1.md',
        'dataset_02/core_edges.json',
        'dataset_02/core_reactions.json',
        'ordered_site_features_01/ordered_features.npz',
        'reaction_center_features_01/center_features.npz',
        'reaction_center_features_01/summary.json',
        'ordered_site_validation_01/INDEPENDENT_QC.json',
        'esm_global_01/global_features.npz',
        'main_baselines_01/fresh_representations.npz',
        'multiaxis_split_01/outer_inner_blocks.json',
        'multiaxis_split_01/axis_assignments.json',
    ]
    write_json(out / 'input_manifest.json', {n: digest_file(root / n) for n in inputs})
    write_json(out / 'progress.json', {'status': evaluation['status'], 'outer_entries': len(receipts), 'updated_utc': now()})
    print(json.dumps({k: v for k, v in evaluation.items() if k != 'aggregate'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--channel', choices=['sequence', 'structure'], required=True)
    run(Path(__file__).resolve().parent, parser.parse_args().channel)
