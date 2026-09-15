"""Independent reconstruction checks for local conditional models with reaction-center profile channel."""

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from scipy.special import logsumexp

from verify_main_baselines import uniform, FIELDS
from run_raw import digest_file, stable, write_json, now


def _center_profile_lookup(root, reaction_ids):
    data = np.load(root / 'reaction_center_features_01/center_features.npz', allow_pickle=False)
    keys = [str(x) for x in data['reaction_keys']]
    vectors = np.asarray(data['feature_vector'], dtype=float)
    if vectors.ndim != 2:
        raise ValueError('Reaction center feature matrix must be 2D')
    index = {k: i for i, k in enumerate(keys)}
    missing = [rid for rid in reaction_ids if rid not in index]
    if missing:
        raise KeyError('Missing reaction center features: ' + ', '.join(missing))
    order = np.array([index[rid] for rid in reaction_ids], dtype=int)
    return vectors[order]


def _sequence_center_profiles(y, center_features):
    counts = y.sum(axis=1, keepdims=True)
    profiles = np.zeros((y.shape[0], center_features.shape[1]), dtype=float)
    available = counts[:, 0] > 0
    profiles[available] = y[available] @ center_features / counts[available]
    return profiles


def _smoothed_center_profiles(query_indices, train_indices, sequence_profiles, homology_matrix):
    if len(query_indices) == 0 or len(train_indices) == 0:
        return np.zeros((len(query_indices), sequence_profiles.shape[1]), dtype=float)
    weights = homology_matrix[np.ix_(query_indices, train_indices)]
    weights = np.maximum(weights, 0.0)
    train_has_profile = np.linalg.norm(sequence_profiles[train_indices], axis=1) > 0
    if not np.any(train_has_profile):
        return np.zeros((len(query_indices), sequence_profiles.shape[1]), dtype=float)
    weights[:, ~train_has_profile] = 0.0
    denom = np.sum(weights, axis=1, keepdims=True)
    valid = denom[:, 0] > 0
    out = np.zeros_like(weights @ sequence_profiles[train_indices], dtype=float)
    if np.any(valid):
        out[valid] = (weights[valid] @ sequence_profiles[train_indices]) / denom[valid]
    return out


def run(root, channel):
    result = root / ('local_conditional_' + channel + '_center_profile_01')
    out = root / ('local_validation_' + channel + '_center_profile_01')
    if out.exists():
        raise FileExistsError(out)
    read = lambda p: json.loads((root / p).read_text())

    failures = []
    inner_count = outer_count = scalar_checks = joint_checks = vectors = 0
    maxdiff = 0.

    for n, expected in read(result / 'input_manifest.json').items():
        if digest_file(root / n) != expected:
            failures.append('source:' + n)

    edges = read('dataset_02/core_edges.json')
    reaction_ids = sorted(read('dataset_02/core_reactions.json'))
    sequences = sorted({e['sequence_sha256'] for e in edges})
    sequence_index = {q: i for i, q in enumerate(sequences)}
    reaction_index = {r: i for i, r in enumerate(reaction_ids)}
    blocks = read('multiaxis_split_01/outer_inner_blocks.json')
    receipts = read(result / 'block_receipts.json')
    evaluation = read(result / 'evaluation.json')
    axes = read('multiaxis_split_01/axis_assignments.json')

    site = np.load(root / 'ordered_site_features_01/ordered_features.npz', allow_pickle=False)
    prefix = channel + '_projected'
    available = {q for q, flag in zip(sequences, site[prefix + '_available']) if flag}

    onehot = site[prefix + '_onehot'].astype(float)
    missing = (~site[prefix + '_present_mask']).astype(float)
    raw_ordered = np.empty((len(sequences), 35, 22), dtype=float)
    raw_ordered[..., :21] = onehot
    raw_ordered[..., 21] = missing
    raw_mask = np.empty((len(sequences), 35, 2), dtype=float)
    raw_mask[..., 0] = onehot[..., 20]
    raw_mask[..., 1] = missing

    esm = np.load(root / 'esm_global_01/global_features.npz', allow_pickle=False)['features'].astype(float)
    rep = np.load(root / 'main_baselines_01/fresh_representations.npz', allow_pickle=False)
    K = rep['chemical_kernel'].astype(float)
    H = rep['homology_bits'].astype(float)
    tie = [stable(r) for r in reaction_ids]

    center_features = _center_profile_lookup(root, reaction_ids)
    features = {
        'global': esm,
        'ordered': raw_ordered.reshape(len(sequences), -1),
        'missing': raw_mask.reshape(len(sequences), -1),
        'composition': onehot.sum(1),
        'center_profile': np.zeros((len(sequences), center_features.shape[1]), dtype=float),
    }
    for f, p in features.items():
        if f == 'center_profile':
            continue
        norm = np.sqrt(np.sum(p * p, axis=1, keepdims=True))
        features[f] = p / np.maximum(norm, 1e-12)

    def construct(block):
        y = np.zeros((len(sequences), len(reaction_ids)))
        truth = defaultdict(set)
        seen = set()
        for edge_index in block['train_edge_indices']:
            edge = edges[edge_index]
            seen.add(reaction_index[edge['reaction_key']])
            if edge['sequence_sha256'] in available:
                y[sequence_index[edge['sequence_sha256']], reaction_index[edge['reaction_key']]] = 1
        for edge_index in block['test_edge_indices']:
            edge = edges[edge_index]
            if edge['sequence_sha256'] in available:
                truth[edge['sequence_sha256']].add(reaction_index[edge['reaction_key']])
        return y, truth, seen

    def reconstruct(y, query_indices, f, model, query_features, train_features):
        train = np.flatnonzero(y.sum(1))
        xtrain = train_features[f][train]
        if not len(train):
            return np.full(len(reaction_ids), np.nan), np.zeros((len(query_indices), len(reaction_ids)), dtype=float)
        target = (y[train] / np.maximum(y[train].sum(1, keepdims=True), 1)) @ K
        prior = target.mean(0)
        prob = np.maximum(prior / prior.sum(), 1e-12)
        prob = prob / prob.sum()
        residual = target - prior
        rms = float(np.sqrt(np.mean(residual ** 2)))
        center = xtrain.mean(0)
        xcentered = xtrain - center
        gram = xcentered @ xcentered.T
        scale = np.trace(gram) / len(train)
        if model['zero_residual']:
            delta = np.zeros((len(query_indices), len(reaction_ids)), dtype=float)
        else:
            coefficient = solve(
                gram / scale + model['alpha'] * np.eye(len(train)),
                residual,
                assume_a='pos',
                check_finite=False
            )
            delta = ((query_features[f] - center) @ xcentered.T / scale) @ coefficient / rms
        return np.log(prob), delta

    def objective(lp, ds, targets, rows, coefficients):
        panels = Counter((r['group'], r['query']) for r in rows)
        group_queries = Counter(g for g, q in panels)
        w = np.array([1 / len(group_queries) / group_queries[r['group']] / panels[(r['group'], r['query'])] for r in rows])
        score = lp.copy()
        for coefficient, d in zip(coefficients, ds):
            score += coefficient * d
        z = logsumexp(score, axis=1)
        probability = np.exp(score - z[:, None])
        target_score = np.array([score[i, t].mean() for i, t in enumerate(targets)])
        gradients = [float(w @ (np.sum(probability * d, axis=1) - np.array([d[i, t].mean() for i, t in enumerate(targets)])))
                     for d in ds]
        return float(w @ (z - target_score)), np.array(gradients)

    metrics = [json.loads(line) for line in (result / 'per_query_metrics.jsonl').read_text().splitlines()]
    by_block = defaultdict(list)
    table = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for row in metrics:
        by_block[row['block_number']].append(row)

    for block_number, block in enumerate(blocks):
        y, truth, seen = construct(block)
        query_ids = sorted(truth)
        training = np.flatnonzero(y.sum(1))
        receipt = receipts[block_number]

        if receipt['query_ids'] != query_ids or receipt['training_sequence_ids'] != [sequences[i] for i in training]:
            failures.append('outer_available_pool')
        if not query_ids or not len(training):
            if receipt['status'] != 'NOT_EVALUABLE' or block_number in by_block:
                failures.append('unevaluable')
            continue

        sequence_center_profiles = _sequence_center_profiles(y, center_features)
        center_train_profiles = _smoothed_center_profiles(training, training, sequence_center_profiles, H)
        center_train_full = np.zeros_like(features['center_profile'])
        center_train_full[training] = center_train_profiles
        center_query_profiles = _smoothed_center_profiles([sequence_index[q] for q in query_ids], training, sequence_center_profiles, H)
        center_query_full = np.zeros_like(features['center_profile'])
        center_query_full[[sequence_index[q] for q in query_ids]] = center_query_profiles

        inner_manifest = read(result / f'block_{block_number:03d}_inner_manifest.json')
        with np.load(result / f'block_{block_number:03d}_inner_predictions.npz', allow_pickle=False) as m:
            inner = {name: m[name] for name in m.files}
        inner_index = defaultdict(list)
        for i, row in enumerate(inner_manifest['rows']):
            inner_index[row['inner_number']].append(i)

        for inner_number, inner_block in enumerate(block['inner_blocks']):
            inner_y, inner_truth, _ = construct(inner_block)
            inner_query_ids = sorted(inner_truth)
            declared = inner_manifest['fits'][inner_number]
            if not inner_query_ids or not np.any(inner_y):
                if declared['status'] != 'NOT_EVALUABLE' or inner_number in inner_index:
                    failures.append('inner_unevaluable')
                continue

            inner_query_indices = [sequence_index[q] for q in inner_query_ids]
            inner_train = np.flatnonzero(inner_y.sum(1))
            inner_sequence_center_profiles = _sequence_center_profiles(inner_y, center_features)
            inner_train_profiles = _smoothed_center_profiles(inner_train, inner_train, inner_sequence_center_profiles, H)
            inner_train_full = np.zeros_like(features['center_profile'])
            inner_train_full[inner_train] = inner_train_profiles
            inner_query_profiles = _smoothed_center_profiles(inner_query_indices, inner_train, inner_sequence_center_profiles, H)

            # Rebuild full feature views exactly as training ran.
            inner_train_features = {k: v.copy() for k, v in features.items()}
            inner_train_features['center_profile'] = inner_train_full

            if declared['query_ids'] != inner_query_ids or declared['training_sequence_ids'] != [sequences[i] for i in inner_train]:
                failures.append('inner_pool')
            if [inner_manifest['rows'][i]['query'] for i in inner_index[inner_number]] != inner_query_ids:
                failures.append('inner_row_order')

            for i, sequence in zip(inner_index[inner_number], inner_query_ids):
                if inner_manifest['targets'][i] != sorted(inner_truth[sequence]) or inner_manifest['rows'][i]['group'] != axes['protein_groups'][sequence]:
                    failures.append('inner_target')

            inner_query_features = {
                'global': features['global'][inner_query_indices],
                'ordered': features['ordered'][inner_query_indices],
                'missing': features['missing'][inner_query_indices],
                'composition': features['composition'][inner_query_indices],
                'center_profile': inner_query_profiles,
            }
            for feature_name in inner_train_features:
                lp, delta = reconstruct(
                    inner_y,
                    inner_query_indices,
                    feature_name,
                    declared['fits'][feature_name],
                    inner_query_features,
                    inner_train_features,
                )
                inner_count += 1
                index = inner_index[inner_number]
                diff = float(np.max(np.abs(delta - inner['delta_' + feature_name][index])))
                maxdiff = max(maxdiff, diff)
                if not np.allclose(lp, inner['log_prior'][index], atol=1e-10, rtol=0):
                    failures.append(f'inner_fit_prior:{block_number}:{inner_number}:{feature_name}')
                if diff > 1e-5:
                    failures.append(f'inner_fit:{block_number}:{inner_number}:{feature_name}')

        if inner_manifest['rows']:
            lp = inner['log_prior']
            targets = inner_manifest['targets']
            records = inner_manifest['rows']
            ds = {name: inner['delta_' + name] for name in features}
            base, _ = objective(lp, [ds['global']], targets, records, [0.])
            for feature_name, choice in inner_manifest['selections'].items():
                value, gradient = objective(lp, [ds[feature_name]], targets, records, [choice['lambda']]); scalar_checks += 1
                if choice['lambda'] < 0 or value > base + 1e-8:
                    failures.append('scalar_loss')
                if choice['reason'] == 'convex_inner_optimum' and abs(gradient[0]) > 1e-6:
                    failures.append('scalar_stationarity')
                if choice['reason'] == 'zero_boundary_optimum' and gradient[0] < -1e-8:
                    failures.append('scalar_zero_boundary')
            choice = inner_manifest['joint_selection']
            coef = np.array(choice['lambdas'])
            value, gradient = objective(lp, [ds['global'], ds['ordered']], targets, records, coef); joint_checks += 1
            single_losses = [base] + [objective(lp, [ds[f]], targets, records, [inner_manifest['selections'][f]['lambda']])[0] for f in ['global', 'ordered']]
            if np.any(coef < 0) or value > min(single_losses) + 1e-8:
                failures.append('joint_loss')
            projected = np.where(coef <= 1e-8, np.minimum(gradient, 0), gradient)
            if choice['joint_optimizer_accepted'] and np.max(np.abs(projected)) > 1e-6:
                failures.append('joint_stationarity')

        with np.load(result / f'block_{block_number:03d}_scores.npz', allow_pickle=False) as m:
            scores = {name: m[name] for name in m.files}

        scores_expected = {}
        ds = {}
        query_indices = [sequence_index[q] for q in query_ids]
        outer_query_features = {
            'global': features['global'][query_indices],
            'ordered': features['ordered'][query_indices],
            'missing': features['missing'][query_indices],
            'composition': features['composition'][query_indices],
            'center_profile': center_query_full[query_indices],
        }
        outer_train_features = {
            'global': features['global'],
            'ordered': features['ordered'],
            'missing': features['missing'],
            'composition': features['composition'],
            'center_profile': center_train_full,
        }
        for feature_name in features:
            lp, delta = reconstruct(
                y,
                query_indices,
                feature_name,
                receipt['fits'][feature_name],
                outer_query_features,
                outer_train_features,
            )
            outer_count += 1
            ds[feature_name] = delta
            scores_expected[feature_name + '_residual'] = lp[None, :] + receipt['selections'][feature_name]['lambda'] * delta
        scores_expected['availability_chemistry_prior'] = lp
        coeff = receipt['joint_selection']['lambdas']
        scores_expected['joint_global_ordered'] = lp[None, :] + coeff[0] * ds['global'] + coeff[1] * ds['ordered']

        h = H[[sequence_index[q] for q in query_ids]]
        balanced = y / np.maximum(y.sum(axis=1, keepdims=True), 1)
        scores_expected['matched_mmseqs_weighted'] = h @ y
        nearest = np.zeros((len(query_ids), len(reaction_ids)))
        for qi, q in enumerate(query_ids):
            candidates = [i for i in training if h[qi, i] > 0]
            if candidates:
                selected = min(candidates, key=lambda i: (-h[qi, i], sequences[i]))
                nearest[qi] = y[selected]
        scores_expected['matched_mmseqs_top1'] = nearest
        scores_expected['matched_homology_chemical_transport'] = ((h[:, training] @ balanced[training]) / np.maximum(h[:, training].sum(axis=1, keepdims=True), 1)) @ K
        for feature_name, method_name in [('global', 'matched_esm'), ('ordered', 'ordered_nearest')]:
            similarity = features[feature_name][[sequence_index[q] for q in query_ids]] @ features[feature_name][training].T
            nearest_profile = np.array([y[training[min(range(len(training)), key=lambda j: (-v[j], sequences[training[j]]))]] for v in similarity])
            scores_expected[method_name + '_chemical_transport'] = (nearest_profile / nearest_profile.sum(axis=1, keepdims=True)) @ K
            if feature_name == 'global':
                scores_expected[method_name + '_top1'] = nearest_profile

        masks = {name: np.ones_like(value, dtype=bool) for name, value in scores_expected.items()}
        for method in ['matched_mmseqs_top1', 'matched_mmseqs_weighted', 'matched_esm_top1', 'matched_homology_chemical_transport']:
            masks[method] = scores_expected[method] > 0 if method != 'matched_homology_chemical_transport' else np.broadcast_to((h[:, training].sum(axis=1, keepdims=True) > 0), (len(query_ids), len(reaction_ids)))

        if list(scores['query_ids']) != query_ids or list(scores['reaction_ids']) != reaction_ids:
            failures.append('candidate_order')
        for name, expected_value in scores_expected.items():
            actual = scores[name]
            vectors += len(query_ids) if expected_value.ndim == 2 else 1
            diff = float(np.max(np.abs(expected_value - actual)))
            maxdiff = max(maxdiff, diff)
            if not np.allclose(expected_value, actual, rtol=1e-7, atol=1e-5):
                failures.append('outer_score:' + name)
            if not np.array_equal(masks[name], scores['domain_' + name]):
                failures.append('outer_domain:' + name)

        rankcache = {}
        for row in by_block[block_number]:
            query = row['sequence_sha256']
            query_idx = query_ids.index(query)
            method = row['method']
            positive = truth[query]
            if row['task'] == 'protein_cold_seen':
                positive = positive & seen
            elif row['task'] == 'protein_cold_unseen':
                positive = positive - seen

            if row['target_reaction_indices'] != sorted(positive) or row['protein_group'] != axes['protein_groups'][query]:
                failures.append('outer_target')
            if method == 'uniform_expectation':
                value = uniform(len(reaction_ids), len(positive))
                covered = len(positive)
                ncovered = len(reaction_ids)
                conditional = value['rr']
            else:
                if (method, query_idx) not in rankcache:
                    score = scores[method]
                    selected = score[query_idx] if score.ndim == 2 else score
                    mask = masks[method][query_idx] if masks[method].ndim == 2 else masks[method]
                    order = sorted(np.flatnonzero(mask), key=lambda i: (-float(selected[i]), tie[i]))
                    rankcache[(method, query_idx)] = {i: j + 1 for j, i in enumerate(order)}
                ranks = rankcache[(method, query_idx)]
                positions = [ranks.get(i, math.inf) for i in positive]
                covered = sum(math.isfinite(i) for i in positions)
                mask = masks[method][query_idx] if masks[method].ndim == 2 else masks[method]
                ncovered = int(np.sum(mask))
                if not positive:
                    continue
                best = min(positions)
                value = {
                    'rr': 1 / best,
                    'mean_positive_rr': math.fsum(1 / i for i in positions) / len(positions),
                }
                conditional = None if math.isinf(best) else 1 / best
                for k in [1, 5, 10]:
                    hits = sum(i <= k for i in positions)
                    value['recall_at_' + str(k)] = hits / len(positions)
                    value['hit_at_' + str(k)] = float(hits > 0)
                if positions and any(i <= 10 for i in positions):
                    value['ndcg_at_10'] = math.fsum(1 / math.log2(i + 1) for i in positions if i <= 10) / math.fsum(1 / math.log2(j + 1) for j in range(1, min(len(positions), 10) + 1))
                else:
                    value['ndcg_at_10'] = 0.0
            for field, v in value.items():
                if abs(v - row[field]) > 1e-10:
                    failures.append('rank:' + field)
            if row['positive_count'] != len(positive) or row['covered_positives'] != covered or row['covered_candidates'] != ncovered or row['candidate_count'] != len(reaction_ids):
                failures.append('coverage_count')
            if (row['conditional_rr'] is None) != (conditional is None):
                failures.append('undefined_conditional')
            elif conditional is not None and abs(row['conditional_rr'] - conditional) > 1e-10:
                failures.append('conditional_rr')
            table[(row['task'], method)][row['protein_group']][query].append(value)
        print(f'{channel} center-profile QA {block_number + 1}/35', flush=True)

    for (task, name), groups in table.items():
        for field in FIELDS:
            value = math.fsum(math.fsum(math.fsum(r[field] for r in panels) / len(panels) for panels in queries.values()) / len(queries) for queries in groups.values()) / len(groups)
            if abs(value - evaluation['aggregate'][task][name]['protein_group_macro_' + field]) > 1e-10:
                failures.append('aggregate:' + field)

    report = {
        'status': 'PASS' if not failures else 'FAIL',
        'created_utc': now(),
        'cohort': channel,
        'failures': failures,
        'inner_feature_models_reconstructed': inner_count,
        'outer_feature_models_reconstructed': outer_count,
        'scalar_coefficient_checks': scalar_checks,
        'joint_coefficient_checks': joint_checks,
        'score_vectors_reconstructed': vectors,
        'metric_rows_recomputed': len(metrics),
        'aggregates_recomputed': len(table),
        'maximum_prediction_difference': maxdiff,
        'independent_biological_validation': False,
        'scope': 'Center-profile augmented local model reconstruction and score/rank parity',
        'verifier_sha256': digest_file(Path(__file__)),
    }

    out.mkdir()
    write_json(out / 'INDEPENDENT_QC.json', report)
    print(json.dumps(report, indent=2))
    if len(metrics) != evaluation['metric_rows'] or len(by_block) != evaluation['eligible_outer_blocks']:
        raise SystemExit(1)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--channel', choices=['sequence', 'structure'], required=True)
    run(Path(__file__).resolve().parent, parser.parse_args().channel)
