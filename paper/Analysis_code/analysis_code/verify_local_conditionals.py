"""Independent pool, linear-solve fit, objective, score and rank reconstruction."""
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


def run(root, channel):
    result = root / ('local_conditional_' + channel + '_01'); out = root / ('local_validation_' + channel + '_01')
    if out.exists(): raise FileExistsError(out)
    read = lambda p: json.loads(p.read_text())
    failures = []; inner_count = outer_count = scalar_checks = joint_checks = vectors = 0; maxdiff = 0.
    for n, h in read(result / 'input_manifest.json').items():
        if digest_file(root / n) != h: failures.append('source:' + n)
    edges = read(root / 'dataset_02/core_edges.json'); rids = sorted(read(root / 'dataset_02/core_reactions.json'))
    seqs = sorted({e['sequence_sha256'] for e in edges}); si = {q: i for i, q in enumerate(seqs)}; ri = {r: i for i, r in enumerate(rids)}
    blocks = read(root / 'multiaxis_split_01/outer_inner_blocks.json'); receipts = read(result / 'block_receipts.json'); ev = read(result / 'evaluation.json')
    axes = read(root / 'multiaxis_split_01/axis_assignments.json')
    site = np.load(root / 'ordered_site_features_01/ordered_features.npz', allow_pickle=False)
    prefix = channel + '_projected'; available = {q for q, flag in zip(seqs, site[prefix + '_available']) if flag}
    onehot = site[prefix + '_onehot'].astype(float); missing = (~site[prefix + '_present_mask']).astype(float)
    raw_ordered = np.empty((len(seqs), 35, 22)); raw_ordered[..., :21] = onehot; raw_ordered[..., 21] = missing
    raw_mask = np.empty((len(seqs), 35, 2)); raw_mask[..., 0] = onehot[..., 20]; raw_mask[..., 1] = missing
    features = {'global': np.load(root / 'esm_global_01/global_features.npz', allow_pickle=False)['features'].astype(float),
                'ordered': raw_ordered.reshape(len(seqs), -1), 'missing': raw_mask.reshape(len(seqs), -1), 'composition': onehot.sum(1)}
    for f, p in features.items(): features[f] = p / np.maximum(np.sqrt(np.sum(p * p, axis=1, keepdims=True)), 1e-12)
    rep = np.load(root / 'main_baselines_01/fresh_representations.npz', allow_pickle=False)
    K = rep['chemical_kernel'].astype(float); H = rep['homology_bits'].astype(float); tie = [stable(r) for r in rids]
    def construct(block):
        y = np.zeros((len(seqs), len(rids))); truth = defaultdict(set); seen = set()
        for i in block['train_edge_indices']:
            e = edges[i]; seen.add(ri[e['reaction_key']])
            if e['sequence_sha256'] in available: y[si[e['sequence_sha256']], ri[e['reaction_key']]] = 1
        for i in block['test_edge_indices']:
            e = edges[i]
            if e['sequence_sha256'] in available: truth[e['sequence_sha256']].add(ri[e['reaction_key']])
        return y, truth, seen
    def reconstruct(y, qs, f, declared):
        train = np.flatnonzero(y.sum(1)); target = (y[train] / y[train].sum(1, keepdims=True)) @ K
        prior = target.mean(0); prob = np.maximum(prior / prior.sum(), 1e-12); prob /= prob.sum()
        residual = target - prior; rms = float(np.sqrt(np.mean(residual ** 2))); x = features[f][train]
        center = x.mean(0); xc = x - center; gram = xc @ xc.T; scale = np.trace(gram) / len(train)
        if declared['zero_residual']: delta = np.zeros((len(qs), len(rids)))
        else:
            coefficient = solve(gram / scale + declared['alpha'] * np.eye(len(train)), residual, assume_a='pos', check_finite=False)
            delta = ((features[f][[si[q] for q in qs]] - center) @ xc.T / scale) @ coefficient / rms
        return np.log(prob), delta
    def objective(lp, ds, targets, records, coefficients):
        panels = Counter((r['group'], r['query']) for r in records); group_queries = Counter(g for g, q in panels)
        w = np.array([1 / len(group_queries) / group_queries[r['group']] / panels[(r['group'], r['query'])] for r in records])
        score = lp.copy()
        for coefficient, d in zip(coefficients, ds): score += coefficient * d
        z = logsumexp(score, axis=1); probability = np.exp(score - z[:, None])
        target_score = np.array([score[i, t].mean() for i, t in enumerate(targets)])
        gradients = [float(w @ (np.sum(probability * d, axis=1) - np.array([d[i, t].mean() for i, t in enumerate(targets)]))) for d in ds]
        return float(w @ (z - target_score)), np.array(gradients)
    metrics = [json.loads(s) for s in (result / 'per_query_metrics.jsonl').read_text().splitlines()]
    byblock = defaultdict(list); table = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for r in metrics: byblock[r['block_number']].append(r)
    for bn, block in enumerate(blocks):
        y, truth, seen = construct(block); qs = sorted(truth); train = np.flatnonzero(y.sum(1)); receipt = receipts[bn]
        if receipt['query_ids'] != qs or receipt['training_sequence_ids'] != [seqs[i] for i in train]: failures.append('outer_available_pool')
        if not qs or not len(train):
            if receipt['status'] != 'NOT_EVALUABLE' or bn in byblock: failures.append('unevaluable')
            continue
        info = read(result / f'block_{bn:03d}_inner_manifest.json')
        with np.load(result / f'block_{bn:03d}_inner_predictions.npz', allow_pickle=False) as a: inner = {n: a[n] for n in a.files}
        indices = defaultdict(list)
        for i, row in enumerate(info['rows']): indices[row['inner_number']].append(i)
        for j, innerblock in enumerate(block['inner_blocks']):
            iy, it, _ = construct(innerblock); iq = sorted(it); declared = info['fits'][j]
            if not iq or not np.any(iy):
                if declared['status'] != 'NOT_EVALUABLE' or j in indices: failures.append('inner_unevaluable')
                continue
            idx = indices[j]
            if declared['query_ids'] != iq or declared['training_sequence_ids'] != [seqs[i] for i in np.flatnonzero(iy.sum(1))]: failures.append('inner_pool')
            if [info['rows'][i]['query'] for i in idx] != iq: failures.append('inner_row_order')
            for i, q in zip(idx, iq):
                if info['targets'][i] != sorted(it[q]) or info['rows'][i]['group'] != axes['protein_groups'][q]: failures.append('inner_target')
            for f in features:
                lp, delta = reconstruct(iy, iq, f, declared['fits'][f]); inner_count += 1
                error = float(np.max(np.abs(delta - inner['delta_' + f][idx]))); maxdiff = max(maxdiff, error)
                if error > 1e-5 or not np.allclose(lp, inner['log_prior'][idx], atol=1e-10, rtol=0): failures.append(f'inner_fit:{bn}:{j}:{f}')
        if info['rows']:
            lp = inner['log_prior']; targets = info['targets']; records = info['rows']; ds = {f: inner['delta_' + f] for f in features}
            base, _ = objective(lp, [ds['global']], targets, records, [0.])
            for f, choice in info['selections'].items():
                value, gradient = objective(lp, [ds[f]], targets, records, [choice['lambda']]); scalar_checks += 1
                if choice['lambda'] < 0 or value > base + 1e-8: failures.append('scalar_loss')
                if choice['reason'] == 'convex_inner_optimum' and abs(gradient[0]) > 1e-6: failures.append('scalar_stationarity')
                if choice['reason'] == 'zero_boundary_optimum' and gradient[0] < -1e-8: failures.append('scalar_zero_boundary')
            choice = info['joint_selection']; coef = np.array(choice['lambdas'])
            value, gradient = objective(lp, [ds['global'], ds['ordered']], targets, records, coef); joint_checks += 1
            single_losses = [base] + [objective(lp, [ds[f]], targets, records, [info['selections'][f]['lambda']])[0] for f in ['global', 'ordered']]
            if np.any(coef < 0) or value > min(single_losses) + 1e-8: failures.append('joint_loss')
            projected = np.where(coef <= 1e-8, np.minimum(gradient, 0), gradient)
            if choice['joint_optimizer_accepted'] and np.max(np.abs(projected)) > 1e-6: failures.append('joint_stationarity')
        with np.load(result / f'block_{bn:03d}_scores.npz', allow_pickle=False) as a: data = {n: a[n] for n in a.files}
        ds = {}; expected = {}
        for f in features:
            lp, delta = reconstruct(y, qs, f, receipt['fits'][f]); ds[f] = delta; outer_count += 1
            expected[f + '_residual'] = lp[None, :] + receipt['selections'][f]['lambda'] * delta
        expected['availability_chemistry_prior'] = lp
        coef = receipt['joint_selection']['lambdas']; expected['joint_global_ordered'] = lp[None, :] + coef[0] * ds['global'] + coef[1] * ds['ordered']
        h = H[[si[q] for q in qs]]; balanced = y / np.maximum(y.sum(1, keepdims=True), 1)
        expected['matched_mmseqs_weighted'] = h @ y; near = np.zeros((len(qs), len(rids)))
        for qi, q in enumerate(qs):
            candidates = [i for i in train if h[qi, i] > 0]
            if candidates: near[qi] = y[min(candidates, key=lambda i: (-h[qi, i], seqs[i]))]
        expected['matched_mmseqs_top1'] = near
        expected['matched_homology_chemical_transport'] = ((h[:, train] @ balanced[train]) / np.maximum(h[:, train].sum(1, keepdims=True), 1)) @ K
        for f, name in [('global', 'matched_esm'), ('ordered', 'ordered_nearest')]:
            similarities = features[f][[si[q] for q in qs]] @ features[f][train].T
            near = np.array([y[train[min(range(len(train)), key=lambda j: (-v[j], seqs[train[j]]))]] for v in similarities])
            expected[name + '_chemical_transport'] = (near / near.sum(1, keepdims=True)) @ K
            if f == 'global': expected[name + '_top1'] = near
        masks = {n: np.ones_like(v, dtype=bool) for n, v in expected.items()}
        for name in ['matched_mmseqs_top1', 'matched_mmseqs_weighted', 'matched_esm_top1']: masks[name] = expected[name] > 0
        masks['matched_homology_chemical_transport'] = np.repeat((h[:, train].sum(1) > 0)[:, None], len(rids), axis=1)
        if list(data['query_ids']) != qs or list(data['reaction_ids']) != rids: failures.append('candidate_order')
        for name, value in expected.items():
            error = float(np.max(np.abs(value - data[name]))); maxdiff = max(maxdiff, error)
            if not np.allclose(value, data[name], rtol=1e-7, atol=1e-5): failures.append('outer_score:' + name)
            if not np.array_equal(masks[name], data['domain_' + name]): failures.append('outer_domain:' + name)
            vectors += len(qs) if value.ndim == 2 else 1
        rankcache = {}
        for row in byblock[bn]:
            q = row['sequence_sha256']; qi = row['query_index_in_block']; name = row['method']; positive = truth[q]
            if row['task'] == 'protein_cold_seen': positive = positive & seen
            elif row['task'] == 'protein_cold_unseen': positive = positive - seen
            if row['target_reaction_indices'] != sorted(positive) or row['protein_group'] != axes['protein_groups'][q]: failures.append('outer_target')
            if name == 'uniform_expectation': value = uniform(len(rids), len(positive)); covered = len(positive); ncovered = len(rids); conditional = value['rr']
            else:
                if (name, qi) not in rankcache:
                    score = data[name]; score = score[qi] if score.ndim == 2 else score
                    mask = masks[name]; mask = mask[qi] if mask.ndim == 2 else mask
                    order = sorted(np.flatnonzero(mask), key=lambda i: (-float(score[i]), tie[i]))
                    rankcache[(name, qi)] = {i: j + 1 for j, i in enumerate(order)}
                ranks = rankcache[(name, qi)]; positions = [ranks.get(i, math.inf) for i in positive]
                covered = sum(math.isfinite(i) for i in positions); ncovered = len(ranks); best = min(positions)
                value = {'rr': 1 / best, 'mean_positive_rr': math.fsum(1 / i for i in positions) / len(positions)}
                conditional = None if math.isinf(best) else 1 / best
                for k in [1, 5, 10]:
                    hits = sum(i <= k for i in positions); value['hit_at_' + str(k)] = float(hits > 0); value['recall_at_' + str(k)] = hits / len(positions)
                value['ndcg_at_10'] = math.fsum(1 / math.log2(i + 1) for i in positions if i <= 10) / math.fsum(1 / math.log2(i + 1) for i in range(1, min(len(positions), 10) + 1))
            for field, v in value.items():
                if abs(v - row[field]) > 1e-10: failures.append('rank:' + field)
            if row['positive_count'] != len(positive) or row['covered_positives'] != covered or row['covered_candidates'] != ncovered or row['candidate_count'] != len(rids): failures.append('coverage_count')
            if (row['conditional_rr'] is None) != (conditional is None): failures.append('undefined_conditional')
            elif conditional is not None and abs(row['conditional_rr'] - conditional) > 1e-10: failures.append('conditional_rr')
            table[(row['task'], name)][row['protein_group']][q].append(value)
        print(f'{channel} independent QA {bn + 1}/35', flush=True)
    for (task, name), groups in table.items():
        for field in FIELDS:
            value = math.fsum(math.fsum(math.fsum(r[field] for r in panels) / len(panels) for panels in queries.values()) / len(queries) for queries in groups.values()) / len(groups)
            if abs(value - ev['aggregate'][task][name]['protein_group_macro_' + field]) > 1e-10: failures.append('aggregate:' + field)
    if len(metrics) != ev['metric_rows'] or len(byblock) != ev['eligible_outer_blocks']: failures.append('overall_counts')
    report = {'status': 'PASS' if not failures else 'FAIL', 'created_utc': now(), 'cohort': channel, 'failures': failures,
              'inner_feature_models_reconstructed': inner_count, 'outer_feature_models_reconstructed': outer_count,
              'scalar_coefficient_checks': scalar_checks, 'joint_coefficient_checks': joint_checks,
              'score_vectors_reconstructed': vectors, 'metric_rows_recomputed': len(metrics), 'aggregates_recomputed': len(table),
              'maximum_prediction_difference': maxdiff, 'independent_biological_validation': False,
              'scope': 'All evaluable pools, fixed-alpha linear-solve fits, coefficient objectives, native masks, ranks and macro metrics',
              'verifier_sha256': digest_file(Path(__file__))}
    out.mkdir(); write_json(out / 'INDEPENDENT_QC.json', report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--channel', choices=['sequence', 'structure'], required=True); args = p.parse_args()
    run(Path(__file__).resolve().parent, args.channel)
