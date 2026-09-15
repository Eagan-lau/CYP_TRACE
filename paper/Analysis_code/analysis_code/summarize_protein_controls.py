"""Audit all fixed-size control replicates and summarize diagnostic contrasts."""
from collections import defaultdict
import json
import math
from pathlib import Path
import numpy as np
from run_raw import write_json, digest_file, stable, now


def run(root):
    out = root / 'protein_control_summary_01'
    if out.exists(): raise FileExistsError(out)
    read = lambda p: json.loads(p.read_text())
    base = root / 'protein_controls_01'; folders = [base / f'replicate_{i:03d}' for i in range(99)]
    if not all((f / 'evaluation.json').exists() for f in folders): raise ValueError('All 99 frozen replicates required')
    blocks = read(root / 'multiaxis_split_01/outer_inner_blocks.json'); edges = read(root / 'dataset_02/core_edges.json')
    seqs = sorted({e['sequence_sha256'] for e in edges}); si = {q: i for i, q in enumerate(seqs)}
    rids = sorted(read(root / 'dataset_02/core_reactions.json')); ri = {r: i for i, r in enumerate(rids)}
    axes = read(root / 'multiaxis_split_01/axis_assignments.json')
    target = {}; pools = {}; fit_pools = {}; seen_byblock = {}
    for bn, b in enumerate(blocks):
        seen_byblock[bn] = {ri[edges[i]['reaction_key']] for i in b['train_edge_indices']}
        pools[bn] = sorted({si[edges[i]['sequence_sha256']] for i in b['train_edge_indices']})
        fit_pools[(bn, 10000)] = pools[bn]
        for j, inner in enumerate(b['inner_blocks']):
            fit_pools[(bn, j + 1)] = sorted({si[edges[i]['sequence_sha256']] for i in inner['train_edge_indices']})
        t = defaultdict(set)
        for i in b['test_edge_indices']: t[edges[i]['sequence_sha256']].add(ri[edges[i]['reaction_key']])
        target[bn] = t
    original = read(root / 'conditional_model_01/evaluation.json')['aggregate']
    pair_map = {'train_profile_permuted_residual': 'esm_conditional_learned', 'query_shuffled_residual': 'esm_conditional_learned',
                'train_profile_permuted_esm_transport': 'esm_cosine_chemical_transport', 'query_shuffled_esm_transport': 'esm_cosine_chemical_transport'}
    failures = []; values = defaultdict(list); outer_fits = inner_fits = rows_checked = 0; maxdiff = 0.; zero_counts = []
    input_hashes = None
    for replicate, folder in enumerate(folders):
        ev = read(folder / 'evaluation.json'); receipts = read(folder / 'block_receipts.json'); manifest = read(folder / 'input_manifest.json')
        if input_hashes is None:
            input_hashes = manifest
            for n, h in manifest.items():
                if digest_file(root / n) != h: failures.append('source:' + n)
        elif manifest != input_hashes: failures.append('replicate_input_manifest:' + str(replicate))
        for bn, r in enumerate(receipts):
            assignments = [(f['inner_number'] + 1, f['assignment']) for f in r['inner_fits']] + [(10000, r['outer_assignment'])]
            for stage, p in assignments:
                expected = fit_pools[(bn, stage)]
                donor = np.random.Generator(np.random.PCG64(np.random.SeedSequence([20260909, replicate, bn, stage]))).permutation(expected).tolist()
                if p['recipient_indices'] != expected or p['donor_indices'] != donor or set(donor) != set(expected): failures.append('training_assignment')
                if p['fixed_points'] != sum(a == b for a, b in zip(expected, donor)): failures.append('fixed_points')
                if p['maximum_chemical_prior_log_difference'] > 1e-10: failures.append('chemical_backbone')
            qn = len(target[bn]); qperm = np.random.Generator(np.random.PCG64(np.random.SeedSequence([20260909, replicate, bn, 20000]))).permutation(qn).tolist()
            if r['query_assignment'] != qperm or r['query_sequences'] != qn: failures.append('query_assignment')
            if r['selection']['lambda'] < 0: failures.append('negative_lambda')
            inner_fits += len(r['inner_fits']); outer_fits += 1
        zero_counts.append(sum(r['selection']['lambda'] == 0 for r in receipts))
        table = defaultdict(lambda: defaultdict(lambda: defaultdict(list))); seenkeys = set(); nrows = 0
        with (folder / 'per_query_metrics.jsonl').open() as h:
            for line in h:
                row = json.loads(line); bn = row['block_number']; q = row['sequence_sha256']; task = row['task']; method = row['method']
                key = (task, bn, q, method)
                if key in seenkeys: failures.append('duplicate_metric')
                seenkeys.add(key)
                expected = target[bn][q]
                if task == 'protein_cold_seen': expected = expected & seen_byblock[bn]
                elif task == 'protein_cold_unseen': expected = expected - seen_byblock[bn]
                if row['target_reaction_indices'] != sorted(expected) or row['positive_count'] != len(expected): failures.append('target')
                if row['candidate_count'] != len(rids) or row['covered_candidates'] != len(rids) or row['covered_positives'] != len(expected): failures.append('coverage')
                if row['protein_group'] != axes['protein_groups'][q] or row['replicate'] != replicate: failures.append('group_or_replicate')
                if not 0 < row['rr'] <= 1 or row['conditional_rr'] != row['rr']: failures.append('rank_range')
                table[(task, method)][row['protein_group']][q].append(row['rr']); nrows += 1
        if nrows != ev['metric_rows'] or len(receipts) != 35 or ev['inner_refits'] != 255: failures.append('replicate_counts')
        for (task, method), groups in table.items():
            actual = math.fsum(math.fsum(math.fsum(panels) / len(panels) for panels in queries.values()) / len(queries) for queries in groups.values()) / len(groups)
            diff = abs(actual - ev['aggregate'][task][method]['protein_group_macro_rr']); maxdiff = max(maxdiff, diff)
            if diff > 1e-12: failures.append('aggregate')
            values[(task, method)].append(actual)
        rows_checked += nrows
        print(f'Control audit {replicate + 1}/99', flush=True)
    summaries = []
    for (task, method), samples in sorted(values.items()):
        obs = original[task][pair_map[method]]['protein_group_macro_rr']; v = np.array(samples)
        summaries.append({'task': task, 'control': method, 'observed_method': pair_map[method], 'observed_mrr': obs,
                          'control_mean_mrr': float(v.mean()), 'control_median_mrr': float(np.median(v)),
                          'empirical_control_2_5_to_97_5_percentiles': list(map(float, np.quantile(v, [.025, .975]))),
                          'observed_minus_control_mean': float(obs - v.mean()),
                          'replicates_at_or_above_observed': int(np.sum(v >= obs)),
                          'tail_count_is_not_a_pvalue': True, 'replicate_mrrs': samples})
    out.mkdir(); report = {'status': 'PASS' if not failures else 'FAIL', 'created_utc': now(), 'failures': failures,
                          'replicates': 99, 'metric_rows_with_verified_targets_and_aggregation': rows_checked,
                          'inner_assignments_verified': inner_fits, 'outer_assignments_verified': outer_fits,
                          'maximum_macro_difference': maxdiff, 'lambda_zero_blocks_per_replicate': zero_counts,
                          'individual_ranks_and_control_refits_independently_recomputed_here': False,
                          'formal_permutation_pvalue_valid': False, 'independent_biological_validation': False}
    write_json(out / 'ASSIGNMENT_AGGREGATION_QC.json', report)
    write_json(out / 'diagnostic_summary.json', {'status': 'DIAGNOSTIC_SUMMARY_INDIVIDUAL_SCORE_QC_PENDING', 'comparisons': summaries,
               'interpretation': 'Diagnostic information ablations, not exchangeability-certified tests or causal protein effects'})
    write_json(out / 'input_manifest.json', {'summarize_protein_controls.py': digest_file(Path(__file__)), **input_hashes})
    print(json.dumps({k: v for k, v in report.items() if k != 'lambda_zero_blocks_per_replicate'}, indent=2))
    if failures: raise SystemExit(1)


if __name__ == '__main__': run(Path(__file__).resolve().parent)
