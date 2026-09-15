"""Independent sorted-list ranks and count-weighted bootstrap quantiles."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import numpy as np
from run_raw import digest_file, stable, write_json, now
from verify_main_baselines import uniform


def run(root, channel):
    result = root / f'local_paired_{channel}_01'; source_dir = root / f'local_conditional_{channel}_01'
    out = root / f'local_paired_validation_{channel}_01'
    if out.exists(): raise FileExistsError(out)
    read = lambda p: json.loads(p.read_text())
    declared = read(result / 'comparison_summary.json'); failures = []; checked = 0; intervals = 0
    for n, h in read(result / 'input_manifest.json').items():
        if digest_file(root / n) != h: failures.append('source:' + n)
    rows = [json.loads(s) for s in (result / 'paired_query_records.jsonl').read_text().splitlines()]
    orig = [json.loads(s) for s in (source_dir / 'per_query_metrics.jsonl').read_text().splitlines()]
    source = {(r['method'], r['task'], r['block_number'], r['sequence_sha256']): r for r in orig}
    rids = sorted(read(root / 'dataset_02/core_reactions.json')); hashes = [stable(r) for r in rids]; byblock = defaultdict(list)
    for r in rows: byblock[r['block_number']].append(r)
    keys = {(r['left_method'], r['right_method'], r['task'], r['block_number'], r['sequence_sha256']) for r in rows}
    if len(keys) != len(rows): failures.append('duplicate_paired_rows')
    for bn, members in sorted(byblock.items()):
        with np.load(source_dir / f'block_{bn:03d}_scores.npz', allow_pickle=False) as saved: d = {k: saved[k] for k in saved.files}
        if list(d['reaction_ids']) != rids: failures.append('reaction_order')
        for row in members:
            values = []; masks = []; targets = []
            for side in ['left', 'right']:
                name = row[side + '_method']; original = source[(name, row['task'], bn, row['sequence_sha256'])]
                if original['rr'] != row[side + '_rr'] or original['covered_positives'] != row[side + '_covered_positives']: failures.append('source_metric')
                if original['protein_group'] != row['protein_group'] or original['positive_count'] != row['target_positives']: failures.append('source_denominator')
                qi = original['query_index_in_block']
                if str(d['query_ids'][qi]) != row['sequence_sha256']: failures.append('query_order')
                if name == 'uniform_expectation': v = None; m = np.ones(len(rids), bool)
                else:
                    v = d[name]; m = d['domain_' + name]
                    v = v[qi] if v.ndim == 2 else v; m = m[qi] if m.ndim == 2 else m
                if int(m.sum()) != row[side + '_candidate_count']: failures.append('native_domain_count')
                values.append(v); masks.append(m); targets.append(original['target_reaction_indices'])
            if targets[0] != targets[1]: failures.append('target_disagreement')
            common = [i for i in range(len(rids)) if masks[0][i] and masks[1][i]]; commonset = set(common)
            truth = [i for i in targets[0] if i in commonset]
            if len(common) != row['common_candidate_count'] or len(truth) != row['common_positive_count'] or len(rids) != row['full_candidate_count']: failures.append('intersection_denominator')
            if bool(truth) != row['common_reranking_defined']: failures.append('undefined_flag')
            for side, value in zip(['left', 'right'], values):
                if not truth:
                    if row['common_' + side + '_rr'] is not None: failures.append('undefined_value')
                    continue
                if value is None: expected = uniform(len(common), len(truth))['rr']
                else:
                    order = sorted(common, key=lambda i: (-float(value[i]), hashes[i])); truthset = set(truth)
                    expected = 1 / next(j + 1 for j, i in enumerate(order) if i in truthset)
                if abs(expected - row['common_' + side + '_rr']) > 1e-12: failures.append('rank')
            checked += 1
        print(f'{channel} independently ranked block {bn}', flush=True)
    rng = np.random.Generator(np.random.PCG64(declared['bootstrap_seed'])); scopes = 0
    for c in declared['comparisons']:
        selected = [r for r in rows if (r['task'], r['left_method'], r['right_method']) == (c['task'], c['left_method'], c['right_method'])]
        for field, key in [('original_target_positive_instances', 'target_positives'), ('left_covered_positive_instances', 'left_covered_positives'), ('right_covered_positive_instances', 'right_covered_positives'), ('common_positive_instances', 'common_positive_count')]:
            if c['coverage'][field] != sum(r[key] for r in selected): failures.append('coverage_summary')
        if c['coverage']['common_undefined_query_panels'] != sum(not r['common_reranking_defined'] for r in selected): failures.append('undefined_summary')
        for scope, a, b, chosen in [('end_to_end', 'left_rr', 'right_rr', selected), ('common_domain_reranking', 'common_left_rr', 'common_right_rr', [r for r in selected if r['common_reranking_defined']])]:
            groups = defaultdict(lambda: defaultdict(list))
            for r in chosen: groups[r['protein_group']][r['sequence_sha256']].append(r)
            values = []
            for group, queries in sorted(groups.items()):
                values.append([math.fsum(math.fsum(r[k] for r in panels) / len(panels) for panels in queries.values()) / len(queries) for k in [a, b]])
            n = len(values); expected = c[scope]
            if n != expected['protein_groups'] or len(chosen) != expected['query_panels'] or len({r['sequence_sha256'] for r in chosen}) != expected['unique_queries']: failures.append('summary_denominator')
            if n:
                v = np.array(values); delta = v[:, 0] - v[:, 1]
                for k, x in [('left_mrr', math.fsum(v[:, 0]) / n), ('right_mrr', math.fsum(v[:, 1]) / n), ('difference', math.fsum(delta) / n)]:
                    if abs(x - expected[k]) > 1e-12: failures.append('summary_mean')
            elif any(expected[k] is not None for k in ['left_mrr', 'right_mrr', 'difference']): failures.append('undefined_mean')
            if n >= 2:
                draws = rng.integers(n, size=(declared['bootstrap_replicates'], n)); counts = np.zeros((len(draws), n), dtype=int)
                np.add.at(counts, (np.repeat(np.arange(len(draws)), n), draws.ravel()), 1)
                boot = np.sort(counts @ delta / n); ci = []
                for p in [.025, .975]:
                    index = p * (len(boot) - 1); low = math.floor(index); high = math.ceil(index)
                    ci.append(boot[low] + (index - low) * (boot[high] - boot[low]))
                if not np.allclose(ci, expected['conditional_protein_bootstrap_95ci'], rtol=0, atol=1e-12): failures.append('interval')
                intervals += 1
            elif expected['conditional_protein_bootstrap_95ci'] is not None: failures.append('undefined_interval')
            scopes += 1
    if checked != declared['paired_query_rows'] or len(declared['comparisons']) != 70: failures.append('total_count')
    out.mkdir(); report = {'status': 'PASS' if not failures else 'FAIL', 'created_utc': now(), 'cohort': channel, 'failures': failures,
                          'paired_rows_recomputed': checked, 'summary_scopes_recomputed': scopes, 'conditional_intervals_recomputed': intervals,
                          'verifier_sha256': digest_file(Path(__file__)), 'independent_biological_validation': False, 'multiplicity_adjusted': False}
    write_json(out / 'INDEPENDENT_QC.json', report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--channel', choices=['sequence', 'structure'], required=True)
    run(Path(__file__).resolve().parent, p.parse_args().channel)
