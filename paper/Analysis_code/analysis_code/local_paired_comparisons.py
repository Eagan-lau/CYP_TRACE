"""All-task exploratory matched-cohort representation contrasts."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import numpy as np
from run_raw import write_json, digest_file, now, stable
from run_main_baselines import rank_result, random_expectations
from paired_expert_comparisons import group_summary

PAIRS = [('ordered_residual', m) for m in ['global_residual', 'missing_residual', 'composition_residual', 'availability_chemistry_prior']]
PAIRS += [('joint_global_ordered', m) for m in ['global_residual', 'ordered_residual', 'availability_chemistry_prior', 'matched_mmseqs_weighted']]
PAIRS += [('ordered_nearest_chemical_transport', m) for m in ['matched_esm_chemical_transport', 'matched_mmseqs_weighted']]
PAIRS += [(m, 'uniform_expectation') for m in ['joint_global_ordered', 'ordered_nearest_chemical_transport', 'matched_esm_chemical_transport', 'global_residual']]


def run(root, channel):
    result = root / f'local_conditional_{channel}_01'; out = root / f'local_paired_{channel}_01'
    if out.exists(): raise FileExistsError(out)
    read = lambda n: json.loads((root / n).read_text())
    qc = f'local_validation_{channel}_01/INDEPENDENT_QC.json'
    if read(qc)['status'] != 'PASS': raise ValueError('Model QC gate')
    original = [json.loads(s) for s in (result / 'per_query_metrics.jsonl').read_text().splitlines()]
    records = {(r['method'], r['task'], r['block_number'], r['sequence_sha256']): r for r in original}
    if len(records) != len(original): raise ValueError('Duplicate records')
    byblock = defaultdict(list)
    for r in original: byblock[r['block_number']].append(r)
    rids = sorted(read('dataset_02/core_reactions.json')); n = len(rids); tie = np.empty(n, int)
    for j, i in enumerate(sorted(range(n), key=lambda i: stable(rids[i]))): tie[i] = j
    rows = []; scorefiles = []
    for bn, blockrows in sorted(byblock.items()):
        file = result / f'block_{bn:03d}_scores.npz'; scorefiles.append(file.relative_to(root).as_posix())
        with np.load(file, allow_pickle=False) as saved: d = {k: saved[k] for k in saved.files}
        if list(d['reaction_ids']) != rids: raise ValueError('Catalogue order')
        def values(name, row):
            qi = row['query_index_in_block']
            if str(d['query_ids'][qi]) != row['sequence_sha256']: raise ValueError('Query order')
            if name == 'uniform_expectation': return None, np.ones(n, bool)
            v = d[name]; m = d['domain_' + name]
            return (v[qi] if v.ndim == 2 else v), (m[qi] if m.ndim == 2 else m)
        for left, right in PAIRS:
            for a in (r for r in blockrows if r['method'] == left):
                b = records[(right, a['task'], bn, a['sequence_sha256'])]
                if a['target_reaction_indices'] != b['target_reaction_indices'] or a['protein_group'] != b['protein_group']: raise ValueError('Pair denominator')
                va, ma = values(left, a); vb, mb = values(right, b); common = ma & mb
                positive = [i for i in a['target_reaction_indices'] if common[i]]; ca = cb = None
                if positive:
                    ca = rank_result(va, common, positive, tie)['rr']
                    cb = random_expectations(int(common.sum()), len(positive))['rr'] if right == 'uniform_expectation' else rank_result(vb, common, positive, tie)['rr']
                rows.append({'cohort': channel, 'task': a['task'], 'block_number': bn, 'sequence_sha256': a['sequence_sha256'], 'protein_group': a['protein_group'],
                             'left_method': left, 'right_method': right, 'left_rr': a['rr'], 'right_rr': b['rr'], 'full_candidate_count': n,
                             'target_positives': a['positive_count'], 'left_covered_positives': a['covered_positives'], 'right_covered_positives': b['covered_positives'],
                             'left_candidate_count': int(ma.sum()), 'right_candidate_count': int(mb.sum()), 'common_candidate_count': int(common.sum()),
                             'common_positive_count': len(positive), 'common_left_rr': ca, 'common_right_rr': cb, 'common_reranking_defined': bool(positive)})
        print(f'{channel} paired block {bn}', flush=True)
    rng = np.random.default_rng(20260909); comparisons = []
    for task in sorted({r['task'] for r in rows}):
        for left, right in PAIRS:
            chosen = [r for r in rows if (r['task'], r['left_method'], r['right_method']) == (task, left, right)]
            common = [r for r in chosen if r['common_reranking_defined']]
            comparisons.append({'task': task, 'left_method': left, 'right_method': right,
                                'end_to_end': group_summary(chosen, 'left_rr', 'right_rr', rng),
                                'common_domain_reranking': group_summary(common, 'common_left_rr', 'common_right_rr', rng),
                                'coverage': {'original_target_positive_instances': sum(r['target_positives'] for r in chosen),
                                             'left_covered_positive_instances': sum(r['left_covered_positives'] for r in chosen),
                                             'right_covered_positive_instances': sum(r['right_covered_positives'] for r in chosen),
                                             'common_positive_instances': sum(r['common_positive_count'] for r in chosen),
                                             'common_undefined_query_panels': len(chosen) - len(common)}})
    out.mkdir()
    with (out / 'paired_query_records.jsonl').open('w') as h:
        for r in rows: h.write(json.dumps(r, allow_nan=False) + '\n')
    write_json(out / 'comparison_summary.json', {'status': 'COMPUTED_QC_PENDING', 'created_utc': now(), 'cohort': channel,
               'comparisons': comparisons, 'paired_query_rows': len(rows), 'bootstrap_seed': 20260909, 'bootstrap_replicates': 5000,
               'comparison_list_recorded_after_point_estimates': True, 'confirmatory_hypothesis_test': False, 'independent_validation': False,
               'interval_scope': 'Conditional fixed-fit protein-group bootstrap; no chemical/publication/refit/external or multiplicity-adjusted inference'})
    names = ['local_paired_comparisons.py', 'paired_expert_comparisons.py', 'LOCAL_PAIRED_ANALYSIS_V1.md',
             result.relative_to(root).as_posix() + '/per_query_metrics.jsonl', qc, 'dataset_02/core_reactions.json'] + scorefiles
    write_json(out / 'input_manifest.json', {n: digest_file(root / n) for n in names})
    print(json.dumps({'cohort': channel, 'paired_rows': len(rows), 'comparisons': len(comparisons)}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--channel', choices=['sequence', 'structure'], required=True)
    run(Path(__file__).resolve().parent, p.parse_args().channel)
