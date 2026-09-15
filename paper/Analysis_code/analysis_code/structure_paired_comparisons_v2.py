"""Paired structural-expert comparisons with identical candidate intersections."""
from collections import defaultdict
import json
from pathlib import Path
import numpy as np
from run_raw import write_json, digest_file, now, stable
from run_main_baselines import rank_result
from paired_expert_comparisons import group_summary

PAIRS = [('foldseek_top1', 'matched_mmseqs_top1'), ('foldseek_top1', 'matched_esm_top1'),
         ('foldseek_chemical_transport', 'matched_mmseqs_chemical_transport'),
         ('foldseek_chemical_transport', 'matched_esm_chemical_transport'),
         ('foldseek_chemical_transport', 'availability_chemical_prior'),
         ('foldseek_chemical_transport', 'deposition_intensity_prior'),
         ('foldseek_chemical_transport', 'training_publication_intensity_prior')]


def run(root):
    result = root / 'structure_matched_02'; out = root / 'structure_paired_02'
    if out.exists(): raise FileExistsError(out)
    if json.loads((root / 'structure_matched_validation_02/INDEPENDENT_QC.json').read_text())['status'] != 'PASS': raise ValueError('Metric gate')
    original = [json.loads(s) for s in (result / 'per_query_metrics.jsonl').read_text().splitlines()]
    records = {(r['method'], r['task'], r['block_number'], r['sequence_sha256']): r for r in original}
    rids = sorted(json.loads((root / 'dataset_02/core_reactions.json').read_text())); tie = np.empty(len(rids), dtype=int)
    for j, i in enumerate(sorted(range(len(rids)), key=lambda i: stable(rids[i]))): tie[i] = j
    rows = []
    for bn in sorted({r['block_number'] for r in original}):
        data = np.load(result / f'block_{bn:03d}_scores.npz', allow_pickle=False)
        for left, right in PAIRS:
            for a in [r for r in original if r['method'] == left and r['block_number'] == bn]:
                b = records[(right, a['task'], bn, a['sequence_sha256'])]
                if a['target_reaction_indices'] != b['target_reaction_indices']: raise ValueError('Unmatched targets')
                qi = a['query_index_in_block']; values = []; masks = []
                for name in [left, right]:
                    value = data[name]; mask = data['domain_' + name]
                    values.append(value[qi] if value.ndim == 2 else value); masks.append(mask[qi] if mask.ndim == 2 else mask)
                common = masks[0] & masks[1]; positive = [i for i in a['target_reaction_indices'] if common[i]]
                ca = cb = None
                if positive:
                    ca, cb = [rank_result(v, common, positive, tie)['rr'] for v in values]
                rows.append({'task': a['task'], 'block_number': bn, 'sequence_sha256': a['sequence_sha256'], 'protein_group': a['protein_group'],
                             'left_method': left, 'right_method': right, 'left_rr': a['rr'], 'right_rr': b['rr'],
                             'target_positives': a['positive_count'], 'left_covered_positives': a['covered_positives'], 'right_covered_positives': b['covered_positives'],
                             'common_candidate_count': int(common.sum()), 'common_positive_count': len(positive),
                             'common_left_rr': ca, 'common_right_rr': cb, 'common_reranking_defined': bool(positive)})
        data.close()
    rng = np.random.default_rng(20260909); comparisons = []
    for task in sorted({r['task'] for r in rows}):
        for left, right in PAIRS:
            chosen = [r for r in rows if (r['task'], r['left_method'], r['right_method']) == (task, left, right)]
            common = [r for r in chosen if r['common_reranking_defined']]
            comparisons.append({'task': task, 'left_method': left, 'right_method': right,
                                'end_to_end': group_summary(chosen, 'left_rr', 'right_rr', rng),
                                'common_domain_reranking': group_summary(common, 'common_left_rr', 'common_right_rr', rng),
                                'common_undefined_query_panels': len(chosen) - len(common)})
    out.mkdir()
    with (out / 'paired_query_records.jsonl').open('w') as h:
        for r in rows: h.write(json.dumps(r, allow_nan=False) + '\n')
    write_json(out / 'comparison_summary.json', {'status': 'COMPUTED_QC_PENDING', 'created_utc': now(), 'comparisons': comparisons,
               'paired_query_rows': len(rows), 'bootstrap_seed': 20260909, 'bootstrap_replicates': 5000,
               'interval_scope': 'Conditional fixed-model protein-group bootstrap only; no chemical/publication/refit/external uncertainty',
               'independent_validation': False, 'confirmatory_hypothesis_test': False})
    names = ['structure_paired_comparisons_v2.py', 'paired_expert_comparisons.py', 'STRUCTURE_MATCHED_ANALYSIS_V1.md',
             'structure_matched_02/per_query_metrics.jsonl', 'structure_matched_validation_02/INDEPENDENT_QC.json']
    write_json(out / 'input_manifest.json', {n: digest_file(root / n) for n in names})
    print(json.dumps({'paired_query_rows': len(rows), 'comparisons': len(comparisons)}))


if __name__ == '__main__': run(Path(__file__).resolve().parent)
