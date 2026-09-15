"""Independent primary-pair, train-score, rank and aggregation reconstruction."""
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import numpy as np
from run_raw import digest_file, write_json, now, stable


def run(root):
    result = root / 'structure_matched_02'; out = root / 'structure_matched_validation_02'
    if out.exists(): raise FileExistsError(out)
    read = lambda n: json.loads((root / n).read_text())
    failures = []; maxdiff = 0.; fits = 0
    for n, sha in read('structure_matched_02/input_manifest.json').items():
        if digest_file(root / n) != sha: failures.append('source:' + n)
    edges = read('dataset_02/core_edges.json'); rids = sorted(read('dataset_02/core_reactions.json'))
    seqs = sorted({e['sequence_sha256'] for e in edges}); si = {q: i for i, q in enumerate(seqs)}; ri = {r: i for i, r in enumerate(rids)}
    selected = read('structure_mapping_audit_02/primary_representatives.json'); available = set(selected)
    census = {r['sequence_sha256']: r for r in read('structure_mapping_audit_02/availability_census.json')}
    blocks = read('multiaxis_split_01/outer_inner_blocks.json'); receipts = read('structure_matched_02/block_receipts.json')
    rep = np.load(root / 'main_baselines_01/fresh_representations.npz', allow_pickle=False)
    esm = np.load(root / 'esm_global_01/global_features.npz', allow_pickle=False)['features'].astype(float)
    esm /= np.linalg.norm(esm, axis=1, keepdims=True)
    K = rep['chemical_kernel'].astype(float); M = rep['homology_bits'].astype(float)
    F = np.zeros(M.shape); pairs = {}
    with (result / 'primary_structure_pair_rows.tsv').open() as h:
        for r in csv.reader(h, delimiter='\t'):
            if (r[0], r[1]) in pairs: failures.append('duplicate_primary_pair')
            pairs[(r[0], r[1])] = r
    for q, qa in selected.items():
        for t, ta in selected.items():
            row = pairs.get((qa['foldseek_identifier'], ta['foldseek_identifier']))
            if row is None: continue
            if float(row[10]) <= .001 and float(row[3]) >= .5 and float(row[4]) >= .5 and float(row[11]) > 0:
                F[si[q], si[t]] = float(row[11])
    saved = np.load(result / 'foldseek_primary_matrix.npz', allow_pickle=False)
    if list(saved['sequence_ids']) != seqs or not np.array_equal(F, saved['homology_bits']): failures.append('foldseek_projection')
    metrics = [json.loads(s) for s in (result / 'per_query_metrics.jsonl').read_text().splitlines()]
    byblock = defaultdict(list)
    for r in metrics: byblock[r['block_number']].append(r)
    fields = ['rr', 'mean_positive_rr', 'recall_at_1', 'recall_at_5', 'recall_at_10', 'hit_at_1', 'hit_at_5', 'hit_at_10', 'ndcg_at_10']
    table = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for bn, block in enumerate(blocks):
        truth = defaultdict(set); labels = defaultdict(set); pubs = defaultdict(set); seen = set()
        for i in block['train_edge_indices']:
            e = edges[i]; seen.add(ri[e['reaction_key']])
            if e['sequence_sha256'] in available:
                labels[e['sequence_sha256']].add(ri[e['reaction_key']]); pubs[e['sequence_sha256']].update(e['publication_ids'])
        for i in block['test_edge_indices']:
            e = edges[i]
            if e['sequence_sha256'] in available: truth[e['sequence_sha256']].add(ri[e['reaction_key']])
        qs = sorted(truth); training = sorted(labels); receipt = receipts[bn]
        if receipt['query_ids'] != qs or receipt['training_sequence_ids'] != training: failures.append('availability:' + str(bn))
        if not qs or not training:
            if receipt['status'] != 'NOT_EVALUABLE' or bn in byblock: failures.append('unevaluable:' + str(bn))
            continue
        data = np.load(result / f'block_{bn:03d}_scores.npz', allow_pickle=False)
        if list(data['query_ids']) != qs or list(data['reaction_ids']) != rids: failures.append('ordering')
        y = np.zeros((len(training), len(rids)))
        for i, t in enumerate(training): y[i, sorted(labels[t])] = 1
        balanced = y / y.sum(1, keepdims=True)
        dw = np.array([math.log1p(census[t]['primary_chain_count']) for t in training])
        pw = np.array([math.log1p(len(pubs[t])) if pubs[t] else 1. for t in training])
        expected = {'availability_chemical_prior': balanced.mean(0) @ K,
                    'deposition_intensity_prior': (dw @ balanced / dw.sum()) @ K,
                    'training_publication_intensity_prior': (pw @ balanced / pw.sum()) @ K}
        masks = {n: np.ones(len(rids), dtype=bool) for n in expected}
        for prefix, matrix in [('foldseek', F), ('matched_mmseqs', M)]:
            near = []; mass = []; transported = []; has = []
            for q in qs:
                w = np.array([matrix[si[q], si[t]] for t in training]); eligible = [i for i, v in enumerate(w) if v > 0]
                best = min(eligible, key=lambda i: (-w[i], training[i])) if eligible else None
                near.append(np.zeros(len(rids)) if best is None else y[best])
                mass.append(sum((w[i] * y[i] for i in range(len(training))), np.zeros(len(rids))))
                transported.append(sum((w[i] * balanced[i] for i in range(len(training))), np.zeros(len(rids))) @ K / max(float(w.sum()), 1.))
                has.append(bool(eligible))
            expected[prefix + '_top1'] = np.array(near); expected[prefix + '_weighted'] = np.array(mass)
            expected[prefix + '_chemical_transport'] = np.array(transported)
            masks[prefix + '_top1'] = np.array(near) > 0; masks[prefix + '_weighted'] = np.array(mass) > 0
            masks[prefix + '_chemical_transport'] = np.repeat(np.array(has)[:, None], len(rids), axis=1)
        near = []
        for q in qs:
            best = min(range(len(training)), key=lambda i: (-float(np.dot(esm[si[q]], esm[si[training[i]]])), training[i]))
            near.append(y[best])
        near = np.array(near); expected['matched_esm_top1'] = near
        expected['matched_esm_chemical_transport'] = (near / near.sum(1, keepdims=True)) @ K
        masks['matched_esm_top1'] = near > 0; masks['matched_esm_chemical_transport'] = np.ones_like(near, dtype=bool)
        for name, value in expected.items():
            diff = float(np.max(np.abs(data[name] - value))); maxdiff = max(maxdiff, diff)
            if not np.allclose(data[name], value, rtol=1e-10, atol=1e-10): failures.append('fit:' + str(bn) + ':' + name)
            if not np.array_equal(data['domain_' + name], masks[name]): failures.append('domain:' + name)
            fits += len(qs) if value.ndim == 2 else 1
        ranks = {}
        for row in byblock[bn]:
            q = row['sequence_sha256']; qi = row['query_index_in_block']; name = row['method']
            positive = truth[q]
            if row['task'] == 'protein_cold_seen': positive = positive & seen
            elif row['task'] == 'protein_cold_unseen': positive = positive - seen
            if set(row['target_reaction_indices']) != positive or qs[qi] != q: failures.append('target_or_query')
            if (name, qi) not in ranks:
                value = expected[name]; value = value[qi] if value.ndim == 2 else value
                mask = masks[name]; mask = mask[qi] if mask.ndim == 2 else mask
                order = sorted(np.flatnonzero(mask), key=lambda i: (-float(value[i]), stable(rids[i])))
                ranks[(name, qi)] = {i: j + 1 for j, i in enumerate(order)}
            ranked = ranks[(name, qi)]; positions = [ranked.get(i, math.inf) for i in sorted(positive)]
            best = min(positions); rr = 1 / best
            metric = {'rr': rr, 'mean_positive_rr': math.fsum(1 / p for p in positions) / len(positive)}
            for k in [1, 5, 10]:
                hits = sum(p <= k for p in positions); metric['hit_at_' + str(k)] = float(hits > 0); metric['recall_at_' + str(k)] = hits / len(positive)
            metric['ndcg_at_10'] = math.fsum(1 / math.log2(p + 1) for p in positions if p <= 10) / math.fsum(1 / math.log2(p + 1) for p in range(1, min(len(positive), 10) + 1))
            for f, actual in metric.items():
                if abs(actual - row[f]) > 1e-12: failures.append('rank:' + f)
            if row['positive_count'] != len(positive) or row['covered_positives'] != sum(math.isfinite(p) for p in positions) or row['candidate_count'] != len(rids) or row['covered_candidates'] != len(ranked): failures.append('denominator')
            if (row['conditional_rr'] is None) != math.isinf(best): failures.append('undefined_conditional')
            table[(row['task'], name)][row['protein_group']][q].append(metric)
    ev = read('structure_matched_02/evaluation.json')
    for (task, name), groups in table.items():
        for f in fields:
            actual = math.fsum(math.fsum(math.fsum(p[f] for p in panels) / len(panels) for panels in queries.values()) / len(queries) for queries in groups.values()) / len(groups)
            if abs(actual - ev['aggregate'][task][name]['protein_group_macro_' + f]) > 1e-12: failures.append('aggregate:' + name)
    if len(metrics) != ev['metric_rows'] or len(byblock) != ev['eligible_outer_blocks']: failures.append('total_counts')
    report = {'status': 'PASS' if not failures else 'FAIL', 'created_utc': now(), 'failures': failures,
              'metric_rows_recomputed': len(metrics), 'score_vectors_recomputed': fits, 'aggregates_recomputed': len(table),
              'eligible_blocks_verified': len(byblock), 'unevaluable_blocks_verified': len(blocks) - len(byblock),
              'primary_sequence_pairs_reconstructed': len(selected) ** 2, 'maximum_score_difference': maxdiff,
              'scope': 'All saved primary-pair projection, fits, masks, targets, ranks and macro metrics; no biological certification',
              'verifier_sha256': digest_file(Path(__file__))}
    out.mkdir(); write_json(out / 'INDEPENDENT_QC.json', report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == '__main__': run(Path(__file__).resolve().parent)
