"""Same-availability Foldseek, sequence and query-independent comparisons."""
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import numpy as np
from Bio import SeqIO
from run_raw import stable, digest_file, write_json, now
from run_main_baselines import rank_result, summarize


def transfer(Y, K, similarities, queries, training, prefix):
    h = similarities[queries][:, training]
    label = Y[training]
    balanced = label / label.sum(1, keepdims=True)
    weighted = h @ label
    near = np.zeros_like(weighted)
    has = h.sum(1) > 0
    if np.any(has): near[has] = label[np.argmax(h[has], axis=1)]
    transported = (h @ balanced / np.maximum(h.sum(1, keepdims=True), 1)) @ K
    scores = {prefix + '_top1': near, prefix + '_weighted': weighted, prefix + '_chemical_transport': transported}
    masks = {prefix + '_top1': near > 0, prefix + '_weighted': weighted > 0,
             prefix + '_chemical_transport': np.broadcast_to(has[:, None], transported.shape).copy()}
    return scores, masks


def build_structure_matrix(root, seqs, selected):
    fsseq = {r.id: str(r.seq) for r in SeqIO.parse(root / 'foldseek_raw_01/structure_sequences.fasta', 'fasta')}
    index = {q: i for i, q in enumerate(seqs)}; aliases = defaultdict(list)
    for q, r in selected.items(): aliases[r['foldseek_identifier']].append(index[q])
    path = root / 'foldseek_raw_01/all_vs_all.tsv'
    expected = (root / 'foldseek_raw_01/output.sha256').read_text().splitlines()
    expected = [s.split()[0] for s in expected if s.split()[-1].endswith('all_vs_all.tsv')]
    if len(expected) != 1 or digest_file(path) != expected[0]: raise ValueError('Foldseek output hash changed')
    H = np.zeros((len(seqs), len(seqs)), dtype=np.float64)
    counts = Counter(); selected_rows = []; self_seen = set(); pair_seen = set()
    with path.open() as handle:
        for row in csv.reader(handle, delimiter='\t'):
            if len(row) != 16: raise ValueError('Unexpected Foldseek columns')
            q, t = row[:2]; ident, qc, tc, length, qs, qe, ts, te, ev, bits, ql, tl, tm, lddt = map(float, row[2:])
            if q not in fsseq or t not in fsseq: raise ValueError('Unknown Foldseek identifier')
            if not all(math.isfinite(v) for v in [ident, qc, tc, length, qs, qe, ts, te, ev, bits, ql, tl]):
                raise ValueError('Non-finite required Foldseek field: ' + repr(row))
            counts['undefined_auxiliary_tm_rows'] += not math.isfinite(tm)
            counts['undefined_auxiliary_lddt_rows'] += not math.isfinite(lddt)
            if not (0 <= ident <= 1 and 0 <= qc <= 1 and 0 <= tc <= 1 and ev >= 0):
                raise ValueError('Foldseek fraction/score range: ' + repr(row))
            if math.isfinite(lddt) and not 0 <= lddt <= 1: raise ValueError('Finite lDDT out of range: ' + repr(row))
            if int(ql) != len(fsseq[q]) or int(tl) != len(fsseq[t]) or not (1 <= qs <= qe <= ql and 1 <= ts <= te <= tl):
                raise ValueError('Foldseek coordinate range')
            counts['alignment_rows'] += 1
            counts['approximate_tm_above_one'] += tm > 1
            counts['nonpositive_bit_score_rows'] += bits <= 0
            if q == t: self_seen.add(q)
            if q not in aliases or t not in aliases: continue
            if (q, t) in pair_seen: raise ValueError('Duplicate primary structure pair')
            pair_seen.add((q, t)); selected_rows.append(row)
            if ev <= .001 and min(qc, tc) >= .5 and bits > 0:
                counts['qualified_primary_structure_pairs'] += 1
                for qi in aliases[q]:
                    for ti in aliases[t]: H[qi, ti] = bits
    counts['database_sequences'] = len(fsseq); counts['self_identifiers_present'] = len(self_seen)
    counts['primary_structure_identifiers'] = len(aliases); counts['primary_pair_rows'] = len(selected_rows)
    counts['database_identifiers_without_self'] = len(set(fsseq) - self_seen)
    if not set(aliases).issubset(self_seen): raise ValueError('Missing primary Foldseek self entries')
    return H, selected_rows, dict(counts)


def run(root):
    out = root / 'structure_matched_02'
    if out.exists(): raise FileExistsError(out)
    read = lambda name: json.loads((root / name).read_text())
    if read('structure_association_validation_02/INDEPENDENT_QC.json')['status'] != 'PASS':
        raise ValueError('Independent structure-association gate failed')
    for file, expected in read('multiaxis_split_01/output_checksums.json').items():
        if digest_file(root / 'multiaxis_split_01' / file) != expected: raise ValueError('Changed frozen split')
    edges = read('dataset_02/core_edges.json'); rids = sorted(read('dataset_02/core_reactions.json'))
    seqs = sorted({e['sequence_sha256'] for e in edges}); si = {q: i for i, q in enumerate(seqs)}; ri = {r: i for i, r in enumerate(rids)}
    selected = read('structure_mapping_audit_02/primary_representatives.json'); available = set(selected)
    census = {r['sequence_sha256']: r for r in read('structure_mapping_audit_02/availability_census.json')}
    blocks = read('multiaxis_split_01/outer_inner_blocks.json'); axes = read('multiaxis_split_01/axis_assignments.json')
    rep = np.load(root / 'main_baselines_01/fresh_representations.npz', allow_pickle=False)
    esm = np.load(root / 'esm_global_01/global_features.npz', allow_pickle=False)
    if list(rep['sequence_ids']) != seqs or list(esm['sequence_ids']) != seqs or list(rep['reaction_ids']) != rids:
        raise ValueError('Representation order mismatch')
    if digest_file(root / 'esm_global_01/global_features.npz') != read('esm_global_01/receipt.json')['features_sha256']:
        raise ValueError('ESM source checksum mismatch')
    K = rep['chemical_kernel'].astype(float); M = rep['homology_bits'].astype(float)
    P = esm['features'].astype(float); P /= np.linalg.norm(P, axis=1, keepdims=True)
    F, fsrows, counts = build_structure_matrix(root, seqs, selected)
    tie = np.empty(len(rids), dtype=int)
    for rank, i in enumerate(sorted(range(len(rids)), key=lambda i: stable(rids[i]))): tie[i] = rank
    out.mkdir(); np.savez_compressed(out / 'foldseek_primary_matrix.npz', sequence_ids=np.array(seqs), homology_bits=F)
    with (out / 'primary_structure_pair_rows.tsv').open('w') as h:
        csv.writer(h, delimiter='\t', lineterminator='\n').writerows(fsrows)
    allrows = []; receipts = []
    for bn, block in enumerate(blocks):
        Y = np.zeros((len(seqs), len(rids))); truth = defaultdict(set); fullseen = set(); pubs = defaultdict(set)
        for i in block['train_edge_indices']:
            e = edges[i]; fullseen.add(ri[e['reaction_key']])
            if e['sequence_sha256'] in available:
                Y[si[e['sequence_sha256']], ri[e['reaction_key']]] = 1
                pubs[e['sequence_sha256']].update(e['publication_ids'])
        for i in block['test_edge_indices']:
            e = edges[i]
            if e['sequence_sha256'] in available: truth[e['sequence_sha256']].add(ri[e['reaction_key']])
        queries = sorted(truth); training = np.flatnonzero(Y.sum(1)); qidx = [si[q] for q in queries]
        receipt = {'block_number': bn, 'block_id': block['block_id'], 'training_sequence_ids': [seqs[i] for i in training],
                   'query_ids': queries, 'training_edges': int(Y.sum()), 'test_edges': sum(map(len, truth.values())),
                   'full_training_seen_labels': len(fullseen), 'matched_training_seen_labels': int(np.sum(Y.sum(0) > 0))}
        if not len(training) or not queries:
            receipts.append({**receipt, 'status': 'NOT_EVALUABLE', 'reason': 'no_primary_training_or_query'}); continue
        scores, masks = transfer(Y, K, F, qidx, training, 'foldseek')
        s, m = transfer(Y, K, M, qidx, training, 'matched_mmseqs'); scores.update(s); masks.update(m)
        best = training[np.argmax(P[qidx] @ P[training].T, axis=1)]
        nearest = Y[best]; balanced = Y[training] / Y[training].sum(1, keepdims=True)
        scores['matched_esm_top1'] = nearest
        scores['matched_esm_chemical_transport'] = (nearest / nearest.sum(1, keepdims=True)) @ K
        scores['availability_chemical_prior'] = balanced.mean(0) @ K
        dw = np.array([math.log1p(census[seqs[i]]['primary_chain_count']) for i in training])
        pw = np.array([math.log1p(len(pubs[seqs[i]])) if pubs[seqs[i]] else 1. for i in training])
        scores['deposition_intensity_prior'] = (dw @ balanced / dw.sum()) @ K
        scores['training_publication_intensity_prior'] = (pw @ balanced / pw.sum()) @ K
        for name in scores:
            if name not in masks: masks[name] = np.ones_like(scores[name], dtype=bool)
        masks['matched_esm_top1'] = nearest > 0
        np.savez_compressed(out / f'block_{bn:03d}_scores.npz', query_ids=np.array(queries), reaction_ids=np.array(rids),
                            **scores, **{'domain_' + k: v for k, v in masks.items()})
        for qi, q in enumerate(queries):
            tasks = {block['task'] + '_all': truth[q]}
            if block['task'] == 'protein_cold':
                tasks.update({'protein_cold_seen': truth[q] & fullseen, 'protein_cold_unseen': truth[q] - fullseen})
            for task, positives in tasks.items():
                if not positives: continue
                for method, value in scores.items():
                    score = value[qi] if value.ndim == 2 else value
                    mask = masks[method][qi] if masks[method].ndim == 2 else masks[method]
                    allrows.append({'task': task, 'block_id': block['block_id'], 'block_number': bn, 'sequence_sha256': q,
                                    'protein_group': axes['protein_groups'][q], 'method': method, 'target_reaction_indices': sorted(positives),
                                    'query_index_in_block': qi, **rank_result(score, mask, positives, tie)})
        receipts.append({**receipt, 'status': 'COMPLETE', 'query_independent_controls': True})
    with (out / 'per_query_metrics.jsonl').open('w') as h:
        for row in allrows: h.write(json.dumps(row, allow_nan=False) + '\n')
    write_json(out / 'block_receipts.json', receipts)
    write_json(out / 'evaluation.json', {'created_utc': now(), 'status': 'COMPLETE_QC_PENDING', 'aggregate': summarize(allrows),
               'metric_rows': len(allrows), 'core_sequences': len(seqs), 'primary_sequences': len(selected),
               'primary_availability_fraction': len(selected) / len(seqs), 'primary_protein_groups': len({axes['protein_groups'][q] for q in selected}),
               'eligible_outer_blocks': sum(r['status'] == 'COMPLETE' for r in receipts), 'foldseek_alignment_audit': counts,
               'independent_biological_validation': False, 'distant_homology_advantage_established': False,
               'limitation': 'Availability-matched development subset; native tool domains differ; simple intensity controls do not eliminate study bias'})
    names = ['run_structure_matched_v2.py', 'STRUCTURE_MATCHED_ANALYSIS_V1.md', 'run_main_baselines.py', 'dataset_02/core_edges.json',
             'dataset_02/core_reactions.json', 'multiaxis_split_01/outer_inner_blocks.json', 'multiaxis_split_01/axis_assignments.json',
             'structure_mapping_audit_02/primary_representatives.json', 'structure_mapping_audit_02/availability_census.json',
             'structure_association_validation_02/INDEPENDENT_QC.json', 'foldseek_raw_01/all_vs_all.tsv',
             'esm_global_01/global_features.npz', 'main_baselines_01/fresh_representations.npz']
    write_json(out / 'input_manifest.json', {n: digest_file(root / n) for n in names})
    print(json.dumps({'status': 'COMPLETE_QC_PENDING', 'metric_rows': len(allrows), 'foldseek_audit': counts}, indent=2))


if __name__ == '__main__':
    run(Path(__file__).resolve().parent)
