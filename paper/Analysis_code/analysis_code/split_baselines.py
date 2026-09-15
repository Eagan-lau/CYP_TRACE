"""Fresh grouping and common-candidate baselines on the rebuilt development core."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from run_raw import stable, write_json, digest_file, now


class Union:
    def __init__(self, ids): self.parent = {key: key for key in ids}
    def root(self, key):
        while self.parent[key] != key:
            self.parent[key] = self.parent[self.parent[key]]
            key = self.parent[key]
        return key
    def join(self, a, b):
        a, b = self.root(a), self.root(b)
        if a != b: self.parent[max(a, b)] = min(a, b)
    def groups(self):
        buckets = defaultdict(list)
        for key in sorted(self.parent): buckets[self.root(key)].append(key)
        return {stable('|'.join(values)): values for values in buckets.values()}


def assign_folds(groups, k):
    counts = [0] * k
    result = {}
    for group, members in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
        fold = min(range(k), key=lambda i: (counts[i], i))
        result[group] = fold
        counts[fold] += len(members)
    return result


def rank_metrics(scores, mask, truth, tie_order, ks=(1, 5, 10)):
    truth = np.asarray(sorted(set(truth)), dtype=int)
    covered = np.where(mask)[0]
    order = covered[np.lexsort((tie_order[covered], -scores[covered]))]
    ranks = np.full(len(mask), np.inf)
    ranks[order] = np.arange(1, len(order) + 1)
    true_ranks = ranks[truth]
    best = float(np.min(true_ranks)) if len(truth) else np.inf
    result = {'rr': 0.0 if not np.isfinite(best) else 1.0 / best,
              'positive_count': int(len(truth)), 'covered_positives': int(np.sum(mask[truth])),
              'covered_candidates': int(len(covered)), 'candidate_count': int(len(mask)),
              'conditional_rr': None if not np.isfinite(best) else 1.0 / best}
    for k in ks:
        result[f'recall_at_{k}'] = float(np.mean(true_ranks <= k)) if len(truth) else None
    return result


def chemical_arrays(chemistry, reaction_ids, spec):
    width = spec['fingerprint']['bits_per_side']
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=spec['fingerprint']['radius'],
        fpSize=width, includeChirality=spec['fingerprint']['include_chirality'])
    matrix = np.zeros((len(reaction_ids), 2 * width), dtype=np.float64)
    scaffolds = {}
    for i, rid in enumerate(reaction_ids):
        keys = set()
        for side_i, side in enumerate(('substrates', 'products')):
            side_smiles = chemistry[rid][side]
            mol = Chem.MolFromSmiles('.'.join(side_smiles))
            if mol is None: raise ValueError(f'Core contains invalid chemistry: {rid}')
            fp = generator.GetFingerprint(mol)
            array = np.zeros(width, dtype=np.uint8)
            DataStructs.ConvertToNumpyArray(fp, array)
            matrix[i, side_i * width:(side_i + 1) * width] = array
            for fragment in Chem.GetMolFrags(mol, asMols=True):
                if not any(atom.GetAtomicNum() == 6 for atom in fragment.GetAtoms()): continue
                scaffold = MurckoScaffold.GetScaffoldForMol(fragment)
                # Do not collapse all acyclic compounds to one empty scaffold.
                key = ('ring:' + Chem.MolToSmiles(scaffold, isomericSmiles=False)) if scaffold.GetNumAtoms() else (
                    'acyclic_connectivity:' + Chem.MolToSmiles(fragment, isomericSmiles=False))
                keys.add(key)
        scaffolds[rid] = sorted(keys)
    intersections = matrix @ matrix.T
    sums = matrix.sum(axis=1)
    union = sums[:, None] + sums[None, :] - intersections
    kernel = np.divide(intersections, union, out=np.zeros_like(intersections), where=union > 0)
    return kernel, scaffolds


def load_similarity(path, sequence_ids, spec, edges):
    index = {key: i for i, key in enumerate(sequence_ids)}
    matrices = {name: np.zeros((len(index), len(index))) for name in ('bits', 'identity', 'bilateral_coverage')}
    identities = [spec['primary_sequence_identity']] + spec['sensitivity_sequence_identities']
    unions = {threshold: Union(sequence_ids) for threshold in identities}
    by_accession = defaultdict(set)
    for edge in edges:
        for acc in edge['accessions']: by_accession[acc].add(edge['sequence_sha256'])
    for values in by_accession.values():
        values = sorted(values)
        for uf in unions.values():
            for value in values[1:]: uf.join(values[0], value)
    rows = 0
    with path.open() as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if len(row) != 8: raise ValueError('MMseqs output schema changed')
            query, target = row[:2]
            identity, qcov, tcov, alnlen, evalue, bits = map(float, row[2:])
            if not (0 <= identity <= 1 and 0 <= qcov <= 1 and 0 <= tcov <= 1):
                raise ValueError('fident/qcov/tcov must be fractions, not percentages')
            if query not in index or target not in index: raise ValueError('Unexpected sequence ID in search result')
            i, j = index[query], index[target]
            rows += 1
            if qcov >= spec['grouping_minimum_bilateral_coverage'] and tcov >= spec['grouping_minimum_bilateral_coverage']:
                for threshold, uf in unions.items():
                    if identity >= threshold: uf.join(query, target)
            if min(qcov, tcov) >= 0.8:
                matrices['identity'][i, j] = max(matrices['identity'][i, j], identity)
            matrices['bilateral_coverage'][i, j] = max(matrices['bilateral_coverage'][i, j], min(qcov, tcov))
            if i != j and evalue <= spec['homology_maximum_evalue'] and min(qcov, tcov) >= spec['homology_minimum_bilateral_coverage']:
                matrices['bits'][i, j] = max(matrices['bits'][i, j], bits)
    return matrices, {threshold: uf.groups() for threshold, uf in unions.items()}, rows


def aggregate(rows, spec):
    results = {}
    rng = np.random.default_rng(spec['bootstrap_seed'])
    for task in sorted({r['task'] for r in rows}):
        results[task] = {}
        for method in spec['methods']:
            selected = [r for r in rows if r['task'] == task and r['method'] == method]
            if not selected: continue
            groups = defaultdict(list)
            for r in selected: groups[r['sequence_group']].append(r)
            result = {'queries': len(selected), 'sequence_groups': len(groups)}
            for metric in ('rr', 'recall_at_1', 'recall_at_5', 'recall_at_10'):
                vals = np.array([np.mean([r[metric] for r in members]) for members in groups.values()])
                result['group_macro_' + metric] = float(vals.mean())
                if metric == 'rr':
                    boot = np.mean(vals[rng.integers(0, len(vals), (spec['bootstrap_replicates'], len(vals)))], axis=1)
                    result['descriptive_rr_interval'] = [float(x) for x in np.quantile(boot, [.025, .975])]
            result['micro_positive_coverage'] = sum(r['covered_positives'] for r in selected) / sum(r['positive_count'] for r in selected)
            result['mean_candidate_coverage'] = float(np.mean([r['covered_candidates'] / r['candidate_count'] for r in selected]))
            results[task][method] = result
    return results


def main(dataset, similarity, output):
    if output.exists(): raise FileExistsError(f'Refusing to overwrite {output}')
    output.mkdir(parents=True)
    spec_path = Path(__file__).with_name('baseline_spec.json')
    spec = json.loads(spec_path.read_text())
    dataset_audit = json.loads((dataset / 'dataset_audit.json').read_text())
    if not dataset_audit['ready_for_split_feasibility']: raise RuntimeError('Dataset audit is not ready')
    edges = json.loads((dataset / 'core_edges.json').read_text())
    chemistry = json.loads((dataset / 'core_reactions.json').read_text())
    sequence_ids = sorted({e['sequence_sha256'] for e in edges})
    reaction_ids = sorted(chemistry)
    si, ri = {s: i for i, s in enumerate(sequence_ids)}, {r: i for i, r in enumerate(reaction_ids)}
    matrices, grouping, search_rows = load_similarity(similarity / 'all_vs_all.tsv', sequence_ids, spec, edges)
    groups = grouping[spec['primary_sequence_identity']]
    write_json(output / 'sequence_groups.json', {str(threshold): mapping for threshold, mapping in grouping.items()})
    group_index = {s: group for group, members in groups.items() for s in members}
    group_count = len(groups)
    if group_count < spec['minimum_sequence_groups']:
        write_json(output / 'status.json', {'status': 'INSUFFICIENT_INDEPENDENT_SEQUENCE_GROUPS', 'groups': group_count,
                   'threshold_not_relaxed': True, 'model_performance_evaluated': False})
        return
    k = min(spec['maximum_outer_folds'], group_count)
    fold_assignment = assign_folds(groups, k)
    protein_fold = {s: fold_assignment[group_index[s]] for s in sequence_ids}
    kernel, scaffold_keys = chemical_arrays(chemistry, reaction_ids, spec)
    cu = Union(reaction_ids)
    first = {}
    for rid in reaction_ids:
        for key in scaffold_keys[rid]:
            if key in first: cu.join(rid, first[key])
            else: first[key] = rid
    chemical_groups = cu.groups()
    chemical_assign = assign_folds(chemical_groups, k)
    chemical_fold = {r: chemical_assign[g] for g, members in chemical_groups.items() for r in members}
    fold_plan = {'protein_folds': protein_fold, 'protein_groups': group_index, 'reaction_scaffold_keys': scaffold_keys,
                 'chemical_component_folds': chemical_fold, 'chemical_component_sizes': sorted([len(v) for v in chemical_groups.values()], reverse=True),
                 'folds': k, 'spec_sha256': digest_file(spec_path), 'dataset_manifest_sha256': digest_file(dataset / 'input_manifest.json')}
    write_json(output / 'split_plan.json', fold_plan)
    np.savez_compressed(output / 'fresh_chemical_kernel.npz', kernel=kernel)
    tie_hashes = sorted(range(len(reaction_ids)), key=lambda i: stable(reaction_ids[i]))
    tie = np.empty(len(reaction_ids), dtype=int)
    tie[tie_hashes] = np.arange(len(reaction_ids))
    alltruth = defaultdict(set)
    for e in edges: alltruth[e['sequence_sha256']].add(ri[e['reaction_key']])
    rows, audits = [], []
    for strict in (False, True):
        for fold in range(k):
            queries = [s for s in sequence_ids if protein_fold[s] == fold]
            query_pubs = {p for e in edges if e['sequence_sha256'] in set(queries) for p in e['publication_ids']}
            train_edges = [e for e in edges if protein_fold[e['sequence_sha256']] != fold
                           and not (set(e['publication_ids']) & query_pubs)
                           and (not strict or chemical_fold[e['reaction_key']] != fold)]
            Y = np.zeros((len(sequence_ids), len(reaction_ids)), dtype=np.float64)
            train_pubs = set()
            for e in train_edges:
                Y[si[e['sequence_sha256']], ri[e['reaction_key']]] = 1.0
                train_pubs.update(e['publication_ids'])
            train_sequences = np.where(Y.sum(axis=1) > 0)[0]
            query_positions = [si[s] for s in queries]
            train_reactions = set(np.where(Y.sum(axis=0) > 0)[0])
            cross = matrices['identity'][np.ix_(query_positions, train_sequences)]
            cross_reverse = matrices['identity'][np.ix_(train_sequences, query_positions)]
            maximum = max(float(cross.max()) if cross.size else 0, float(cross_reverse.max()) if cross_reverse.size else 0)
            if maximum >= spec['primary_sequence_identity']:
                raise AssertionError('Cross-fold sequence similarity threshold violated')
            if train_pubs & query_pubs: raise AssertionError('Recorded publication leakage')
            audits.append({'track': 'strict_double_cold' if strict else 'protein_cold', 'fold': fold,
                'query_sequences': len(queries), 'training_sequences': len(train_sequences),
                'training_edges': len(train_edges), 'training_reactions': len(train_reactions),
                'max_search_detected_bilateral_identity': maximum, 'recorded_publication_overlap': 0,
                'doi_pmid_alias_completeness_certified': False})
            if not train_edges:
                audits[-1]['status'] = 'NO_TRAINING_EDGES_AFTER_PURGE'
            frequency = Y.sum(axis=0)
            prior = frequency @ kernel
            H = matrices['bits'][query_positions] @ Y
            transport = H @ kernel
            track = 'strict_double_cold' if strict else 'protein_cold'
            np.savez_compressed(output / f'scores_{track}_{fold}.npz',
                query_ids=np.asarray(queries), reaction_ids=np.asarray(reaction_ids),
                frequency=frequency, chemistry_prior=prior, homology=H, transport=transport,
                training_edge_count=np.asarray([len(train_edges)]))
            for qn, s in enumerate(queries):
                truth = alltruth[s]
                if strict:
                    tasks = {'strict_double_cold': [r for r in truth if chemical_fold[reaction_ids[r]] == fold]}
                else:
                    tasks = {'protein_cold_all': sorted(truth),
                             'protein_cold_seen_reaction': sorted(truth & train_reactions),
                             'protein_cold_unseen_reaction': sorted(truth - train_reactions)}
                score_by_method = {'uniform': np.ones(len(reaction_ids)), 'training_frequency': frequency,
                                   'chemistry_prior': prior, 'mmseqs_transfer': H[qn],
                                   'homology_chemical_transport': transport[qn]}
                for task, task_truth in tasks.items():
                    if not task_truth: continue
                    for method, scores in score_by_method.items():
                        if not train_edges and method != 'uniform': mask = np.zeros(len(reaction_ids), dtype=bool)
                        elif method == 'mmseqs_transfer': mask = H[qn] > 0
                        elif method == 'homology_chemical_transport': mask = np.full(len(reaction_ids), bool(H[qn].sum() > 0))
                        else: mask = np.ones(len(reaction_ids), dtype=bool)
                        metrics = rank_metrics(scores, mask, task_truth, tie)
                        rows.append({'task': task, 'fold': fold, 'sequence_sha256': s, 'sequence_group': group_index[s],
                                     'method': method, **metrics})
            if train_edges: audits[-1]['status'] = 'EVALUATED'
            print(f'{audits[-1]["track"]} fold {fold}: {len(queries)} queries; {len(train_edges)} training edges', flush=True)
    with (output / 'per_query_metrics.jsonl').open('w') as stream:
        for row in rows: stream.write(json.dumps(row, allow_nan=False) + '\n')
    write_json(output / 'fold_audits.json', audits)
    write_json(output / 'evaluation.json', {'created_utc': now(), 'aggregate': aggregate(rows, spec),
        'sequences': len(sequence_ids), 'reactions': len(reaction_ids), 'positive_edges': len(edges),
        'sequence_groups': group_count, 'chemical_components': len(chemical_groups), 'search_rows': search_rows,
        'initial_core_only': True, 'independent_validation': False,
        'selected_joint_model_evaluated': False, 'paper_complete': False,
        'spec': spec, 'spec_sha256': digest_file(spec_path), 'script_sha256': digest_file(Path(__file__)),
        'input_hashes': {str(p): digest_file(p) for p in [dataset / 'core_edges.json', dataset / 'core_reactions.json', similarity / 'all_vs_all.tsv']},
        'limitations': ['Source-limited development cohort, not complete integrated paper dataset',
                        'Catalogue is supplied retrospectively; out-of-catalogue discovery not assessed',
                        'Recorded publication purge does not certify all DOI/PMID aliases',
                        'No selected or calibrated joint model in this stage',
                        'Acyclic coldness uses full connectivity, not a generalized acyclic scaffold family',
                        'Search-detected similarity isolation is not exhaustive alignment certification']})
    write_json(output / 'status.json', {'status': 'DEVELOPMENT_BASELINES_COMPLETE', 'paper_complete': False,
                                       'selected_joint_model_trained': False, 'independent_validation': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', type=Path, required=True)
    parser.add_argument('--similarity', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.dataset, args.similarity, args.output)
