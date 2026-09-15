"""Independent metric reconstruction from fresh score archives (not scorer imports)."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import numpy as np


root = Path(__file__).resolve().parent
folder = root / 'baselines_01'
evaluation = json.loads((folder / 'evaluation.json').read_text())
plan = json.loads((folder / 'split_plan.json').read_text())
edges = json.loads((root / 'dataset_01/core_edges.json').read_text())
rows = [json.loads(line) for line in (folder / 'per_query_metrics.jsonl').read_text().splitlines()]
by_sequence = defaultdict(set)
for edge in edges: by_sequence[edge['sequence_sha256']].add(edge['reaction_key'])
archives = {}
failures = []
groups = defaultdict(lambda: defaultdict(list))
max_error = 0.
for row in rows:
    strict = row['task'] == 'strict_double_cold'
    track = 'strict_double_cold' if strict else 'protein_cold'
    key = (track, row['fold'])
    if key not in archives:
        archives[key] = np.load(folder / f'scores_{track}_{row["fold"]}.npz', allow_pickle=False)
    arrays = archives[key]
    candidates = arrays['reaction_ids'].tolist()
    queries = arrays['query_ids'].tolist()
    qi = queries.index(row['sequence_sha256'])
    seen = {candidates[i] for i, value in enumerate(arrays['frequency']) if value > 0}
    truth = by_sequence[row['sequence_sha256']]
    if row['task'].endswith('_seen_reaction'): truth = truth & seen
    if row['task'].endswith('_unseen_reaction'): truth = truth - seen
    if strict: truth = {r for r in truth if plan['chemical_component_folds'][r] == row['fold']}
    method = row['method']
    has_train = arrays['training_edge_count'][0] > 0
    if method == 'uniform': scores, present = np.ones(len(candidates)), set(candidates)
    elif method == 'training_frequency': scores, present = arrays['frequency'], set(candidates) if has_train else set()
    elif method == 'chemistry_prior': scores, present = arrays['chemistry_prior'], set(candidates) if has_train else set()
    elif method == 'mmseqs_transfer':
        scores = arrays['homology'][qi]
        present = {candidates[i] for i, value in enumerate(scores) if value > 0}
    else:
        scores = arrays['transport'][qi]
        present = set(candidates) if sum(arrays['homology'][qi]) > 0 else set()
    index = {r: i for i, r in enumerate(candidates)}
    ordered = sorted(present, key=lambda r: (-float(scores[index[r]]), hashlib.sha256(r.encode()).hexdigest()))
    positions = {r: rank + 1 for rank, r in enumerate(ordered)}
    recovered = [positions[r] for r in truth if r in positions]
    expected = {'rr': 1 / min(recovered) if recovered else 0,
                'positive_count': len(truth), 'covered_positives': len(recovered),
                'covered_candidates': len(present), 'candidate_count': len(candidates)}
    for k in (1, 5, 10): expected[f'recall_at_{k}'] = sum(p <= k for p in recovered) / len(truth)
    for name, value in expected.items():
        error = abs(row[name] - value)
        max_error = max(max_error, error)
        if error > 1e-12: failures.append((row['task'], row['sequence_sha256'], method, name))
    groups[(row['task'], method)][row['sequence_group']].append(row['rr'])
for (task, method), values in groups.items():
    rebuilt = sum(sum(v) / len(v) for v in values.values()) / len(values)
    reported = evaluation['aggregate'][task][method]['group_macro_rr']
    if abs(rebuilt - reported) > 1e-12: failures.append(('aggregate', task, method))
audit = {'per_query_method_rows_recomputed': len(rows), 'group_macro_mrrs_recomputed': len(groups),
         'maximum_numeric_difference': max_error, 'failures': failures,
         'status': 'PASS' if not failures else 'FAIL', 'independent_implementation': True,
         'independent_biological_validation': False, 'whole_paper_complete': False,
         'verification_scope': 'Metric arithmetic, coverage masks and group aggregation; not proof of source identity, model validity or full leakage elimination'}
(folder / 'INDEPENDENT_METRIC_QC.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps(audit, indent=2))
if failures: raise SystemExit(1)
