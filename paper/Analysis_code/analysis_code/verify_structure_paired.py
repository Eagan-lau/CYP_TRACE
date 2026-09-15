"""Independent sorted-list ranks and bootstrap multiplicity reconstruction."""
from collections import defaultdict
import json
import math
from pathlib import Path
import numpy as np
from run_raw import digest_file, stable, write_json, now

root = Path(__file__).resolve().parent; result = root / 'structure_paired_01'; out = root / 'structure_paired_validation_01'
if out.exists(): raise FileExistsError(out)
read = lambda p: json.loads(p.read_text())
declared = read(result / 'comparison_summary.json'); failures = []; checked = 0; intervals = 0
for n, h in read(result / 'input_manifest.json').items():
    if digest_file(root / n) != h: failures.append('source:' + n)
rows = [json.loads(s) for s in (result / 'paired_query_records.jsonl').read_text().splitlines()]
orig = [json.loads(s) for s in (root / 'structure_matched_01/per_query_metrics.jsonl').read_text().splitlines()]
source = {(r['method'], r['task'], r['block_number'], r['sequence_sha256']): r for r in orig}
rids = sorted(read(root / 'dataset_02/core_reactions.json')); byblock = defaultdict(list)
for r in rows: byblock[r['block_number']].append(r)
for bn, members in byblock.items():
    d = np.load(root / 'structure_matched_01' / f'block_{bn:03d}_scores.npz', allow_pickle=False)
    for row in members:
        values = []; masks = []; targets = []
        for side in ['left', 'right']:
            name = row[side + '_method']; original = source[(name, row['task'], bn, row['sequence_sha256'])]
            if original['rr'] != row[side + '_rr'] or original['covered_positives'] != row[side + '_covered_positives']: failures.append('end_to_end')
            qi = original['query_index_in_block']; v = d[name]; m = d['domain_' + name]
            values.append(v[qi] if v.ndim == 2 else v); masks.append(m[qi] if m.ndim == 2 else m); targets.append(original['target_reaction_indices'])
        if targets[0] != targets[1]: failures.append('unmatched_target')
        common = [i for i in range(len(rids)) if masks[0][i] and masks[1][i]]; commonset = set(common)
        truth = [i for i in targets[0] if i in commonset]
        if len(common) != row['common_candidate_count'] or len(truth) != row['common_positive_count']: failures.append('common_denominator')
        if bool(truth) != row['common_reranking_defined']: failures.append('undefined_flag')
        for side, value in zip(['left', 'right'], values):
            if not truth:
                if row['common_' + side + '_rr'] is not None: failures.append('undefined_value')
            else:
                order = sorted(common, key=lambda i: (-float(value[i]), stable(rids[i])))
                ranks = {i: j + 1 for j, i in enumerate(order)}
                if abs(1 / min(ranks[i] for i in truth) - row['common_' + side + '_rr']) > 1e-12: failures.append('common_ranking')
        checked += 1
    d.close()
rng = np.random.Generator(np.random.PCG64(declared['bootstrap_seed'])); scopes = 0
for c in declared['comparisons']:
    selected = [r for r in rows if (r['task'], r['left_method'], r['right_method']) == (c['task'], c['left_method'], c['right_method'])]
    for scope, a, b, chosen in [('end_to_end', 'left_rr', 'right_rr', selected),
                               ('common_domain_reranking', 'common_left_rr', 'common_right_rr', [r for r in selected if r['common_reranking_defined']])]:
        groups = defaultdict(lambda: defaultdict(list))
        for r in chosen: groups[r['protein_group']][r['sequence_sha256']].append(r)
        values = []
        for group, queries in sorted(groups.items()):
            values.append([math.fsum(math.fsum(r[k] for r in panels) / len(panels) for panels in queries.values()) / len(queries) for k in [a, b]])
        n = len(values); expected = c[scope]
        if n != expected['protein_groups'] or len(chosen) != expected['query_panels']: failures.append('summary_denominator')
        if n:
            v = np.array(values); delta = v[:, 0] - v[:, 1]
            for k, x in [('left_mrr', math.fsum(v[:, 0]) / n), ('right_mrr', math.fsum(v[:, 1]) / n), ('difference', math.fsum(delta) / n)]:
                if abs(x - expected[k]) > 1e-12: failures.append('summary_mean')
        if n >= 2:
            draws = rng.integers(n, size=(declared['bootstrap_replicates'], n)); counts = np.zeros((len(draws), n), dtype=int)
            np.add.at(counts, (np.repeat(np.arange(len(draws)), n), draws.ravel()), 1)
            boot = np.sort(counts @ delta / n); ci = []
            for p in [.025, .975]:
                i = p * (len(boot) - 1); low = math.floor(i); high = math.ceil(i)
                ci.append(boot[low] + (i - low) * (boot[high] - boot[low]))
            if not np.allclose(ci, expected['conditional_protein_bootstrap_95ci'], rtol=0, atol=1e-12): failures.append('interval')
            intervals += 1
        elif expected['conditional_protein_bootstrap_95ci'] is not None: failures.append('undefined_interval')
        scopes += 1
out.mkdir(); report = {'status': 'PASS' if not failures else 'FAIL', 'created_utc': now(), 'failures': failures,
                      'paired_rows_recomputed': checked, 'summary_scopes_recomputed': scopes, 'conditional_intervals_recomputed': intervals,
                      'verifier_sha256': digest_file(Path(__file__)), 'external_or_full_dependence_inference_certified': False}
write_json(out / 'INDEPENDENT_QC.json', report); print(json.dumps(report, indent=2))
if failures: raise SystemExit(1)
