"""Separate reconciliation of P02 artifacts, without calling either dataset builder."""
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
from run_raw import now, write_json, digest_file

root = Path(__file__).resolve().parent
output = root / 'p02_validation_01'
if output.exists(): raise FileExistsError(f'Refusing to replace {output}')
output.mkdir()


def read(name): return json.loads((root / name).read_text(encoding='utf-8'))
def lines(name): return [json.loads(x) for x in (root / name).read_text(encoding='utf-8').splitlines()]


checks = {}
for folder in ['dataset_02', 'named_support_01']:
    frozen = read(folder + '/output_checksums.json')
    checks[folder + '_output_hashes'] = all(digest_file(root / folder / f) == sha for f, sha in frozen.items())
    inputs = read(folder + '/input_manifest.json')
    checks[folder + '_input_hashes'] = all(digest_file(Path(f)) == sha for f, sha in inputs.items())
core = read('dataset_02/dataset_audit.json')
edges = read('dataset_02/core_edges.json')
chemistry = read('dataset_02/core_reactions.json')
membership = lines('dataset_02/membership.jsonl')
projections = read('dataset_02/full_reaction_projections.json')
named = lines('named_support_01/named_context_resolution.jsonl')
named_audit = read('named_support_01/named_support_audit.json')
members = {r['assertion_id']: r for r in membership}
edge_keys = {(r['sequence_sha256'], r['reaction_key']) for r in edges}
checks['unique_main_edges_and_membership_ids'] = len(edge_keys) == len(edges) and len(members) == len(membership)
checks['edge_membership_links_reconcile'] = Counter(a for e in edges for a in e['assertion_ids']) == Counter(
    m['assertion_id'] for m in membership if m['eligible'])
checks['edge_provenance_and_publications_reconcile'] = all(
    e['publication_ids'] == sorted({p for a in e['assertion_ids'] for p in members[a]['admissible_publication_ids']})
    and all(members[a]['sequence_sha256'] == e['sequence_sha256'] and members[a]['main_reaction_keys'] == [e['reaction_key']]
            for a in e['assertion_ids']) for e in edges)
checks['main_keys_recompute_without_builder'] = all(k == 'MAIN:' + hashlib.sha256(json.dumps(
    {'substrate': v['substrates'], 'product': v['products']}, sort_keys=True).encode()).hexdigest() for k, v in chemistry.items())
checks['full_component_multiset_independently_reconciled'] = all(
    Counter(p['full_reaction'][side + 's']) == Counter(p[side + 's']) + Counter(
        x['smiles'] for x in p['removed_participants'] if x['side'] == side)
    for p in projections.values() for side in ['substrate', 'product'])
fasta = (root / 'dataset_02/core_sequences.fasta').read_text().splitlines()
sequence_map = dict(zip((h[1:] for h in fasta[::2]), fasta[1::2]))
checks['fasta_hashes_exact_and_all_edges_covered'] = all(hashlib.sha256(s.encode()).hexdigest() == key for key, s in sequence_map.items()) and set(sequence_map) == {e['sequence_sha256'] for e in edges}
checks['reported_counts_reconciled'] = (core['unique_edges'] == len(edges) and core['labelled_sequences'] == len(sequence_map)
    and core['directed_main_transformations'] == len(chemistry)
    and core['single_pair_transformations'] == sum(len(r['substrates']) == len(r['products']) == 1 for r in chemistry.values()))
dbpath = root / 'run_03/raw_rebuild.sqlite'
db = sqlite3.connect(dbpath.resolve().as_uri() + '?mode=ro', uri=True)
db.row_factory = sqlite3.Row
checks['raw_database_integrity'] = db.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
actual_counts = Counter(r['source'] for r in membership)
checks['all_source_denominators_match_raw_database'] = all(db.execute('SELECT COUNT(*) FROM assertions WHERE source=?', (s,)).fetchone()[0] == n for s, n in actual_counts.items())
family = {r[0] for r in db.execute("SELECT sequence_sha256 FROM protein_assertions WHERE family_evidence='PF00067'")}
checks['all_model_sequences_raw_family_confirmed'] = set(sequence_map) <= family
conflicts = lines('dataset_02/quarantined_component_name_conflicts.jsonl')
bad_locators = {(r['source'], r['locator']) for r in conflicts}
checks['source_name_structure_conflicts_quarantined'] = len(conflicts) == core['quarantined_helper_name_structure_components'] and all(not r['eligible'] for r in membership if (r['source'], r['source_locator']) in bad_locators)
initial = {r['assertion_id']: r for r in lines('dataset_01/membership.jsonl')}
checks['P450_initial_admission_not_relaxed'] = all(initial[r['assertion_id']]['eligible'] for r in membership if r['source'] == 'P450Rdb' and r['eligible'])
checks['no_independence_or_construct_claim'] = core['independent_test_records_certified'] == 0 and all(not r['assayed_construct_verified'] for r in edges)
source_full, source_main = defaultdict(set), defaultdict(set)
for r in membership:
    if r['eligible']:
        for key in r['full_reaction_keys']: source_full[r['source']].add((r['sequence_sha256'], key))
        for key in r['main_reaction_keys']: source_main[r['source']].add((r['sequence_sha256'], key))
overlap = {}
for source in ['P450Rdb', 'SwissProt', 'UniProt_API']:
    overlap[source] = {'eligible_assertions': sum(r['eligible'] for r in membership if r['source'] == source),
                       'full_component_unique_edges': len(source_full[source]), 'main_unique_edges': len(source_main[source])}
cross = {'admitted_assertions_identical_for_both_representations': True,
    'full_component_shared_P450_UniProt_edges': len(source_full['P450Rdb'] & source_full['UniProt_API']),
    'main_shared_P450_UniProt_edges': len(source_main['P450Rdb'] & source_main['UniProt_API']),
    'source_denominators': overlap,
    'main_shared_fraction_of_P450_edges': len(source_main['P450Rdb'] & source_main['UniProt_API']) / len(source_main['P450Rdb']),
    'main_shared_fraction_of_UniProt_edges': len(source_main['P450Rdb'] & source_main['UniProt_API']) / len(source_main['UniProt_API']),
    'this_is_representation_linkage_not_model_improvement_or_new_independent_experiments': True}
checks['uniprot_exports_not_independent'] = source_main['SwissProt'] == source_main['UniProt_API'] and all(
    e['source_lineages'] == sorted({'P450Rdb' if s == 'P450Rdb' else 'UniProt' for s in e['sources']}) for e in edges)
checks['named_support_never_adds_labels'] = all(not r['new_model_label_generated'] for r in named) and named_audit['exact_model_positives_added'] == 0
checks['named_qualified_matches_link_existing_edges'] = all((r['matched_core_edge']['sequence_sha256'], r['matched_core_edge']['reaction_key']) in edge_keys for r in named if r['qualified_named_context_match'])
checks['named_source_counts_recomputed'] = dict(Counter(r['source'] for r in named)) == {s: v['source_assertions'] for s, v in named_audit['source_counts'].items()}
checks['named_qualified_counts_recomputed'] = all(sum(r['qualified_named_context_match'] for r in named if r['source'] == s) == v['qualified_named_context_matches_not_new_labels'] for s, v in named_audit['source_counts'].items())
db.close()
tests = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(root), '-p', 'test_*.py', '-v'],
    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
(output / 'test_output.log').write_bytes(tests.stdout)
checks['all_regression_tests_pass'] = tests.returncode == 0
receipt = {'status': 'P02_VIEW_CHECKS_PASS' if all(checks.values()) else 'FAIL', 'created_utc': now(),
    'checks': checks, 'cross_source_representation_audit': cross,
    'scope': 'Executed source-linked MAIN development view and conservative named-context audit',
    'BRENDA_SABIO_semantic_integration_complete': False, 'exact_construct_truth_set_complete': False,
    'new_grouping_or_model_scores_computed': False, 'cluster_sync_completed': False,
    'selected_joint_predictor_complete': False, 'independent_validation_complete': False,
    'whole_paper_complete': False,
    'source_code_sha256': {p.name: digest_file(p) for p in sorted(root.glob('*.py'))},
    'raw_store_sha256': digest_file(dbpath)}
write_json(output / 'P02_RECEIPT.json', receipt)
print(json.dumps(receipt, indent=2))
if not all(checks.values()): raise SystemExit(1)
