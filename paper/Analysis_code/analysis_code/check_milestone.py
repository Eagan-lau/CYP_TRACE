"""Executable receipt for this completed foundation, not whole-paper acceptance."""
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys
from run_raw import digest_file, write_json, now

root = Path(__file__).resolve().parent
def read(name): return json.loads((root / name).read_text(encoding='utf-8'))
checks = {}
core = read('dataset_01/dataset_audit.json')
checks['initial_core_membership_checks'] = all(core['checks'].values())
metrics = read('baselines_01/INDEPENDENT_METRIC_QC.json')
checks['independent_metric_arithmetic'] = metrics['status'] == 'PASS' and not metrics['failures']
human = read('human_substrate_01/human_substrate_audit.json')
checks['human_label_mechanical_checks'] = all(human['checks'].values())
geometry = read('structure_assets_01/structure_audit.json')
checks['geometry_no_parse_errors'] = not geometry['errors']
checks['geometry_rows_reconcile'] = sum(1 for _ in (root / 'structure_assets_01/heme_contact_residues.jsonl').open()) == geometry['counts']['heme_contact_residue_rows']
checks['geometry_chains_reconcile'] = sum(1 for _ in (root / 'structure_assets_01/chain_assets.jsonl').open()) == geometry['counts']['resolved_protein_chains']
brenda = read('brenda_01/brenda_audit.json')
rows = [json.loads(line) for line in (root / 'brenda_01/brenda_assertions.jsonl').read_text(encoding='utf-8').splitlines()]
checks['brenda_rows_reconcile'] = len(rows) == brenda['counts']['context_assertions_retained']
checks['brenda_scopes_reconcile'] = dict(Counter(r['scope'] for r in rows)) == brenda['scopes']
checks['brenda_no_broken_protein_or_literature_ids'] = all(not r['missing_protein_ids'] and not r['missing_reference_ids'] for r in rows)
checks['brenda_no_automatic_exact_positives'] = not any(r['eligible_exact_reaction_positive'] for r in rows)
annotation = read('annotation_02/annotation_resolution_audit.json')
checks['database_name_not_evidence_independence'] = all(not link['cross_database_evidence_independence_certified'] and link['independent_sources'] is not True
                                                      for link in annotation['exact_sequence_cross_source_links'])
tests = subprocess.run([sys.executable, '-m', 'unittest', 'discover', '-s', str(root), '-p', 'test_*.py'],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
(root / 'milestone_test_output.log').write_bytes(tests.stdout)
checks['automated_tests_exit_zero'] = tests.returncode == 0
files = ['run_03/raw_rebuild.sqlite', 'dataset_01/core_edges.json', 'dataset_01/core_reactions.json',
         'baselines_01/evaluation.json', 'baselines_01/split_plan.json', 'baselines_01/INDEPENDENT_METRIC_QC.json',
         'human_substrate_01/normalized_labels.json', 'structure_assets_01/structure_audit.json',
         'brenda_01/brenda_audit.json', 'annotation_02/annotation_resolution_audit.json']
receipt = {'created_utc': now(), 'status': 'FOUNDATION_CHECKS_PASS' if all(checks.values()) else 'FOUNDATION_CHECKS_FAIL',
           'checks': checks, 'artifact_sha256': {f: digest_file(root / f) for f in files},
           'receipt_script_sha256': digest_file(Path(__file__)),
           'acceptance_scope': 'Raw-derived foundation, initial source-limited dataset and baseline arithmetic only',
           'all_sources_semantically_integrated': False, 'selected_joint_model_complete': False,
           'independent_biological_validation': False, 'whole_paper_complete': False}
write_json(root / 'MILESTONE_RECEIPT.json', receipt)
print(json.dumps(receipt, indent=2))
if not all(checks.values()): raise SystemExit(1)
