"""Adjudicate initial raw assertions and prepare a new, inspectable core.

No old integrated table, label, split or model is an input. This first exact
reaction cohort is source-limited to P450Rdb while other sources are adjudicated.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sqlite3

from run_raw import digest_file, stable, write_json, now


def build(raw_run: Path, output: Path):
    if output.exists():
        raise FileExistsError(f'Refusing to replace {output}')
    status = json.loads((raw_run / 'status.json').read_text())
    if status['status'] != 'RAW_EXTRACTION_COMPLETE_DOWNSTREAM_PENDING':
        raise RuntimeError('Raw extraction did not complete successfully')
    output.mkdir(parents=True)
    db = sqlite3.connect((raw_run / 'raw_rebuild.sqlite').resolve().as_uri() + '?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    accessions = defaultdict(set)
    exact_family = set()
    accession_taxa = defaultdict(set)
    for row in db.execute('SELECT accession,sequence_sha256,taxid,family_evidence FROM protein_assertions'):
        if not row['sequence_sha256']: continue
        for acc in row['accession'].split('|'):
            if acc:
                accessions[acc].add(row['sequence_sha256'])
                if row['taxid']: accession_taxa[acc].add(row['taxid'])
        if row['family_evidence'] == 'PF00067': exact_family.add(row['sequence_sha256'])
    conflicting_accessions = {acc: sorted(values) for acc, values in accessions.items() if len(values) > 1}
    write_json(output / 'accession_sequence_disagreements.json', conflicting_accessions)
    edge_assertions = defaultdict(list)
    exclusions = Counter()
    membership = []
    for row in db.execute("SELECT * FROM assertions WHERE source='P450Rdb' ORDER BY assertion_id"):
        reasons = json.loads(row['exclusion_reasons'])
        if row['accession'] in conflicting_accessions:
            reasons.append('accession_sequence_disagreement_requires_adjudication')
        if len(accession_taxa.get(row['accession'], set())) > 1:
            reasons.append('accession_taxonomy_disagreement_requires_adjudication')
        if row['sequence_sha256'] not in exact_family:
            reasons.append('exact_sequence_PF00067_confirmation_unavailable')
        reasons = sorted(set(reasons))
        exclusions.update(reasons)
        membership.append({'assertion_id': row['assertion_id'], 'eligible': not reasons, 'reasons': reasons,
                           'sequence_sha256': row['sequence_sha256'], 'reaction_key': row['reaction_key']})
        if not reasons: edge_assertions[(row['sequence_sha256'], row['reaction_key'])].append(dict(row))
    with (output / 'membership.jsonl').open('w', encoding='utf-8') as stream:
        for row in membership: stream.write(json.dumps(row, sort_keys=True) + '\n')
    edges = []
    for (seq, reaction), rows in sorted(edge_assertions.items()):
        pubs = sorted({pub for row in rows for pub in json.loads(row['publication_ids'])})
        edges.append({'sequence_sha256': seq, 'reaction_key': reaction,
                      'publication_ids': pubs, 'assertion_ids': sorted(row['assertion_id'] for row in rows),
                      'accessions': sorted({r['accession'] for r in rows if r['accession']}),
                      'taxids': sorted({r['taxid'] for r in rows}),
                      'assayed_construct_verified': False, 'source_scope': 'P450Rdb_development_core',
                      'historical_exposure': 'DEVELOPMENT_EXPOSED_OR_UNVERIFIED'})
    write_json(output / 'core_edges.json', edges)
    sequences = sorted({e['sequence_sha256'] for e in edges})
    reactions = sorted({e['reaction_key'] for e in edges})
    proteins = {row['sequence_sha256']: row['sequence'] for row in db.execute('SELECT * FROM proteins')
                if row['sequence_sha256'] in set(sequences)}
    with (output / 'core_sequences.fasta').open('w', encoding='ascii') as stream:
        for seq in sequences: stream.write(f'>{seq}\n{proteins[seq]}\n')
    chemistry = {row['reaction_key']: {'substrates': json.loads(row['substrates']), 'products': json.loads(row['products'])}
                 for row in db.execute('SELECT * FROM reactions') if row['reaction_key'] in set(reactions)}
    write_json(output / 'core_reactions.json', chemistry)
    assertions_by_source = dict(db.execute('SELECT source,COUNT(*) FROM assertions GROUP BY source'))
    checks = {
        'unique_edge_keys': len(edges) == len({(e['sequence_sha256'], e['reaction_key']) for e in edges}),
        'all_edge_sequences_present': all(e['sequence_sha256'] in proteins for e in edges),
        'all_edge_reactions_present': all(e['reaction_key'] in chemistry for e in edges),
        'all_edges_publication_traceable': all(e['publication_ids'] for e in edges),
        'all_edges_exact_sequence_family_confirmed': all(seq in exact_family for seq in sequences),
        'membership_denominator_preserved': len(membership) == assertions_by_source.get('P450Rdb', 0),
        'eligible_assertions_reconcile': sum(len(e['assertion_ids']) for e in edges) == sum(m['eligible'] for m in membership),
        'no_independent_test_claim': True,
    }
    summary = {'created_utc': now(), 'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks,
        'source_assertion_counts': assertions_by_source, 'raw_P450Rdb_assertions': len(membership),
        'initial_P450Rdb_eligible_assertions': sum(m['eligible'] for m in membership),
        'unique_edges': len(edges), 'labelled_sequences': len(sequences), 'directed_reaction_catalogue': len(reactions),
        'nonexclusive_exclusion_counts': dict(exclusions), 'accessions_with_multiple_sequences': len(conflicting_accessions),
        'independent_test_records_certified': 0,
        'cohort_scope': 'initial source-limited family-confirmed development core; other sources pending adjudication',
        'new_raw_build_sha256': digest_file(raw_run / 'raw_rebuild.sqlite'),
        'script_sha256': digest_file(Path(__file__)), 'protocol_sha256': digest_file(Path(__file__).with_name('PROTOCOL.md')),
        'ready_for_split_feasibility': all(checks.values()) and len(sequences) >= 2,
        'model_performance_evaluated': False, 'paper_complete': False}
    write_json(output / 'dataset_audit.json', summary)
    write_json(output / 'input_manifest.json', {p.name: digest_file(p) for p in
        [raw_run / 'raw_rebuild.sqlite', raw_run / 'consumed_raw_inputs.json', Path(__file__).with_name('PROTOCOL.md')]})
    db.close()
    print(json.dumps(summary, indent=2), flush=True)
    if not all(checks.values()): raise RuntimeError('Dataset audit failed')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--raw-run', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    build(a.raw_run, a.output)
