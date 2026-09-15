"""Fresh source-level resolution and conservative cross-source identity audit."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import sqlite3

from run_raw import digest_file, stable, write_json, now, normalize_structure


def rhea_links(source, raw):
    """Return declared crossrefs; never assume that a master means LR."""
    if source == 'UniProt_API':
        refs = [x.get('id', '') for x in raw.get('reaction', {}).get('reactionCrossReferences', [])
                if x.get('database') == 'Rhea']
        physiological = [x.get('reactionCrossReference', {}).get('id', '')
                         for x in raw.get('physiologicalReactions', [])
                         if x.get('reactionCrossReference', {}).get('database') == 'Rhea']
        clean = lambda values: sorted({v.split(':')[1] for v in values if re.fullmatch(r'RHEA:\d+', v)})
        return clean(refs), clean(physiological)
    if source == 'SwissProt':
        text = raw['comment']
        primary = text.split('PhysiologicalDirection=', 1)[0]
        physiological = text.split('PhysiologicalDirection=', 1)[1] if 'PhysiologicalDirection=' in text else ''
        return sorted(set(re.findall(r'Rhea:RHEA:(\d+)', primary))), sorted(set(re.findall(r'Rhea:RHEA:(\d+)', physiological)))
    return [], []


def exact_reaction(smiles):
    try: left, agents, right = smiles.split('>')
    except ValueError: return {'status': 'invalid_reaction_syntax', 'reaction_key': None}
    if agents: return {'status': 'agents_present_requires_adjudication', 'reaction_key': None}
    sides = {}
    for side, text in [('substrate', left), ('product', right)]:
        if not text: return {'status': 'empty_side', 'reaction_key': None}
        parsed = [normalize_structure(s) for s in text.split('.')]
        statuses = [r['structure_status'] for r in parsed if not r['canonical_smiles']]
        if statuses: return {'status': 'unresolved_component:' + '|'.join(sorted(set(statuses))), 'reaction_key': None}
        sides[side] = sorted(r['canonical_smiles'] for r in parsed)
    return {'status': 'all_components_parseable_not_main_pair_adjudicated',
            'reaction_key': 'RXN:' + stable(json.dumps(sides, sort_keys=True)), 'sides': sides}


def build(raw_run, dataset, output):
    if output.exists(): raise FileExistsError(f'Refusing to overwrite {output}')
    output.mkdir(parents=True)
    dbpath = raw_run / 'raw_rebuild.sqlite'
    db = sqlite3.connect(dbpath.resolve().as_uri() + '?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    sequences_by_accession, source_sequences = defaultdict(set), defaultdict(set)
    family_sequences = set()
    for row in db.execute('SELECT source,accession,sequence_sha256,family_evidence FROM protein_assertions'):
        if not row['sequence_sha256']: continue
        source_sequences[row['source']].add(row['sequence_sha256'])
        for accession in row['accession'].split('|'):
            if accession: sequences_by_accession[accession].add(row['sequence_sha256'])
        if row['family_evidence'] == 'PF00067': family_sequences.add(row['sequence_sha256'])
    rhea = {r['rhea_id']: dict(r) for r in db.execute('SELECT * FROM rhea_structures')}
    masters = {r['master_id'] for r in rhea.values()}
    core = json.loads((dataset / 'core_edges.json').read_text())
    core_keys = {(e['sequence_sha256'], e['reaction_key']) for e in core}
    counts, rows, chemistry_cache = defaultdict(Counter), [], {}
    uniprot_semantics = defaultdict(set)
    for record in db.execute('SELECT * FROM assertions ORDER BY source,assertion_id'):
        source, raw = record['source'], json.loads(record['raw_json'])
        c = counts[source]
        c['source_assertions'] += 1
        accessions = [a for a in record['accession'].split('|') if a]
        identity_candidates = sorted({s for a in accessions for s in sequences_by_accession.get(a, set())})
        sequence = record['sequence_sha256']
        # A unique accession link is reported, not substituted as an assayed construct.
        linked = sequence or (identity_candidates[0] if len(accessions) == 1 and len(identity_candidates) == 1 else '')
        primary, physiological = rhea_links(source, raw)
        declared = sorted(set(primary + physiological))
        directed = [rid for rid in declared if rid in rhea]
        keys = set()
        for rid in directed:
            if rid not in chemistry_cache: chemistry_cache[rid] = exact_reaction(rhea[rid]['reaction_smiles'])
            if chemistry_cache[rid]['reaction_key']: keys.add(chemistry_cache[rid]['reaction_key'])
        if source == 'P450Rdb' and record['reaction_key']: keys.add(record['reaction_key'])
        pubs = json.loads(record['publication_ids'])
        flags = {'source_has_explicit_sequence': bool(sequence),
            'unique_sequence_link_available': bool(linked) and len(identity_candidates) <= 1 and len(accessions) <= 1,
            'linked_sequence_has_PF00067': linked in family_sequences,
            'has_recorded_publication': bool(pubs),
            'declared_Rhea_link_resolves': any(r in rhea or r in masters for r in declared),
            'declared_directional_Rhea_structure_available': bool(directed),
            'complete_component_set_parseable': bool(keys),
            'exact_component_set_edge_matches_initial_core': any((linked, key) in core_keys for key in keys),
            'multiple_accessions': len(accessions) > 1,
            'accession_has_multiple_sequence_versions': len(identity_candidates) > 1}
        c.update({name: int(value) for name, value in flags.items()})
        if source in {'UniProt_API', 'SwissProt'} and sequence:
            semantic_rhea = sorted({rhea[r]['master_id'] if r in rhea else r for r in declared if r in rhea or r in masters})
            for rid in semantic_rhea: uniprot_semantics[source].add((sequence, rid))
        row = {'assertion_id': record['assertion_id'], 'source': source, 'flags': flags,
            'source_sequence_sha256': sequence, 'linked_sequence_sha256': linked,
            'identity_candidates': identity_candidates, 'accessions': accessions,
            'declared_primary_rhea_ids': primary, 'declared_physiological_rhea_ids': physiological,
            'declared_directional_rhea_ids': directed, 'normalized_full_component_reaction_keys': sorted(keys),
            'publication_ids': pubs, 'source_evidence_type': record['evidence_type'],
            'source_locator': record['locator'], 'source_path': record['source_path'], 'source_sha256': record['source_sha256'],
            'exact_assayed_construct_verified': False, 'main_pair_adjudicated': False,
            'new_model_label_generated': False}
        rows.append(row)
    with (output / 'assertion_resolution.jsonl').open('w', encoding='utf-8') as stream:
        for row in rows: stream.write(json.dumps(row, ensure_ascii=False) + '\n')
    for rid, record in chemistry_cache.items():
        record['declared_direction'] = rhea[rid]['direction']
        record['master_id'] = rhea[rid]['master_id']
    write_json(output / 'declared_rhea_chemistry.json', chemistry_cache)
    rates = {s: {k: {'numerator': value, 'denominator': c['source_assertions'], 'fraction': value / c['source_assertions']}
                 for k, value in c.items() if k != 'source_assertions'} for s, c in counts.items()}
    links = []
    for a in sorted(source_sequences):
        for b in sorted(source_sequences):
            if a >= b: continue
            common = source_sequences[a] & source_sequences[b]
            links.append({'source_a': a, 'source_b': b, 'unique_exact_sequences_shared': len(common),
                          'source_a_sequence_denominator': len(source_sequences[a]),
                          'source_b_sequence_denominator': len(source_sequences[b]),
                          'independent_sources': False if {a, b} <= {'UniProt_API', 'SwissProt'} else None,
                          'cross_database_evidence_independence_certified': False})
    api, swiss = uniprot_semantics['UniProt_API'], uniprot_semantics['SwissProt']
    summary = {'created_utc': now(), 'source_counts': {s: dict(c) for s, c in counts.items()}, 'source_assertion_rates': rates,
        'exact_sequence_cross_source_links': links,
        'uniprot_source_lineage_overlap': {'API_sequence_master_Rhea_pairs': len(api), 'SwissProt_sequence_master_Rhea_pairs': len(swiss),
            'shared_pairs': len(api & swiss), 'union_pairs': len(api | swiss), 'direction_not_distinguished_in_this_overlap_count': True,
            'these_exports_are_not_independent_databases': True},
        'status': 'SOURCE_RESOLUTION_AUDIT_COMPLETE_ADJUDICATION_PENDING',
        'fresh_raw_sha256': digest_file(dbpath), 'script_sha256': digest_file(Path(__file__)),
        'independent_validation': False, 'paper_complete': False,
        'limitations': ['No metadata link is proof of catalytic activity, physiological use or exact assayed construct',
            'All-component exact matching is conservative; cofactor, protonation and tautomer differences remain',
            'Absent explicit reaction direction is not imputed from a master Rhea ID',
            'Protein-source inventory includes candidate/name-only entries; PF00067 is reported separately',
            'BRENDA typed assertions and named SABIO chemistry still require separate chemical adjudication',
            'Per-source assertion rates are not independent experimental sample sizes']}
    write_json(output / 'annotation_resolution_audit.json', summary)
    print(json.dumps({'source_counts': summary['source_counts'], 'uniprot_source_lineage_overlap': summary['uniprot_source_lineage_overlap']}, indent=2))
    db.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--raw-run', type=Path, required=True)
    p.add_argument('--dataset', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    build(a.raw_run, a.dataset, a.output)
