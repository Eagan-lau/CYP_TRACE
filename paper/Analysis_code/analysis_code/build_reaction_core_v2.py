"""P02: source-linked main transformations from the fresh raw store, not old labels.

This is a development cohort, not a construct-verified or independent truth set.
Original full-component reactions are preserved alongside each MAIN projection.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import sqlite3

from rdkit import rdBase
from audit_annotations import exact_reaction, rhea_links
from reaction_roles_v1 import RoleRegistry, MACRO_ROLES
from run_raw import digest_file, now, write_json

SOURCES = ('P450Rdb', 'SwissProt', 'UniProt_API')
VARIANT = re.compile(r'\b(mutant|mutated|mutation|variant|engineered|chimeric|chimera)\b', re.I)


def experimental_publications(source, raw):
    """Only the reaction evidence; a physiological-direction citation is insufficient."""
    if source == 'UniProt_API':
        return sorted({'PMID:' + e['id'] for e in raw.get('reaction', {}).get('evidences', [])
            if e.get('evidenceCode') == 'ECO:0000269' and e.get('source') == 'PubMed'
            and re.fullmatch(r'\d+', e.get('id', ''))})
    text = raw.get('comment', '').split('PhysiologicalDirection=', 1)[0]
    text = re.sub(r'\s+', '', text)
    return sorted({'PMID:' + p for p in re.findall(r'ECO:0000269\|PubMed:(\d+)', text)})


def macro_ids(source, raw):
    if source == 'UniProt_API':
        refs = [x.get('id', '') for x in raw.get('reaction', {}).get('reactionCrossReferences', [])]
        return sorted({x.split(':')[1] for x in refs if re.fullmatch(r'RHEA-COMP:\d+', x)})
    return sorted(set(re.findall(r'RHEA-COMP:(\d+)', re.sub(r'\s+', '', raw.get('comment', '')))))


def declared_directions(source, raw):
    if source == 'UniProt_API':
        return {x.get('reactionCrossReference', {}).get('id', '').removeprefix('RHEA:'):
                x.get('directionType', '') for x in raw.get('physiologicalReactions', [])}
    text = re.sub(r'\s+', '', raw.get('comment', ''))
    return {rid: direction for direction, rid in re.findall(
        r'PhysiologicalDirection=([^;]+);Xref=Rhea:RHEA:(\d+)', text)}


def json_lines(path, rows):
    with path.open('w', encoding='utf-8') as out:
        for row in rows: out.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + '\n')


def build(raw_run, initial, output):
    if output.exists(): raise FileExistsError(f'Refusing to overwrite {output}')
    root = Path(__file__).resolve().parent
    dbpath = raw_run / 'raw_rebuild.sqlite'
    initial_audit = json.loads((initial / 'dataset_audit.json').read_text())
    raw_hash = digest_file(dbpath)
    if initial_audit['new_raw_build_sha256'] != raw_hash: raise ValueError('Raw input hash changed')
    if initial_audit['status'] != 'PASS': raise ValueError('Initial membership did not pass')
    output.mkdir(parents=True)
    db = sqlite3.connect(dbpath.resolve().as_uri() + '?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    registry = RoleRegistry(db)
    write_json(output / 'participant_role_manifest.json', registry.manifest())
    accession_seqs, accession_taxa, family = defaultdict(set), defaultdict(set), set()
    for r in db.execute('SELECT accession,sequence_sha256,taxid,family_evidence FROM protein_assertions'):
        if r['family_evidence'] == 'PF00067': family.add(r['sequence_sha256'])
        for a in r['accession'].split('|'):
            if a and r['sequence_sha256']: accession_seqs[a].add(r['sequence_sha256'])
            if a and r['taxid']: accession_taxa[a].add(r['taxid'])
    proteins = {r['sequence_sha256']: dict(r) for r in db.execute('SELECT * FROM proteins')}
    initial_membership = {r['assertion_id']: r for r in
        map(json.loads, (initial / 'membership.jsonl').read_text().splitlines())}
    # Keep row-level source faults even when the source's formula field agrees.
    conflict_keys = {(r['name'], r['canonical_smiles']) for r in registry.conflicts}
    bad_components = [dict(r) for r in db.execute('SELECT * FROM components')
                      if (r['name'], r['canonical_smiles']) in conflict_keys]
    bad_locators = {(r['source'], r['locator']) for r in bad_components}
    json_lines(output / 'quarantined_component_name_conflicts.jsonl', bad_components)
    unavailable_names = {name for role, names in __import__('reaction_roles_v1').ROLE_NAMES.items()
                         if role in registry.unavailable for name in names}
    unavailable_locators = {(r['source'], r['locator']) for r in db.execute('SELECT source,locator,name FROM components')
                            if r['name'] in unavailable_names}
    rhea = {r['rhea_id']: dict(r) for r in db.execute('SELECT * FROM rhea_structures')}
    masters = {r['master_id'] for r in rhea.values()}
    full = {r['reaction_key']: {'substrates': json.loads(r['substrates']),
                              'products': json.loads(r['products'])} for r in db.execute('SELECT * FROM reactions')}
    rhea_cache, projections = {}, {}

    def project(key):
        if key not in projections:
            projections[key] = registry.project(full[key]['substrates'], full[key]['products'])
        return projections[key]

    def parse_rhea(rid):
        if rid not in rhea_cache:
            parsed = exact_reaction(rhea[rid]['reaction_smiles'])
            if parsed['reaction_key']:
                full[parsed['reaction_key']] = {'substrates': parsed['sides']['substrate'],
                                                'products': parsed['sides']['product']}
            rhea_cache[rid] = dict(parsed, master_id=rhea[rid]['master_id'], direction=rhea[rid]['direction'])
        return rhea_cache[rid]

    membership, admitted, counts, exclusions = [], defaultdict(list), defaultdict(Counter), defaultdict(Counter)
    for record in db.execute('SELECT * FROM assertions ORDER BY source,assertion_id'):
        source = record['source']
        if source not in SOURCES: continue
        r, raw = dict(record), json.loads(record['raw_json'])
        seq, acc, aid = r['sequence_sha256'], r['accession'], r['assertion_id']
        reasons, full_keys, main_keys = [], [], []
        primary, physiological, directed, macros = [], [], [], []
        pubs = json.loads(r['publication_ids'])
        if not proteins.get(seq, {}).get('model_sequence_valid'): reasons.append('no_valid_explicit_sequence')
        if seq not in family: reasons.append('exact_sequence_PF00067_confirmation_unavailable')
        if not acc or '|' in acc or accession_seqs[acc] != {seq}:
            reasons.append('accession_sequence_ambiguous_or_unavailable')
        if not r['taxid'] or accession_taxa[acc] != {r['taxid']}:
            reasons.append('accession_taxonomy_ambiguous_or_unavailable')
        if source == 'P450Rdb':
            reasons.extend(initial_membership[aid]['reasons'])
            if (source, r['locator']) in bad_locators: reasons.append('helper_name_structure_conflict')
            if (source, r['locator']) in unavailable_locators: reasons.append('unavailable_helper_identity')
            if r['reaction_key'] in full: full_keys.append(r['reaction_key'])
        else:
            pubs = experimental_publications(source, raw)
            if not pubs: reasons.append('reaction_level_experimental_PubMed_evidence_unavailable')
            text = raw.get('reaction', {}).get('name', '') if source == 'UniProt_API' else raw.get('comment', '')
            if VARIANT.search(text): reasons.append('explicit_variant_or_engineered_assertion_requires_adjudication')
            primary, physiological = rhea_links(source, raw)
            declared = sorted(set(primary + physiological))
            directed = [rid for rid in declared if rid in rhea]
            if not directed: reasons.append('explicit_directional_Rhea_structure_unavailable')
            if any(rid not in rhea and rid not in masters for rid in declared):
                reasons.append('declared_Rhea_identifier_unresolved_in_frozen_release')
            primary_masters = {rhea[rid]['master_id'] if rid in rhea else rid for rid in primary}
            if physiological and (not primary_masters or any(
                    rhea[rid]['master_id'] not in primary_masters for rid in directed if rid in physiological)):
                reasons.append('physiological_direction_primary_reaction_mismatch')
            macros = macro_ids(source, raw)
            if any(m not in MACRO_ROLES for m in macros): reasons.append('unadjudicated_macromolecule_participant')
            direction_text = declared_directions(source, raw)
            for rid in directed:
                expected = {'left-to-right': 'LR', 'right-to-left': 'RL'}.get(direction_text.get(rid, ''))
                if expected and expected != rhea[rid]['direction']: reasons.append('direction_text_identifier_conflict')
                parsed = parse_rhea(rid)
                if not parsed['reaction_key']: reasons.append('Rhea_chemistry_' + parsed['status'])
                else: full_keys.append(parsed['reaction_key'])
        if not full_keys: reasons.append('no_parseable_complete_reaction')
        for key in sorted(set(full_keys)):
            projection = project(key)
            reasons.extend(projection['problems'])
            removed_roles = {p['role'] for p in projection['removed_participants']}
            if any(m in MACRO_ROLES and MACRO_ROLES[m] not in removed_roles for m in macros):
                reasons.append('declared_macro_reactive_part_not_resolved_as_auxiliary')
            main_keys.append(projection['reaction_key'])
        main_keys = sorted(set(main_keys))
        if len(main_keys) > 1: reasons.append('multiple_main_directions_or_transformations_in_one_assertion')
        reasons = sorted(set(reasons))
        row = {'assertion_id': aid, 'source': source, 'source_lineage': 'UniProt' if source != 'P450Rdb' else source,
            'eligible': not reasons, 'reasons': reasons, 'sequence_sha256': seq, 'accession': acc, 'taxid': r['taxid'],
            'full_reaction_keys': sorted(set(full_keys)), 'main_reaction_keys': main_keys,
            'declared_primary_rhea_ids': primary, 'declared_physiological_rhea_ids': physiological,
            'directional_rhea_ids': directed, 'macro_component_ids': macros,
            'admissible_publication_ids': pubs, 'all_source_publication_ids': json.loads(r['publication_ids']),
            'source_locator': r['locator'], 'source_path': r['source_path'], 'source_sha256': r['source_sha256'],
            'exact_assayed_construct_verified': False, 'CYP_domain_activity_attribution_verified': False,
            'historical_exposure': 'DEVELOPMENT_EXPOSED_OR_UNVERIFIED'}
        membership.append(row)
        counts[source]['source_assertions'] += 1
        counts[source]['admitted_assertions'] += int(row['eligible'])
        counts[source]['parseable_complete_chemistry_assertions'] += bool(full_keys)
        counts[source]['candidate_main_transformation_assertions'] += bool(main_keys)
        exclusions[source].update(reasons)
        if row['eligible']: admitted[(seq, main_keys[0])].append(row)
    edges = []
    for (seq, key), rows in sorted(admitted.items()):
        edges.append({'sequence_sha256': seq, 'reaction_key': key,
            'publication_ids': sorted({p for r in rows for p in r['admissible_publication_ids']}),
            'assertion_ids': sorted(r['assertion_id'] for r in rows),
            'accessions': sorted({r['accession'] for r in rows}), 'taxids': sorted({r['taxid'] for r in rows}),
            'sources': sorted({r['source'] for r in rows}), 'source_lineages': sorted({r['source_lineage'] for r in rows}),
            'full_reaction_keys': sorted({k for r in rows for k in r['full_reaction_keys']}),
            'assayed_construct_verified': False, 'CYP_domain_activity_attribution_verified': False,
            'source_scope': 'multi_source_reference_sequence_development_core',
            'historical_exposure': 'DEVELOPMENT_EXPOSED_OR_UNVERIFIED'})
    sequences, reactions = {e['sequence_sha256'] for e in edges}, {e['reaction_key'] for e in edges}
    chemistry = {p['reaction_key']: {k: p[k] for k in ['substrates', 'products', 'single_pair']}
                 for p in projections.values() if p['reaction_key'] in reactions}
    json_lines(output / 'membership.jsonl', membership)
    write_json(output / 'core_edges.json', edges)
    write_json(output / 'core_reactions.json', chemistry)
    write_json(output / 'full_reaction_projections.json', {k: dict(p, full_reaction=full[k]) for k, p in projections.items()})
    write_json(output / 'declared_rhea_resolution.json', rhea_cache)
    with (output / 'core_sequences.fasta').open('w', encoding='ascii') as out:
        for seq in sorted(sequences): out.write('>' + seq + '\n' + proteins[seq]['sequence'] + '\n')
    edge_sets = {source: {(e['sequence_sha256'], e['reaction_key']) for e in edges if source in e['sources']} for source in SOURCES}
    pairwise = []
    for i, a in enumerate(SOURCES):
        for b in SOURCES[i+1:]:
            pairwise.append({'source_a': a, 'source_b': b, 'shared_unique_edges': len(edge_sets[a] & edge_sets[b]),
                'source_a_edge_denominator': len(edge_sets[a]), 'source_b_edge_denominator': len(edge_sets[b]),
                'independent_evidence_certified': False, 'same_source_lineage': a != 'P450Rdb' and b != 'P450Rdb'})
    conservation = all(Counter(full[k][side + 's']) == Counter(p[side + 's']) + Counter(
        x['smiles'] for x in p['removed_participants'] if x['side'] == side)
        for k, p in projections.items() for side in ['substrate', 'product'])
    checks = {'raw_store_hash_unchanged': digest_file(dbpath) == raw_hash,
        'membership_denominators_preserved': all(counts[s]['source_assertions'] == db.execute(
            'SELECT COUNT(*) FROM assertions WHERE source=?', (s,)).fetchone()[0] for s in SOURCES),
        'source_assertion_admission_reconciles': sum(r['eligible'] for r in membership) == sum(len(e['assertion_ids']) for e in edges),
        'unique_sequence_main_edges': len(edges) == len(admitted),
        'all_main_chemistry_and_sequences_present': all(e['reaction_key'] in chemistry and e['sequence_sha256'] in proteins for e in edges),
        'all_admitted_assertions_have_one_main_label': all(len(r['main_reaction_keys']) == 1 for r in membership if r['eligible']),
        'all_edges_have_admissible_publications': all(e['publication_ids'] for e in edges),
        'all_sequences_exact_family_confirmed': sequences <= family,
        'full_component_multiset_conserved': conservation,
        'conflicting_helper_source_rows_not_admitted': all(not r['eligible'] for r in membership if (r['source'], r['source_locator']) in bad_locators),
        'construct_or_independence_not_claimed': all(not e['assayed_construct_verified'] for e in edges)}
    audit = {'status': 'PASS' if all(checks.values()) else 'FAIL', 'created_utc': now(), 'checks': checks,
        'source_counts': {k: dict(v) for k, v in counts.items()},
        'nonexclusive_exclusion_counts': {k: dict(v) for k, v in exclusions.items()},
        'labelled_sequences': len(sequences), 'unique_edges': len(edges), 'directed_main_transformations': len(chemistry),
        'single_pair_transformations': sum(c['single_pair'] for c in chemistry.values()),
        'multicomponent_transformations': sum(not c['single_pair'] for c in chemistry.values()),
        'unique_edges_by_export': {s: len(x) for s, x in edge_sets.items()},
        'cross_source_edge_overlap': pairwise, 'source_lineage_counts_are_not_independent_evidence_counts': True,
        'quarantined_helper_name_structure_components': len(bad_components),
        'quarantined_helper_name_structure_parent_assertions': len(bad_locators),
        'reaction_scope': 'source-linked main transformation, not exact assayed construct certified',
        'independent_test_records_certified': 0, 'model_performance_evaluated': False,
        'paper_complete': False, 'rdkit_version': rdBase.rdkitVersion, 'fresh_raw_sha256': raw_hash,
        'limitations': ['Historical data remain development-exposed',
            'UniProt exports share a lineage; distinct databases are not independent experiments',
            'Exact reference sequence is not necessarily the experimental construct',
            'CYP domain activity in multifunctional proteins is not certified',
            'Main compound protonation and tautomer differences are deliberately not merged',
            'Incomplete source stoichiometry and unidentified auxiliary participants may remain',
            'This changed label estimand cannot be compared directly with full-component MRR']}
    write_json(output / 'dataset_audit.json', audit)
    inputs = [dbpath, initial / 'membership.jsonl', initial / 'dataset_audit.json',
        root / 'REACTION_RULES_V1.md', root / 'reaction_roles_v1.py', root / 'audit_annotations.py',
        root / 'run_raw.py', Path(__file__)]
    write_json(output / 'input_manifest.json', {str(p.resolve()): digest_file(p) for p in inputs})
    write_json(output / 'output_checksums.json', {p.name: digest_file(p) for p in sorted(output.iterdir()) if p.is_file()})
    db.close()
    print(json.dumps(audit, indent=2))
    if not all(checks.values()): raise RuntimeError('Core v2 audit failed')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--raw-run', type=Path, required=True)
    p.add_argument('--initial-dataset', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    build(a.raw_run, a.initial_dataset, a.output)
