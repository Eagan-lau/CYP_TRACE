"""Conservative BRENDA/SABIO name-dictionary links; adds zero training labels."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import sqlite3
import unicodedata

from build_reaction_core_v2 import json_lines
from build_brenda import VARIANT
from reaction_roles_v1 import RoleRegistry, ROLE_NAMES
from run_raw import digest_file, now, write_json

REACTION_FIELDS = {'reaction', 'natural_substrates_products', 'substrates_products'}
EXPLICIT_INORGANIC_NAMES = {'H2O': 'O', 'H+': '[H+]', 'O2': 'O=O', 'H2O2': 'OO'}


def name_key(name):
    # No case-folding: R/r, D/d and other chemical descriptors can be meaningful.
    return ' '.join(unicodedata.normalize('NFC', name).split())


def brenda_equation(value):
    suffix = re.search(r'\s*\{(r|ir|\?)\}\s*$', value)
    reversibility = suffix.group(1) if suffix else 'not_recorded'
    text = value[:suffix.start()] if suffix else value
    if text.count('=') != 1: return None, reversibility, 'not_one_explicit_equation'
    left, right = text.split('=')
    return ([x.strip() for x in re.split(r'\s+\+\s+', left.strip())],
            [x.strip() for x in re.split(r'\s+\+\s+', right.strip())]), reversibility, None


def resolve_names(names, dictionary, ambiguous):
    structures, unresolved, resolved = [], [], []
    for raw in names:
        key, count = name_key(raw), 1
        coefficient = re.fullmatch(r'([1-9]\d?)\s+(.+)', key)
        if coefficient: count, key = int(coefficient.group(1)), coefficient.group(2)
        if key in dictionary:
            smiles = dictionary[key]['smiles']
            structures.extend([smiles] * count)
            resolved.append({'source_name': raw, 'dictionary_name': key, 'coefficient': count,
                             'smiles': smiles, 'basis': dictionary[key]['basis']})
        else:
            unresolved.append({'source_name': raw, 'dictionary_name': key,
                               'reason': 'ambiguous_source_name' if key in ambiguous else 'name_not_in_closed_dictionary'})
    return {'structures': structures, 'resolved_terms': resolved, 'unresolved_terms': unresolved}


def build(raw_run, core, brenda, output):
    if output.exists(): raise FileExistsError(f'Refusing to replace {output}')
    audit = json.loads((core / 'dataset_audit.json').read_text())
    dbpath = raw_run / 'raw_rebuild.sqlite'
    if audit['status'] != 'PASS' or audit['fresh_raw_sha256'] != digest_file(dbpath):
        raise ValueError('Core/raw provenance did not match')
    output.mkdir(parents=True)
    db = sqlite3.connect(dbpath.resolve().as_uri() + '?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    registry = RoleRegistry(db)
    names = defaultdict(lambda: defaultdict(list))
    for r in db.execute("SELECT * FROM components WHERE structure_status='valid' AND formula_status='match'"):
        if r['name']:
            names[name_key(r['name'])][r['canonical_smiles']].append({'source': r['source'], 'locator': r['locator'],
                'side': r['role'], 'slot': r['slot']})
    dictionary = {name: {'smiles': next(iter(structures)), 'basis': 'unique_exact_case_source_name_structure',
                         'component_locators': next(iter(structures.values()))}
                  for name, structures in names.items() if len(structures) == 1}
    ambiguous = {name: dict(values) for name, values in names.items() if len(values) > 1}
    for role, aliases in ROLE_NAMES.items():
        if role not in registry.records: continue
        for alias in aliases:
            dictionary[name_key(alias)] = {'smiles': registry.records[role][0]['canonical_smiles'],
                'basis': 'predeclared_helper_role:' + role}
    for name, smiles in EXPLICIT_INORGANIC_NAMES.items():
        dictionary[name] = {'smiles': smiles, 'basis': 'explicit_inorganic_formula'}
    write_json(output / 'name_structure_dictionary.json', dictionary)
    write_json(output / 'ambiguous_source_names.json', ambiguous)
    seqs, taxa, family = defaultdict(set), defaultdict(set), set()
    for r in db.execute('SELECT accession,sequence_sha256,taxid,family_evidence FROM protein_assertions'):
        for acc in r['accession'].split('|'):
            if r['sequence_sha256']: seqs[acc].add(r['sequence_sha256'])
            if r['taxid']: taxa[acc].add(r['taxid'])
        if r['family_evidence'] == 'PF00067': family.add(r['sequence_sha256'])
    core_edges = json.loads((core / 'core_edges.json').read_text())
    edge_index = {(r['sequence_sha256'], r['reaction_key']): r for r in core_edges}
    brenda_rows = [json.loads(line) for line in (brenda / 'brenda_assertions.jsonl').read_text(encoding='utf-8').splitlines()]
    source_proteins = {(r['ec'], r['protein_id']): r['raw'] for r in
                       json.loads((brenda / 'source_proteins.json').read_text(encoding='utf-8'))}
    brenda_audit = json.loads((brenda / 'brenda_audit.json').read_text())
    normalized = []
    for r in brenda_rows:
        pids = r['protein_ids']
        accessions = []
        if len(pids) == 1:
            protein = source_proteins.get((r['ec'], pids[0]), {})
            if str(protein.get('source', '')).lower() in {'uniprot', 'swissprot', 'trembl'}:
                accessions = sorted(set(protein.get('accessions', [])))
        terms, reversibility, problem = (None, 'not_applicable', None)
        if r['field'] in REACTION_FIELDS:
            terms, reversibility, problem = brenda_equation(r['raw'].get('value', ''))
        normalized.append({'source': 'BRENDA', 'assertion_id': r['assertion_id'], 'endpoint': r['field'],
            'accessions': accessions, 'source_protein_context_unambiguous': len(pids) == 1 and len(accessions) == 1,
            'taxid': '', 'terms': terms, 'parse_problem': problem, 'reversibility': reversibility,
            'variant_flag': r['variant_text_flag'], 'publication_ids': r['publication_ids'],
            'original_value': r['raw'].get('value', ''),
            'source_locator': {'ec': r['ec'], 'field': r['field'], 'index': r['source_array_index_0based'], 'protein_ids': pids},
            'source_path': brenda_audit['source_path'], 'source_sha256': brenda_audit['source_sha256']})
    sabio_count = 0
    for record in db.execute("SELECT * FROM assertions WHERE source='SABIO_RK' ORDER BY assertion_id"):
        sabio_count += 1
        r, raw = dict(record), json.loads(record['raw_json'])
        accessions = [x for x in r['accession'].split('|') if x]
        terms = (raw.get('Substrates', '').split('|'), raw.get('Products', '').split('|'))
        text = ' '.join(str(raw.get(k, '')) for k in ['Comment', 'EnzymeType', 'Specification'])
        normalized.append({'source': 'SABIO_RK', 'assertion_id': r['assertion_id'], 'endpoint': 'kinetic_reaction_context',
            'accessions': accessions, 'source_protein_context_unambiguous': len(accessions) == 1,
            'taxid': r['taxid'], 'terms': terms, 'parse_problem': None, 'reversibility': 'not_certified',
            'variant_flag': bool(VARIANT.search(text)), 'publication_ids': json.loads(r['publication_ids']),
            'original_value': raw.get('Reaction', ''), 'source_locator': r['locator'],
            'source_path': r['source_path'], 'source_sha256': r['source_sha256']})
    counts, endpoints, problems, outputs, matched_edges = defaultdict(Counter), defaultdict(Counter), defaultdict(Counter), [], defaultdict(set)
    for r in normalized:
        reasons = []
        candidates = sorted({s for a in r['accessions'] for s in seqs[a]})
        linked = candidates[0] if r['source_protein_context_unambiguous'] and len(candidates) == 1 else ''
        if not linked: reasons.append('source_protein_or_sequence_not_unique')
        if linked not in family: reasons.append('exact_sequence_PF00067_confirmation_unavailable')
        if linked and any(len(taxa[a]) != 1 for a in r['accessions']): reasons.append('reference_accession_taxonomy_not_unique')
        if linked and r['taxid'] and any(taxa[a] != {r['taxid']} for a in r['accessions']): reasons.append('source_taxonomy_disagreement')
        if not r['publication_ids']: reasons.append('no_recorded_publication')
        if r['variant_flag']: reasons.append('explicit_variant_text_requires_adjudication')
        left, right, projection = None, None, None
        if r['parse_problem']: reasons.append(r['parse_problem'])
        if r['terms'] is not None:
            left, right = [resolve_names(terms, dictionary, ambiguous) for terms in r['terms']]
            if left['unresolved_terms'] or right['unresolved_terms']: reasons.append('some_named_participants_unresolved')
            else:
                projection = registry.project(left['structures'], right['structures'])
                reasons.extend(projection['problems'])
        else: reasons.append('endpoint_not_an_explicit_substrate_product_equation')
        key = projection['reaction_key'] if projection and not projection['problems'] else None
        matched = edge_index.get((linked, key)) if linked and key else None
        if not matched: reasons.append('no_exact_reference_sequence_main_core_edge_match')
        qualified = not reasons
        counts[r['source']]['source_assertions'] += 1
        endpoints[r['source']][r['endpoint']] += 1
        counts[r['source']]['complete_name_dictionary_resolution'] += left is not None and not (left['unresolved_terms'] or right['unresolved_terms'])
        counts[r['source']]['valid_main_projection'] += bool(key)
        counts[r['source']]['same_order_core_edge_matches_before_context_gates'] += matched is not None
        counts[r['source']]['qualified_named_context_matches_not_new_labels'] += qualified
        if qualified: matched_edges[r['source']].add((linked, key))
        problems[r['source']].update(reasons)
        row = {k: v for k, v in r.items() if k != 'terms'}
        row.update({'linked_reference_sequence_sha256': linked, 'sequence_candidates': candidates,
            'resolved_substrates': left, 'resolved_products': right, 'main_projection': projection,
            'matched_core_edge': {'sequence_sha256': linked, 'reaction_key': key} if matched else None,
            'same_publication_ids_as_core': sorted(set(r['publication_ids']) & set(matched['publication_ids'])) if matched else [],
            'qualified_named_context_match': qualified, 'reasons': sorted(set(reasons)),
            'source_order_is_not_certified_physiological_direction': True,
            'SABIO_list_multiplicity_not_certified': r['source'] == 'SABIO_RK',
            'name_dictionary_match_is_not_external_compound_identity_certification': True,
            'new_model_label_generated': False, 'assayed_construct_verified': False,
            'independent_evidence_certified': False})
        outputs.append(row)
    json_lines(output / 'named_context_resolution.jsonl', outputs)
    checks = {'BRENDA_denominator_preserved': counts['BRENDA']['source_assertions'] == len(brenda_rows),
        'SABIO_denominator_preserved': counts['SABIO_RK']['source_assertions'] == sabio_count,
        'one_output_per_source_assertion': len(outputs) == len({(r['source'], r['assertion_id']) for r in outputs}),
        'no_new_labels_or_independent_claims': all(not r['new_model_label_generated'] and not r['independent_evidence_certified'] for r in outputs),
        'qualified_links_have_existing_core_edge': all(r['matched_core_edge'] for r in outputs if r['qualified_named_context_match']),
        'qualified_BRENDA_links_not_inhibitor_or_kinetic_endpoints': all(r['endpoint'] in REACTION_FIELDS for r in outputs
             if r['source'] == 'BRENDA' and r['qualified_named_context_match']),
        'raw_store_hash_unchanged': audit['fresh_raw_sha256'] == digest_file(dbpath)}
    summary = {'status': 'PASS' if all(checks.values()) else 'FAIL', 'created_utc': now(), 'checks': checks,
        'source_counts': {s: dict(c) for s, c in counts.items()}, 'endpoint_counts': {s: dict(c) for s, c in endpoints.items()},
        'nonexclusive_resolution_reasons': {s: dict(c) for s, c in problems.items()},
        'unique_core_edges_with_qualified_named_context': {s: len(matched_edges[s]) for s in counts},
        'exact_case_dictionary_entries': len(dictionary), 'ambiguous_source_name_entries': len(ambiguous),
        'exact_model_positives_added': 0, 'independent_validation': False,
        'limitation': 'Name matching is an adjudication aid, not certified compound identity, direction, construct or independent evidence.',
        'BRENDA_reversibility_source': 'https://brenda-enzymes.org/schemas/docs/2.0.0/brenda.schema.html',
        'license_sensitive_records_remain_private_workspace_outputs': True}
    write_json(output / 'named_support_audit.json', summary)
    inputs = [dbpath, core / 'core_edges.json', core / 'dataset_audit.json', brenda / 'brenda_assertions.jsonl',
        brenda / 'source_proteins.json', brenda / 'brenda_audit.json', Path(__file__),
        Path(__file__).with_name('reaction_roles_v1.py'), Path(__file__).with_name('build_brenda.py')]
    write_json(output / 'input_manifest.json', {str(p.resolve()): digest_file(p) for p in inputs})
    write_json(output / 'output_checksums.json', {p.name: digest_file(p) for p in sorted(output.iterdir()) if p.is_file()})
    db.close()
    print(json.dumps(summary, indent=2))
    if not all(checks.values()): raise RuntimeError('Named support audit failed')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--raw-run', type=Path, required=True)
    p.add_argument('--core', type=Path, required=True)
    p.add_argument('--brenda', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    build(a.raw_run, a.core, a.brenda, a.output)
