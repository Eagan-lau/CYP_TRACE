"""Extract typed CYP-related source assertions from the original BRENDA JSON.

The EC block is a context boundary, never a protein--reaction Cartesian product.
No molecular structure or exact assayed sequence is invented from a name.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import sqlite3
import tarfile

from run_raw import digest_file, stable, now, write_json

FIELDS = ('reaction', 'natural_substrates_products', 'substrates_products', 'km_value',
          'turnover_number', 'kcat_km', 'inhibitor', 'ki_value', 'ic50_value')
CYP_TEXT = re.compile(r'\b(?:P-?450|CYP\d+[A-Z0-9]*)\b', re.I)
VARIANT = re.compile(r'\b(?:mutant|mutagenesis|engineered|chimera|chimeric|variant)\b|\b[A-Z]\d{1,4}[A-Z]\b', re.I)


def accession_links(protein, known):
    if str(protein.get('source', '')).lower() not in {'uniprot', 'swissprot', 'trembl'}:
        return []
    return sorted(set(protein.get('accessions', [])) & known)


def typed_row(ec, field, number, row, proteins, references, known):
    pids = [str(p) for p in row.get('proteins', [])]
    links = {pid: accession_links(proteins.get(pid, {}), known) for pid in pids}
    all_accessions = sorted({a for values in links.values() for a in values})
    reference_ids = [str(r) for r in row.get('references', [])]
    pmids = sorted({'PMID:' + str(references[r]['pmid']) for r in reference_ids
                    if r in references and references[r].get('pmid')})
    missing_pids = [p for p in pids if p not in proteins]
    missing_refs = [r for r in reference_ids if r not in references]
    if not pids: scope = 'EC_CONTEXT_NO_PROTEIN_LINK'
    elif not all_accessions: scope = 'SOURCE_PROTEIN_CONTEXT_NO_LINK_TO_CURRENT_CYP_SEQUENCE_STORE'
    elif len(pids) == 1 and len(all_accessions) == 1: scope = 'SINGLE_SOURCE_PROTEIN_ACCESSION_LINK_NOT_CONSTRUCT_VERIFIED'
    else: scope = 'MULTIPLE_SOURCE_PROTEIN_LINKS_NOT_EXPANDED'
    return {'assertion_id': stable(f'BRENDA:2026.1:{ec}:{field}:{number}'), 'ec': ec,
            'field': field, 'source_array_index_0based': number, 'protein_ids': pids,
            'accession_links_by_protein_id': links, 'linked_cyp_accessions': all_accessions,
            'publication_ids': pmids, 'reference_ids': reference_ids,
            'missing_protein_ids': missing_pids, 'missing_reference_ids': missing_refs,
            'scope': scope, 'variant_text_flag': bool(VARIANT.search(row.get('comment', ''))),
            'variant_flag_absence_does_not_certify_wild_type': True,
            'chemical_structure_resolved': False, 'eligible_exact_reaction_positive': False,
            'historical_exposure': 'DEVELOPMENT_EXPOSED_OR_UNVERIFIED', 'raw': row}


def build(source, raw_run, output):
    if output.exists(): raise FileExistsError(f'Refusing to overwrite {output}')
    output.mkdir(parents=True)
    dbpath = raw_run / 'raw_rebuild.sqlite'
    db = sqlite3.connect(dbpath.resolve().as_uri() + '?mode=ro', uri=True)
    known = {acc for (value,) in db.execute('SELECT DISTINCT accession FROM protein_assertions')
             for acc in value.split('|') if acc}
    db.close()
    print('Reading original BRENDA JSON (not a historical extract)', flush=True)
    with tarfile.open(source, 'r:gz') as archive:
        members = [m for m in archive.getmembers() if m.isfile() and m.name.endswith('.json')]
        if len(members) != 1: raise ValueError('Expected exactly one JSON archive member')
        member = members[0]
        with archive.extractfile(member) as stream: payload = json.load(stream)
    if payload.get('release') != '2026.1' or not isinstance(payload.get('data'), dict):
        raise ValueError('Unexpected raw BRENDA release/schema')
    counts, per_field, scopes = Counter(), Counter(), Counter()
    contexts, proteins_out, references_out = [], [], []
    with (output / 'brenda_assertions.jsonl').open('w', encoding='utf-8') as stream:
        for ec, entry in payload['data'].items():
            counts['all_source_ec_blocks_scanned'] += 1
            proteins = entry.get('protein', {})
            references = entry.get('reference', {})
            if not isinstance(proteins, dict) or not isinstance(references, dict):
                raise ValueError(f'Unexpected protein/reference schema at {ec}')
            direct = {p: accession_links(record, known) for p, record in proteins.items()}
            direct = {p: a for p, a in direct.items() if a}
            names = ' '.join([str(entry.get('recommended_name', '')), str(entry.get('systematic_name', ''))]
                             + [str(x.get('value', '')) for x in entry.get('synonyms', [])])
            name_context = bool(CYP_TEXT.search(names))
            if not direct and not name_context: continue
            counts['selected_context_ec_blocks'] += 1
            counts['ec_blocks_with_current_CYP_accession_link'] += bool(direct)
            contexts.append({'ec': ec, 'name_context_match': name_context,
                             'protein_accession_links': direct,
                             'recommended_name': entry.get('recommended_name', ''),
                             'context_selection_not_CYP_truth_for_every_record': True})
            for pid, record in proteins.items():
                proteins_out.append({'ec': ec, 'protein_id': pid, 'raw': record,
                                     'current_cyp_accessions': direct.get(pid, [])})
            for refid, record in references.items():
                references_out.append({'ec': ec, 'reference_id': refid, 'raw': record})
            for field in FIELDS:
                rows = entry.get(field, [])
                if not isinstance(rows, list): raise ValueError(f'Unexpected {ec}/{field} schema')
                for number, row in enumerate(rows):
                    out = typed_row(ec, field, number, row, proteins, references, known)
                    stream.write(json.dumps(out, ensure_ascii=False, allow_nan=False) + '\n')
                    counts['context_assertions_retained'] += 1
                    counts['assertions_explicitly_linking_current_CYP_accession'] += bool(out['linked_cyp_accessions'])
                    counts['assertions_with_PMID'] += bool(out['publication_ids'])
                    counts['variant_text_flagged_assertions'] += out['variant_text_flag']
                    counts['assertions_with_broken_protein_reference'] += bool(out['missing_protein_ids'])
                    counts['assertions_with_broken_literature_reference'] += bool(out['missing_reference_ids'])
                    per_field[field] += 1
                    scopes[out['scope']] += 1
    write_json(output / 'ec_contexts.json', contexts)
    write_json(output / 'source_proteins.json', proteins_out)
    write_json(output / 'source_references.json', references_out)
    summary = {'created_utc': now(), 'counts': dict(counts), 'per_field': dict(per_field), 'scopes': dict(scopes),
        'status': 'TYPED_SOURCE_EXTRACTION_COMPLETE_CHEMICAL_ADJUDICATION_PENDING',
        'source_path': str(source.resolve()), 'source_sha256': digest_file(source),
        'member': member.name, 'release': payload['release'], 'format_version': payload.get('version'),
        'fresh_raw_store_sha256': digest_file(dbpath), 'script_sha256': digest_file(Path(__file__)),
        'exact_model_positives_added': 0, 'source_assertions_not_independent_experiments': True,
        'not_a_complete_CYP_census': True, 'independent_validation': False, 'paper_complete': False,
        'limitations': ['Context includes mixed enzyme classes and organism-level protein records',
            'Known accession linkage is not exact assayed construct verification',
            'Text variant flags require adjudication; absence does not prove wild type',
            'Named substrate/product terms are not silently converted to precise chemical structures',
            'Inhibitor/Ki/IC50 records remain separate from substrate positives',
            'License-sensitive source records are private workspace outputs, not a public release']}
    write_json(output / 'brenda_audit.json', summary)
    print(json.dumps(summary, ensure_ascii=True, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--raw-run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    build(a.source, a.raw_run, a.output)
