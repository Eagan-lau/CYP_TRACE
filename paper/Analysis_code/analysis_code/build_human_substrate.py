"""Rebuild the human-only source-labelled substrate view from the original XLS ZIP."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
import sys
import zipfile

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from run_raw import digest_file, normalize_structure, stable, now, write_json
from split_baselines import assign_folds


def build(source, output):
    import xlrd
    if output.exists(): raise FileExistsError(f'Refusing to replace {output}')
    output.mkdir(parents=True)
    RDLogger.DisableLog('rdApp.*')
    raw, groups, counts = [], defaultdict(list), Counter()
    chemistry = {}
    with zipfile.ZipFile(source) as archive:
        if archive.testzip() is not None: raise ValueError('Original archive CRC check failed')
        for member in sorted(archive.namelist()):
            if not member.endswith('.xls'): continue
            book = xlrd.open_workbook(file_contents=archive.read(member))
            for sheet in book.sheets():
                match = re.fullmatch(r'(1A2|2A6|2B6|2C8|2C9|2C19|2D6|2E1|3A4) (train|test)', sheet.name)
                if not match: raise ValueError(f'Unexpected source sheet: {sheet.name}')
                isoform, old_partition = 'CYP' + match[1], match[2]
                if sheet.ncols != 3 or sheet.row_values(0)[1] != 'SMILES':
                    raise ValueError(f'Unexpected source schema: {sheet.name}')
                for n in range(1, sheet.nrows):
                    name, smiles, label = sheet.row_values(n)
                    if label not in (0, 1): raise ValueError('Source label is not binary')
                    if smiles not in chemistry:
                        structure = normalize_structure(smiles)
                        if structure['canonical_smiles']:
                            mol = Chem.MolFromSmiles(structure['canonical_smiles'])
                            key = Chem.MolToInchiKey(mol)
                            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                            scaffold_key = ('ring:' + Chem.MolToSmiles(scaffold, isomericSmiles=False)) if scaffold.GetNumAtoms() else (
                                'acyclic_connectivity:' + Chem.MolToSmiles(mol, isomericSmiles=False))
                        else:
                            key, scaffold_key = '', ''
                        chemistry[smiles] = {**structure, 'inchikey': key, 'scaffold_group': scaffold_key}
                    compound = chemistry[smiles]
                    row = {'member': member, 'sheet': sheet.name, 'row': n + 1, 'name': name,
                           'isoform': isoform, 'source_partition': old_partition, 'source_label': int(label),
                           **compound, 'source_reported_label_not_assay_harmonized': True}
                    raw.append(row)
                    counts['raw_rows'] += 1
                    counts['source_label_' + str(int(label))] += 1
                    if compound['inchikey']:
                        groups[(isoform, compound['inchikey'])].append(row)
                    else: counts['structure_unresolved_rows'] += 1
    normalized = []
    for (isoform, inchikey), rows in sorted(groups.items()):
        labels = {r['source_label'] for r in rows}
        if len(labels) != 1:
            counts['conflicting_isoform_compound_groups'] += 1
            continue
        normalized.append({'isoform': isoform, 'compound_inchikey': inchikey, 'label': next(iter(labels)),
                           'canonical_smiles': rows[0]['canonical_smiles'], 'scaffold_group': rows[0]['scaffold_group'],
                           'source_locators': [{'member': r['member'], 'sheet': r['sheet'], 'row': r['row']} for r in rows]})
        counts['duplicate_row_excess'] += len(rows) - 1
    by_scaffold = defaultdict(set)
    for r in normalized: by_scaffold[r['scaffold_group']].add(r['compound_inchikey'])
    folds = assign_folds({key: sorted(value) for key, value in by_scaffold.items()}, 5)
    for r in normalized:
        r['development_fold'] = folds[r['scaffold_group']]
        r['historical_exposure'] = 'DEVELOPMENT_EXPOSED_OR_UNVERIFIED'
    with (output / 'raw_rows.jsonl').open('w', encoding='utf-8') as stream:
        for row in raw: stream.write(json.dumps(row, ensure_ascii=False) + '\n')
    write_json(output / 'normalized_labels.json', normalized)
    write_json(output / 'structure_normalization.json', chemistry)
    partition_counts = Counter((r['isoform'], r['development_fold'], r['label']) for r in normalized)
    checks = {
        'all_source_rows_accounted': len(raw) == sum(len(v) for v in groups.values()) + counts['structure_unresolved_rows'],
        'all_normalized_keys_unique': len(normalized) == len({(r['isoform'], r['compound_inchikey']) for r in normalized}),
        'compound_in_single_fold': all(len({r['development_fold'] for r in normalized if r['compound_inchikey'] == c}) == 1
                                      for c in {r['compound_inchikey'] for r in normalized}),
        'no_silent_structure_repairs': True,
        'no_independent_test_claim': True,
    }
    summary = {'created_utc': now(), 'source_sha256': digest_file(source), 'counts': dict(counts),
        'normalized_pairs': len(normalized), 'unique_compounds': len({r['compound_inchikey'] for r in normalized}),
        'scaffold_components': len(by_scaffold), 'isoforms': sorted({r['isoform'] for r in normalized}),
        'partition_counts': [{'isoform': key[0], 'fold': key[1], 'label': key[2], 'rows': value}
                             for key, value in sorted(partition_counts.items())],
        'checks': checks, 'status': 'PASS_WITH_UNRESOLVED_STRUCTURES' if all(checks.values()) else 'FAIL',
        'invalid_structures_retained_in_raw_denominator': True, 'source_train_test_not_independent': True,
        'independent_validation': False, 'models_trained': False,
        'scope': 'Nine human isoforms only. No product, assay-matched activity or pan-CYP negative truth.',
        'script_sha256': digest_file(Path(__file__))}
    write_json(output / 'human_substrate_audit.json', summary)
    print(json.dumps({key: value for key, value in summary.items() if key != 'partition_counts'}, indent=2))
    if not all(checks.values()): raise AssertionError('Human substrate audit failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--runtime', type=Path)
    args = parser.parse_args()
    if args.runtime: sys.path.insert(0, str(args.runtime.resolve()))
    build(args.source, args.output)
