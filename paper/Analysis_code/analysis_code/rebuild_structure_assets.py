"""Re-extract physical heme neighbourhoods from raw CIFs, without legacy features."""
from __future__ import annotations
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree
from run_raw import digest_file, stable, now, write_json


def build(raw_folder, output):
    if output.exists(): raise FileExistsError(f'Refusing to overwrite {output}')
    output.mkdir(parents=True)
    parser = MMCIFParser(QUIET=True)
    counts, errors, inputs = Counter(), [], []
    files = sorted(raw_folder.glob('*.cif.gz'))
    if not files: raise ValueError('No original mmCIF files')
    with (output / 'heme_contact_residues.jsonl').open('w') as contacts, (output / 'chain_assets.jsonl').open('w') as chains:
        for number, path in enumerate(files, 1):
            sha = digest_file(path)
            inputs.append({'file': path.name, 'sha256': sha, 'bytes': path.stat().st_size})
            counts['files_attempted'] += 1
            try:
                with gzip.open(path, 'rt') as handle: structure = parser.get_structure(path.name, handle)
                model = next(structure.get_models())
                hemes = []
                for chain in model:
                    for residue in chain:
                        if residue.resname not in {'HEM', 'HEC', 'HEA'}: continue
                        atoms = [a for a in residue.get_atoms() if a.element not in {'H', 'D'}]
                        irons = [a for a in atoms if a.element.upper() == 'FE']
                        if not irons: continue
                        hemes.append((chain.id, residue.id, residue.resname, atoms, irons[0]))
                if not hemes:
                    counts['files_without_recognized_iron_heme'] += 1
                    continue
                counts['files_with_recognized_iron_heme'] += 1
                counts['physical_heme_instances'] += len(hemes)
                for chain in model:
                    residues = [r for r in chain if is_aa(r, standard=False)]
                    if not residues: continue
                    sequence = ''.join(seq1(r.resname, custom_map={'MSE': 'M'}, undef_code='X') for r in residues)
                    chain_record = {'file': path.name, 'file_sha256': sha, 'model': model.id, 'chain': chain.id,
                        'sequence': sequence, 'sequence_sha256': stable(sequence), 'resolved_residues': len(residues),
                        'sequence_scope': 'resolved_coordinate_residues_not_full_construct',
                        'catalysis_label_generated': False}
                    chains.write(json.dumps(chain_record) + '\n')
                    counts['resolved_protein_chains'] += 1
                    for heme_chain, heme_id, heme_name, atoms, iron in hemes:
                        tree = cKDTree(np.asarray([a.coord for a in atoms], dtype=float))
                        for position, residue in enumerate(residues):
                            heavy = [a for a in residue.get_atoms() if a.element not in {'H', 'D'}]
                            if not heavy: continue
                            coords = np.asarray([a.coord for a in heavy], dtype=float)
                            distance = float(np.min(tree.query(coords)[0]))
                            if distance > 5.0: continue
                            row = {'file': path.name, 'file_sha256': sha, 'model': model.id, 'protein_chain': chain.id,
                                'coordinate_sequence_position_0based': position, 'auth_residue_number': residue.id[1],
                                'insertion_code': residue.id[2].strip(), 'residue_name': residue.resname,
                                'aa': seq1(residue.resname, custom_map={'MSE': 'M'}, undef_code='X'),
                                'heme_chain': heme_chain, 'heme_residue_id': list(heme_id), 'heme_component': heme_name,
                                'minimum_heme_heavy_atom_distance_A': distance,
                                'minimum_iron_distance_A': float(np.min(np.linalg.norm(coords - iron.coord, axis=1))),
                                'cross_chain_contact': chain.id != heme_chain, 'label': 'GEOMETRY_ONLY_NOT_CATALYTIC_TRUTH'}
                            contacts.write(json.dumps(row) + '\n')
                            counts['heme_contact_residue_rows'] += 1
                            counts['cross_chain_contact_rows'] += int(chain.id != heme_chain)
            except Exception as error:
                errors.append({'file': path.name, 'error': str(error), 'type': type(error).__name__})
            if number % 100 == 0: print(f'{number}/{len(files)} raw structures processed', flush=True)
    write_json(output / 'source_manifest.json', inputs)
    write_json(output / 'structure_audit.json', {'created_utc': now(), 'counts': dict(counts), 'errors': errors,
        'status': 'RAW_GEOMETRY_COMPLETE' if not errors else 'RAW_GEOMETRY_PARTIAL_WITH_ERRORS',
        'definition': 'First deposited model; Bio.PDB selected altlocs; protein heavy atoms within 5 A of HEM/HEC/HEA with iron; cross-chain contacts retained and flagged',
        'not_implemented': ['other heme chemical components', 'query-specific substrate contacts', 'channel dynamics', 'canonical-query transfer', 'a trained structural expert'],
        'independent_validation': False, 'legacy_features_read': False,
        'script_sha256': digest_file(Path(__file__)), 'paper_complete': False})
    print(json.dumps(dict(counts)), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--raw-folder', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    build(a.raw_folder, a.output)
