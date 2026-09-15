"""Raw polymer-aware positional reference selection without functional scores."""
from collections import defaultdict, Counter
import gzip
import json
from pathlib import Path
import re
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from run_raw import digest_file, write_json, stable, now


def run(root):
    out = root / 'ordered_anchor_01'
    if out.exists(): raise FileExistsError(out)
    read = lambda n: json.loads((root / n).read_text())
    if read('structure_association_validation_02/INDEPENDENT_QC.json')['status'] != 'PASS': raise ValueError('Corrected identity gate')
    selected = read('structure_mapping_audit_02/primary_representatives.json')
    bychain = defaultdict(list)
    for q, r in selected.items(): bychain[(r['file'], r['chain'])].append((q, r))
    manifest = {r['file']: r['sha256'] for r in read('structure_assets_01/source_manifest.json')}
    parser = MMCIFParser(QUIET=True); census = []; contacts = []; eligible = []
    for (file, chainid), members in sorted(bychain.items()):
        path = root.parent / 'data/raw/rcsb_cyp_ligand_templates_v1/mmcif' / file
        if digest_file(path) != manifest[file]: raise ValueError('Changed raw geometry source')
        with gzip.open(path, 'rt') as h: d = MMCIF2Dict(h)
        with gzip.open(path, 'rt') as h: model = next(parser.get_structure(file, h).get_models())
        fields = [d['_atom_site.' + k] for k in ['auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code', 'auth_comp_id', 'label_seq_id', 'pdbx_PDB_model_num']]
        polymer = {(c, int(n), '' if ins in ['.', '?'] else ins, aa)
                   for c, n, ins, aa, label, m in zip(*fields) if label not in ['.', '?'] and int(m) == model.serial_num}
        raw_res = [r for r in model[chainid] if is_aa(r, standard=False)]
        residues = [r for r in raw_res if (chainid, r.id[1], r.id[2].strip(), r.resname) in polymer]
        sequence = ''.join(seq1(r.resname, custom_map={'MSE': 'M'}, undef_code='X') for r in residues)
        expected = {r['target_sequence_sha256'] for q, r in members}
        record = {'file': file, 'chain': chainid, 'source_sha256': manifest[file], 'sequence': sequence,
                  'sequence_sha256': stable(sequence), 'associated_core_sequences': [q for q, r in members],
                  'foldseek_identifier': members[0][1]['foldseek_identifier'], 'polymer_residues': len(residues),
                  'nonpolymer_amino_acids_removed': len(raw_res) - len(residues), 'eligible_reference': False,
                  'coordinate_policy': 'First model; Bio.PDB selected alternative atoms; polymer identity from raw label_seq_id'}
        if expected != {stable(sequence)}:
            census.append({**record, 'reason': 'polymer_sequence_differs_from_association'}); continue
        hemes = [(c.id, r) for c in model for r in c if r.resname in {'HEM', 'HEC', 'HEA'} and any(a.element.upper() == 'FE' for a in r.get_atoms())]
        pairs = []
        for pos, res in enumerate(residues):
            if res.resname != 'CYS' or 'SG' not in res: continue
            for hc, heme in hemes:
                for iron in [a for a in heme.get_atoms() if a.element.upper() == 'FE']:
                    distance = float(np.linalg.norm(res['SG'].coord - iron.coord))
                    if distance <= 3.0: pairs.append((distance, hc, str(heme.id), pos, heme, iron))
        if not pairs:
            census.append({**record, 'reason': 'no_qualified_axial_cysteine_iron_pair'}); continue
        distance, hc, hid, axial, heme, iron = min(pairs, key=lambda r: (r[0], r[1], r[2], r[3]))
        heme_xyz = np.array([a.coord for a in heme.get_atoms() if a.element.upper() not in {'H', 'D'}], dtype=float)
        local = []
        for pos, res in enumerate(residues):
            heavy = [a for a in res.get_atoms() if a.element.upper() not in {'H', 'D'}]
            if not heavy: continue
            xyz = np.array([a.coord for a in heavy], dtype=float)
            dmin = float(np.linalg.norm(xyz[:, None, :] - heme_xyz[None, :, :], axis=2).min())
            if dmin <= 5.0:
                local.append({'file': file, 'chain': chainid, 'position_0based': pos, 'auth_residue_id': list(res.id),
                              'aa': sequence[pos], 'residue_name': res.resname, 'minimum_heme_distance_A': dmin,
                              'minimum_iron_distance_A': float(np.linalg.norm(xyz - iron.coord, axis=1).min()),
                              'CA_present': 'CA' in res, 'CA_altloc': res['CA'].altloc if 'CA' in res else None,
                              'polymer_verified': True})
        exxr = [m.start() for m in re.finditer(r'(?=E..R)', sequence)]
        perf = [m.start() for m in re.finditer(r'PERF', sequence)]
        triples = [(e, p, axial) for e in exxr for p in perf if e + 3 < p and p + 3 < axial]
        record.update({'axial_cysteine_position_0based': axial, 'heme_chain': hc, 'heme_residue_id': list(heme.id),
                       'heme_component': heme.resname, 'sulfur_iron_distance_A': distance, 'geometric_anchor_pair_count': len(pairs),
                       'heme_contact_positions': [r['position_0based'] for r in local], 'contact_count': len(local),
                       'ordered_motif_triples': triples})
        contacts.extend(local)
        if len(triples) != 1: record['reason'] = 'no_unique_exact_ordered_motif_triple'
        elif not local or not all(r['CA_present'] for r in local): record['reason'] = 'missing_contact_CA'
        else:
            record.update({'eligible_reference': True, 'reason': 'eligible'}); eligible.append((record, local))
        census.append(record)
    out.mkdir(); write_json(out / 'geometry_census.json', census); write_json(out / 'polymer_heme_contacts.json', contacts)
    if not eligible:
        write_json(out / 'audit.json', {'status': 'NO_QUALIFIED_REFERENCE', 'created_utc': now(), 'counts': dict(Counter(r['reason'] for r in census))})
        raise SystemExit(2)
    reference, local = sorted(eligible, key=lambda x: (x[0]['file'], x[0]['chain']))[0]
    reference = {**reference, 'ordered_anchor_positions': [r['position_0based'] for r in local], 'ordered_contacts': local,
                 'selection': 'Lexicographically first geometrically/motif-qualified file-chain; no functional scores read'}
    write_json(out / 'reference.json', reference)
    with (out / 'reference.fasta').open('w') as h: h.write('>' + reference['foldseek_identifier'] + '\n' + reference['sequence'] + '\n')
    keys = []
    for line in (root / 'foldseek_raw_01/structures.lookup').read_text().splitlines():
        cells = line.split('\t')
        if cells[1] == reference['foldseek_identifier']: keys.append(cells[0])
    if len(keys) != 1: raise ValueError('Reference Foldseek lookup key ambiguous')
    with (out / 'reference_foldseek_key.txt').open('w') as h: h.write(keys[0] + '\n')
    summary = {'status': 'REFERENCE_DEFINED_MAPPING_NOT_YET_EVALUATED', 'created_utc': now(), 'raw_primary_chains': len(census),
               'eligible_reference_chains': len(eligible), 'counts': dict(Counter(r['reason'] for r in census)),
               'reference_file': reference['file'], 'reference_chain': reference['chain'], 'ordered_positions': len(local),
               'SRS_or_channel_boundaries_certified': False, 'functional_model_trained': False}
    write_json(out / 'audit.json', summary)
    names = ['prepare_ordered_heme_anchor.py', 'ORDERED_HEME_ANCHOR_V1.md', 'structure_mapping_audit_02/primary_representatives.json',
             'structure_association_validation_02/INDEPENDENT_QC.json', 'structure_assets_01/source_manifest.json', 'foldseek_raw_01/structures.lookup']
    write_json(out / 'input_manifest.json', {n: digest_file(root / n) for n in names}); print(json.dumps(summary, indent=2))


if __name__ == '__main__': run(Path(__file__).resolve().parent)
