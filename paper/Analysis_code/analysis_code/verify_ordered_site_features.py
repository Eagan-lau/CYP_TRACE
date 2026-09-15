"""Independent alignment-column indexing and raw anchor-distance verification."""
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from run_raw import write_json, digest_file, now


def column_mapping(query, target, row, columns):
    qstart, qend, tstart, tend, qa, ta = [row[i] for i in columns]
    qstart, qend, tstart, tend = map(int, [qstart, qend, tstart, tend])
    if qa.replace('-', '') != query[qstart - 1:qend] or ta.replace('-', '') != target[tstart - 1:tend]: raise ValueError('Source sequence mismatch')
    qpresent = np.array(list(qa)) != '-'; tpresent = np.array(list(ta)) != '-'
    qindex = np.cumsum(qpresent) + qstart - 2; tindex = np.cumsum(tpresent) + tstart - 2
    return {int(tindex[i]): int(qindex[i]) for i in np.flatnonzero(qpresent & tpresent)}


def run(root):
    out = root / 'ordered_site_validation_01'
    if out.exists(): raise FileExistsError(out)
    read = lambda n: json.loads((root / n).read_text())
    failures = []
    for n, h in read('ordered_site_features_01/input_manifest.json').items():
        if digest_file(root / n) != h: failures.append('source:' + n)
    ref = read('ordered_anchor_01/reference.json'); anchor = ref['ordered_anchor_positions']; reference = ref['sequence']
    qseq = {r.id: str(r.seq) for r in SeqIO.parse(root / 'dataset_02/core_sequences.fasta', 'fasta')}; seqs = sorted(qseq)
    fsseq = {r.id: str(r.seq) for r in SeqIO.parse(root / 'foldseek_raw_01/structure_sequences.fasta', 'fasta')}
    selected = read('structure_mapping_audit_02/primary_representatives.json')
    records = {(r['channel'], r['sequence_sha256']): r for r in read('ordered_site_features_01/projection_records.json')}
    matrices = np.load(root / 'ordered_site_features_01/ordered_features.npz', allow_pickle=False)
    if list(matrices['sequence_ids']) != seqs or list(matrices['reference_positions_0based']) != anchor: failures.append('feature_order')
    maps = {}; alignment_qualified = {}; rawmaps = {}; backtraces = 0
    with (root / 'ordered_anchor_01/structure_to_reference.tsv').open() as h:
        for r in csv.reader(h, delimiter='\t'):
            if r[1] != ref['foldseek_identifier']: raise ValueError('Wrong reference target')
            rawmaps[r[0]] = (column_mapping(fsseq[r[0]], reference, r, [2, 3, 4, 5, 12, 13]), r)
            backtraces += 1
    for q, s in selected.items():
        if s['foldseek_identifier'] not in rawmaps: continue
        mapping, r = rawmaps[s['foldseek_identifier']]; inverse = {t: p for p, t, same in s['aligned_position_map_0based'] if same}
        key = ('structure_projected', q)
        maps[key] = {p: inverse[t] for p, t in mapping.items() if t in inverse}
        alignment_qualified[key] = float(r[10]) <= .001 and float(r[8]) >= .5 and float(r[9]) >= .5
    with (root / 'structure_blast_01/core_to_resolved_chains.tsv').open() as h:
        for r in csv.reader(h, delimiter='\t'):
            if r[1] != ref['sequence_sha256']: continue
            key = ('sequence_projected', r[0]); maps[key] = column_mapping(qseq[r[0]], reference, r, [4, 5, 6, 7, 12, 13])
            alignment_qualified[key] = float(r[8]) <= .001 and (int(r[7]) - int(r[6]) + 1) / int(r[11]) >= .5
            backtraces += 1
    alphabet = 'ACDEFGHIKLMNPQRSTVWYX'; e, p, c = ref['ordered_motif_triples'][0]
    motif_positions = list(range(e, e + 4)) + list(range(p, p + 4)) + [c]; available = defaultdict(int); positions_checked = 0
    for channel in ['sequence_projected', 'structure_projected']:
        if matrices[channel + '_onehot'].shape != (600, len(anchor), 21): failures.append('onehot_shape')
        for qi, q in enumerate(seqs):
            key = (channel, q); r = records[key]['projection']; mapping = maps.get(key)
            expected_codes = []; expected_onehot = np.zeros((len(anchor), 21)); qualified = False
            if mapping is None:
                if r is not None: failures.append('unexpected_projection')
                aa = '-' * len(anchor)
            else:
                indices = [mapping.get(j, -1) for j in anchor]; aa = ''.join(qseq[q][i] if i >= 0 else '-' for i in indices)
                if r is None or r['query_positions_0based'] != indices or r['ordered_residues'] != aa: failures.append('position_projection')
                mapped = all(j in mapping for j in motif_positions)
                conserved = mapped and qseq[q][mapping[e]] == 'E' and qseq[q][mapping[e + 3]] == 'R' and qseq[q][mapping[c]] == 'C'
                known = sum(a in alphabet[:-1] for a in aa) / len(anchor)
                qualified = bool(alignment_qualified[key] and conserved and known >= .8)
                if r['motif_query_positions'] != [mapping.get(j, -1) for j in motif_positions]: failures.append('motif_mapping')
            for j, residue in enumerate(aa):
                code = -1 if residue == '-' else alphabet.index(residue) if residue in alphabet else 20
                expected_codes.append(code)
                if code >= 0: expected_onehot[j, code] = 1
            if not np.array_equal(matrices[channel + '_residue_codes'][qi], expected_codes) or not np.array_equal(matrices[channel + '_onehot'][qi], expected_onehot): failures.append('positional_encoding')
            if not np.array_equal(matrices[channel + '_present_mask'][qi], np.array(expected_codes) >= 0): failures.append('missing_mask')
            if matrices[channel + '_available'][qi] != qualified or records[key]['available'] != qualified: failures.append('availability')
            available[channel] += int(qualified); positions_checked += len(anchor)
    path = root.parent / 'data/raw/rcsb_cyp_ligand_templates_v1/mmcif' / ref['file']
    if digest_file(path) != ref['source_sha256']: failures.append('raw_reference_hash')
    with gzip.open(path, 'rt') as h: d = MMCIF2Dict(h)
    with gzip.open(path, 'rt') as h: model = next(MMCIFParser(QUIET=True).get_structure(ref['file'], h).get_models())
    fields = [d['_atom_site.' + k] for k in ['auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code', 'auth_comp_id', 'label_seq_id', 'pdbx_PDB_model_num']]
    polymer = {(chain, int(number), '' if ins in ['.', '?'] else ins, aa) for chain, number, ins, aa, label, m in zip(*fields)
               if label not in ['.', '?'] and int(m) == model.serial_num}
    residues = [r for r in model[ref['chain']] if is_aa(r, standard=False) and (ref['chain'], r.id[1], r.id[2].strip(), r.resname) in polymer]
    if ''.join(seq1(r.resname, custom_map={'MSE': 'M'}, undef_code='X') for r in residues) != reference: failures.append('polymer_reference_sequence')
    heme = model[ref['heme_chain']][tuple(ref['heme_residue_id'])]
    xyz = np.array([a.coord for a in heme.get_atoms() if a.element.upper() not in ['H', 'D']], dtype=float)
    measured = []; geometry_error = 0.
    for j, res in enumerate(residues):
        heavy = np.array([a.coord for a in res.get_atoms() if a.element.upper() not in ['H', 'D']], dtype=float)
        if not len(heavy): continue
        distance = float(cdist(heavy, xyz).min())
        if distance <= 5:
            measured.append(j)
            match = next((r for r in ref['ordered_contacts'] if r['position_0based'] == j), None)
            if match is None: failures.append('reference_contact_missing')
            else: geometry_error = max(geometry_error, abs(distance - match['minimum_heme_distance_A']))
    if measured != anchor or geometry_error > 1e-10: failures.append('raw_geometry_reconstruction')
    out.mkdir(); report = {'status': 'PASS' if not failures else 'FAIL', 'created_utc': now(), 'failures': failures,
                          'backtraces_independently_indexed': backtraces, 'positional_cells_verified': positions_checked,
                          'available_sequences': dict(available), 'raw_reference_contacts_recomputed': len(measured),
                          'maximum_geometry_difference_A': geometry_error, 'functional_accuracy_certified': False,
                          'universal_SRS_or_channel_identity_certified': False, 'verifier_sha256': digest_file(Path(__file__))}
    write_json(out / 'INDEPENDENT_QC.json', report); print(json.dumps(report, indent=2))
    if failures: raise SystemExit(1)


if __name__ == '__main__': run(Path(__file__).resolve().parent)
