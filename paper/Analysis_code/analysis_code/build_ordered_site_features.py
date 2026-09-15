"""Positional residue identities from reference geometry, with no residue pooling."""
from collections import Counter
import csv
import json
from pathlib import Path
import numpy as np
from Bio import SeqIO
from run_raw import stable, write_json, digest_file, now


def position_map(query, target, qstart, qend, tstart, tend, qa, ta):
    if len(qa) != len(ta) or query[qstart - 1:qend] != qa.replace('-', '') or target[tstart - 1:tend] != ta.replace('-', ''):
        raise ValueError('Aligned residues do not reconstruct source')
    q = qstart - 1; t = tstart - 1; mapping = {}
    for a, b in zip(qa, ta):
        if a != '-' and b != '-': mapping[t] = q
        q += a != '-'; t += b != '-'
    return mapping


def project(query, mapping, anchor, motifs):
    positions = [mapping.get(i, -1) for i in anchor]
    aa = [query[i] if i >= 0 else '-' for i in positions]
    known = sum(a in 'ACDEFGHIKLMNPQRSTVWY' for a in aa) / len(anchor)
    e, p, c = motifs
    motif_positions = list(range(e, e + 4)) + list(range(p, p + 4)) + [c]
    mapped_motifs = all(i in mapping for i in motif_positions)
    conserved = mapped_motifs and query[mapping[e]] == 'E' and query[mapping[e + 3]] == 'R' and query[mapping[c]] == 'C'
    return {'query_positions_0based': positions, 'ordered_residues': ''.join(aa), 'known_anchor_fraction': known,
            'all_motif_positions_mapped': mapped_motifs, 'EXXR_and_axial_C_conserved': bool(conserved),
            'motif_query_positions': [mapping.get(i, -1) for i in motif_positions],
            'motif_query_residues': ''.join(query[mapping[i]] if i in mapping else '-' for i in motif_positions),
            'positional_qualification': bool(known >= .8 and conserved)}


def run(root):
    out = root / 'ordered_site_features_01'
    if out.exists(): raise FileExistsError(out)
    read = lambda p: json.loads((root / p).read_text())
    ref = read('ordered_anchor_01/reference.json'); anchor = ref['ordered_anchor_positions']; motifs = ref['ordered_motif_triples'][0]
    qseq = {r.id: str(r.seq) for r in SeqIO.parse(root / 'dataset_02/core_sequences.fasta', 'fasta')}; seqs = sorted(qseq)
    fsseq = {r.id: str(r.seq) for r in SeqIO.parse(root / 'foldseek_raw_01/structure_sequences.fasta', 'fasta')}
    selected = read('structure_mapping_audit_02/primary_representatives.json')
    refseq = ref['sequence']; t = ref['sequence_sha256']; fsid = ref['foldseek_identifier']
    sequence_maps = {}; structure_maps = {}; rawrows = {}; counts = Counter()
    alignment = root / 'ordered_anchor_01/structure_to_reference.tsv'
    if digest_file(alignment) != (root / 'ordered_anchor_01/alignment.sha256').read_text().split()[0]: raise ValueError('Structure alignment checksum')
    with alignment.open() as h:
        for r in csv.reader(h, delimiter='\t'):
            if len(r) != 14 or r[1] != fsid or r[0] in rawrows: raise ValueError('Unexpected reference structural alignment')
            if int(r[6]) != len(fsseq[r[0]]) or int(r[7]) != len(refseq): raise ValueError('Structure source length')
            mapping = position_map(fsseq[r[0]], refseq, *map(int, r[2:6]), r[12], r[13])
            rawrows[r[0]] = (r, mapping); counts['structural_backtraces_reconstructed'] += 1
    for q, s in selected.items():
        if s['foldseek_identifier'] not in rawrows: continue
        r, mapping = rawrows[s['foldseek_identifier']]
        inverse = {chainpos: qpos for qpos, chainpos, identical in s['aligned_position_map_0based'] if identical}
        coremap = {refpos: inverse[chainpos] for refpos, chainpos in mapping.items() if chainpos in inverse}
        for rp, qp in coremap.items():
            if fsseq[r[0]][mapping[rp]] != qseq[q][qp]: raise ValueError('Primary core-chain mapping mismatch')
        fields = project(qseq[q], coremap, anchor, motifs)
        eligible = float(r[10]) <= .001 and min(float(r[8]), float(r[9])) >= .5
        structure_maps[q] = {**fields, 'alignment_qualified': eligible, 'available': bool(eligible and fields['positional_qualification']),
                             'reference_coverage': float(r[9]), 'coordinate_chain_coverage': float(r[8]), 'evalue': float(r[10]),
                             'file': s['file'], 'chain': s['chain'], 'foldseek_identifier': s['foldseek_identifier'],
                             'source_type': 'direct_structure_supported_reference_position_projection_not_measured_substrate_contacts'}
    with (root / 'structure_blast_01/core_to_resolved_chains.tsv').open() as h:
        for r in csv.reader(h, delimiter='\t'):
            if r[1] != t: continue
            q = r[0]
            if q in sequence_maps: raise ValueError('Multiple reference BLAST HSPs')
            mapping = position_map(qseq[q], refseq, int(r[4]), int(r[5]), int(r[6]), int(r[7]), r[12], r[13])
            fields = project(qseq[q], mapping, anchor, motifs)
            tcov = (int(r[7]) - int(r[6]) + 1) / int(r[11]); qcov = (int(r[5]) - int(r[4]) + 1) / int(r[10])
            eligible = float(r[8]) <= .001 and tcov >= .5
            sequence_maps[q] = {**fields, 'alignment_qualified': eligible, 'available': bool(eligible and fields['positional_qualification']),
                                'reference_coverage': tcov, 'query_coverage': qcov, 'evalue': float(r[8]),
                                'source_type': 'sequence_alignment_inferred_positions_not_query_geometry'}
            counts['sequence_backtraces_reconstructed'] += 1
    alphabet = 'ACDEFGHIKLMNPQRSTVWYX'; arrays = {'sequence_ids': np.array(seqs), 'reference_positions_0based': np.array(anchor), 'amino_acid_alphabet': np.array(alphabet)}
    rows = []
    for name, maps in [('sequence_projected', sequence_maps), ('structure_projected', structure_maps)]:
        codes = np.full((len(seqs), len(anchor)), -1, dtype=np.int16); available = np.zeros(len(seqs), dtype=bool)
        onehot = np.zeros((len(seqs), len(anchor), len(alphabet)), dtype=np.float32)
        for qi, q in enumerate(seqs):
            r = maps.get(q)
            if r:
                available[qi] = r['available']
                for j, aa in enumerate(r['ordered_residues']):
                    if aa != '-':
                        code = alphabet.index(aa) if aa in alphabet else 20; codes[qi, j] = code; onehot[qi, j, code] = 1
            rows.append({'sequence_sha256': q, 'channel': name, 'projection': r, 'available': bool(r and r['available'])})
        arrays[name + '_residue_codes'] = codes; arrays[name + '_present_mask'] = codes >= 0
        arrays[name + '_available'] = available; arrays[name + '_onehot'] = onehot
        counts[name + '_available_sequences'] = int(available.sum()); counts[name + '_mapped_sequences'] = len(maps)
    out.mkdir(); np.savez_compressed(out / 'ordered_features.npz', **arrays); write_json(out / 'projection_records.json', rows)
    audit = {'status': 'FEATURES_COMPUTED_INDEPENDENT_MAPPING_QC_PENDING', 'created_utc': now(), 'sequences': len(seqs),
             'ordered_positions': len(anchor), 'onehot_channels_per_position': len(alphabet), 'counts': dict(counts),
             'reference_file': ref['file'], 'reference_chain': ref['chain'], 'averaged_across_positions': False,
             'functional_model_trained': False, 'actual_query_catalytic_contacts_certified': False, 'independent_validation': False}
    write_json(out / 'audit.json', audit)
    files = ['build_ordered_site_features.py', 'ORDERED_SITE_PROJECTION_V1.md', 'ordered_anchor_01/reference.json',
             'ordered_anchor_01/structure_to_reference.tsv', 'structure_blast_01/core_to_resolved_chains.tsv',
             'structure_mapping_audit_02/primary_representatives.json', 'dataset_02/core_sequences.fasta', 'foldseek_raw_01/structure_sequences.fasta']
    write_json(out / 'input_manifest.json', {f: digest_file(root / f) for f in files}); print(json.dumps(audit, indent=2))


if __name__ == '__main__': run(Path(__file__).resolve().parent)
