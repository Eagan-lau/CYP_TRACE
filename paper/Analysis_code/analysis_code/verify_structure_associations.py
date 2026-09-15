"""Explain excluded Foldseek entries and independently verify retained mappings.

No excluded sequence is silently rescued. Comparisons use the frozen primary
representatives until any alternate atom-selection rule is separately audited.
"""
from collections import Counter, defaultdict
from difflib import SequenceMatcher
import gzip
import json
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from run_raw import digest_file, write_json, now, stable


def run(root):
    out = root / 'structure_association_validation_01'
    if out.exists():
        raise FileExistsError(out)
    read = lambda name: json.loads((root / name).read_text())
    chains = [json.loads(s) for s in (root / 'structure_assets_01/chain_assets.jsonl').read_text().splitlines()]
    byfile = defaultdict(list)
    for c in chains:
        byfile[c['file']].append(c)
    fs = {r.id: str(r.seq) for r in SeqIO.parse(root / 'foldseek_raw_01/structure_sequences.fasta', 'fasta')}
    core = {r.id: str(r.seq) for r in SeqIO.parse(root / 'dataset_02/core_sequences.fasta', 'fasta')}
    seq = {r.id: str(r.seq) for r in SeqIO.parse(root / 'structure_mapping_01/resolved_chain_sequences.fasta', 'fasta')}
    manifest = {r['file']: r['sha256'] for r in read('structure_assets_01/source_manifest.json')}
    excluded = read('structure_mapping_audit_01/foldseek_chain_reconciliation.json')['unmatched']
    parser = MMCIFParser(QUIET=True)
    details = []; reasons = Counter(); cache = {}; problems = []
    for item in excluded:
        name = item['foldseek_identifier']; file = name[:4] + '.cif.gz'; expected = fs[name]
        if file not in cache:
            path = root.parent / 'data/raw/rcsb_cyp_ligand_templates_v1/mmcif' / file
            if digest_file(path) != manifest[file]:
                raise ValueError('Changed raw CIF: ' + file)
            with gzip.open(path, 'rt') as h:
                model = next(parser.get_structure(file, h).get_models())
            cache[file] = model
        model = cache[file]
        candidates = [c for c in model if name == file[:4] or name == file[:4] + '_' + c.id]
        tests = []
        for c in candidates:
            residues = [r for r in c if is_aa(r, standard=False)]
            variants = {
                'all_recognized_amino_acids': residues,
                'CA_present': [r for r in residues if 'CA' in r],
                'complete_N_CA_C_backbone': [r for r in residues if all(a in r for a in ['N', 'CA', 'C'])],
                'standard_amino_acids_only': [r for r in residues if is_aa(r, standard=True)],
            }
            sequences = {k: ''.join(seq1(r.resname, custom_map={'MSE': 'M'}, undef_code='X') for r in rs) for k, rs in variants.items()}
            matches = [k for k, s in sequences.items() if s == expected]
            full = sequences['all_recognized_amino_acids']
            diffs = [{'operation': op, 'raw_range_0based': [a, b], 'foldseek_range_0based': [x, y],
                      'raw': full[a:b], 'foldseek': expected[x:y]}
                     for op, a, b, x, y in SequenceMatcher(None, full, expected, autojunk=False).get_opcodes() if op != 'equal']
            tests.append({'chain': c.id, 'raw_length': len(full), 'exact_variant_matches': matches,
                          'differences': diffs, 'residues_without_CA': sum('CA' not in r for r in residues),
                          'nonstandard_residue_names': dict(Counter(r.resname for r in residues if not is_aa(r, standard=True)))})
        matches = [(r['chain'], v) for r in tests for v in r['exact_variant_matches']]
        if file not in byfile:
            reason = 'file_excluded_from_heme_assets'
        elif any(v == 'CA_present' for _, v in matches):
            reason = 'exact_after_CA_presence_filter'
        elif any(v == 'complete_N_CA_C_backbone' for _, v in matches):
            reason = 'exact_after_complete_backbone_filter'
        elif any(v == 'standard_amino_acids_only' for _, v in matches):
            reason = 'exact_after_standard_residue_filter'
        else:
            reason = 'unresolved_parser_difference_excluded'
        reasons[reason] += 1
        details.append({**item, 'file': file, 'reason': reason, 'raw_chain_tests': tests,
                        'admitted_to_primary_analysis': False})
    associations = read('structure_mapping_audit_01/qualified_chain_associations.json')
    selected = read('structure_mapping_audit_01/primary_representatives.json')
    nmap = 0; naa = 0
    for q, members in associations.items():
        primary = []
        for r in members:
            t = r['target_sequence_sha256']; mapping = r['aligned_position_map_0based']
            if stable(seq[t]) != t or fs[r['foldseek_identifier']] != seq[t]:
                problems.append('sequence_identity:' + q)
            positions_q = [a for a, b, match in mapping]; positions_t = [b for a, b, match in mapping]
            if positions_q != sorted(set(positions_q)) or positions_t != sorted(set(positions_t)):
                problems.append('non_monotone_mapping:' + q)
            for a, b, match in mapping:
                if (core[q][a] == seq[t][b]) != match:
                    problems.append('residue_mapping:' + q)
                naa += 1
            if r['primary_exact_aligned_identity']:
                if not all(flag for a, b, flag in mapping) or r['identity'] != 1 or r['query_coverage'] < .8 or r['resolved_chain_coverage'] < .95 or r['evalue'] > .001:
                    problems.append('primary_qualification:' + q)
                primary.append(r)
            nmap += 1
        if not primary:
            if q in selected: problems.append('unexpected_selected:' + q)
            continue
        # Independent staged elimination rather than importing the selector.
        primary = [r for r in primary if r['query_coverage'] == max(v['query_coverage'] for v in primary)]
        primary = [r for r in primary if r['resolved_chain_coverage'] == max(v['resolved_chain_coverage'] for v in primary)]
        meta = [r['experimental_metadata'] for r in primary]
        if all(m['experimental_methods'] == meta[0]['experimental_methods'] for m in meta) and all(m['resolution_A'] is not None for m in meta):
            best = min(m['resolution_A'] for m in meta)
            primary = [r for r in primary if r['experimental_metadata']['resolution_A'] == best]
        answer = sorted(primary, key=lambda r: (r['file'], r['chain']))[0]
        if q not in selected or selected[q] != answer: problems.append('representative:' + q)
    census = read('structure_mapping_audit_01/availability_census.json')
    if len(census) != len(core) or {r['sequence_sha256'] for r in census if r['primary_available']} != set(selected):
        problems.append('availability_census')
    report = {'status': 'PASS' if not problems else 'FAIL', 'created_utc': now(), 'problems': problems,
              'qualified_associations_verified': nmap, 'aligned_residue_pairs_verified': naa,
              'primary_representatives_verified': len(selected), 'excluded_entries': len(details),
              'excluded_reason_counts': dict(reasons), 'excluded_entries_rescued': 0,
              'scope': 'Sequence and representative-selection arithmetic; geometry and catalytic-function validity not certified'}
    out.mkdir(); write_json(out / 'excluded_entry_adjudication.json', details)
    write_json(out / 'INDEPENDENT_QC.json', report)
    names = ['verify_structure_associations.py', 'structure_assets_01/source_manifest.json',
             'structure_mapping_audit_01/qualified_chain_associations.json', 'structure_mapping_audit_01/primary_representatives.json',
             'foldseek_raw_01/structure_sequences.fasta', 'structure_mapping_01/resolved_chain_sequences.fasta', 'dataset_02/core_sequences.fasta']
    write_json(out / 'input_manifest.json', {n: digest_file(root / n) for n in names})
    print(json.dumps(report, indent=2))
    if problems: raise SystemExit(1)


if __name__ == '__main__':
    run(Path(__file__).resolve().parent)
