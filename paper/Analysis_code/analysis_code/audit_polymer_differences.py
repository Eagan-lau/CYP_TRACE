"""Second-stage source adjudication: amino-acid ligands are not chain residues."""
from collections import Counter
import gzip
import json
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from run_raw import digest_file, write_json, now


def run(root):
    out = root / 'structure_polymer_adjudication_01'
    if out.exists(): raise FileExistsError(out)
    read = lambda f: json.loads((root / f).read_text())
    entries = read('structure_association_validation_01/excluded_entry_adjudication.json')
    fs = {r.id: str(r.seq) for r in SeqIO.parse(root / 'foldseek_raw_01/structure_sequences.fasta', 'fasta')}
    manifest = {r['file']: r['sha256'] for r in read('structure_assets_01/source_manifest.json')}
    parser = MMCIFParser(QUIET=True); cache = {}; rows = []; counts = Counter()
    for entry in entries:
        file = entry['file']; name = entry['foldseek_identifier']
        if file not in cache:
            path = root.parent / 'data/raw/rcsb_cyp_ligand_templates_v1/mmcif' / file
            if digest_file(path) != manifest[file]: raise ValueError('Raw source changed')
            with gzip.open(path, 'rt') as h: d = MMCIF2Dict(h)
            with gzip.open(path, 'rt') as h: model = next(parser.get_structure(file, h).get_models())
            fields = [d['_atom_site.' + k] for k in ['auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code', 'auth_comp_id', 'label_seq_id', 'pdbx_PDB_model_num']]
            polymer = {(c, int(n), '' if ins in ['.', '?'] else ins, aa)
                       for c, n, ins, aa, label, m in zip(*fields)
                       if label not in ['.', '?'] and int(m) == int(model.serial_num)}
            cache[file] = model, polymer
        model, polymer = cache[file]; tests = []
        for chain in model:
            if not (name == file[:4] or name == file[:4] + '_' + chain.id): continue
            allaa = [r for r in chain if is_aa(r, standard=False)]
            isaa = lambda r: (chain.id, r.id[1], r.id[2].strip(), r.resname) in polymer
            aa = [r for r in allaa if isaa(r)]
            variants = {'polymer_residues': aa, 'CA_present_polymer_residues': [r for r in aa if 'CA' in r]}
            exact = [k for k, residues in variants.items()
                     if ''.join(seq1(r.resname, custom_map={'MSE': 'M'}, undef_code='X') for r in residues) == fs[name]]
            tests.append({'chain': chain.id, 'exact_variant_matches': exact,
                          'nonpolymer_amino_acids': [{'residue_id': list(r.id), 'resname': r.resname} for r in allaa if not isaa(r)]})
        matched = any(t['exact_variant_matches'] for t in tests)
        nonpoly = any(t['nonpolymer_amino_acids'] for t in tests)
        if entry['reason'] == 'file_excluded_from_heme_assets': reason = 'no_recognized_heme_file'
        elif matched and nonpoly: reason = 'nonpolymer_amino_acid_ligands_with_optional_CA_filter'
        elif matched: reason = 'CA_presence_filter'
        else: reason = 'unresolved_excluded'
        counts[reason] += 1
        rows.append({'foldseek_identifier': name, 'first_audit_reason': entry['reason'], 'reason': reason, 'tests': tests,
                     'raw_polymer_definition': 'atom_site.label_seq_id not . or ? in first deposited model', 'admitted_to_primary': False})
    out.mkdir(); write_json(out / 'entries.json', rows)
    report = {'status': 'SOURCE_ADJUDICATION_COMPLETE', 'created_utc': now(), 'counts': dict(counts), 'entries': len(rows),
              'primary_set_unchanged': True, 'rescued_entries': 0,
              'geometry_warning': 'is_aa alone includes amino-acid ligands. Future local features must use polymer label_seq_id and retain nonpolymer ligand contacts separately; no old contact row is automatically a protein site.'}
    write_json(out / 'audit.json', report)
    write_json(out / 'input_manifest.json', {n: digest_file(root / n) for n in ['audit_polymer_differences.py',
               'structure_association_validation_01/excluded_entry_adjudication.json', 'structure_assets_01/source_manifest.json',
               'foldseek_raw_01/structure_sequences.fasta']})
    print(json.dumps(report, indent=2))


if __name__ == '__main__': run(Path(__file__).resolve().parent)
