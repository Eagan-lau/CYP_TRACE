"""Fresh raw-source extraction. No imports or inputs from legacy analyses."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import sqlite3
import sys
import tarfile
import time
import zipfile

MISSING = {'', '/', '-', 'na', 'n/a', 'none', 'null', 'nan', 'unknown'}
AA = set('ACDEFGHIKLMNPQRSTVWY')
ACCESSION = re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)(?:-\d+)?')


def now():
    return datetime.now(timezone.utc).isoformat()


def clean(x):
    s = str(x or '').strip()
    return '' if s.casefold() in MISSING else s


def digest_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def stable(value):
    return hashlib.sha256(value.encode('utf-8')).hexdigest()


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def formula_heavy(formula):
    text = clean(formula)
    if not text:
        return None
    # In C2H3O2- the 2 is an oxygen count, not a charge magnitude.
    # A single-element Fe3+ notation is ambiguous without a charge delimiter.
    if re.fullmatch(r'[A-Z][a-z]?\d+[+-]', text):
        return None
    text = re.sub(r'\^\d*[+-]$', '', text)
    text = re.sub(r'[+-]\d*$', '', text)
    pieces = re.findall(r'([A-Z][a-z]?)(\d*)', text)
    if not pieces or ''.join(a + b for a, b in pieces) != text:
        return None
    counts = Counter()
    from rdkit import Chem
    elements = {Chem.GetPeriodicTable().GetElementSymbol(n) for n in range(1, 119)}
    for element, count in pieces:
        if element not in elements:
            return None
        if element != 'H':
            counts[element] += int(count or 1)
    return dict(counts)


def normalize_structure(smiles, formula=''):
    from rdkit import Chem
    raw = clean(smiles)
    result = {'raw_smiles': raw, 'raw_formula': clean(formula), 'canonical_smiles': None,
              'structure_status': 'missing', 'formula_status': 'unavailable'}
    if not raw:
        return result
    mol = Chem.MolFromSmiles(raw)
    if mol is None:
        result['structure_status'] = 'invalid'
        return result
    if any(atom.GetAtomicNum() == 0 or atom.HasQuery() for atom in mol.GetAtoms()):
        result['structure_status'] = 'generic_query_not_exact_molecule'
        return result
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    result['canonical_smiles'] = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    result['structure_status'] = 'valid'
    expected = formula_heavy(formula)
    if expected is not None:
        actual = dict(Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1))
        result['formula_status'] = 'match' if expected == actual else 'heavy_atom_conflict'
    elif clean(formula):
        result['formula_status'] = 'unparsed'
    return result


def sequence_info(sequence):
    seq = re.sub(r'\s+', '', clean(sequence)).upper()
    return seq, stable(seq) if seq else '', bool(len(seq) >= 100 and set(seq) <= AA)


def accession_tokens(text):
    return sorted(set(ACCESSION.findall(clean(text).upper())))


def publications(text, doi=''):
    out = ['PMID:' + x for x in re.findall(r'(?<!\d)\d{5,9}(?!\d)', clean(text))]
    for value in re.findall(r'10\.\d{4,9}/[^\s;|]+', clean(doi), flags=re.I):
        out.append('DOI:' + value.lower().rstrip('.,'))
    return sorted(set(out))


def csv_encoding(path):
    # Legacy database CSVs contain single-byte characters. Detect once before
    # iteration so a late decoder failure cannot leave a partly imported file.
    data = Path(path).read_bytes()
    for encoding in ('utf-8-sig', 'cp1252', 'latin1'):
        try:
            data.decode(encoding, errors='strict')
            return encoding
        except UnicodeDecodeError:
            continue
    raise UnicodeError(f'No lossless supported decoding for {path}')


class RawBuild:
    def __init__(self, project, output):
        self.project = project.resolve()
        self.raw = self.project / 'data/raw'
        self.output = output.resolve()
        self.output.mkdir(parents=True, exist_ok=True)
        self.inputs = []
        self.summary = {'started_utc': now(), 'status': 'RUNNING', 'stages': {},
                        'independent_validation': False, 'legacy_scores_read': False,
                        'historical_exposure': 'DEVELOPMENT_EXPOSED_OR_UNVERIFIED',
                        'paper_complete': False}
        self.db_path = self.output / 'raw_rebuild.sqlite'
        if self.db_path.exists():
            raise FileExistsError(f'Refusing to overwrite {self.db_path}; select a new run directory')
        self.db = sqlite3.connect(self.db_path)
        self.db.executescript('''
        CREATE TABLE proteins(sequence_sha256 TEXT PRIMARY KEY, sequence TEXT, length INTEGER, model_sequence_valid INTEGER);
        CREATE TABLE protein_assertions(source TEXT, locator TEXT, source_path TEXT, source_sha256 TEXT,
            accession TEXT, sequence_sha256 TEXT, taxid TEXT, symbol TEXT, family_evidence TEXT, raw_json TEXT);
        CREATE TABLE assertions(assertion_id TEXT PRIMARY KEY, source TEXT, locator TEXT, source_path TEXT,
            source_sha256 TEXT, sequence_sha256 TEXT, accession TEXT, taxid TEXT, reaction_key TEXT,
            publication_ids TEXT, evidence_type TEXT, model_eligible INTEGER, exclusion_reasons TEXT, raw_json TEXT);
        CREATE TABLE reactions(reaction_key TEXT PRIMARY KEY, substrates TEXT, products TEXT, directed INTEGER);
        CREATE TABLE components(source TEXT, locator TEXT, role TEXT, slot INTEGER, name TEXT,
            raw_smiles TEXT, canonical_smiles TEXT, raw_formula TEXT, structure_status TEXT, formula_status TEXT);
        CREATE TABLE rhea_structures(rhea_id TEXT PRIMARY KEY, reaction_smiles TEXT, direction TEXT, master_id TEXT);
        CREATE TABLE cypstrate_rows(source_path TEXT, member TEXT, sheet TEXT, row_number INTEGER, row_json TEXT);
        ''')

    def log(self, message):
        print(f'[{now()}] {message}', flush=True)

    def source(self, path):
        path = Path(path).resolve()
        allowed = [self.raw.resolve(), (self.output.parent / 'raw_additions').resolve()]
        if not any(path.is_relative_to(root) for root in allowed):
            raise ValueError(f'Input is not inside an approved raw root: {path}')
        rel = str(path.relative_to(self.project)) if path.is_relative_to(self.project) else str(path)
        sha = digest_file(path)
        self.inputs.append({'path': rel, 'bytes': path.stat().st_size, 'sha256': sha,
                            'historical_exposure': 'DEVELOPMENT_EXPOSED_OR_UNVERIFIED'})
        return rel, sha

    def checkpoint(self, stage, result):
        self.db.commit()
        self.summary['stages'][stage] = result
        self.summary['updated_utc'] = now()
        write_json(self.output / 'status.json', self.summary)
        write_json(self.output / 'consumed_raw_inputs.json', self.inputs)
        self.log(f'{stage}: {json.dumps(result, ensure_ascii=False)}')

    def inventory(self):
        totals = defaultdict(lambda: {'files': 0, 'bytes': 0, 'empty_files': 0})
        errors = []
        def onerror(error):
            errors.append(str(error))
        with (self.output / 'raw_inventory.jsonl').open('w', encoding='utf-8') as stream:
            for directory, _, files in os.walk(self.raw, onerror=onerror, followlinks=False):
                for name in sorted(files):
                    path = Path(directory) / name
                    try:
                        stat = path.stat()
                        rel = path.relative_to(self.raw)
                        record = {'path': rel.as_posix(), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
                                  'classification': 'metadata_sidecar' if name.endswith('.curlmeta') else 'raw_tree_file_not_yet_adjudicated'}
                        stream.write(json.dumps(record) + '\n')
                        item = totals[rel.parts[0]]
                        item['files'] += 1
                        item['bytes'] += stat.st_size
                        item['empty_files'] += int(stat.st_size == 0)
                    except OSError as error:
                        errors.append(str(error))
        self.checkpoint('raw_inventory', {'resources': dict(sorted(totals.items())), 'errors': errors,
                                          'all_files_content_verified': False})

    def protein(self, source, locator, path, sha, accession, seq, taxid, symbol, family, raw):
        seq, seqhash, valid = sequence_info(seq)
        if seqhash:
            self.db.execute('INSERT OR IGNORE INTO proteins VALUES(?,?,?,?)', (seqhash, seq, len(seq), int(valid)))
        self.db.execute('INSERT INTO protein_assertions VALUES(?,?,?,?,?,?,?,?,?,?)',
                        (source, locator, path, sha, accession, seqhash, str(taxid or ''), symbol, family, json.dumps(raw)))
        return seqhash, valid

    def assertion(self, source, locator, path, sha, seqhash, accession, taxid, reaction_key,
                  pubs, evidence_type, reasons, raw):
        aid = stable('|'.join((source, path, locator)))
        self.db.execute('INSERT INTO assertions VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                        (aid, source, locator, path, sha, seqhash, accession, str(taxid or ''), reaction_key,
                         json.dumps(pubs), evidence_type, int(not reasons), json.dumps(sorted(set(reasons))), json.dumps(raw)))

    def p450rdb(self):
        folder = self.raw / 'p450rdb'
        candidates = list(folder.glob('*'))
        counts = Counter()
        for prefix, mode in [('2.', 'protein'), ('1.', 'reaction')]:
            matches = [p for p in candidates if p.name.startswith(prefix) and p.suffix.lower() == '.csv']
            if len(matches) != 1:
                raise ValueError(f'Expected one raw P450Rdb {mode} file, found {matches}')
            path = matches[0]
            rel, sha = self.source(path)
            encoding = csv_encoding(path)
            self.inputs[-1]['text_encoding'] = encoding
            with path.open(encoding=encoding, newline='') as stream:
                for number, row in enumerate(csv.DictReader(stream), 2):
                    counts[mode + '_rows'] += 1
                    locator = f'csv_row:{number}'
                    accessions = accession_tokens(row.get('Uniprot ID'))
                    acc = '|'.join(accessions)
                    taxid = clean(row.get('Txid'))
                    seqhash, seqvalid = self.protein('P450Rdb', mode + ':' + locator, rel, sha, acc,
                        row.get('Sequence' if mode == 'protein' else 'sequence'), taxid, clean(row.get('Symbol')),
                        'source_curated_P450_not_independent_domain_confirmation', row)
                    if mode == 'protein':
                        continue
                    reasons = []
                    if not seqvalid: reasons.append('missing_or_noncanonical_sequence')
                    if not taxid.isdigit(): reasons.append('missing_or_invalid_taxonomy')
                    if len(accessions) > 1: reasons.append('multiple_accessions_ambiguous_assignment')
                    pubs = publications(row.get('PMID'), row.get('DOI'))
                    if not pubs: reasons.append('missing_traceable_publication')
                    sides = {}
                    for role, prefix2, nmax in [('substrate', 'sub', 4), ('product', 'pro', 6)]:
                        values = []
                        for n in range(1, nmax + 1):
                            name = clean(row.get(('Substrate' if role == 'substrate' else 'Product') + str(n)))
                            smiles = clean(row.get(f'{prefix2}_Smiles{n}'))
                            cid = clean(row.get(f'{prefix2}_CID{n}'))
                            if not any((name, smiles, cid)):
                                continue
                            item = normalize_structure(smiles, row.get(f'{prefix2}_formula{n}'))
                            self.db.execute('INSERT INTO components VALUES(?,?,?,?,?,?,?,?,?,?)',
                                ('P450Rdb', locator, role, n, name, item['raw_smiles'], item['canonical_smiles'],
                                 item['raw_formula'], item['structure_status'], item['formula_status']))
                            counts['structure_' + item['structure_status']] += 1
                            counts['formula_' + item['formula_status']] += 1
                            if not item['canonical_smiles']:
                                reasons.append(role + '_incomplete_structure')
                            else:
                                values.append(item['canonical_smiles'])
                            if item['formula_status'] == 'heavy_atom_conflict':
                                reasons.append('heavy_atom_formula_conflict')
                        if not values: reasons.append(role + '_no_valid_structure')
                        sides[role] = sorted(values)
                    reaction_key = ''
                    if sides['substrate'] and sides['product'] and not any('incomplete_structure' in r for r in reasons):
                        reaction_key = 'RXN:' + stable(json.dumps(sides, sort_keys=True))
                        self.db.execute('INSERT OR IGNORE INTO reactions VALUES(?,?,?,1)',
                            (reaction_key, json.dumps(sides['substrate']), json.dumps(sides['product'])))
                        if sides['substrate'] == sides['product']: reasons.append('identity_reaction')
                    self.assertion('P450Rdb', locator, rel, sha, seqhash, acc, taxid, reaction_key, pubs,
                                   'source_curated_publication_reaction_construct_unverified', reasons, row)
                    counts['initial_model_eligible_rows'] += int(not reasons)
        self.checkpoint('p450rdb', dict(counts))

    def rhea(self):
        folder = self.raw / 'rhea/tsv'
        direction_path = folder / 'rhea-directions.tsv'
        self.source(direction_path)
        directions = {}
        with direction_path.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                master = row['RHEA_ID_MASTER']
                for column, direction in [('RHEA_ID_MASTER', 'unspecified'), ('RHEA_ID_LR', 'LR'),
                                          ('RHEA_ID_RL', 'RL'), ('RHEA_ID_BI', 'bidirectional')]:
                    directions[row[column]] = (master, direction)
        path = folder / 'rhea-reaction-smiles.tsv'
        self.source(path)
        count = 0
        with path.open() as stream:
            for row in csv.reader(stream, delimiter='\t'):
                if len(row) != 2 or not row[0].isdigit():
                    raise ValueError(f'Unexpected raw Rhea SMILES row: {row[:1]}')
                master, direction = directions.get(row[0], ('', 'unresolved'))
                self.db.execute('INSERT INTO rhea_structures VALUES(?,?,?,?)', (row[0], row[1], direction, master))
                count += 1
        self.checkpoint('rhea', {'direction_entries': len(directions), 'raw_directional_structures': count,
                               'master_direction_not_assumed': True, 'protein_truth_generated': False})

    def uniprot_api(self):
        counts = Counter()
        for path in sorted((self.raw / 'uniprot/api').glob('cyp_batch_*.json')):
            rel, sha = self.source(path)
            payload = json.loads(path.read_text())
            if not isinstance(payload.get('results'), list):
                raise ValueError(f'Not a UniProt results response: {path}')
            for entry in payload['results']:
                counts['entries'] += 1
                acc = entry['primaryAccession']
                refs = entry.get('uniProtKBCrossReferences', [])
                family = 'PF00067' if any(x.get('database') == 'Pfam' and x.get('id') == 'PF00067' for x in refs) else 'not_verified_PF00067'
                taxid = entry.get('organism', {}).get('taxonId', '')
                seqhash, _ = self.protein('UniProt_API', acc, rel, sha, acc, entry.get('sequence', {}).get('value'),
                    taxid, entry.get('uniProtkbId', ''), family, {'entryAudit': entry.get('entryAudit'), 'organism': entry.get('organism')})
                for n, comment in enumerate(entry.get('comments', [])):
                    if comment.get('commentType') != 'CATALYTIC ACTIVITY': continue
                    counts['catalytic_assertions'] += 1
                    reaction = comment.get('reaction', {})
                    evidence = reaction.get('evidences', [])
                    pubs = sorted({'PMID:' + str(e['id']) for e in evidence if e.get('source') == 'PubMed' and e.get('id')})
                    codes = sorted({e.get('evidenceCode', '') for e in evidence})
                    self.assertion('UniProt_API', f'{acc}:comment:{n}', rel, sha, seqhash, acc, taxid, '', pubs,
                        'experimental_code' if 'ECO:0000269' in codes else 'annotation_not_experimental_code',
                        ['pending_direction_and_exact_construct_adjudication'], comment)
        self.checkpoint('uniprot_api', dict(counts))

    def swissprot(self):
        from Bio import SwissProt
        candidates = [self.raw / 'uniprot/uniprot_sprot.dat.gz', self.raw / 'uniprot_swissprot/uniprot_sprot.dat.gz']
        path = next((p for p in candidates if p.is_file()), None)
        if path is None:
            self.checkpoint('swissprot', {'status': 'MISSING_RAW_FILE'})
            return
        rel, sha = self.source(path)
        counts = Counter()
        with gzip.open(path, 'rt', encoding='utf-8') as stream:
            for entry in SwissProt.parse(stream):
                counts['all_records_scanned'] += 1
                domain = any(x[0] == 'Pfam' and x[1] == 'PF00067' for x in entry.cross_references)
                name_match = bool(re.search(r'cytochrome\s+P450|cytochrome\s+P-450', entry.description, flags=re.I))
                if not domain and not name_match: continue
                counts['cyp_records'] += 1
                counts['pfam_confirmed'] += int(domain)
                acc = entry.accessions[0]
                taxid = str(entry.taxonomy_id[0]) if len(entry.taxonomy_id) == 1 else ''
                family = 'PF00067' if domain else 'name_only_family_unverified'
                seqhash, _ = self.protein('SwissProt', acc, rel, sha, '|'.join(entry.accessions), entry.sequence,
                    taxid, entry.entry_name, family, {'description': entry.description, 'organism': entry.organism,
                    'lineage': entry.organism_classification, 'sequence_update': entry.sequence_update,
                    'comments': entry.comments, 'cross_references': entry.cross_references})
                for n, comment in enumerate(entry.comments):
                    if not comment.startswith('CATALYTIC ACTIVITY:'): continue
                    counts['catalytic_assertions'] += 1
                    pubs = sorted({'PMID:' + x for x in re.findall(r'PubMed:(\d+)', comment)})
                    self.assertion('SwissProt', f'{acc}:comment:{n}', rel, sha, seqhash, acc, taxid, '', pubs,
                        'experimental_code' if 'ECO:0000269' in comment else 'annotation_not_experimental_code',
                        ['pending_direction_and_exact_construct_adjudication'], {'comment': comment})
                if counts['cyp_records'] % 500 == 0:
                    self.log(f'SwissProt scan: {counts["all_records_scanned"]} records; {counts["cyp_records"]} CYPs')
        self.checkpoint('swissprot', dict(counts))

    def sabio(self):
        accessions = set()
        for (value,) in self.db.execute('SELECT DISTINCT accession FROM protein_assertions'):
            accessions.update(value.split('|'))
        accessions.discard('')
        counts = Counter()
        seen = set()
        for path in sorted((self.raw / 'sabio_rk/export_csv').glob('*.csv')):
            rel, sha = self.source(path)
            encoding = csv_encoding(path)
            self.inputs[-1]['text_encoding'] = encoding
            with path.open(encoding=encoding, newline='') as stream:
                for n, row in enumerate(csv.DictReader(stream), 2):
                    counts['all_rows'] += 1
                    ids = accession_tokens(row.get('UniprotIDs'))
                    matches = sorted(set(ids) & accessions)
                    if not matches: continue
                    counts['cyp_linked_rows'] += 1
                    signature = stable(json.dumps(row, sort_keys=True))
                    if signature in seen:
                        counts['repeated_export_rows'] += 1
                        continue
                    seen.add(signature)
                    reasons = ['pending_chemical_identity_and_exact_construct_adjudication']
                    if len(ids) > 1: reasons.append('ambiguous_multiple_proteins')
                    self.assertion('SABIO_RK', f'row:{n}:entry:{row.get("EntryID", "")}', rel, sha, '', '|'.join(ids),
                        row.get('NCBITaxonomyID', ''), '', publications(row.get('PubMedID')),
                        'source_kinetic_record_not_condition_matched', reasons, row)
        self.checkpoint('sabio_rk', dict(counts))

    def brenda(self):
        folder = self.raw / 'brenda/2026_1'
        files = sorted(folder.glob('*.tar.gz'))
        archives = []
        for path in files:
            rel, sha = self.source(path)
            with tarfile.open(path, 'r:gz') as archive:
                members = [{'name': m.name, 'bytes': m.size} for m in archive.getmembers() if m.isfile()]
            archives.append({'path': rel, 'sha256': sha, 'members': members})
        self.checkpoint('brenda', {'archives': archives, 'semantic_extraction': 'PENDING',
                                  'ec_support_not_upgraded_to_protein_truth': True})

    def cypstrate(self):
        candidates = list((self.output.parent / 'raw_additions').glob('*.zip'))
        candidates += list((self.raw / 'cypstrate').glob('**/*.zip'))
        if not candidates:
            self.checkpoint('cypstrate', {'status': 'MISSING_RAW_ZIP'})
            return
        path = candidates[0]
        rel, sha = self.source(path)
        with zipfile.ZipFile(path) as archive:
            bad = archive.testzip()
            if bad: raise ValueError(f'Corrupt ZIP member: {bad}')
            members = archive.namelist()
            if importlib.util.find_spec('xlrd') is None:
                self.checkpoint('cypstrate', {'status': 'RAW_VERIFIED_PARSER_DEPENDENCY_MISSING', 'sha256': sha,
                    'members': members, 'missing_dependency': 'xlrd', 'labels_rebuilt': False})
                return
            import xlrd
            counts = Counter()
            schemas = []
            for member in members:
                if not member.lower().endswith('.xls'): continue
                book = xlrd.open_workbook(file_contents=archive.read(member))
                for sheet in book.sheets():
                    schemas.append({'member': member, 'sheet': sheet.name, 'rows': sheet.nrows,
                                    'columns': sheet.ncols, 'header': sheet.row_values(0)})
                    for rownum in range(1, sheet.nrows):
                        self.db.execute('INSERT INTO cypstrate_rows VALUES(?,?,?,?,?)',
                            (rel, member, sheet.name, rownum + 1, json.dumps(sheet.row_values(rownum))))
                        counts['raw_rows_extracted'] += 1
            self.checkpoint('cypstrate', {'status': 'RAW_ROWS_EXTRACTED_SEMANTIC_NORMALIZATION_PENDING',
                                        'counts': dict(counts), 'schemas': schemas})

    def finish(self):
        cur = self.db
        queries = {
            'unique_sequences': 'SELECT COUNT(*) FROM proteins',
            'source_protein_assertions': 'SELECT COUNT(*) FROM protein_assertions',
            'source_reaction_assertions': 'SELECT COUNT(*) FROM assertions',
            'p450rdb_initial_eligible_assertions': 'SELECT COUNT(*) FROM assertions WHERE model_eligible=1',
            'p450rdb_initial_unique_sequence_reactions': 'SELECT COUNT(*) FROM (SELECT DISTINCT sequence_sha256,reaction_key FROM assertions WHERE model_eligible=1)',
            'p450rdb_initial_labelled_sequences': 'SELECT COUNT(DISTINCT sequence_sha256) FROM assertions WHERE model_eligible=1',
        }
        result = {name: cur.execute(sql).fetchone()[0] for name, sql in queries.items()}
        result['assertions_by_source'] = dict(cur.execute('SELECT source,COUNT(*) FROM assertions GROUP BY source'))
        core = cur.execute('SELECT DISTINCT p.sequence_sha256,p.sequence FROM proteins p JOIN assertions a ON a.sequence_sha256=p.sequence_sha256 WHERE a.model_eligible=1 ORDER BY p.sequence_sha256').fetchall()
        with (self.output / 'initial_core_sequences.fasta').open('w', encoding='ascii') as stream:
            for sha, seq in core:
                stream.write(f'>{sha}\n{seq}\n')
        cur.executescript('CREATE INDEX assertion_seq ON assertions(sequence_sha256); CREATE INDEX protein_acc ON protein_assertions(accession);')
        cur.commit()
        result['sqlite_integrity'] = cur.execute('PRAGMA integrity_check').fetchone()[0]
        self.summary['status'] = 'RAW_EXTRACTION_COMPLETE_DOWNSTREAM_PENDING'
        self.checkpoint('initial_raw_rebuild_summary', result)
        cur.close()
        write_json(self.output / 'artifact_checksums.json', {p.name: digest_file(p) for p in
            [self.db_path, self.output / 'initial_core_sequences.fasta', self.output / 'consumed_raw_inputs.json']})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--skip-swissprot', action='store_true', help='Parser test only; not the complete raw build')
    args = ap.parse_args()
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.warning')
    RDLogger.DisableLog('rdApp.error')
    run = RawBuild(args.project, args.output)
    runtime = {'python': sys.version, 'executable': sys.executable, 'platform': platform.platform(),
               'hostname': platform.node(), 'argv': sys.argv, 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
               'source_code_sha256': digest_file(Path(__file__)), 'packages': {}}
    from importlib.metadata import version, PackageNotFoundError
    for package in ['rdkit', 'biopython', 'pandas', 'numpy', 'scipy', 'scikit-learn', 'xlrd', 'nbclient']:
        try: runtime['packages'][package] = version(package)
        except PackageNotFoundError: runtime['packages'][package] = None
    write_json(args.output / 'runtime.json', runtime)
    try:
        run.inventory()
        run.p450rdb()
        run.rhea()
        run.uniprot_api()
        if not args.skip_swissprot: run.swissprot()
        run.sabio()
        run.brenda()
        run.cypstrate()
        run.finish()
    except Exception as error:
        run.summary.update(status='FAILED', error_type=type(error).__name__, error=str(error), updated_utc=now())
        write_json(args.output / 'status.json', run.summary)
        write_json(args.output / 'consumed_raw_inputs.json', run.inputs)
        raise


if __name__ == '__main__':
    main()
