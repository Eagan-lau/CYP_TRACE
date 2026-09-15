"""Post-evaluation label/identity audits and fixed-score sensitivity analyses.

All inputs are local copies of the frozen reference and predictions. No fitting
occurs. Results are source-label evaluations, not assay-level truth adjudication.
"""
from collections import Counter, defaultdict
from functools import lru_cache
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize as std
from rdkit.Chem import rdMolDescriptors
from sklearn.metrics import average_precision_score
from revise_comparisons import ap_template, weighted_ap, write_tsv

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'audit_inputs'
OUT = ROOT / 'tables'
SEED = 20260915
METHODS = ['isoform_specific_knn', 'pooled_chemical_knn', 'isoform_logistic']
RDLogger.DisableLog('rdApp.warning')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((INPUT/name).read_text(encoding='utf-8'))


def jsonl(name):
    return [json.loads(s) for s in (INPUT/name).read_text(encoding='utf-8').splitlines() if s]


@lru_cache(maxsize=None)
def parent_identity(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {'parent_status': 'parse_failure', 'parent_smiles': '', 'parent_key': ''}
    mol = std.Cleanup(mol)
    mol = std.FragmentParent(mol)
    mol = std.Uncharger().uncharge(mol)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    Chem.RemoveStereochemistry(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    enumerator = std.TautomerEnumerator()
    result = enumerator.Enumerate(mol)
    status = str(result.status)
    if status != 'Completed':
        return {'parent_status': 'tautomer_' + status, 'parent_smiles': '', 'parent_key': '', 'parent_formula': formula}
    mol = enumerator.PickCanonical(result)
    Chem.RemoveStereochemistry(mol)
    notation = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    return {'parent_status': 'ok', 'parent_smiles': notation,
            'parent_key': Chem.MolToInchiKey(mol), 'parent_formula': formula}


@lru_cache(maxsize=None)
def steroid_topology(smiles):
    """Broad structural sensitivity flag: fused all-carbon 6-6-6-5 tetracycle.

    Bond orders are deliberately unrestricted. This flags a ring topology, not
    a pharmacological class, a label derivation or a biochemical error.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()
             if len(r) in (5, 6) and all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in r)]
    for group in itertools.combinations(rings, 4):
        if sorted(map(len, group)) != [5, 6, 6, 6] or len(set.union(*group)) != 17:
            continue
        edges = [(i,j) for i in range(4) for j in range(i) if len(group[i] & group[j]) == 2]
        if len(edges) == 3 and sorted(Counter(i for edge in edges for i in edge).values()) == [1,1,2,2]:
            return True
    return False


def expanded_table_rows(table):
    pending = {}
    result = []
    for tr in table.findall('./tbody/tr'):
        current = {}
        for column, (value, remaining) in list(pending.items()):
            current[column] = value
            if remaining == 1:
                del pending[column]
            else:
                pending[column] = (value, remaining-1)
        column = 0
        for cell in tr.findall('td'):
            while column in current:
                column += 1
            value = ''.join(cell.itertext()).strip()
            current[column] = value
            span = int(cell.get('rowspan', 1))
            if span > 1:
                pending[column] = (value, span-1)
            column += 1
        result.append([current.get(i, '') for i in range(6)])
    return result


def metrics(rows):
    isoforms = sorted({r['isoform'] for r in rows})
    values = {}
    counts = []
    for iso in isoforms:
        local = [r for r in rows if r['isoform'] == iso]
        y = [r['label'] for r in local]
        counts.append(dict(isoform=iso, labels=len(y), positives=sum(y)))
        if not sum(y):
            raise ValueError('No positives in ' + iso)
    for method in METHODS:
        if not all(method in r for r in rows):
            continue
        values[method] = float(np.mean([
            average_precision_score([r['label'] for r in rows if r['isoform'] == iso],
                                    [r[method] for r in rows if r['isoform'] == iso])
            for iso in isoforms]))
    return {'labels': len(rows), 'compounds': len({r['compound_inchikey'] for r in rows}),
            'scaffolds': len({r['scaffold_group'] for r in rows}), 'isoforms': len(isoforms),
            'macro_AP': values, 'per_isoform_counts': counts}


def ap_interval(rows, bootstrap):
    scaffolds = sorted({r['scaffold_group'] for r in rows})
    index = {s:i for i,s in enumerate(scaffolds)}
    templates = []
    for iso in sorted({r['isoform'] for r in rows}):
        local = [r for r in rows if r['isoform'] == iso]
        y = np.array([r['label'] for r in local])
        templates.append((ap_template(y,np.array([r[METHODS[0]] for r in local])),
                          ap_template(y,np.array([r[METHODS[1]] for r in local])),
                          np.array([index[r['scaffold_group']] for r in local])))
    rng = np.random.default_rng(SEED)
    deltas = []
    for _ in range(bootstrap):
        counts = rng.multinomial(len(scaffolds), np.full(len(scaffolds), 1/len(scaffolds))).astype(float)
        delta = [weighted_ap(a,counts[ix])-weighted_ap(b,counts[ix]) for a,b,ix in templates]
        if np.isfinite(delta).all():
            deltas.append(float(np.mean(delta)))
    return {'specific_minus_pooled_AP_95_interval': np.quantile(deltas,[.025,.975]).tolist(),
            'bootstrap_valid': len(deltas), 'bootstrap_requested': bootstrap}


def selection(rows, method, fraction):
    mask = np.zeros(len(rows), dtype=bool)
    for iso in sorted({r['isoform'] for r in rows}):
        indices = [i for i,r in enumerate(rows) if r['isoform'] == iso]
        order = sorted(indices,key=lambda i:(-rows[i][method],hashlib.sha256(rows[i]['compound_inchikey'].encode()).hexdigest()))
        mask[order[:math.ceil(len(order)*fraction)]] = True
    return mask


def release_values(rows, mask):
    out = []
    for iso in sorted({r['isoform'] for r in rows}):
        ix = [i for i,r in enumerate(rows) if r['isoform'] == iso]
        chosen = [i for i in ix if mask[i]]
        positives = sum(rows[i]['label'] for i in ix)
        hits = sum(rows[i]['label'] for i in chosen)
        out.append(dict(isoform=iso, labels=len(ix), positives=positives,
                        released=len(chosen), reported_positive_hits=hits,
                        precision=hits/len(chosen), recall=hits/positives))
    return out


def main(bootstrap):
    development = read('development_labels.json')
    external = read('external_labels.json')
    predictions = jsonl('external_predictions.jsonl')
    raw = jsonl('external_raw_rows.jsonl')
    strict = [dict(r) for r in predictions if r['strata']['scaffold_and_lineage_eligible']]
    assert len(strict) == 3035
    with (OUT/'strict_logistic_paired_predictions.tsv').open(encoding='utf-8') as f:
        logistic = {(r['isoform'],r['compound_inchikey']):r for r in csv.DictReader(f,delimiter='\t')}
    for r in strict:
        pair = logistic[(r['isoform'],r['compound_inchikey'])]
        assert r['label'] == int(pair['label'])
        assert abs(r['isoform_specific_knn']-float(pair['isoform_specific_knn'])) < 1e-14
        r['isoform_logistic'] = float(pair['isoform_logistic'])
    external_index = {(r['isoform'],r['compound_inchikey']):r for r in external}
    train_compounds = {r['compound_inchikey']:r['canonical_smiles'] for r in development}
    external_compounds = {r['compound_inchikey']:r['canonical_smiles'] for r in external}
    print('Auditing parent identities on both sides',flush=True)
    parents = {}
    parent_train = defaultdict(set)
    identities = []
    cache_path = OUT/'parent_identity_cache.json'
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    for kind, compounds in [('development',train_compounds),('external',external_compounds)]:
        for key,smiles in sorted(compounds.items()):
            if smiles not in cache:
                cache[smiles] = parent_identity(smiles)
                if len(cache) % 100 == 0:
                    cache_path.write_text(json.dumps(cache),encoding='utf-8')
            result = cache[smiles]
            parents[(kind,key)] = result
            if kind == 'development' and result['parent_status'] == 'ok':
                parent_train[result['parent_smiles']].add(key)
            identities.append(dict(dataset=kind, compound_inchikey=key, canonical_smiles=smiles,**result))
    cache_path.write_text(json.dumps(cache),encoding='utf-8')
    write_tsv(OUT/'chemical_parent_identity_map.tsv',identities)
    unresolved_training = [r for r in identities if r['dataset']=='development' and r['parent_status']!='ok']
    assert all(r.get('parent_formula') for r in unresolved_training), 'Training parse failure requires manual review'
    unresolved_training_formulas = {r['parent_formula'] for r in unresolved_training}
    matched = set()
    unresolved = set()
    for row in identities:
        if row['dataset'] != 'external':
            continue
        if row['parent_status'] != 'ok' or row.get('parent_formula') in unresolved_training_formulas:
            unresolved.add(row['compound_inchikey'])
        elif row['parent_smiles'] in parent_train:
            matched.add(row['compound_inchikey'])
    audit_rows=[]
    for r in external:
        key=r['compound_inchikey']; parent=parents[('external',key)]
        audit_rows.append(dict(isoform=r['isoform'],compound_inchikey=key,label=r['label'],
            exact_training_overlap=key in train_compounds, parent_training_overlap=key in matched,
            parent_status=parent['parent_status'], strict_external=r['strata']['scaffold_and_lineage_eligible'],
            conservative_identity_exclusion=key in unresolved,
            primary_external=r['strata']['lineage_eligible'],
            matched_training_ids=';'.join(sorted(parent_train.get(parent['parent_smiles'],set()))),
            names=';'.join(sorted({x['name'] for x in r['source_locators']})),
            sources=';'.join(r['sources'])))
    write_tsv(OUT/'external_identity_overlap_audit.tsv',audit_rows)

    table = ET.parse(INPUT/'ni_2025_fulltext.xml').find(".//table-wrap[@id='Tab2']/table")
    table_rows = expanded_table_rows(table)
    write_tsv(OUT/'source_table2_steroid_classification.tsv',[
        dict(zip(['class','name','source','drug_name','dose','interaction'],r)) for r in table_rows])
    explicit_indirect = {(r[1].casefold(),r[5].split()[0]) for r in table_rows
                         if 'inhibitor' in r[5].lower() or 'inducer' in r[5].lower()}
    semantics=[]
    strict_keys={(r['isoform'],r['compound_inchikey']) for r in strict}
    for r in external:
        names={x['name'].casefold().strip() for x in r['source_locators']}
        indirect = r['label']==0 and any((name,r['isoform']) in explicit_indirect for name in names)
        structural_flag = r['label']==0 and steroid_topology(r['canonical_smiles'])
        semantics.append(dict(isoform=r['isoform'],compound_inchikey=r['compound_inchikey'],label=r['label'],
            strict_external=(r['isoform'],r['compound_inchikey']) in strict_keys,
            evidence_category='source_reported_classification',
            assay_level_verification='unresolved_from_available_row_metadata',
            explicit_table2_indirect_role_match=indirect,steroid_topology_negative_sensitivity_flag=structural_flag,
            names=';'.join(sorted(names)),sources=';'.join(r['sources']),
            source_locators=json.dumps(r['source_locators'],ensure_ascii=False)))
    write_tsv(OUT/'external_label_semantics_audit.tsv',semantics)
    sem_index={(r['isoform'],r['compound_inchikey']):r for r in semantics}
    strict_indirect={k for k in strict_keys if sem_index[k]['explicit_table2_indirect_role_match']}
    strict_steroid={k for k in strict_keys if sem_index[k]['steroid_topology_negative_sensitivity_flag']}
    keep_parent=lambda r:r['compound_inchikey'] not in matched|unresolved
    collections={
        'strict_original':strict,
        'strict_parent_disjoint':[r for r in strict if keep_parent(r)],
        'strict_exclude_table2_indirect_role':[r for r in strict if (r['isoform'],r['compound_inchikey']) not in strict_indirect],
        'strict_exclude_steroid_topology_negatives':[r for r in strict if (r['isoform'],r['compound_inchikey']) not in strict_steroid],
        'strict_parent_disjoint_exclude_steroid_topology_negatives':[r for r in strict if keep_parent(r) and (r['isoform'],r['compound_inchikey']) not in strict_steroid],
        'primary_parent_disjoint':[r for r in predictions if r['strata']['lineage_eligible'] and keep_parent(r)],
    }
    sensitivities={}
    for name,rows in collections.items():
        print(f'Sensitivity {name}: {len(rows)} labels',flush=True)
        result=metrics(rows)
        result['removed_labels']=(3035 if name.startswith('strict') else 4301)-len(rows)
        if result['removed_labels'] or name=='strict_original':
            result.update(ap_interval(rows,bootstrap))
        sensitivities[name]=result
    releases=[]; per_isoform=[]
    for fraction in (.1,.25):
        for method in METHODS:
            mask=selection(strict,method,fraction)
            changed=[{**r,'label':1-r['label']} for r in strict]
            assert np.array_equal(mask,selection(changed,method,fraction)), 'Selection used labels'
            local=release_values(strict,mask)
            releases.append(dict(method=method,nominal_release=fraction,released=int(mask.sum()),labels=len(strict),
                micro_release=float(mask.mean()),reported_positive_hits=sum(r['reported_positive_hits'] for r in local),
                isoform_macro_precision=float(np.mean([r['precision'] for r in local])),
                isoform_macro_recall=float(np.mean([r['recall'] for r in local])),
                status='exploratory_fixed_score_comparison'))
            per_isoform.extend([dict(method=method,nominal_release=fraction,**r) for r in local])
    write_tsv(OUT/'early_release_model_comparison.tsv',releases)
    write_tsv(OUT/'early_release_per_isoform.tsv',per_isoform)
    for expected,observed in [('strict_original',sensitivities['strict_original']['macro_AP']['isoform_specific_knn'])]:
        assert abs(observed-.5428801073616204)<1e-12
    final=dict(status='COMPLETED_WITH_SOURCE_SEMANTICS_LIMIT',analysis_status='post_evaluation_fixed_scores',
        seed=SEED,rdkit_version=rdBase.rdkitVersion,models_refitted=False,original_strata_modified=False,
        development_compounds=len(train_compounds),external_compounds=len(external_compounds),
        identity_status_counts=dict(Counter(r['dataset']+':'+r['parent_status'] for r in identities)),
        exact_overlap_compounds=len(set(train_compounds)&set(external_compounds)),
        parent_overlap_compounds=len(matched),additional_parent_overlap_compounds=len(matched-set(train_compounds)),
        strict_additional_parent_overlap_labels=sum(r['compound_inchikey'] in matched for r in strict),
        strict_parent_unresolved_labels=sum(r['compound_inchikey'] in unresolved for r in strict),
        unresolved_training_compounds=len(unresolved_training),
        identity_uncertainty_rule='Exclude unresolved external identities and any external parent formula matching an unresolved training parent formula. Formula screening is conservative, not an identity match.',
        strict_label_counts=dict(Counter(str(r['label']) for r in strict)),
        assay_level_labels_verified=0,
        available_row_fields=sorted({k for r in raw for k in r}),
        explicit_table2_indirect_roles=sorted(explicit_indirect),
        strict_table2_indirect_role_matches=len(strict_indirect),
        strict_steroid_topology_negative_flags=len(strict_steroid),
        semantics_note='No assay-method/result fields support individual experimental negative certification. Table 2 role matches are specific cues; topology flags are conservative sensitivity probes, not label-error findings.',
        sensitivities=sensitivities,early_release=releases,
        input_sha256={p.name:digest(p) for p in sorted(INPUT.iterdir()) if p.is_file()},
        prediction_table_sha256=digest(OUT/'strict_logistic_paired_predictions.tsv'),
        code_sha256=digest(Path(__file__)))
    (OUT/'external_evidence_audit.json').write_text(json.dumps(final,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in final.items() if k not in ('input_sha256','available_row_fields','sensitivities')},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--bootstrap',type=int,default=5000)
    args=parser.parse_args()
    assert parent_identity('CCO')['parent_smiles']==parent_identity('CCO.[Na+]')['parent_smiles']
    assert parent_identity('C[C@H](O)F')['parent_smiles']==parent_identity('C[C@@H](O)F')['parent_smiles']
    assert parent_identity('O=c1cccc[nH]1')['parent_smiles']==parent_identity('Oc1ccccn1')['parent_smiles']
    main(args.bootstrap)
