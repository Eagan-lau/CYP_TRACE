"""Conservative auxiliary-participant separation; not a learned label selector."""
from collections import Counter, defaultdict
from functools import lru_cache
import json

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from run_raw import stable

# These are biochemical definitions, not hyperparameters selected with scores.
ROLE_NAMES = {
    'FMNH2': ['FMNH2'], 'FMN': ['FMN'],
    'NADPH': ['NADPH(4-)', 'NADPH'], 'NADP': ['NADP+'],
    'NADH': ['NADH'], 'NAD': ['NAD+'],
    'FeS_reduced': ['[2Fe-2S]1+'], 'FeS_oxidized': ['[2Fe-2S]2+'],
    'heme_reduced': ['heme b'], 'heme_oxidized': ['Fe(III)-heme b'],
}
PAIRS = [('FMNH2', 'FMN'), ('NADPH', 'NADP'), ('NADH', 'NAD'),
         ('FeS_reduced', 'FeS_oxidized'), ('heme_reduced', 'heme_oxidized')]
INORGANIC = {'O': 'water', '[H+]': 'proton', 'O=O': 'dioxygen', 'OO': 'hydrogen_peroxide'}
APPROVED_MACRO_IDS = {'11964', '11965', '9998', '9999', '10000', '10001'}
MACRO_ROLES = {'11964': 'FMNH2', '11965': 'FMN', '9998': 'FeS_reduced',
               '9999': 'FeS_oxidized', '10000': 'FeS_reduced', '10001': 'FeS_oxidized'}
# Independently documented reactive-part identities resolve a source name error.
# These are not selected by majority frequency or by downstream performance.
REFERENCE_ANCHORS = {
    'FMNH2': {
        'smiles': 'Cc1cc2Nc3c([nH]c(=O)[nH]c3=O)N(C[C@H](O)[C@H](O)[C@H](O)COP([O-])([O-])=O)c2cc1C',
        'chebi_id': 'CHEBI:57618', 'formula': 'C17H21N4O9P', 'charge': -2,
        'source_url': 'https://www.rhea-db.org/rhea/46308',
        'identity_url': 'https://www.ebi.ac.uk/chebi/CHEBI:57618', 'consulted_date': '2026-09-09'},
    'FMN': {
        'smiles': 'C12=NC([N-]C(C1=NC=3C(N2C[C@@H]([C@@H]([C@@H](COP(=O)([O-])[O-])O)O)O)=CC(=C(C3)C)C)=O)=O',
        'chebi_id': 'CHEBI:58210', 'formula': 'C17H18N4O9P', 'charge': -3,
        'source_url': 'https://www.rhea-db.org/rhea/46308', 'consulted_date': '2026-09-09'},
}


@lru_cache(maxsize=30000)
def heavy_signature(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError('Invalid structure passed to participant classifier')
    return tuple(sorted(Counter(a.GetAtomicNum() for a in mol.GetAtoms() if a.GetAtomicNum() > 1).items()))


@lru_cache(maxsize=30000)
def helper_equivalence_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError('Invalid helper SMILES')
    if any(a.GetAtomicNum() in (26,) for a in mol.GetAtoms()):
        return 'metal_exact:' + Chem.MolToSmiles(mol, isomericSmiles=True)
    # Only for recognition of a predeclared helper. Never applied to main labels.
    uncharged = rdMolStandardize.Uncharger().uncharge(mol)
    standard = rdMolStandardize.TautomerEnumerator().Canonicalize(uncharged)
    return 'helper_only:' + Chem.MolToSmiles(standard, isomericSmiles=True)


class RoleRegistry:
    def __init__(self, db):
        self.exact, self.normalized, self.signatures = {}, {}, set()
        self.records, self.unavailable, self.conflicts = {}, {}, []
        for role, names in ROLE_NAMES.items():
            records = [dict(r) for name in names for r in db.execute(
                "SELECT DISTINCT name,canonical_smiles,formula_status FROM components WHERE name=? AND structure_status='valid' AND formula_status='match'", (name,))]
            keys = {helper_equivalence_key(r['canonical_smiles']) for r in records}
            if role in REFERENCE_ANCHORS:
                key = helper_equivalence_key(REFERENCE_ANCHORS[role]['smiles'])
                self.conflicts.extend(dict(r, expected_role=role,
                    reason='source_helper_name_structure_disagrees_with_reference_anchor')
                    for r in records if helper_equivalence_key(r['canonical_smiles']) != key)
                records = [r for r in records if helper_equivalence_key(r['canonical_smiles']) == key]
                records.append({'name': role, 'canonical_smiles': Chem.MolToSmiles(
                    Chem.MolFromSmiles(REFERENCE_ANCHORS[role]['smiles']), isomericSmiles=True),
                    'formula_status': 'authoritative_reference_anchor'})
                keys = {key}
            if len(keys) != 1:
                self.unavailable[role] = {'reason': 'zero_or_ambiguous_source_structures', 'equivalence_keys': sorted(keys)}
                continue
            self.records[role] = records
            key = next(iter(keys))
            if key in self.normalized and self.normalized[key] != role: raise ValueError('Helper role collision')
            self.normalized[key] = role
            for record in records:
                self.exact[record['canonical_smiles']] = role
                self.signatures.add(heavy_signature(record['canonical_smiles']))

    def classify(self, smiles):
        if smiles in INORGANIC: return INORGANIC[smiles]
        if smiles in self.exact: return self.exact[smiles]
        if heavy_signature(smiles) not in self.signatures: return None
        return self.normalized.get(helper_equivalence_key(smiles))

    def project(self, substrates, products):
        remaining = {'substrate': list(substrates), 'product': list(products)}
        removed, problems = [], []
        for side in remaining:
            for smiles in list(remaining[side]):
                if smiles in INORGANIC:
                    remaining[side].remove(smiles)
                    removed.append({'side': side, 'smiles': smiles, 'role': INORGANIC[smiles], 'rule': 'closed_inorganic_list'})
        for reduced, oxidized in PAIRS:
            for left_role, right_role in [(reduced, oxidized), (oxidized, reduced)]:
                left = [s for s in remaining['substrate'] if self.classify(s) == left_role]
                right = [s for s in remaining['product'] if self.classify(s) == right_role]
                # Do not silently erase unmatched or unequal carrier stoichiometry.
                if not left or not right: continue
                if len(left) != len(right):
                    problems.append('unequal_paired_carrier_stoichiometry')
                    continue
                for side, items, role in [('substrate', left, left_role), ('product', right, right_role)]:
                    for smiles in items:
                        remaining[side].remove(smiles)
                        removed.append({'side': side, 'smiles': smiles, 'role': role, 'rule': 'opposite_redox_carrier_pair'})
        for side, values in remaining.items():
            if not values: problems.append('empty_' + side + '_after_auxiliary_separation')
            if any(self.classify(s) for s in values): problems.append('unpaired_known_auxiliary_remains')
            if not any(any(a.GetAtomicNum() == 6 for a in Chem.MolFromSmiles(s).GetAtoms()) for s in values):
                problems.append('no_carbon_containing_' + side)
        if sorted(remaining['substrate']) == sorted(remaining['product']): problems.append('identity_main_transformation')
        sides = {side: sorted(values) for side, values in remaining.items()}
        return {'substrates': sides['substrate'], 'products': sides['product'],
                'reaction_key': 'MAIN:' + stable(json.dumps(sides, sort_keys=True)),
                'single_pair': len(sides['substrate']) == len(sides['product']) == 1,
                'removed_participants': removed, 'problems': sorted(set(problems)),
                'main_charge_stereochemistry_and_multiplicity_preserved': True}

    def manifest(self):
        return {'roles': self.records, 'unavailable_roles': self.unavailable, 'redox_pairs': PAIRS,
                'inorganic_exact_smiles': INORGANIC, 'recognized_macro_reactive_parts': sorted(APPROVED_MACRO_IDS),
                'main_molecules_never_uncharged_or_tautomer_standardized': True,
                'registry_source': 'formula-consistent original P450Rdb component records in fresh raw store',
                'authoritative_anchor_overrides': REFERENCE_ANCHORS,
                'source_name_structure_conflicts': self.conflicts,
                'registry_selected_using_model_scores': False}
