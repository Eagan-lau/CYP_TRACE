"""Frozen fixed-human CYP substrate-ranking route."""
from __future__ import annotations

import gzip
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold


def _stable(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _scaffold(molecule: Chem.Mol) -> str:
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    if scaffold.GetNumAtoms():
        return "ring:" + Chem.MolToSmiles(scaffold, isomericSmiles=False)
    return "acyclic_connectivity:" + Chem.MolToSmiles(molecule, isomericSmiles=False)


class HumanSubstrateModel:
    """Similarity-weighted kNN model frozen before external evaluation."""

    def __init__(self, bundle: dict):
        self.bundle = bundle
        algorithm = bundle["algorithm"]
        self.isoform_k = int(algorithm["isoform_k"])
        self.pooled_k = int(algorithm["pooled_k"])
        self.isoforms = tuple(bundle["scope"]["development_isoforms"])
        self.validated_isoforms = frozenset(bundle["scope"]["externally_validated_isoforms"])
        self.compounds = bundle["compounds"]
        self.ids = [row["compound_inchikey"] for row in self.compounds]
        self.smiles = [row["canonical_smiles"] for row in self.compounds]
        self.scaffolds = frozenset(row["scaffold_group"] for row in self.compounds)
        self.id_set = frozenset(self.ids)
        self._generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=int(algorithm["morgan_radius"]), fpSize=int(algorithm["morgan_bits"])
        )
        self._fingerprints = []
        labels = np.full((len(self.compounds), len(self.isoforms)), np.nan, dtype=np.float64)
        isoform_index = {name: index for index, name in enumerate(self.isoforms)}
        for row_number, row in enumerate(self.compounds):
            molecule = Chem.MolFromSmiles(row["canonical_smiles"])
            if molecule is None:
                raise ValueError(f"model bundle contains an invalid molecule at row {row_number}")
            self._fingerprints.append(self._generator.GetFingerprint(molecule))
            for isoform, label in row["labels"].items():
                labels[row_number, isoform_index[isoform]] = int(label)
        self.labels = labels
        self.pooled_labels = np.nanmean(labels, axis=1)
        self.pooled_fallback = float(np.nanmean(self.pooled_labels))
        self.isoform_fallback = np.nanmean(labels, axis=0)
        self.isoform_index = isoform_index
        self.tie = np.empty(len(self.ids), dtype=np.int64)
        for position, index in enumerate(sorted(range(len(self.ids)), key=lambda value: _stable(self.ids[value]))):
            self.tie[index] = position

    @lru_cache(maxsize=64)
    def _similarity_order(self, canonical_smiles: str) -> tuple[np.ndarray, np.ndarray]:
        molecule = Chem.MolFromSmiles(canonical_smiles)
        query_fp = self._generator.GetFingerprint(molecule)
        similarity = np.asarray(DataStructs.BulkTanimotoSimilarity(query_fp, self._fingerprints), dtype=np.float32)
        return similarity, np.lexsort((self.tie, -similarity))

    @classmethod
    def load(cls, path: Path | None = None) -> "HumanSubstrateModel":
        if path is None:
            path = Path(__file__).resolve().parent / "data" / "human_model_v1.json.gz"
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            bundle = json.load(stream)
        if bundle.get("schema_version") != "cyptrace-human-model-v1":
            raise ValueError("unsupported human model bundle schema")
        return cls(bundle)

    def score(self, smiles: str, isoform: str) -> dict:
        isoform = isoform.upper()
        if isoform not in self.isoform_index:
            return {
                "isoform": isoform,
                "input_smiles": smiles,
                "route": "ABSTAIN_UNSEEN_PROTEIN",
                "abstained": True,
                "reason": "The frozen model contains only nine fixed human CYP isoforms; a name match is required.",
                "score": None,
            }
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return {
                "isoform": isoform,
                "input_smiles": smiles,
                "route": "ABSTAIN_INVALID_CHEMISTRY",
                "abstained": True,
                "reason": "RDKit could not parse the supplied SMILES.",
                "score": None,
            }
        canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
        inchikey = Chem.MolToInchiKey(molecule)
        scaffold = _scaffold(molecule)
        similarity, order = self._similarity_order(canonical)

        pooled_chosen = order[: min(self.pooled_k, len(order))]
        pooled_weight = similarity[pooled_chosen].astype(np.float64)
        pooled_score = (
            float(np.dot(pooled_weight, self.pooled_labels[pooled_chosen]) / pooled_weight.sum())
            if pooled_weight.sum()
            else self.pooled_fallback
        )

        iso_index = self.isoform_index[isoform]
        valid = order[~np.isnan(self.labels[order, iso_index])]
        chosen = valid[: min(self.isoform_k, len(valid))]
        weight = similarity[chosen].astype(np.float64)
        score = (
            float(np.dot(weight, self.labels[chosen, iso_index]) / weight.sum())
            if weight.sum()
            else float(self.isoform_fallback[iso_index])
        )
        neighbours = [
            {
                "compound_inchikey": self.ids[index],
                "similarity": float(similarity[index]),
                "source_label": int(self.labels[index, iso_index]),
            }
            for index in chosen[:5]
        ]
        externally_validated = isoform in self.validated_isoforms
        return {
            "isoform": isoform,
            "input_smiles": smiles,
            "canonical_smiles": canonical,
            "compound_inchikey": inchikey,
            "scaffold_group": scaffold,
            "route": "FIXED_HUMAN_EXTERNALLY_VALIDATED" if externally_validated else "FIXED_HUMAN_DEVELOPMENT_ONLY",
            "abstained": False,
            "score": score,
            "score_kind": "uncalibrated_similarity_weighted_substrate_ranking_score",
            "pooled_chemistry_score": pooled_score,
            "isoform_conditioning_delta": score - pooled_score,
            "domain": {
                "fixed_human_isoform": True,
                "isoform_externally_validated": externally_validated,
                "exact_training_compound": inchikey in self.id_set,
                "training_scaffold_seen": scaffold in self.scaffolds,
                "maximum_training_tanimoto": float(similarity.max()),
                "new_protein_supported": False,
                "product_or_reaction_centre_supported": False,
            },
            "nearest_labelled_neighbours": neighbours,
            "interpretation": "Use for within-isoform ranking only; this score is not an assay-harmonized probability or a clinical decision.",
        }


RDLogger.DisableLog("rdApp.*")
