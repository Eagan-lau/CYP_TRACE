# Frozen protocol: nine-isoform human substrate endpoint v1

Status: frozen before any model performance from this endpoint was calculated.
The previously built source audit and scaffold assignments were visible; no
classifier, neighbour or ensemble result was available when this protocol and
its acceptance rules were fixed.

## Scientific question and claim boundary

For compounds whose Bemis--Murcko scaffold is absent from model training, does
an isoform-conditioned chemical-neighbour strategy improve the classification
of source-reported substrate/non-substrate labels for nine fixed human CYP
isoforms? A second question is whether shrinkage toward a chemistry-only model
improves upon the isoform-specific model alone.

The source is the CYPstrate supplementary collection (Molecules 2021, 26,
4678), which compiled binary labels for CYP1A2, CYP2A6, CYP2B6, CYP2C8, CYP2C9,
CYP2C19, CYP2D6, CYP2E1 and CYP3A4 from earlier collections. These labels are
used as reported. They are not reinterpreted as assay-harmonized catalytic
rates, products, reaction centres, pan-CYP negatives or evidence for a novel
protein. The original source train/test assignment is historical exposure and
is not an independent test here.

## Fixed data and denominators

- Input: `human_substrate_01/normalized_labels.json` made directly from the
  checksum-locked original XLS archive.
- All 14,955 unambiguous isoform--compound rows with a normalized structure are
  in the model denominator. The 549 raw rows without a resolvable structure
  remain in source accounting and are not silently repaired or scored.
- A compound and all its isoform labels occur in exactly one development fold.
  The existing five folds are assigned by scaffold component without labels or
  model performance. They are not changed after this protocol is frozen.
- Missing isoform labels are unknown and are not converted to class zero.

## Chemical representation and three prespecified methods

Each compound is represented by a radius-2, 2,048-bit Morgan fingerprint from
its normalized canonical SMILES. Pairwise chemical similarity is Tanimoto.
Neighbour ties are resolved only by a stable hash of the compound InChIKey.

For a query compound `c` and isoform `p`:

1. `pooled_chemical_knn` assigns every training compound the mean of all its
   available isoform labels, then takes the similarity-weighted mean among its
   `k_pool` closest training compounds. This score has chemistry but no query-
   isoform identity.
2. `isoform_specific_knn` uses only the training labels reported for isoform
   `p`, taking the similarity-weighted mean among its `k_iso` closest compounds.
3. `data_driven_shrinkage` is
   `(1 - w) * pooled_chemical_knn + w * isoform_specific_knn`.
   `w = 0` is a legal solution and means that the data do not support an
   isoform-specific correction.

If the selected neighbours have zero total Tanimoto similarity, the component
falls back to its corresponding outer-training prevalence. This is an explicit
defined-domain fallback, not a missing score treated as negative evidence.

## Nested parameter selection

The outer loop evaluates each of the five existing scaffold folds exactly once.
Within an outer-training set, each of the remaining four folds is used once as
inner validation and the other three as inner training. The grids are fixed:

- `k_pool` and `k_iso`: 1, 3, 5, 11, 25, 51 or 101;
- `w`: 0.0, 0.1, ..., 1.0.

The pooled and isoform-specific baselines select their own `k` from inner
out-of-fold predictions. The shrinkage strategy selects `k_pool`, `k_iso` and
`w` jointly. Selection maximizes the unweighted mean of average precision over
the nine isoforms. Exact ties prefer higher macro AUROC, then fewer total
neighbours, then smaller `w`. Outer labels are not inspected during selection.

## Metrics, uncertainty and falsification

All five outer predictions are pooled after each row has been scored exactly
once. Average precision is primary; ROC AUC and Brier score are secondary.
Metrics are reported per isoform and as an unweighted isoform macro average.
Five thousand bootstrap replicates resample scaffold components, not individual
rows or labels. Differences use paired resamples and the same scored rows.

To test whether alignment between an isoform and its conditional score matters,
the final shrinkage scores are permuted among the available isoforms within
each compound. Ninety-nine deterministic permutations retain the compound,
labels, score multiset and chemistry-only score. The plus-one upper-tail Monte
Carlo p value is reported; it is not called an exact p value and the model is
not refitted in this score-level falsification test.

Two conclusions have separate conjunctive gates:

- `isoform_condition_signal` requires a positive 95% scaffold-bootstrap lower
  bound for macro-average-precision difference versus `pooled_chemical_knn`,
  an upper-tail permutation p value no greater than 0.05, complete scoring of
  the 14,955 normalized rows, and at least 20 positives and 20 negatives for
  every isoform.
- `shrinkage_improves_both_components` additionally requires a positive 95%
  lower bound versus `isoform_specific_knn`.

Failure of either gate is retained as a result. No threshold, accuracy claim,
clinical interpretation or independent biological validation is inferred from
this development endpoint.

## Reproducibility gate

An independent program must verify input hashes, regenerate all fingerprints
and similarities, confirm compound/scaffold separation, reconstruct every outer
score and metric, repeat nested choices, and recompute the bootstrap and
permutation summaries before any number enters the manuscript.
