# Candidate-specific ordered-site × reaction-centre model

This development protocol is frozen after the lineage-control results and before
candidate-specific interaction performance is calculated. It replaces the
incorrect interpretation of the centre-profile 01 branch. It remains a
development analysis on historically exposed data, not an independent test.

## Biological estimand

Ask whether the ordered protein microenvironment provides candidate-specific
information about a directed substrate-to-product reaction centre beyond the
training-only chemical prior, global sequence descriptor, missingness,
composition and taxonomic lineage controls. The estimand is ranking improvement
for documented protein–reaction relationships. It is not catalytic probability,
negative-activity classification, residue causality or a mechanistic transition
state model.

## Label-free reaction representation

For each of the 1,341 directed MAIN candidates, use its frozen RXNMapper output
only after map validation. Reconstruct the changed atoms and bonds from mapped
substrate and product sides. Encode rooted circular atom environments at radii
0, 1 and 2 around changed atoms separately for substrate and product, retaining
element, formal charge, isotope, aromaticity, chirality and bond order. Hash the
two sides into separate fixed 512-bin count vectors with a documented stable
hash, and append directed added, removed and bond-order-change counts. Record
collisions, centre size, mapping warning and an explicit centre-availability
mask. Do not use protein labels, CYP names, publications or model performance to
construct this representation.

Candidates without a validated nonempty centre retain the chemical prior and
zero interaction contribution; they are not assigned a negative score. The
primary analysis retains them in the full 1,341-candidate denominator and reports
centre coverage. A sensitivity analysis excludes mapping-warning candidates
without changing model selection.

## Separable bilinear residual

Let `X` be the position-preserving ordered protein representation and `Z` the
candidate reaction-centre representation. Both are normalized and centered;
`X` is centered within the current training pool and `Z` with a fixed
label-free catalogue transform. The interaction residual is

`Delta_int(p,r) = X_p W Z_r^T`.

Fit `W` to the same equal-protein, chemical-kernel-smoothed training residual
used by the existing conditional model. Use the exact separable ridge solution
from the singular-value decompositions of `X` and `Z`; do not choose a hidden
embedding rank from outer performance. Select the nonnegative ridge penalty by
training-only generalized cross-validation over the full Kronecker spectrum,
with the zero-interaction solution included. Numerical rank uses a machine-
precision tolerance recorded in every receipt, not a performance threshold.

Outer scores are

`log P_chem(r) + lambda_global Delta_global(p,r) + lambda_int Delta_int(p,r)`.

The two nonnegative contribution coefficients are selected jointly using the
same group-weighted inner catalogue-choice log loss; zero is legal for either
coefficient. Also retain chemistry-only, global-only, ordered protein-only,
lineage-only, interaction-only and global-plus-ordered comparators. No outer
target, final test record or candidate-specific manual rule selects a penalty,
coefficient, feature channel or fallback.

## Controls and evaluation

Run the sequence-projected and structure-supported cohorts separately with the
existing outer/inner partitions, publication purges, candidates, positives and
availability definitions. The primary interaction contrast is joint
global-plus-interaction versus global alone; interaction-only versus the
protein-only ordered residual tests whether candidate centres add information.
Report protein-cold seen and unseen reactions separately, with chemical-cold and
double-cold panels retained. All methods use the same candidate and positive
denominators before any centre-availability sensitivity.

Required controls are: missingness-only and unordered site composition; the
completed six-rank lineage model; protein-position permutation preserving each
protein's residue composition and missingness; and reaction-centre row
permutation within centre-size/mapping-status strata. Every permutation is
confined to its training fit, keeps outer targets unchanged, refits the ridge
penalty and reselects coefficients. These are descriptive information ablations
unless exchangeability is separately justified.

Save feature vocabularies, transforms, singular spectra, GCV curves/optima,
coefficients, inner predictions, outer score arrays, masks, ranks, aggregation
and source hashes. An independent implementation must reconstruct every feature
cell, separable-ridge prediction, selection objective, score, denominator and
reported interval before manuscript promotion.

## Acceptance boundary

Success requires a positive, denominator-matched interaction increment whose
direction is not reproduced by missingness, composition, lineage or frozen
permutations, with interpretable coverage in the intended novelty domain. A
signal confined to the 23-sequence structure cohort or seen reactions is reported
as such. Failure, a selected zero coefficient or unresolved cross-runtime
near-tie sensitivity is a valid outcome and triggers abstention rather than a
forced predictor claim. Independent biological generalization remains a later,
one-time public-evidence gate.
