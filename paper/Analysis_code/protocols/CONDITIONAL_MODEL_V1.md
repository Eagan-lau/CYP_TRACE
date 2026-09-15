# First fresh protein-conditioned model: prescoring specification

This model is a diagnostic component of the eventual coverage-domain strategy,
not the final calibrated/abstaining tool. All 35 outer blocks are retrospective
development evaluations. Existing MAIN and BLAST baseline results have already
been viewed; this model is not advertised as preregistered independent validation.

## Information and fitting

Use fresh ESM2-650M layer-33 full-residue mean features from `esm_features_01`.
L2-normalize a sequence vector and center it using training proteins only.
No PCA or feature scaling is fitted to the complete query collection.
Long sequences use the explicitly recorded overlapping-window representation.

For each labelled training protein, distribute unit mass over its retained MAIN
labels and transport that mass through the frozen two-sided Morgan Tanimoto
kernel. The mean transported vector is the chemical backbone. Its nonnegative
values, normalized over the supplied catalogue, define P_chem; only a 1e-12
numerical floor prevents undefined logarithms. This is a catalogue-relative
retrieval prior, not a biochemical positive probability.

Fit the centered, protein-specific transported chemical vector by linear kernel
ridge regression on ESM features. A continuous generalized-cross-validation
criterion chooses ridge regularization using only the current fitting rows.
Compare its objective explicitly against the infinite-regularization (zero
residual) limit. GCV is a training regularizer heuristic, not a homology-independent
validation claim. Group/publication-purged inner predictions determine whether
the residual survives.

The candidate score is `log P_chem(c) + lambda * Delta(p,c)`. Divide predicted
residuals by the training residual RMS so the coefficient has a recorded scale.
Estimate nonnegative lambda from the frozen inner out-of-fold predictions by
minimizing catalogue-choice log loss, with equal mass per protein group, equal
mass per protein within a group, and equal mass per chemical panel within a
protein. Each query's observed target labels share unit target mass. The other
candidates are retrieval alternatives, not experimentally demonstrated negatives.
The resulting softmax must not be described as calibrated catalytic probability.

Lambda=0 is an explicit boundary solution. Choose it when the derivative at zero
does not favor protein information, no residual is learned, fewer than two inner
protein groups are represented, or numerical optimization is not identifiable.
This minimum is a mechanical identifiability rule, not a power guarantee.
Use a convex derivative root with an adaptively expanded bracket, not a search
over dozens of ranking policies or outer-fold MRR. Boundary/numerical conditions
are recorded, never silently treated as successful selection.

## Evaluation and controls

Fit each inner model exclusively on that inner block's training edges; fit the
final outer residual only on the outer training edges. Preserve all memberships,
publication purges, 1,341 candidates, target sets and group aggregation rules.
Save chemical-only scores, the learned conditional scores, lambda=1 diagnostic
scores, and fresh ESM nearest-label/nearest-chemical transfer. The lambda=1 method
is an ablation, not a selected deployed model. All-domain diagnostic scoring is
not evidence that double-cold queries should be accepted.

Save every inner target, score, learned regularizer and selection coefficient.
Run an explicit zero-feature test, no-target-leakage tests and independent rank
arithmetic verification. Repeated training-permutation/lineage controls,
dependence-aware confidence/power, branch-wise risk calibration and the ordered
local-site interaction are subsequent required acceptance stages. They cannot be
substituted with one fitted coefficient or one seed.

## Numerical settings and sources

Normalize the training protein kernel to mean diagonal one. Search log ridge in
[-14,14] with bounded scalar minimization, also comparing the zero-model limit;
these are numerical bounds, not biological choices. Flag bound-adjacent optima.
Expand the lambda derivative bracket from 1 by doubling up to 1024; if no finite
root is bracketed, choose zero and record non-identifiability rather than release
an arbitrary large coefficient. Numerical tolerance is 1e-8 for the scalar root.

- [Golub, Heath and Wahba, generalized cross-validation (1979)](https://pages.stat.wisc.edu/~wahba/stat860/pdf1/golub.heath.wahba.pdf)
- [ESM official model and extraction documentation](https://github.com/facebookresearch/esm)
- [SciPy scalar optimization implementation](https://github.com/scipy/scipy/blob/main/scipy/optimize/_minimize.py)

Those references motivate algorithms and software, not claims of efficacy on CYP.
