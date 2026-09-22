# Candidate-interaction permutation controls

This negative-control specification is frozen after the primary interaction
readout and before either control is evaluated.  The controls are descriptive
falsification analyses, not confirmatory tests and not an independent dataset.

## Protein-position permutation

For every protein independently, reorder the 35 anchored positions by the
SHA-256 lexical order of
`protein_position_permutation_v1|sequence_sha256|position_index`.  Apply the
same position order to residue codes, presence masks and the full 21-state
one-hot array.  This exactly preserves each protein's residue composition,
noncanonical-state count and missing-position count while destroying shared
positional correspondence.  Availability and every non-ordered input remain
unchanged.  The transformation uses no labels, split identities or scores.

## Reaction-centre row permutation

Stratify all 1,341 candidates by exact changed-centre atom count and mapping
warning status.  Within each stratum, sort reaction identifiers by SHA-256 of
`reaction_center_row_permutation_v1|reaction_key`, then cyclically assign each
candidate the next candidate's complete reaction-centre row.  Strata of size
one remain fixed and are reported.  The standardized reaction matrix is row
permuted exactly; singular values and right singular vectors remain unchanged,
while left singular vectors and the chemical-kernel-times-left-basis cache are
recomputed.  Candidate identities, labels, chemical prior, warning status and
outer targets do not move.

## Fitting and interpretation

Each control reruns every evaluable inner and outer fit for both sequence and
structure cohorts.  The ridge penalty is reselected by training-only GCV and
the nonnegative interaction/global coefficients are reselected on the same
inner records.  All candidate and positive denominators remain unchanged.  A
primary gain is biologically informative only if its direction is not matched
by either control under the same task and aggregation.  Because the control
details were finalized after viewing the primary interaction result, all
control comparisons remain development evidence.
