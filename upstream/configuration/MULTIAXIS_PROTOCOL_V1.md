# MAIN-cohort multi-axis development evaluation, version 1

This definition is recorded before any MAIN-cohort model scores are computed.
It changes neither the original source labels nor the frozen protein groups.
The whole project still includes structural/local modelling, human substrate
prediction, reverse retrieval, independent evidence and manuscript writing.

## Chemical definition and publication identity

Use the complete retained MAIN substrate and product sets. For every carbon-
containing connected molecular fragment, use an achiral RDKit Bemis–Murcko ring
scaffold; for an acyclic fragment use its achiral canonical connectivity rather
than one empty key. Preserve charge in this key. This is scaffold-key equality,
not a claim that all structurally similar molecules are excluded across folds.
Full MAIN labels still retain stereochemistry, multiplicity and charge.

A reaction belongs to the connected component formed by sharing **any** such
key on either side. Keep bridges, carbon coproducts and large components. Report
their contribution to component size; never remove a common fragment just to
obtain a convenient split. Charge/tautomer-invariant grouping can be a separately
declared sensitivity analysis, not a silent replacement of this definition.

Use only raw P450Rdb records with exactly one DOI and one PMID as explicit alias
links. Preserve the raw IDs, link locators and ambiguity. If an alias component
contains multiple distinct PMIDs, do not assert one paper: conservatively purge
the entire connected provenance unit and report the conflict. This improves
known alias handling but does not certify complete literature independence.

## Outer tasks

- Protein-cold: retain the five existing protein-fold assignments. Test every
  edge of the query proteins; train on other protein folds after publication
  purge. Report known and unseen reaction identities relative to that training set.
- Chemical-cold: hold out each chemical fold in turn, allowing known proteins;
  test its edges, train on other chemical folds after publication purge.
- Double-cold: use **all** protein-fold × chemical-fold combinations. Test edges
  in that rectangle; train outside both held axes. Purge publications associated
  with any edge of the outer query proteins, not just the selected rectangle.
  An edge is tested once across the complete set of rectangles. A protein can
  occur in several rectangles; those evaluations are not independent replicates.

All tasks retain the same supplied retrospective MAIN catalogue as the primary
candidate universe. Known positives outside a target subset are not negatives.
For chemical and double-cold tasks, a secondary novel-candidate-only analysis may
use the held chemical fold's catalogue, explicitly reporting its different
candidate denominator. Never compare these two denominators as model improvement.

Empty rectangles and publication-purged empty training sets remain visible and
are not replaced by easier folds. Group counts and retained edges are reported
before any model fitting. The one-group protein fold is preserved; aggregate
inference is by groups, not treating folds or sequences as independent replicates.

## Inner selection

Within each outer training set, deterministically assign up to three inner
protein folds (and, where the task requires chemical novelty, up to three inner
chemical folds). Do not split an existing component. Repeat the same axis and
publication exclusions using only outer-training records. Do not let outer test
labels enter fitting, selection, calibration, candidate-frequency estimates or
regularization choice. Supplied candidate structures are allowed as queries;
their test-label prevalence is not a fitting input.

For chemical-cold inner selection hold chemical groups; for double-cold selection
use all inner protein × chemical combinations. Record non-evaluable inner blocks
and their reason. Sparse inner evidence must not trigger a search over many
strategies: the zero-protein correction remains the explicit default and a legal
selected solution. Model grids and the exact selection objective are frozen in
the model specification before scoring.

## Scientific limits

An exact reference-sequence association is not proof of the assayed construct
or CYP-domain activity in a fusion protein. Search/scaffold/paper disjointness
does not establish that historical data are independently unused. The full
candidate catalogue is retrospective, not a prospective product generator.
Unknown pairs cannot support claims of biochemical specificity or true-activity
probability. These are development assessments; independent discovery requires
the separate qualification gate in WHOLE_PROJECT_ACCEPTANCE.md.
