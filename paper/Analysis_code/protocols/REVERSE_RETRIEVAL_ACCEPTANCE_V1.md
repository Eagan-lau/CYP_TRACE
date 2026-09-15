# Reverse-retrieval operational acceptance addendum

This addendum is frozen after the reaction-blind candidate census and before any
candidate-to-core alignment or reverse-ranking result is read. It makes the
qualitative acceptance boundary in `REVERSE_RETRIEVAL_V1.md` executable.

The primary panel is the 1,803-sequence source-observed candidate set. Alignment
qualification reuses the forward baseline thresholds without tuning: E-value at
most 0.001 and at least 0.5 coverage on both sequences. MMseqs2 weighted
similarity to training proteins carrying the query reaction is the prespecified
seen-reaction method. MMseqs2 homology-weighted chemical transport is the
prespecified unseen-reaction method. Corresponding nearest-positive and BLASTp
scores are secondary comparisons; no best-method selection is performed.

The non-protein comparator ranks candidates by the number of retained training
reactions on that exact sequence, exposing annotation/study intensity. The
analytic uniform expectation is also required. A cell is accepted only if its
prespecified method has taxid-clustered 95% interval lower bounds above zero
against both comparators, spans at least 20 taxids, covers at least 80% of the
documented positive instances, and has an exact within-taxid protein-permutation
P value at most 0.05 using 99 frozen permutations plus the observed statistic.
All gates are conjunctive. The candidate and positive denominator is unchanged
for every comparison, and uncovered positives remain zero.

Cells are defined by seen versus unseen query reaction and represented versus
new documented-positive protein. Acceptance in a represented-protein cell does
not transfer to a new-protein cell. The source-observed panel remains an
incomplete-proteome limitation even if a cell passes.

