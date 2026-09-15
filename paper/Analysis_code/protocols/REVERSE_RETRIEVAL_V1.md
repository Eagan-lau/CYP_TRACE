# Species-matched reverse-retrieval protocol

This protocol is frozen after the CLEAN pruning readout and before the reverse-
retrieval candidate census or model performance is calculated. It is a
development analysis on historically exposed public records, not independent
biological validation.

## Biological question

Given a reaction of interest and a set of cytochrome P450 sequences from the
same organism, can the strategy place a documented catalyst ahead of the other
candidate CYPs? This is a reaction-to-enzyme retrieval task. It is distinct from
forward protein-to-reaction annotation and cannot be inferred by transposing a
forward MRR table.

## Candidate and positive units

The primary candidate panel is the source-observed exact-family panel. A
candidate must have a valid explicit sequence, an exact PF00067 family record in
the fresh raw store or membership in the independently family-confirmed 600-
sequence core, and one unambiguous NCBI taxid across its retained provenance.
Sequences with conflicting taxids are quarantined. The panel is not called a
complete proteome because the downloaded sources were reaction-led rather than
an exhaustive proteome census. A core-only panel is retained as a sensitivity,
not substituted for the larger primary denominator.

Within every frozen outer block, held edges are grouped by exact taxid and
reaction identity. A query panel requires at least one documented held positive
and at least two candidate CYP sequences. All candidates are retained whether
or not they have a documented reaction; undocumented candidates are unknown,
not negatives. Repeated source exports of one sequence--reaction relationship
do not create additional positives. Candidate membership is fixed without
using the held reaction label.

## Applicability and experts

Seen and unseen reactions are explicit domains. For a seen reaction, MMseqs2 and
BLASTp may transfer that reaction only from labelled training proteins through
qualifying alignments. For an unseen reaction their direct-transfer score is
undefined, not zero. Homology-weighted chemical transport may score an unseen
reaction through similarities to training reactions, but only for candidate
proteins with a qualifying labelled homologue. Global protein conditioning and
candidate-specific interaction may be added only after source-observed candidate
features are rebuilt with the same frozen encoders; unavailable features have
zero routing weight rather than a numerical zero score. CLEAN remains contextual
EC evidence and does not prune candidates after failing its exact-unexposed gate.

The route is also stratified by whether a documented positive sequence carries
any retained training edge. This distinguishes a new reaction for a represented
protein from retrieval of a new protein. Scores from these strata or from
seen/unseen reactions are not mixed before evaluation.

## Metrics, controls and uncertainty

All methods use the same taxid-specific candidates, held positives and original
denominator. Report first-positive reciprocal rank, recall and hit at 1, 5 and
10, candidate coverage, documented-positive coverage and panel-size strata.
An uncovered documented positive contributes zero. The analytic uniform-ranking
expectation is computed for each candidate/positive count.

The primary uncertainty unit is taxid: all reactions, proteins and repeated
outer appearances for a taxid stay together in each bootstrap draw. Chemical-
component and protein-group clustered sensitivities are reported separately;
outer blocks are not treated as independent replicates. Required falsification
controls permute candidate protein scores within exact-taxid panels and permute
protein-conditioned residuals while retaining panel size, score distribution
and positive denominator. Any data-driven coefficient or selective cutoff is
chosen only from the corresponding inner blocks; no outer reaction or rank may
select it.

## Acceptance boundary

A route is promoted only when it improves over the analytic uniform expectation
and the appropriate non-protein comparator on identical candidates, has a
positive taxid-clustered conditional interval, is not reproduced by the protein
permutation, retains documented-positive coverage, and is supported across at
least 20 taxids. Represented-protein results cannot establish new-protein
retrieval. Source-observed panels cannot establish proteome-wide screening.
Failure, insufficient taxid support and abstention are valid outcomes.

