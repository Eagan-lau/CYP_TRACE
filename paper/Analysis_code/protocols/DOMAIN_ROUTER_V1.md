# Coverage-layer CYP routing and abstention protocol

This development protocol is frozen after the candidate-interaction and its
permutation-control readouts, but before any routed-panel result is calculated.
It therefore evaluates a development-derived strategy on historically exposed
data; it is not the later one-time independent public-evidence test.

## Biological and operational question

The router asks which information source is defined for a particular query
protein and candidate reaction.  It does not create a universal cross-domain
score.  The intended output is either an evidence record, a ranked list within
one explicitly defined applicability cell, or abstention.  Unknown
protein--reaction pairs remain unknown rather than biochemical negatives.

## Observable applicability variables

All variables are reconstructed separately from each frozen outer training
pool, without using its outer targets.

* `exact_protein_represented`: the exact sequence SHA256 has at least one
  retained training edge.
* `candidate_chemistry_seen`: the candidate reaction identity has at least one
  retained training edge.  Its complement is an unseen candidate; a numerical
  score of zero is never used to encode this state.
* `labelled_homologue_available`: the query has at least one MMseqs2 match to a
  retained labelled training protein at E-value <= 0.001 and bilateral coverage
  >= 0.5, using the already frozen all-versus-all matrix.
* `ordered_site_available`: the sequence-projected 35-position representation
  passes its existing availability rules.
* `documented_training_edge`: the exact protein--reaction pair occurs in the
  retained training evidence.

These are availability facts, not learned gates.  Candidate reactions are
partitioned into seen and unseen panels before ranking.  Scores from different
panels are never compared or added.

## Frozen routing table

1. A documented training edge returns the underlying evidence and is not
   presented as a new prediction.
2. An exact represented protein with an unseen candidate and an available
   ordered-site projection uses `joint_global_interaction`.  It is ranked only
   against the unseen-candidate panel.  This is the candidate-interaction
   route supported by the sequence chemical-cold development result.
3. A new protein with a seen candidate and at least one labelled homologue uses
   MMseqs2 label transfer.  `mmseqs_weighted` is the operational expert because
   it assigns evidence only to labels carried by qualifying homologues and
   aggregates all such evidence.  `mmseqs_top1` remains a reported sparse,
   precision-oriented comparator; the choice is not reselected from outer
   targets.
4. A new protein with an unseen candidate and a labelled homologue evaluates
   `homology_chemical_transport` as a provisional transport expert, only within
   the unseen-candidate panel.  It is promoted to an accepted route only if its
   paired development interval is positive against both the chemistry prior
   and the analytic uniform expectation on the same candidates and positives.
   Otherwise the router abstains.
5. A protein-specific request without a labelled homologue abstains.  A
   chemistry-only prior may be shown separately as non-protein-conditioned
   context, but it is not a CYP functional prediction.
6. A previously undocumented pair in the represented-protein/seen-chemistry
   cell is not evaluated by the current novelty splits.  The development router
   abstains rather than inferring support from absence in a positive-only
   catalogue.

Structure-subset interaction results are excluded from routing because their
apparent improvements were reproduced or exceeded by frozen permutation
controls.

## Denominators and estimands

The routed panels are:

* `represented_protein_unseen_chemistry`: chemical-cold target panels for which
  the exact query protein occurs in the corresponding outer training pool and
  the ordered sequence representation is available; candidates are all and
  only outer-unseen reactions.
* `new_protein_seen_chemistry_with_homologue`: protein-cold seen-target panels
  with a qualifying labelled homologue; candidates are all and only outer-seen
  reactions.  Sparse MMseqs2 domain masks are intersected with that candidate
  panel, and uncovered documented positives contribute reciprocal rank zero.
* `new_protein_unseen_chemistry_with_homologue`: the unrepresented-query subset
  of chemical-cold, protein-cold unseen-target and double-cold panels with a
  qualifying labelled homologue; candidates are all and only outer-unseen
  reactions.  The protein-cold and double-cold panels remain the stricter
  protein-novelty evidence; the chemical-cold subset is retained for accounting
  completeness and is not substituted for them.
* corresponding no-homologue and unavailable-site cells are counted as
  abstentions and receive no fabricated rank.

Within every reported comparison, methods share the identical query panels,
candidate set, positives and tie rule.  Report query-panel count, unique
proteins, protein groups, positive instances, candidate counts, candidate
coverage and positive coverage.  Collapse repeated chemical panels within
query before giving equal weight to protein groups.  Conditional 95% intervals
use the frozen 5,000-replicate protein-group bootstrap and are descriptive:
they exclude split, chemical-component, publication, fitting and multiplicity
uncertainty.

## Selective-risk description

For accepted routes, confidence is a label-free ordering statistic, not a
probability: top-1 minus top-2 routed-score margin for the interaction route and
maximum qualifying MMseqs2 bit score for the homology route.  At cumulative
coverage cut points 10%, 25%, 50%, 75% and 100%, report documented-positive
retrieval failure (`1 - hit@10`) and MRR.  No threshold is selected from the
outer readout, and the curves must not be called precision, catalytic
probability or calibration.  A deployable cutoff requires the later frozen
external evidence set.

## Acceptance and claim boundary

The router succeeds if it eliminates cross-domain score mixing, preserves all
denominators, and exposes where a defined expert retrieves documented evidence
better than its matched alternatives.  It does not establish biochemical
specificity, prospective precision or independent generalization.  Any cell
without supported protein-conditioned signal is an explicit abstention.  An
independent implementation must reconstruct every availability variable,
candidate mask, rank, coverage quantity, aggregate and interval before a routed
claim enters the manuscript.
