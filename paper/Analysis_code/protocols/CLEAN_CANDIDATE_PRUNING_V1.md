# CLEAN candidate-pruning protocol

This protocol is frozen after the zero-acceptance coverage-router readout and
before the new CLEAN feature, mapping or pruning results are calculated.  It is
a nested development analysis on historically exposed data, not an independent
or temporally matched evaluation.

## Intended role

CLEAN is used as a coarse catalytic-class gate, never as an exact
substrate--product ranking score.  It may remove candidate reactions whose
source-supported EC class is incompatible with the query's CLEAN calls.  It
does not contribute a zero or numerical score outside its domain and is never
added to MMseqs2, chemical, interaction or structure scores.

## Official model and exposure boundary

Use the official research-only CLEAN v1.0.0 `split100` 128-dimensional head over
ESM-1b mean sequence embeddings, its supplied checkpoint, 241,025 training
assignment embeddings and 5,242 EC labels.  The official maximum-separation
rule over the ten nearest EC centres supplies one to six EC calls.  Sequences
longer than 1,022 residues are unavailable and are not truncated.  Existing
ESM-1b rows from the checksum-verified pre-cutoff CLEAN run may be reused by
exact sequence SHA256; all remaining supported sequences are recomputed from
the same frozen model.  Exact split100 accession and sequence overlap are
reported.  Neither exact-overlap nor ESM-1b pretraining exposure is treated as
independent validation.

CLEAN source, training data, weights and derived model assets remain in the
licensed private workspace.  They are not part of a public redistribution
bundle; required attribution and publication-notification terms remain in
force.

## Reaction-to-EC bridge

Build the bridge anew from the fresh `run_03` raw store and the admitted
`dataset_02` records.  Parse an exact four-level numeric EC only from the same
source assertion that states the reaction: UniProt API `reaction.ecNumber`,
Swiss-Prot catalytic-activity `EC=` text or the P450Rdb `EC number` field.  The
old integrated hierarchy table is not an input.  BRENDA accession-linked ECs
are retained as protein-level context and coverage evidence but do not assign
an EC to a precise candidate reaction without a structure-resolved reaction
join.

For each candidate EC, retain its supporting assertion, source and sequence.
When evaluating a query, remove EC support contributed by that exact query
sequence before declaring candidate compatibility.  This leave-query-out rule
prevents a held protein--reaction record from supplying its own pruning label.
Partial ECs, unmapped candidates and conflicts remain explicit rather than
being silently promoted to exact mappings.

## Frozen pruning modes

The no-pruning mode is always legal.  Six CLEAN modes cross EC resolution
`L2`, `L3` or `L4` with either:

* `open`: retain every candidate without a leave-query-out exact EC mapping and
  remove only mapped candidates incompatible with all CLEAN calls;
* `strict`: retain only mapped candidates compatible with at least one CLEAN
  call.

Compatibility is equality of the first two, first three or all four EC fields,
respectively.  No distance cutoff, hand-set score weight or result-dependent
reaction rule is added.

## Nested selection and evaluation

For each outer block and routed expert separately, construct all seven masks on
its inner out-of-fold query panels.  Select the mode with the highest inner
protein-group-macro first-positive reciprocal rank; uncovered positives have
reciprocal rank zero.  Exact ties prefer the less destructive mode in this
order: no pruning, open L2, open L3, open L4, strict L2, strict L3, strict L4.
The outer target, outer score, outer rank and other outer blocks cannot select
the mode.  Zero candidate coverage or no evaluable inner panel selects no
pruning.

Apply the selected mask to the already frozen expert score within its existing
coverage cell.  Primary cells are MMseqs2-weighted transfer for new proteins
with seen chemistry, homology chemical transport for new-protein/unseen-
chemistry panels, and global-plus-interaction for the exact-represented/
unseen-chemistry panel.  Report unpruned and pruned MRR, hit/recall, candidate
coverage, positive coverage, selected-mode frequencies and exact-CLEAN-exposed
versus unexposed strata.  All paired comparisons retain the same query,
positive and pre-pruning candidate denominator.  A removed documented positive
contributes zero rather than disappearing.

## Acceptance boundary

CLEAN pruning is useful only if nestedly selected outer performance improves
over the same unpruned expert without an undisclosed denominator change and the
effect is not confined to exact CLEAN training overlap.  Failure or universal
selection of no pruning is a valid result.  Even a positive result supports
only coarse EC-compatible candidate reduction; it does not show that CLEAN
predicts a CYP substrate, product, reaction centre or exact reaction.
