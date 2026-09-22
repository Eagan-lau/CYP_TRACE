# Exploratory paired local-representation comparisons

This comparison list is recorded after the local aggregate point estimates were
viewed. It is an exploratory development analysis, not preregistration or a
confirmatory test. Report the entire list in both availability cohorts and all
five task panels, including negative and undefined comparisons.

Ordered residual versus global residual, missingness residual, composition
residual and chemistry prior tests the added information of the positional
description. Joint global/ordered versus global, ordered, chemistry prior and
MMseqs2 weighted transfer tests complementarity and the native coverage tradeoff.
Ordered nearest chemical transport versus ESM nearest chemical transport and
MMseqs2 weighted transfer tests simple retrieval without fitting a residual.
Compare joint, ordered transport, ESM transport and global residual separately
against analytic uniform expectation. This is not a method-selection rule.

Use exactly the same queries, target positives and 1,341-candidate catalogue
within each pair. Native masks remain explicit. Compute end-to-end RR and
reranking restricted to the intersection of both masks and target positives.
Without a common positive, reranking is undefined, never zero. Uniform RR in
each eligible domain is the exact random-ranking expectation, not a seeded
random ranking. Average panels within query, queries within protein group and
groups equally. Keep original full-training seen/unseen labels and the separate
289-sequence and 23-structure availability cohorts.

The 5,000 paired protein-group bootstrap draws (seed 20260909, stable task/pair
order) condition on the existing fits, partitions, catalogue and source evidence.
They do not account for chemical/publication dependence, fitting uncertainty or
multiple exploratory comparisons. Intervals are descriptive; no p-values,
discovery-set designation or universal superiority claims follow from them.

Require independent fitting/rank QC before comparison; preserve all score files.
Independently reconstruct domain intersections, sorted-list ranks, denominators,
group summaries and bootstrap quantiles before reporting numerical contrasts.
