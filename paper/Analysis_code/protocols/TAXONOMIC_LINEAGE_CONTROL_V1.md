# Taxonomic-lineage control for ordered-site conditional models

This protocol is frozen after inspection of the availability-matched local-model
point estimates and before calculation of any taxonomic-lineage conditional
score. It is an exploratory post hoc sensitivity analysis on the development
data, not a confirmatory test or independent discovery evidence.

## Fixed taxonomy representation

Use the frozen NCBI `new_taxdump.tar.gz` in the project raw-data registry. The
600 core sequence hashes and their taxids come only from
`dataset_02/core_edges.json`. Redirect obsolete identifiers through
`merged.dmp`; deleted or unresolved identifiers are errors, not missing zeros.
For an identical sequence linked to more than one taxid, use the lowest common
ancestor (LCA) of the resolved taxids and report every such case.

Represent each sequence by one label-free binary ancestor indicator at each of
six ranks: NCBI `domain` (the current-taxdump replacement for the historical
`superkingdom` rank), phylum, class, order, family and genus. If the LCA path
does not contain a rank, use one explicit rank-specific missing category. Do not
include species, accessions, CYP names, reaction labels, publications or source
names. Construct the vocabulary over all 600 sequences before fitting; this is
permitted transductive use of outcome-free taxonomy and must be disclosed.

## Matched conditional comparisons

Run the sequence-projected (289-sequence) and structure-projected (23-sequence)
cohorts separately. Retain the frozen outer/inner edge partitions, publication
purges, candidate catalogue, availability restriction, protein-group weighting,
chemical prior, GCV ridge fitting and inner-only nonnegative coefficient
selection used by `LOCAL_CONDITIONAL_MODEL_V1.md`.

Fit four single protein descriptions on identical pools: global ESM, ordered
sites, missingness-only, and the fixed taxonomy representation. Also fit

`log P_chem + lambda_lineage * Delta_lineage + lambda_ordered * Delta_ordered`,

where both nonnegative coefficients are selected jointly from the same inner
catalogue-choice log loss and either coefficient may equal zero. Outer labels
must never select a representation or coefficient. Assert that the independently
rerun chemical, global, ordered and missingness scores equal their existing
availability-matched counterparts before interpreting the lineage extension.

The primary diagnostic contrasts are ordered versus lineage, ordered versus
missingness, joint lineage-plus-ordered versus lineage, and joint versus ordered.
Report protein-group macro MRR, the same candidate/positive denominators, and
conditional protein-group bootstrap intervals. These intervals remain
conditional on the fixed data, partitions, features and refits and are not
multiplicity-adjusted confidence intervals for a universal biological effect.

## Interpretation limits

A positive ordered-versus-lineage contrast is evidence that the ordered
representation contains ranking information not reproduced by this six-rank
taxonomy control on the matched cohort. It does not isolate causal active-site
residues, remove fine-scale sequence ancestry, correct structure-availability
bias, establish catalytic activity, or validate a novel public discovery.
Failure to improve is a valid result. The 23-sequence structure cohort and all
seen-reaction subgroups must retain their exact small denominators in prose.

Every taxonomy mapping, feature cell, fit, inner prediction, outer score, rank,
aggregate and input hash requires independent reconstruction before promotion to
the claim ledger or manuscript.

The first construction attempt found 600/600 missing historical
`superkingdom` values because this NCBI snapshot uses `domain`. That attempt is
retained as `taxonomy_lineages_01_before_domain_fix` and is excluded. The rank
name correction was made before any lineage-model performance was calculated.
