# CYP-TRACE raw-data rebuild, 9 September 2026

## Status and scientific scope

This is a new computational analysis, not a new independent experiment. Old
models, feature matrices, processed labels, split assignments and scores are
not inputs. Historical exposure is retained: all rebuilt records default to
`DEVELOPMENT_EXPOSED_OR_UNVERIFIED`. No independent test set is certified by
renaming files, rerunning parsers, changing a seed or freezing this protocol.
No new wet-lab measurements, public upload or manuscript submission is involved.

The questions are (1) what functional resolution public CYP evidence supports,
(2) whether complementary protein and chemical information improves recovery
of supported reactions at the same candidate budget, and (3) when the workflow
must abstain. Failure and null results remain publishable observations, not
reasons to select another test partition.

## Inputs and isolation

Read only original source exports under the existing project's `data/raw`.
The locally downloaded CYPstrate ZIP can be copied unchanged into this run's
`raw_additions`; its checksum must be recorded. The original source files and
all previous analyses remain unchanged. A metadata inventory of every raw file
is distinct from a checksum manifest of files actually consumed. Metadata or
directory presence alone never establishes a successful download.

All output goes under this directory. Every parsed assertion retains a source
relative path, file SHA-256, and row/accession/comment locator. Original values
are retained separately from normalization. Exact sequence hashes define
sequence entities; accession-to-sequence disagreements are explicit, not
resolved by keeping the last sequence. Multiple accession links are not
expanded into multiple independent observations.

## Dataset views

* D1: all source assertions, including unresolved and non-model-compatible rows.
* D2: traceable protein--directed reaction associations. The initial P450Rdb
  core requires an explicit valid sequence, taxonomy, publication identifier,
  complete parseable structures on both sides, no heavy-atom formula conflict,
  no exact identity reaction, and no ambiguous multiple-accession assignment.
  The initial primary core additionally requires PF00067 evidence for the exact
  sequence in raw UniProt/Swiss-Prot and quarantines accession--sequence or
  accession--taxonomy disagreements. Identical sequences in different species
  alone are not an identity conflict. Retain the broader source-curated view
  separately; do not change this gate after model scoring.
  This is a source-linked development cohort, not proof of the exact assayed
  construct, activity, or physiological function. UniProt/Rhea direction and
  sequence-version adjudication must pass before adding their records to D2.
* D3: source-reported human isoform--compound substrate labels. No transfer of
  non-substrate labels to other proteins, organisms, products or assay endpoints.
* D4: structures and local-residue features, not catalytic truth. Missing
  structure is retained; structural ascertainment is controlled in comparisons.
* D5: application candidates. A species-enzyme retrieval claim requires an
  audited complete candidate census; otherwise the task is panel retrieval.
* D6: independent evidence reserve. Empty until historical-use and provenance
  audits establish eligibility. Downloaded does not imply unused or independent.

These are related views, not six independent replicates. Within and across tasks,
shared proteins, publications, reaction chemistry and database lineage are audited.

## Normalization before scoring

Preserve reaction direction, molecular stereochemistry, charge and component
multiplicity. Use RDKit canonical isomeric SMILES, removing atom-map numbering
only; do not invent a structure from a name or silently repair invalid SMILES.
Compare source formulas on heavy-atom counts, not protonation or charge.
Unresolved formula syntax is recorded as unverified rather than a match.
No EC-level record is upgraded to an exact protein-reaction positive. An
unrecorded pair is UNKNOWN, never an assayed negative. Publicly curated and
experimentally supported assertions remain separate evidence types.

## Evaluation design

The first computational milestone is raw parsing and split feasibility, not
an immediate performance claim. Freeze data membership and candidate rules
before model scoring; record subsequent amendments and their reasons.

Protein novelty uses exact-sequence/accession grouping plus search-detected
similarity components at 40% identity and 80% coverage of both sequences.
Prespecified 30% and 50% sensitivity analyses assess definition dependence,
not opportunities to select a favorable result. Report actual maximum
cross-partition identity and coverage; a clustering label alone is insufficient.
Publication-overlapping training assertions are purged for each outer/inner
assessment. Missing publication identifiers cannot establish independence.

Within the development pool, use group-aware outer assessment with inner model
selection. Choose the number of folds from eligible independent groups before
scoring; five is the maximum, not a condition to weaken grouping. Probability
calibration and abstention thresholds must not be fitted on an evaluated fold.

Evaluate new protein/seen pair, new protein/unseen pair, strict two-sided chemical
scaffold novelty, and same-species reaction-to-enzyme recovery separately. Seen
means present in the relevant training fold. Partition chemical graphs using
both substrate and product; report quarantined bridges and coverage losses.

Every method uses identical query/candidate/positive sets. Distinguish scoring
coverage, within-domain ranking, and full-budget retrieval. Treat out-of-domain
scores as unavailable. Report MRR, Recall@K and covered-query/group counts;
uncertainty resamples the declared independent grouping unit. Report both pooled
and lineage-balanced estimands. Multiplicity and development exposure remain
visible. A nominal interval crossing zero is neither confirmed benefit nor
equivalence. Formal biochemical risk cannot be inferred from unrecorded pairs.

## Model sequence

1. Recompute homology (MMseqs2 and BLAST where installed), frequency and chemical
   similarity controls from training data only.
2. Fit a small chemical-prior plus protein-conditional correction. Zero protein
   contribution is a valid solution. Keep homology-supported and exploration
   domains explicit; test same-budget allocation rather than summing raw scales.
3. Rebuild ordered local residue features from raw structures and evaluate them
   against global, missingness-only, lineage and repeated-permutation controls.
4. Consider additional experts only through the declared inner selection
   procedure. CLEAN pruning must report candidate loss and positive recall.
5. Reproduce the selected procedure in a new environment before independent
   verification or general new-sequence tool claims.

## Whole-paper work plan

P01 raw inventory/checksums and fresh source parsing.
P02 protein/reaction/compound identity reconciliation and source-level audit.
P03 dataset membership, grouping, leakage/power feasibility and frozen splits.
P04 baseline training and common-denominator evaluation.
P05 conditional modelling, local-structure features and controls.
P06 independent-evidence eligibility, application utility and abstention.
P07 figures, claim ledger and manuscript computed from this run only.

No stage is marked complete by creating its directory or passing a unit test.
Runtime receipts distinguish implemented, executed, scientifically admissible
and independently validated. Existing model scores are never used to populate
missing results.
