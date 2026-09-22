# Frozen protocol: external six-isoform human CYP substrate audit v1

Status: frozen on 2026-09-13 before downloading or opening the Figshare data
files. The source publication and repository metadata were read only to define
the acquisition and eligibility rules below. No row labels, structures, source
fields or model performance from this resource had been inspected.

## Scientific question and claim boundary

Can the already locked CYPstrate-trained chemical-neighbour models transfer to
newly curated substrate/non-substrate labels for the six shared human isoforms
CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP2E1 and CYP3A4? This is a fixed-protein,
new-chemistry endpoint. It cannot test new proteins, exact products, reaction
centres, catalytic rates or pan-CYP negatives.

The candidate external resource is Figshare record 26630515, version 4
(2025-03-28), associated with Scientific Data 12, 1427 (2025). Its labels were
compiled from DrugBank, SuperCYP, Cytochrome P450 Knowledgebase, interaction
tables and literature. They are source-reported classifications, not assay-
harmonized measurements. Publication after model development does not by itself
establish record-level independence.

## Acquisition and immutability

- Acquire every file listed by the official Figshare v4 API and record file ID,
  name, download URL, byte count, supplied MD5, UTC acquisition time and local
  SHA-256 in a manifest.
- Preserve the archive/files byte-for-byte. Parsing writes only derived files.
- The Figshare data are CC BY 4.0; attribution and the source DOI are retained.
- A failed or partial acquisition is an auditable outcome and supplies no test.

## Locked predictor

The predictor is fixed before the external labels are opened. Molecules use the
same RDKit radius-2, 2,048-bit Morgan fingerprints and Tanimoto similarity as
`HUMAN_SUBSTRATE_ENDPOINT_V1.md`. For each shared isoform, the score is the
similarity-weighted mean label of the 25 closest CYPstrate compounds carrying a
label for that isoform. `k=25` was selected independently in all five nested
outer-development folds. All 14,955 normalized CYPstrate rows are used for the
production fit. Zero total similarity falls back to that isoform's training
prevalence. No external label, source field, train/test assignment or metric may
change representation, neighbour count, fallback or routing.

The previously evaluated pooled-chemical model is retained only as a locked
control. Its production neighbour count is 51, the modal nested choice (three
of five folds); ties are resolved without external performance. The primary
comparison is isoform-specific kNN minus pooled-chemical kNN on identical
eligible rows.

## Identity, overlap and evidence-lineage audit

Structures from the external files are normalized with the same parent-
fragment and RDKit conventions used for CYPstrate. Every row receives:

- parse status, canonical SMILES and standard InChIKey;
- exact-InChIKey overlap with any CYPstrate compound;
- Bemis--Murcko scaffold and scaffold overlap with CYPstrate;
- source text as supplied, plus a normalized source-lineage category;
- agreement or conflict where an exact isoform--compound label is present in
  CYPstrate; and
- original Figshare file, split and row locator.

The source-provided training/test split is descriptive and is never treated as
our independent split. Exact duplicate rows, stereochemical variants sharing a
standard identity and conflicting duplicates are retained in the audit but not
double-counted in evaluation.

## Prespecified evaluation strata

1. `all_parseable_unique`: one deterministic record per isoform--InChIKey,
   regardless of CYPstrate overlap. This is a transfer/compatibility analysis,
   not independent validation.
2. `exact_novel`: no exact InChIKey anywhere in CYPstrate. This is the primary
   new-compound stratum.
3. `scaffold_novel`: `exact_novel` and no Bemis--Murcko scaffold in CYPstrate.
   This is the stringent chemical-extrapolation stratum.
4. `lineage_eligible`: `exact_novel` plus a supplied source lineage that is not
   CYPstrate, Tian, Hunt or an unverifiable aggregate reference. This stratum is
   required for an independent-evidence claim.
5. `scaffold_and_lineage_eligible`: intersection of strata 3 and 4.

Rows with an unparseable structure, unknown/conflicting binary label, ambiguous
isoform, exact identity overlap, or unverifiable source lineage remain in the
fixed acquisition denominator with an exclusion reason. A database copy of the
same underlying assertion is not an independent biological replicate.

## Metrics, uncertainty and falsification

Average precision is primary; ROC AUC and Brier score are secondary. Report per
isoform values and an unweighted isoform macro average. Both classes and at
least 20 positives and 20 negatives must be present for an isoform to enter a
macro value. Five thousand paired bootstrap replicates resample external
Bemis--Murcko scaffolds. The primary effect is the isoform-specific minus pooled
macro average precision with its percentile 95% interval.

As a score-level falsification control, isoform-specific scores are permuted
among the available shared isoforms within each compound for 99 deterministic
permutations. The plus-one upper-tail Monte Carlo p value is reported. Neither
bootstrap nor permutation changes the model.

## Acceptance and fail-closed rules

An `independent_external_transfer` conclusion requires all of the following:

- a complete checksum-locked v4 acquisition;
- at least four evaluable isoforms and, in each, at least 20 positives and 20
  negatives in `lineage_eligible`;
- at least 100 unique compounds and 20 scaffold groups in that stratum;
- a positive 95% scaffold-bootstrap lower bound for macro-average-precision
  difference versus the locked pooled control; and
- a permutation p value no greater than 0.05.

A `scaffold_external_transfer` conclusion additionally requires the same class
and group minima in `scaffold_and_lineage_eligible`. If source lineage cannot be
verified, the data may support only a compatibility analysis. Insufficient
novel compounds, class collapse, source overlap, conflicting labels or no
increment over the pooled model are valid results and must not be replaced by a
weaker post hoc subset.

## Reproducibility gate

An independent program must verify acquisition and input hashes, regenerate
identities/fingerprints/scaffolds, reconstruct every prediction and exclusion,
repeat all summaries, bootstraps and permutations, and confirm that this
protocol and model implementation predate acquisition before any external
number enters the manuscript.
