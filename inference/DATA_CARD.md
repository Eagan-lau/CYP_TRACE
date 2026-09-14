# Data card: CYP-TRACE strategy v1

## Fixed-human model data

The packaged human model derives from the 2021 CYPstrate article supplement,
licensed CC BY 4.0 with attribution. The source contains binary substrate/non-
substrate labels for nine human isoforms. After retaining all source rows in the
audit denominator, 549 structurally unresolved rows were excluded from modelling;
14,955 unique, non-conflicting isoform–compound labels and 1,768 compounds remain.
The original train/test labels were treated as development-exposed and were not
claimed as independent partitions.

The external evaluation uses the 26-file Figshare record 26630515 version 4
(DOI: 10.6084/m9.figshare.26630515), associated with Ni et al., *Scientific
Data* (2025), DOI: 10.1038/s41597-025-05753-8. It is not included in the model
bundle and is never read during inference.

## General evidence index

`build-evidence-index` consumes three qualified `dataset_02` products:
`core_sequences.fasta`, `core_reactions.json`, and `core_edges.json`. The current
development view contains 600 exact sequences, 1,341 MAIN reaction
representations and 2,304 sequence–reaction relationships. Exact retrieval can
return publication identifiers, accessions, taxids, source databases and source
assertion identifiers.

These records combine sources with different redistribution terms. UniProt and
Rhea are CC BY 4.0. P450Rdb and SABIO-RK did not have an explicit redistribution
licence established in the project audit. BRENDA 2026.1 was downloaded after
licence acceptance and is not bundled raw. Therefore the mixed-source exact-
evidence index is local-only and carries `public_redistribution_ready: false`.
It must be regenerated from an authorized local project copy and must not be
placed in a public release without source-specific filtering or permission.

## Known quality constraints

- Reaction labels are positive–unknown, not positive–negative.
- The MAIN representation removes declared helpers and may differ from the full
  biochemical equation; all representation rules are versioned upstream.
- Exact sequence, construct, CYP-domain attribution and publication traceability
  are audited separately. Evidence retrieval does not repair incomplete source
  assertions.
- The fixed-human labels are not assay-condition harmonized.
- Structure availability, publication density, phylogeny and database exposure
  can be confounded.

Every generated model or evidence bundle records SHA256 hashes of its direct
inputs. The public release must retain source citations and the project licence
matrix.

