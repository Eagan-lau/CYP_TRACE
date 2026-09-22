# Frozen protocol: prospective public CYP evidence refresh v1

Status: frozen before downloading or inspecting the current P450Rdb export.
The July 2026 Mahood CYP benchmark, the July 2026 OpenADMET/Octant assay, the
previous P450Rdb v2 export and all resources listed in the project exposure
register are development-exposed and cannot become independent by regrouping.

## Acquisition

Acquire all five files linked by the current official P450Rdb v2 download page:
reactions, P450 annotations, compounds, sequences and reaction cascades. Record
the page URL, resolved file URLs, UTC acquisition time, HTTP response metadata,
byte counts and SHA-256 hashes. The downloaded snapshot is read-only after its
manifest is written. Raw files are not placed in a public redistribution bundle
unless an explicit redistribution licence is documented.

## Difference and exposure gate

First compare every acquired hash with every earlier P450Rdb file in the raw
manifest and exposure register. If the reaction, enzyme and sequence exports are
identical to an exposed snapshot, the eligible independent set is exactly zero.

If any required export differs, parse the old and current snapshots without
model scores. A candidate record must be absent as the same protein--reaction
association from all development and external-exposure resources. Exact
sequence, accession, publication, full-reaction and main substrate--product-pair
overlaps are recorded separately; lack of one exact match is not itself proof
of independence.

## Record qualification

Before predictions are opened, an eligible exact-reaction record must have:

- an unambiguous full amino-acid sequence and taxid;
- parseable, directional substrate and product structures with component
  multiplicity retained;
- a source-stated protein--reaction association and primary-publication
  identifier;
- no exact protein--reaction association in any development-exposed source;
- no primary publication in the development or external-exposure register;
- a documented wild-type/variant/construct state, or an explicit unresolved
  construct flag that prevents a sequence-exact biological claim.

Records failing a field remain in the acquisition denominator with a reason.
Named chemistry is not silently structure-resolved from a different database.
Database copies of one publication are one evidence lineage, not independent
replicates.

## Sealed evaluation eligibility

The refresh can supply an independent model test only if at least 20 qualified
protein groups remain and their labels were not viewed before model and routing
artifacts were hash-locked. Fewer records may be reported as a feasibility or
case audit but cannot pass the independent generalization gate.

For a species-matched reverse task, a versioned complete species proteome and a
reaction-blind CYP family rule must define the candidate denominator before
labels are opened. Source-observed CYP lists are insufficient for this gate.

The one-time primary comparison, if eligible, is the frozen domain-routed score
against its domain-matched MMseqs2/chemical control on identical candidates and
positives. Positive--unknown panels report candidate coverage and rank, not
precision or false-positive rate. Protein, chemical, publication and model-
pretraining exposure are reported independently.

## Fail-closed outcome

Identical files, zero eligible records, fewer than 20 protein groups, an
incomplete candidate proteome or an unavailable primary-publication lineage are
valid outcomes. In those cases the manuscript states that independent public-
evidence validation remains unavailable; no old benchmark is substituted.
