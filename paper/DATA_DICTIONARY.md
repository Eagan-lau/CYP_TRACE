# Data dictionary

JSON is UTF-8. `.jsonl.gz` files contain one object per line under gzip.
TSV headers are machine-readable; scholarly titles are in the Table index.
Hashes are SHA-256 unless named MD5. Undefined scores are not negatives.

## Open biological core

| File | Row unit and principal fields |
|---|---|
| `qualification.jsonl.gz` | One source assertion: assertion ID, source, eligibility/reasons, sequence hash, accession/taxid, Rhea IDs, full/main reaction keys, admissible publications and source checksums. Decisions are frozen. |
| `uniprot_assertions.jsonl.gz` | One extracted assertion, including original `raw_json` reaction evidence, identifiers and source locator; admitted and excluded records are retained. |
| `proteins.jsonl.gz` | One exact sequence with SHA-256, amino-acid string, length and validity flag. |
| `protein_provenance.jsonl.gz` | One source/protein association with accession, taxid, sequence hash and family evidence. Multiple rows may support one sequence. |
| `rhea_reactions.jsonl.gz` | One directional Rhea ID with reaction SMILES, direction and master ID. |
| `rhea_main_projections.json` | One full-component reaction mapped to main substrates/products, removed auxiliaries with sides/roles, problems and main-reaction hash. |
| `expected_core.json` | Expected counts and the sorted edge-identity digest. |

Rebuilt edges deduplicate source assertions by exact sequence and main
reaction while retaining publication and assertion IDs. An edge is not
necessarily an independent experiment or construct-verified observation.

## General metrics

`identity_assertion_audit.jsonl.gz` contains lineage, admission decisions,
exclusion reasons and non-invertible identity keys. It omits source reaction
text, chemical structures and protein sequences; it supports identity-overlap
recalculation rather than substituting for a raw annotation database.

Pruning rows identify the panel, group, training exposure, chosen policy,
candidate/positive counts and selected/unpruned reciprocal ranks. Selection
records preserve development choices; inner selection predictions are not
included. Interaction and structure rows identify task, paired methods,
query/group, reciprocal ranks and denominators. Group-macro estimates first
average panels within query and then queries within group.

Reverse rows identify panel/method, reaction key, taxid, candidate/positive
counts, reciprocal rank and covered counts. MRR averages within taxid and
then across taxids. Positive coverage is total covered positives divided by
total documented positives, not a taxid-macro average. One panel may occur
under several methods or positive-availability regimes.

## Human data and software

Root-level `data/development_labels.json` and `data/external_labels.json`
contain compound/isoform labels and normalized structures.
`evaluation/strict_external_scores.tsv` fixes the 3,035-label specific,
pooled and exploratory logistic comparison. Macro precision/recall give
equal weight to six isoforms, unlike pooled hit counts. Labels are source
classifications, not harmonized negatives for all CYPs.

Exact-evidence matches return source records without a prediction score.
Unknown sequence/reaction pairs return empty evidence and a null score.
Human-substrate scores are within-isoform ranking scores, not probabilities.
