# Fresh structural expert and availability controls

The raw PDB coordinate collection is a geometry resource, not a new activity
dataset. All resolved chains remain in the search reference, including repeated
structures, modified constructs, chains without local heme contacts and unresolved
residues. Deduplicate identical resolved sequences for sequence-to-chain search;
preserve the many-to-many chain/source mapping.

Run fresh BLASTp from the 600 core sequences to this chain reference. Primary
direct structure association requires 100% identity over the reported aligned
region, >=95% resolved-chain coverage and >=80% full-query coverage. This denotes
sequence-compatible partial experimental geometry, not certification of the full
assayed construct. Explicitly report a stricter complete-coverage subset and
identity-relaxed sensitivities at 95%; do not change the primary rule after scores.
Every aligned residue and heme-contact index must reconstruct the raw sequence.

Run Foldseek 10-941cd33 with 3Di+AA alignment and exhaustive pairwise search on
the raw coordinate collection. Preserve database sequences, structural scores and
alignment coordinates. Approximate alignment TM score is identified as such.
Reconcile Foldseek chain sequences with the independently parsed chain assets
before transferring scores to reference-sequence identities.

When more than one compatible structure exists, choose a single representative
without activity labels or retrieval scores: highest query coverage, highest
resolved-chain coverage, experimental resolution when comparable, then stable
source/chain identifier. Do not maximize a query's final ranking across all of its
deposited structures. Retain the complete association and missingness census.

Evaluate sequence, ESM and structural experts on exactly the same primary
structure-available query and training subsets as well as reporting full-cohort
coverage. Record structure method, resolution, construct mismatch, number of
available structures, and annotation/publication density. The availability-only
and study-intensity controls are required before claiming distant-homology gain.

The ordered local-site model is a separate next stage: structural anchor mapping
and positional residue identity, not the average of all heme-contact residues.
A BLAST search alone does not establish SRS boundaries; heme proximity alone does
not prove a substrate-access channel or reaction-center interaction.
