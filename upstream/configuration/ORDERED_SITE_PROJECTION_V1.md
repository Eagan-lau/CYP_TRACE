# Frozen ordered-site projection rules

Use the 35 heme-contact positions in ordered_anchor_01 (3c6g chain A), in their
reference order. This reference was selected without functional scores. Preserve
sequence-only and measured-structure-supported projection as separate channels;
never call sequence-inferred contacts a query's measured active site.

For sequence projection use the previously computed raw core-to-chain BLAST
alignment to the exact reference-chain sequence. For structure projection use
the existing Foldseek amino-acid backtrace from a corrected primary structure
to the reference, then invert the exact-aligned sequence-to-chain position map
to obtain core-query positions. Reconstruct all non-gap strings against their
source sequences, enforce monotone one-to-one indices and retain every missing
reference position explicitly. Do not fill absent residues from a close homolog.

An available local representation requires E <=0.001, reference coverage >=0.5,
at least 80% of the 35 contact positions mapped to non-unknown query residues,
mapping of all positions in the reference EXXR/PERF motifs and the axial cysteine,
and conserved E/R at EXXR and C at the axial position. Full query coverage is
reported but is not a qualification requirement, because fusion proteins can
contain an otherwise covered CYP domain. Structure projection additionally
requires >=0.5 coverage of its query coordinate chain. Exact PERF sequence in
every query is not required; residue identities and mismatches remain visible.
These a priori geometry/mapping criteria are not evidence of substrate or
reaction specificity. Do not optimize them on retrieval results.

At each fixed position retain the query index, amino-acid identity and a missing
flag. One-hot residue channels have 20 standard identities plus unknown (X);
unmapped positions have all-zero identity channels with a separate mask. Flatten
position-major only for models requiring vectors, never average across positions.
Export all 600 sequences including unavailable rows. Model comparisons must use
explicit availability masks and missingness-only controls. A single-template
projection and motif-conservation filter have lineage/structural-selection bias;
their performance does not establish universal SRS or reaction-centre mapping.
