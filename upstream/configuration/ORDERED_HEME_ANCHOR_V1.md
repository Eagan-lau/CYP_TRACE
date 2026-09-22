# Ordered, geometry-defined heme neighbourhood

Use corrected structure_mapping_audit_02, not superseded 01 assignments. This
stage defines a positional reference, not a trained model or proven substrate
recognition site/channel. It never reads activity scores to select a structure.

Reparse the corrected primary representative CIFs. A chain residue must belong
to a polymer in atom_site.label_seq_id, and its author chain/residue/insertion
identifier must reconcile with Bio.PDB's first-model residues. Exclude free
amino-acid ligands from protein positions. Use selected alternative coordinates
consistently and record all source identifiers. Require the complete extracted
polymer sequence to agree with the sequence used by the primary association;
otherwise flag that candidate rather than silently shifting indices.

Identify an axial CYS sulfur within 3.0 Angstrom of an iron in HEM/HEC/HEA. When
several geometric pairs qualify, use shortest sulfur-iron distance, then stable
heme/residue identifier. This is an anchor qualification criterion, not proof
of catalysis. Collect all polymer residues with any non-H/D heavy atom <=5.0
Angstrom from that heme's heavy atoms, in sequence order. Do not force 20-40
positions, discard inconvenient residues or pool them into a mean vector.

A reference candidate additionally requires one unambiguous ordered EXXR,
PERF, axial-CYS motif triple and C-alpha availability at every contact position.
Select the lexicographically first eligible file/chain, independent of activity
labels, annotation density and retrieval performance. Keep the full geometry
census, unsuccessful motif qualifications and source-coordinate distances.
If no candidate qualifies, report no reference; do not relax after performance.

Extract all existing Foldseek alignments to this single reference with their
amino-acid backtraces. Subsequent mapping must reconstruct raw sequences and
preserve the ordered positions and missingness. Sequence-only projection to
other CYPs is a distinct inferred alignment, not measured geometry. Call this a
heme-reference-aligned positional representation, not universal SRS boundaries
or an experimentally verified substrate-access channel. Exact PERF absence may
exclude legitimate CYPs; this is a limitation of reference qualification, not a
biological absence claim. Quantify representation coverage before modeling.
