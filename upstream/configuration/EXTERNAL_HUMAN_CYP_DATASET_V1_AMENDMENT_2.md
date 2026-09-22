# Protocol clarification 2: cross-file tautomer representation

Recorded on 2026-09-13 after overlap parsing but before any external model score
was calculated. The first evaluation attempt stopped during input construction
because 22 standard InChIKeys had more than one RDKit canonical SMILES across
isoform files. The variants were tautomeric, proton-location or charge-form
representations that map to the same standard InChIKey.

To prevent isoform-file-specific SMILES notation from becoming an artificial
conditional signal, every external query is deterministically represented by
the RDKit molecule obtained after a standard-InChI round trip, followed by
canonical isomeric SMILES. The same round trip is used on both CYPstrate and
external molecules only when defining scaffold-overlap strata. It does not
alter the already frozen CYPstrate training fingerprints or labels. Standard
InChIKey remains the exact-identity key. Source SMILES and their first-pass
canonical forms remain in the row audit.

This resolution uses no label, metric or model score. The failed attempt and
its pre-score exception are retained in the execution record.
