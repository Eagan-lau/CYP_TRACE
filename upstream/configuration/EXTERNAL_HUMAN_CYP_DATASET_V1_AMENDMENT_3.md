# Protocol clarification 3: failed InChI reconstruction

Recorded on 2026-09-13 before any external prediction was calculated. RDKit
could generate a standard InChI and InChIKey for tirapazamine but could not
reconstruct a molecule from that InChI. This stopped the second pre-score input
attempt.

For this defined failure mode, retain the component-preserving canonical SMILES
as the query representation. If more than one retained representation shares
an InChIKey, choose the lexicographically first canonical SMILES for every row
with that identity. This deterministic fallback uses no isoform, label or model
score. Its use is counted explicitly. The same key-level resolution is applied
when constructing the training-side scaffold-overlap vocabulary, but the
already locked training fingerprints remain unchanged.
