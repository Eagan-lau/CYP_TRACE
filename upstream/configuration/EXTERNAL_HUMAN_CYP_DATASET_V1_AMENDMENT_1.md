# Protocol clarification 1: structure components

Recorded on 2026-09-13 after acquisition and schema/source-vocabulary
inspection, but before structure overlap, fingerprints, predictions or external
performance were calculated.

The phrase "parent-fragment ... conventions used for CYPstrate" in
`EXTERNAL_HUMAN_CYP_DATASET_V1.md` is inaccurate: the locked CYPstrate builder
did not select or repair a parent fragment. To preserve predictor identity,
external SMILES are processed by the exact component-preserving
`run_raw.normalize_structure` function used by `build_human_substrate.py`.
Multi-component structures remain multi-component and receive an audit flag;
they are not silently desalted. This clarification introduces no data-driven
choice and makes the external representation match the frozen training model.

The original protocol is retained unchanged. Both its original SHA-256 and the
SHA-256 of this clarification must be recorded in downstream manifests.
