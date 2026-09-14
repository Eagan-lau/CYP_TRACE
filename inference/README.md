# CYP-TRACE frozen inference package 0.1.0

From the **repository root**, install with `python -m pip install ./inference`.
Then run `python -m cyptrace_pipeline doctor` or
`python -m cyptrace_pipeline human-substrate --smiles CCO --isoform CYP3A4 --output result.json`.

The bundled model uses similarity-weighted 25-neighbour within-isoform ranking
and a 51-neighbour pooled control. It stores all 14,955 normalized development
labels over 1,768 compounds. Scores are not calibrated catalytic probabilities.
Six named human isoforms have external-evaluation status; three are development-
only. This status applies to the task, not proof about an individual query.

`exact_training_compound` checks exact InChIKey and `training_scaffold_seen`
checks the original scaffold definition. Neither performs a parent-identity
novelty audit. The historical model schema and prediction values are unchanged.

General `reaction-screen` retrieves exact documented evidence from a separately
authorized local index and otherwise abstains. That mixed-source index is not
bundled. The synthetic evidence test in the repository tests routing, not biology.

See the repository-root README, `test_inference.py`, `rebuild_human_bundle.py`,
`REPRODUCIBILITY_MATRIX.md` and `ATTRIBUTION.md` for complete execution and source
information. Code licensing and source-data terms are distinct.
