# CYP-TRACE evaluation and fixed-human inference

Evaluation and frozen-inference package for https://github.com/Eagan-lau/CYP_TRACE.
The companion in [`paper/`](paper/) provides assertion/query-level reproduction
for Figures 1, 3 and 4, original human label inputs, and an open UniProt–Rhea
subset (375 sequences, 639 reactions, 1,232 edges) with tested reconstruction
and biological evidence lookup. Start with the
[reproduction guide](paper/REPRODUCTION_GUIDE.md).
The project owner confirmed MIT for project-owned code and authorized publication
of this scoped package on 2026-09-15. A persistent archival release and its author
metadata are not yet completed. See `RELEASE_METADATA.json` and the reproduction
matrix for what this package does and does not reconstruct.

The V20 reproduction supplement adds development-only logistic parameter
selection, external refitting/evaluation, original-fold and fingerprint checks,
and a single acceptance command. It preserves the V19 manuscript results and
the frozen inference model. Read the
[restricted-data/access note](paper/RESTRICTED_DATA_AND_REVIEWER_ACCESS.md)
before interpreting a passing acceptance result as full upstream reconstruction.

## Complete declared-scope acceptance

Extract or clone a fixed version into a new directory. With Python 3.12 in an
isolated environment, run from the repository root:

```sh
python -m venv .venv
# Activate .venv with the command appropriate to your operating system.
python -m pip install -r requirements-reproduction.txt
python -m pip install --no-deps ./inference
python run_acceptance.py --output acceptance_run
```

The output directory must not already exist. Tests read only package-relative
inputs and place generated files in that directory. No network connection is
used after dependency installation. `ACCEPTANCE.json` and per-command logs
record exactly what ran. `paper/Provenance/V20_ACCEPTANCE.json` is the release
author's clean-copy receipt; rerun the command to produce your own receipt.

The paper studies separation of evidence reachability from ranking, and
fixed-human-isoform conditioning. It does not claim a strongest general CYP
predictor. The score is not an assay-harmonized probability.

## Quick reproduction

From this directory, Python 3.10 or newer:

```sh
python verify_package.py
python -m unittest discover -s evaluation -p test_evaluator.py -v
python evaluation/evaluate_candidates.py evaluation/figure2_candidates.tsv.gz --output figure2_fixed.json
python evaluation/evaluate_candidates.py evaluation/figure2_candidates.tsv.gz --ties average --output figure2_average.json
python evaluation/check_external.py
```

These commands use the standard library. Figure 2 includes 254,790 rows,
95 query panels and 27 protein groups. Common-domain MRR is 0.7500 under the
original fixed tie rule and 0.6591667 under exact uniform tie-order averaging
for BOTH methods. Equality follows from constant within-domain scores, not
equivalent discriminatory ability. External AP checks use all 3,035 labels
retained after the original identity, source-lineage and scaffold exclusions.

## Install and test actual inference

Use a fresh virtual environment, then:

```sh
python -m pip install ./inference
python -m cyptrace_pipeline doctor
python -m cyptrace_pipeline human-substrate --smiles CCO --isoform CYP3A4 --output human_example.json
python test_inference.py
python rebuild_human_bundle.py --output rebuilt_human_model.json.gz
```

The rebuild reconstructs the frozen kNN bundle from all 14,955 normalized
development labels (1,768 compounds); it does not repeat model selection.
The inference check recalculates every retained external label's two kNN scores
and compares them with the frozen values, and tests invalid-input, unseen-
protein and synthetic exact-evidence/abstention paths. The synthetic evidence
fixture checks software behavior only, not biological validity. The additional
`paper/test_open_evidence.py` checks all 1,232 biological edges in the rebuilt
open subset without treating retrieval as novel prediction.

The historical model schema uses `FIXED_HUMAN_EXTERNALLY_VALIDATED` to identify
six isoforms covered by the original external evaluation. It does not certify
an individual query. `exact_training_compound` is an exact InChIKey flag;
`training_scaffold_seen` uses the original scaffold definition. Neither checks
parent-identity novelty. Three other isoforms are development-only. The model
and score values are unchanged in this revision.

Tested scientific versions are recorded in `TESTED_ENVIRONMENT.json`; broader
package dependency ranges are not a guarantee of identical chemistry behavior.

## Scope and source rights

`REPRODUCIBILITY_MATRIX.md` states precisely what is included and missing.
`ATTRIBUTION.md` identifies the human data sources and transformations.
`LICENSE_SCOPE.md` distinguishes the existing inference-code notice from
source-data terms. No BRENDA export,
CLEAN weights, or mixed-source general sequence/reaction index is included.
The paper companion now contains the scientific analysis and earlier plotting
code. Final author-edited figure PDFs and the manuscript remain submission
files, not a pixel-identical plotting-code output. Full mixed-source raw/model
refitting is still outside the executable public scope.

`SHA256SUMS.txt` and `verify_package.py` check all packaged payload hashes.
Archive metadata must be finalized by the authors before a persistent release
is created. No archival DOI or author identity has been invented.
