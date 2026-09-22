# CYP-TRACE

Code and data for **CYP-TRACE: evidence-aware evaluation of cytochrome P450
function prediction**.

CYP-TRACE evaluates candidate reachability and ranking, prioritizes substrates
within fixed human CYP isoforms, and retrieves documented sequence–reaction
evidence. The package contains analysis code, source tables, human labels,
a frozen inference model and an open UniProt–Rhea reconstruction subset.

## Install and reproduce

Use a fixed commit and a fresh Python 3.12 environment. From the repository root:

```sh
python -m venv .venv
# Activate .venv for your operating system.
python -m pip install -r requirements-reproduction.txt
python -m pip install --no-deps ./inference
python run_acceptance.py --output acceptance_run
```

Choose an output directory that does not exist. The 19 checks evaluate file
integrity, saved-prediction metrics, open-evidence reconstruction, human kNN
inference and development-only logistic selection followed by external
evaluation, plus reconstruction of human labels from the original public
supplements. Results, commands and software versions are written to
`acceptance_run/ACCEPTANCE.json`. After dependency installation, the checks
run offline with package-relative inputs.

## Use the software

```sh
python -m cyptrace_pipeline doctor
python -m cyptrace_pipeline human-substrate --smiles CCO --isoform CYP3A4 --output human_example.json
python paper/rebuild_open_core.py --output rebuilt_open_core
python -m cyptrace_pipeline build-evidence-index --dataset-dir rebuilt_open_core --output open_evidence.json.gz
```

Human-substrate scores provide within-isoform rankings with chemical neighbours
and exposure diagnostics. Six isoforms have external evaluation; three are
development-only. Scores are uncalibrated. Exact sequence–reaction lookup
returns documented evidence, with unresolved outputs for unmatched queries.

## Data and analysis

- [Reproduction guide](paper/REPRODUCTION_GUIDE.md): commands and data locations.
- [Reproduction matrix](REPRODUCIBILITY_MATRIX.md): executable scope for each result.
- [Result map](paper/RESULT_REPRODUCTION_MAP.tsv): figure-to-input mapping.
- [Source reconstruction](upstream/README.md): complete-core and human raw-data entry points, source hashes, environments and downstream commands.
- [Inference guide](inference/README.md): interface and output fields.
- [Evaluator guide](evaluation/README.md): candidate-table format and ranking conventions.

Human development data comprise 14,955 labels for 1,768 compounds. The strict
external evaluation has 3,035 labels across six isoforms. The open biological
subset contains 375 sequences, 639 reactions and 1,232 relationships. It is
the UniProt–Rhea component of the study's mixed-source development core.

The complete 600-sequence / 1,341-reaction / 2,304-edge core was reconstructed
from the original source files using the supplied portable entry point.
Those 27 input files are identified in `upstream/core_input_manifest.json`.
The P450Rdb v2.0 reaction and protein files were downloaded on **29 July 2026**.
Their filenames, sizes and SHA-256 checksums identify the study inputs;
the website-reference access date of 22 September 2026 records a later link check.
Exact-snapshot acquisition and separately obtained pretrained assets remain
requirements for a complete upstream refit. They are described in
[the data-access statement](paper/RESTRICTED_DATA_AND_REVIEWER_ACCESS.md).

## Licence and citation

Project-owned code is MIT-licensed. Data attribution and source-specific terms
are in [ATTRIBUTION.md](ATTRIBUTION.md) and [LICENSE_SCOPE.md](LICENSE_SCOPE.md).
Use the exact Git commit when citing or reproducing this package.
`RELEASE_METADATA.json` records the repository and archive identifiers.
Archive DOI: <https://doi.org/10.5281/zenodo.22899238>.
The archive corresponds to tag `submission-20260922` and includes the
portable source-reconstruction entry points and their input manifests.
