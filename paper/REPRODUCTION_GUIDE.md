# CYP-TRACE reproduction guide

## Obtain a fixed copy

Download or clone https://github.com/Eagan-lau/CYP_TRACE and use the exact commit
cited in the manuscript. Software is at the root and paper data are in
`paper/`. The unified submission ZIP preserves this same layout. Extract it
into a new directory and run the commands below from that directory.
No cluster, VPN or credentials are needed.
Archive DOI: https://doi.org/10.5281/zenodo.22899238, corresponding to tag
`submission-20260922`. It includes portable raw-source reconstruction.
See the root release metadata and scope matrix.

## Install and run

The software supports Python 3.10 or newer; the complete acceptance suite
was tested on Python 3.12.14. `requirements-reproduction.txt` pins its scientific
dependencies. Install these in a fresh environment:

```sh
python -m venv .venv
# Activate .venv using the command for your operating system.
python -m pip install -r requirements-reproduction.txt
python -m pip install --no-deps ./inference
python run_acceptance.py --output acceptance_run
```

Run from the repository root, choosing unused output paths:

```sh
python verify_package.py
python paper/reviewer_checks.py --software-root .
python paper/reproduce_general_metrics.py
python paper/rebuild_open_core.py --output rebuilt_open_core
python paper/test_open_evidence.py --core rebuilt_open_core
python test_inference.py
python rebuild_human_bundle.py --output rebuilt_human_model.json.gz
```

| Command | Recomputed inputs and outputs | Scope |
|---|---|---|
| `reproduce_general_metrics.py` | Figure 1 overlap from 11,593 assertion audit rows; 15 pruning strata and their 5,000-replicate intervals; 40 interaction and 35 structure comparisons; 14 reverse method/regime summaries from 20,000 panel-metric rows | Fixed predictions, not refitting. Interaction and structure intervals are supplied but not resampled by this command. |
| `evaluation/evaluate_candidates.py` | Figure 2 metrics from 254,790 candidate rows, 95 panels and 27 groups, including fixed and exact averaged tie handling | Saved scores, labels and applicability masks |
| `evaluation/check_external.py` | Strict external AP and release comparisons | Same 3,035 labels and frozen scores |
| `reproduce_external_intervals.py` | Recompute the original paired scaffold intervals in both S7 cohorts | Fixed scores; no refitting |
| `rebuild_open_core.py` | Exact sequences, reaction-level PubMed links, 642 Rhea reactions, participant conservation and 1,232 edge identities | Open subset; frozen admission decisions |
| `test_open_evidence.py` | Retrieval of all 1,232 biological edges and unresolved queries | Evidence fidelity, not new prediction accuracy |
| `test_inference.py` | Both kNN scores for every strict external label and invalid/unknown input paths | Frozen inference |
| `rebuild_human_bundle.py` | kNN asset from 14,955 normalized development labels | No hyperparameter reselection |
| `reproduce_logistic.py tune` then `evaluate` | Original five-fold development selection, six final fits, 3,035 external predictions and paired scaffold interval | Exploratory comparison; starts from normalized public structures/labels and frozen folds |

The full suite executes logistic tuning with a data root containing only
development labels. External records and the saved parameter table are absent
from that tuning data root. Evaluation subsequently compares the selected
parameters, fold APs, regenerated fingerprints, predictions and interval with
the historical records. For separate execution:

```sh
python paper/reproduce_logistic.py tune --output logistic_selection
python paper/reproduce_logistic.py evaluate --selection logistic_selection/logistic_selection.json --output logistic_evaluation
```

`Provenance/LOGISTIC_INPUT_PROVENANCE.json` records hashes of the original
arrays; `Provenance/human_development_folds.tsv` makes the original fold mapping
inspectable. `Provenance/original_logistic_functions.py.txt` preserves the
historical functions, and `reproduce_logistic.py` is their portable implementation.

The integrity checker verifies file hashes; the other commands perform the
listed numerical and biological operations.

## Open biological inputs

`Source_data/open_biological_core/` contains sequences, UniProt provenance,
original extracted catalytic assertions, Rhea structures, participant
projections and qualification decisions. The reconstruction script produces
`core_sequences.fasta`, `core_reactions.json` and `core_edges.json`. Build a
biological lookup index with:

```sh
python -m cyptrace_pipeline build-evidence-index --dataset-dir rebuilt_open_core --output open_evidence.json.gz
```

The 375-sequence / 639-reaction / 1,232-edge subset preserves the UniProt
component of the development core. API and Swiss-Prot are one annotation
lineage. Admission decisions remain frozen; reconstruction verifies the
admitted evidence and chemistry. Full mixed-source qualification additionally
requires cross-source ambiguity screening. See
`DATA_DICTIONARY.md` and `ATTRIBUTION_OPEN_CORE.md`.

## Original human data and downloads

`Source_data/raw_development/molecules-26-04678-s001.zip` is the original
CYPstrate supplement. `Source_data/raw_external/` contains the 12 original
Figshare v4 label files used. The 14 unused supplier fingerprint files are
linked with hashes in S1b but are not bundled.

```sh
python paper/fetch_external_sources.py --output fresh_external_sources
```

This downloads the 12 used files and checks byte count, official MD5 and
historical SHA-256. Existing files are checked rather than overwritten.
`--include-unused` additionally retrieves supplier fingerprints without
adding evidence to the paper. Normalized human inputs are in root-level
`data/`. Source normalization/model-selection code is under
`Analysis_code/analysis_code/`. The repository-root command
`python upstream/rebuild_human_sources.py --output human_source_run` performs
tested raw normalization. `upstream/prepare_workspace.py` supplies scripts
and configurations in a new directory; `upstream/README.md` covers branches.

## Remaining requirements for a full mixed-source refit

The complete 600-sequence / 2,304-edge core can be reconstructed using
`upstream/rebuild_core.py` and its 27-file source manifest. Its acceptance
receipt records successful reconstruction from the original inputs.
General model refitting additionally requires search/feature assets and
separately obtained CLEAN resources.
`RESULT_REPRODUCTION_MAP.tsv` records executable scope; the figure map links
supplied data to figures. Public fixed-prediction checks do not replace these
upstream requirements.

P450Rdb redistribution permission has not been established. The provider files
checked on 15 September 2026 differ from the frozen snapshot: 3,849 versus
3,821 reaction records and 1,015 versus 1,012 protein records. File hashes and
acquisition details are in `Provenance/p450_current_download_audit.json`.

Closing full mixed-source raw reproduction requires a provider-approved
reviewer route for those exact historical files or redistribution permission.
Exact reconstruction requires the recorded historical snapshot. BRENDA raw
archives and SABIO-RK exports are not bundled; source-specific access is recorded in S1a.
CLEAN research-use weights must be acquired from their source.

The resource-by-resource access routes, distinction between licensing and
technical omissions, and outstanding reviewer-access decision are documented
in `RESTRICTED_DATA_AND_REVIEWER_ACCESS.md`.

## Tables, figures and documentation

S1a lists 15 used resource/model entries. The 42-entry historical inventory
is technical provenance, not a supplementary results table. S22b contains
380 method/domain records, including panels without a reachable positive;
only three common-domain panels qualify for conditional positive retrieval.
Table dimensions and file names are in `Tables/Table_index.tsv`, not the
scholarly table legends.

Analysis and earlier plotting code are supplied. Author-edited final figure
PDFs are in the submission folder; the earlier plotting code does not
regenerate their final layout pixel for pixel. Manuscript/SI files are
submitted separately, not publicly posted by this code/data release.
AI assistance is documented in `Provenance/AI_ASSISTANCE.md`.
