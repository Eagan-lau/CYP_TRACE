# Reconstruction from source files

All paths below are relative to this repository or a newly chosen output
directory. No author account, cluster directory or pre-existing feature cache
is required by the two reconstruction entry points.

## 1. Human source labels: included inputs

Install the root `requirements-reproduction.txt`, then run:

```sh
python upstream/rebuild_human_sources.py --output human_source_run
```

The original CYPstrate supplementary ZIP and 12 used Figshare v4 label CSVs
are included. This command rebuilds 14,955 development labels and 14,526
normalized external records, including the 3,035-label primary evaluation
cohort. It verifies every normalized record against the original computation,
not only the counts. `HUMAN_SOURCE_ACCEPTANCE.json` reports the result.
To retrieve the external files again, use
`python paper/fetch_external_sources.py --output external_download`.

## 2. Complete biological core: separately acquired inputs

`core_input_manifest.json` lists all 27 original consumed files, provider
URLs, byte counts and SHA-256 hashes. Arrange acquired inputs under `RAW/`
with its `p450rdb/`, `rhea/` and `uniprot/` subdirectories exactly as listed.
The P450Rdb v2.0 reaction and protein files were downloaded on **29 July 2026**
from the [official download page](https://www.cellknowledge.com.cn/p450rdb_v2/download.html).
Their download date, byte counts and hashes are recorded in the manifest.
Install `requirements-core.txt`, then run:

```sh
python upstream/rebuild_core.py --raw-root RAW --output core_run --check-inputs-only
python upstream/rebuild_core.py --raw-root RAW --output core_run
```

The input check does not create `core_run`. Reconstruction refuses mismatched
files and occupied output directories. It parses source records, applies
qualification and normalization, and constructs the complete 600-sequence,
1,341-reaction, 2,304-edge core. `CORE_ACCEPTANCE.json` checks all biological
edge fields, reaction identities and sequences. Path-derived assertion IDs
and FASTA line endings are excluded from semantic comparison.

This entry point was executed successfully on Linux with Python 3.11.5;
the machine-readable receipt is `core_acceptance.json`. Original input bytes
were supplied for that test. Public provider URLs are acquisition routes,
not a guarantee that a changing download still serves the same snapshot.
P450Rdb source files are acquired from the provider, not redistributed here.
Use the dated inputs identified by the manifest for exact reconstruction. The bundled UniProt–Rhea
subset has its own offline reconstruction in `paper/rebuild_open_core.py`.

## 3. General analyses and model fitting

Create a workspace with `python upstream/prepare_workspace.py --output WORK`.
It copies scientific scripts and protocol/configuration files to
`WORK/analysis`. Put the complete reconstruction's `run_03`, `dataset_01`
and `dataset_02` directories there before running downstream general analyses.
Do not substitute the smaller open subset for the complete core in published
general-model comparisons. Human normalization needs no complete core.

`ANALYSIS_INPUTS.tsv` records 1,287 original stage/input/hash relationships.
It is an input inventory, not a claim that every listed artifact is bundled.
`commands/` contains the original computation commands with cluster-specific
launch settings removed. Invoke each launcher with the **absolute** path to
`WORK/analysis`; install the named external binaries first. Optional
`CYPTRACE_PYTHON`, `CYPTRACE_MMSEQS` and `CYPTRACE_FOLDSEEK` select executables.
Launcher names preserve links to the original stage provenance.
Array launchers require a second positional argument giving the task index.
Their headers state the complete index range; execute every index, for example
`bash esm_features.sh /absolute/path/to/WORK/analysis 0` through index 3.

The main dependency order is:

| Branch | Required order and inputs |
|---|---|
| Sequence baselines and splits | Complete core → `p02_similarity.sh` → protein splitting with `split_baselines.py` → `build_multiaxis_splits.py` → main and BLAST baselines |
| Global conditional models | Core and splits → ESM2 weights/download receipt → ESM residue features and collection → conditional models and paired comparisons |
| CLEAN pruning | Core, splits and main baseline scores → official CLEAN/ESM1b assets → `build_clean_core_features.py` and verification → EC bridge → pruning and verification |
| Structure and ordered sites | Original PDB mmCIF files → geometry extraction → chain search and polymer-identity audit → matched retrieval → ordered heme anchor/site vectors → local/interaction models |
| Same-species retrieval | Complete raw database/core and splits → reverse candidate census → MMseqs2/BLAST alignments → alignment normalization → reverse baselines |
| Human model development | Source normalization → `evaluate_human_substrate.py` → model verification → `evaluate_external_human_cyp.py` → selective-behaviour evaluation |

Python entry points with command-line arguments display them with `--help`.
Inspect each stage's input inventory before execution; a later stage cannot
replace a missing upstream model or alignment with a saved result silently.
The general refit branches require additional PDB, taxonomy, pretrained-model
and binary assets, and were **not** rerun end-to-end during package acceptance.
They are not advertised as a tested single-command full-paper refit.

The original search versions were MMseqs2 18-8cc5c, BLAST+ 2.14.0 and Foldseek
10-941cd33. Neural feature generation used a separate Linux/PyTorch environment;
the root Windows reproduction environment is not its substitute. Original
protocols and stage manifests specify numerical settings.

## CLEAN and ESM assets

Obtain CLEAN v1.0.0 from <https://github.com/tttianhao/CLEAN>, commit
`0cf2cac2cac71cd626f36ec2ddd6f2c223ef7054`; cite the published Science paper,
<https://doi.org/10.1126/science.adf2465>. Follow the provider's checkpoint
download instructions and research-use terms. Place the source checkout at
`WORK/tools/external/CLEAN_v1_0_0` and its pretrained files at
`WORK/data/restricted/clean_v1_0_0/pretrained/`. Obtain ESM1b from the official
ESM project (<https://github.com/facebookresearch/esm>) and place its checkpoint
at `WORK/data/restricted/esm1b_v1/esm1b_t33_650M_UR50S.pt`.
`build_clean_core_features.py` contains and enforces all three checkpoint hashes.

The old precomputed ESM1b cache is optional: absent caches trigger extraction
from the supplied core FASTA. Partial caches fail explicitly. This portability
change is unit-tested; fresh inference of the entire CLEAN feature set has
not been rerun. CLEAN is an external published dependency, not part of the
project's MIT grant. Its original source and weights are not redistributed.

## Figure inputs

`paper/Source_data/figure_inputs` and `paper/Source_data/analysis_tables` hold
the plotted data. The plotting scripts now resolve these package-relative
locations. They require Matplotlib and an installed Arial font. Final
author-edited PDF layouts are supplied with the manuscript; earlier plotting
code is not claimed to recreate subsequent manual layout edits exactly.
