# Data access and reconstruction scope

This statement accompanies *CYP-TRACE: evidence-aware evaluation of cytochrome
P450 function prediction*. Paths are relative to the code repository.

## Public data and software

The package provides original and normalized human substrate labels, a frozen
human kNN model, candidate- and query-level evaluation records, supplementary
tables and analysis code. `paper/RESULT_REPRODUCTION_MAP.tsv` links each main
figure to its inputs and executable checks. `paper/Tables/Table_index.tsv`
locates all 30 supplementary table files, including S22. The files can be
downloaded without cluster credentials or a VPN.

The open UniProt–Rhea subset contains 375 sequences, 639 reactions and 1,232
relationships. It preserves the admitted UniProt component of the full
mixed-source development core, which contains 600 sequences, 1,341 reactions
and 2,304 relationships. Offline reconstruction and lookup tests cover the
open subset. Separately, `upstream/rebuild_core.py` was executed against all
27 original raw inputs and reproduced the complete core, including biological
edge fields, sequences and reaction identities. The input manifest, hashes and
acceptance receipt are supplied under `upstream/`. This does not replace the
requirement for readers to acquire the identical raw files.

Human-source normalization was rerun from the original public supplements;
all 14,955 development and 14,526 normalized external records matched.
Human analyses include all 3,035 strict external labels, saved scores,
reconstruction of the frozen kNN asset and recalculation of its predictions.
The logistic implementation repeats parameter selection using development
folds, then fits the selected models and evaluates the external labels.
This model-class comparison was introduced after external inspection and
remains exploratory.

## Separately acquired inputs

| Input | Analytical role | Access route and distribution status |
|---|---|---|
| Historical P450Rdb reaction and protein files | Mixed-source evidence qualification and general-CYP development core | Obtain from the [P450Rdb provider](https://www.cellknowledge.com.cn/p450rdb_v2/download.html). Redistribution and a provider-approved reviewer route for the historical snapshot have not been established. Raw files and the full mixed-source index are excluded. Recorded sizes and hashes are in `paper/Provenance/p450_current_download_audit.json`; derived assertion decisions and query metrics are supplied. |
| CLEAN training and model assets | Pretrained EC context, exposure stratification and reaction pruning | Obtain under the terms and asset instructions in the [CLEAN repository](https://github.com/tttianhao/CLEAN). The package supplies analysis code and saved pruning comparisons; running these comparisons from saved predictions does not rerun CLEAN. |
| Historical BRENDA bulk archive | Biochemical context and accession-linked EC audit | Obtain through [BRENDA downloads](https://www.brenda-enzymes.org/download.php). Project acquisition records identify CC BY 4.0 with provider acceptance/DSI notices. The bulk archive is omitted; audit summaries, source identifiers and acquisition references are supplied. BRENDA contributes no additional exact core edges. |
| Historical SABIO-RK exports | Annotation-resolution audit | Obtain through [SABIO-RK](https://sabiork.h-its.org/) under the provider's terms. Redistribution terms for the historical exports were not certified. Audit categories and counts are supplied; no SABIO-RK assertions enter the exact core. |

The P450Rdb files checked on 15 September 2026 differ from the study snapshot:
3,849 versus 3,821 reaction records and 1,015 versus 1,012 protein records,
with different hashes. Exact reconstruction therefore requires the historical
files. Permission and historical-copy enquiries must be addressed to the
provider; no approved application procedure or response time is known.

Complete general-CYP upstream refitting also requires the source inputs,
sequence-search tools, feature generation and model assets identified in
Table S1a. PDB coordinates and ESM resources are obtained through their
respective providers. These technical dependencies are distinct from the
unresolved historical P450Rdb permission.

## Reproduction and reuse

`python run_acceptance.py --output acceptance_run` checks payload hashes,
rebuilds human labels from original source files, recalculates saved-prediction
metrics, reconstructs and queries the open
evidence subset, rebuilds the frozen human asset and repeats kNN inference
and the exploratory logistic analysis. Commands, environments and results are
recorded in `ACCEPTANCE.json` and per-command logs.

The supported operations are specified in `REPRODUCIBILITY_MATRIX.md`.
They exclude a complete mixed-source upstream refit and pixel-identical
recreation of final author-edited figure layouts. Project-owned code is MIT;
third-party data retain their source terms and attribution.

A provider-approved historical P450Rdb access route or an arrangement agreed
with the journal is still required for delivery of the exact original files.
Archive DOI: https://doi.org/10.5281/zenodo.22899238, corresponding to tag
`submission-20260922`, including portable source reconstruction. CLEAN is cited
as an external published tool with its source, version and checkpoints in
`upstream/README.md`; its research-use terms are not the project's MIT licence.
