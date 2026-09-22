# Computational reproduction scope

Detailed commands, data dictionary and access boundaries are in
[`paper/REPRODUCTION_GUIDE.md`](paper/REPRODUCTION_GUIDE.md).

| Paper component | Publicly executable operation | Remaining upstream requirement |
|---|---|---|
| Figure 1 identity and core | Recalculate source unions and full/MAIN overlap from 11,593 audit records; rebuild the complete core with `upstream/rebuild_core.py` (tested) | Original 27 source files for complete-core reconstruction; see hashes and provider routes |
| Figure 2 | Full saved-score coverage/ranking evaluation and exact tie averaging | Original searches, features and score reconstruction |
| Figure 3 | Pruning-stratum estimates/intervals; paired interaction and structure point estimates from query records | Full protein models, CLEAN assets, structure/features and refitting |
| Figure 4 | Taxid-macro MRR and coverage from 20,000 panel-metric rows | Census/search/transport score generation |
| Figures 5–6 | Raw-source label reconstruction, saved-score AP/release checks, all-label kNN inference and model-asset reconstruction | Original nested kNN selection is supplied but not rerun in acceptance |
| Exploratory logistic comparison (S14b and paired comparison) | Regenerate fingerprints, repeat 180 development-fold fits to select six models, refit them and recompute all 3,035 external scores and the 5,000-replicate paired interval | Starts from supplied normalized labels and original folds; does not repeat source normalization or establish a new independent test |
| General biological evidence lookup | Rebuild 375 sequences, 639 reactions, 1,232 UniProt/Rhea edges and test every exact match | Admission decisions are frozen; not a new test or complete mixed-source atlas |
| Human source data | Original CYPstrate supplement, 12 used Figshare v4 label files, all normalized development/external labels | Original selection workflow is distinct from frozen inference |
| Supplementary analyses | All 30 table files, figure objects, external identity audits and scientific source scripts in `paper/` | See per-figure dependencies in the guide; not every model refit was rerun |
| Final figure artwork | Submission PDFs retained with the manuscript | Earlier plotting code does not reproduce final author layout edits pixel for pixel |
| Persistent archive | DOI 10.5281/zenodo.22899238; tag submission-20260922 | Includes portable source-reconstruction entry points, manifests and ordered commands |

The P450Rdb v2.0 study files were downloaded on 29 July 2026 and are identified
by their sizes and hashes in `upstream/core_input_manifest.json`. Raw files
are acquired from the provider rather than redistributed in this package.
CLEAN weights and BRENDA/SABIO-RK bulk exports are also omitted. The open
biological subset and human datasets are supplied with their source identifiers
and transformations.

Run `python run_acceptance.py --output acceptance_run` after installing
`requirements-reproduction.txt` and `./inference`. The detailed access note is
[`paper/RESTRICTED_DATA_AND_REVIEWER_ACCESS.md`](paper/RESTRICTED_DATA_AND_REVIEWER_ACCESS.md).
The complete-core entry point and human-source reconstruction are documented
in [`upstream/README.md`](upstream/README.md). Exact reconstruction uses the
inputs identified by the manifest; provider download contents may change.
