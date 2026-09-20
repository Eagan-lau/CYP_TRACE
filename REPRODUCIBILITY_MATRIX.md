# Computational reproduction scope

Detailed commands, data dictionary and access boundaries are in
[`paper/REPRODUCTION_GUIDE.md`](paper/REPRODUCTION_GUIDE.md).

| Paper component | Publicly executable operation | Remaining upstream requirement |
|---|---|---|
| Figure 1 identity | Recalculate source unions and full/MAIN overlap from 11,593 assertion audit records | Historical raw mixed-source qualification |
| Figure 2 | Full saved-score coverage/ranking evaluation and exact tie averaging | Original searches, features and score reconstruction |
| Figure 3 | Pruning-stratum estimates/intervals; paired interaction and structure point estimates from query records | Full protein models, CLEAN assets, structure/features and refitting |
| Figure 4 | Taxid-macro MRR and coverage from 20,000 panel-metric rows | Census/search/transport score generation |
| Figures 5–6 | Saved-score AP/release checks, all-label kNN inference and frozen model-asset reconstruction | Original normalization and nested selection are preserved as scientific code, not a turnkey full refit |
| Exploratory logistic comparison (S14b and paired comparison) | Regenerate fingerprints, repeat 180 development-fold fits to select six models, refit them and recompute all 3,035 external scores and the 5,000-replicate paired interval | Starts from supplied normalized labels and original folds; does not repeat source normalization or establish a new independent test |
| General biological evidence lookup | Rebuild 375 sequences, 639 reactions, 1,232 UniProt/Rhea edges and test every exact match | Admission decisions are frozen; not a new test or complete mixed-source atlas |
| Human source data | Original CYPstrate supplement, 12 used Figshare v4 label files, all normalized development/external labels | Original selection workflow is distinct from frozen inference |
| Supplementary analyses | All 30 table files, figure objects, external identity audits and scientific source scripts in `paper/` | See per-figure dependencies in the guide; not every model refit was rerun |
| Final figure artwork | Submission PDFs retained with the manuscript | Earlier plotting code does not reproduce final author layout edits pixel for pixel |
| Persistent archive | GitHub commit provides review access | Archive deposition and DOI not yet completed |

P450Rdb historical files have no verified redistribution permission, and the
current official downloads differ from the frozen input. They were not uploaded.
CLEAN weights and BRENDA/SABIO-RK bulk exports are also omitted. The open
biological subset and human datasets are supplied with their source identifiers
and transformations.

Run `python run_acceptance.py --output acceptance_run` after installing
`requirements-reproduction.txt` and `./inference`. The detailed access note is
[`paper/RESTRICTED_DATA_AND_REVIEWER_ACCESS.md`](paper/RESTRICTED_DATA_AND_REVIEWER_ACCESS.md).
The historical third-party access arrangement remains unresolved.
