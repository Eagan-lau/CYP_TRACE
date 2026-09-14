# Reproduction scope, not a blanket completeness claim

| Operation / evidence | This candidate | Remaining requirement |
|---|---|---|
| Figure 2 fixed and tie-averaged retrieval metrics | Complete saved-score fixture, groups, labels, masks, evaluator, tests, expected outputs | Upstream reference search and general biological atlas are not reconstructed |
| Figure 5 original 3,035-label AP comparison | Complete saved scores/labels for specific kNN, pooled kNN and logistic; evaluator | Model-selection and logistic fitting workflows are in the separate full local capsule |
| Figure 6 10%/25% release comparison | Complete saved scores, stable tie keys, checks | Retrospective release fractions are not prospective calibrated thresholds |
| Fixed-human kNN inference | Inference code, frozen model, all normalized development labels, full external labels/predictions, rebuild and inference tests | Original raw-source normalization and nested selection are separate steps, not re-executed here |
| Parent-identity sensitivity (2,767 labels) | Row-level flags, normalized structures and reported summaries | Full audit computation is in the separate local revision workspace; no new independent test is claimed |
| Figures 1, 3 and 4 and general raw reconstruction | Not provided by this subset | Full local capsule exists; identify a redistributable representative biological dataset and supply its reconstruction recipe, or establish permissions for mixed-source inputs |
| All six figure drawing sources and manuscript | Available in the separate local v5 workspace | File-level publication review and inclusion in repository/archive are not yet completed |
| Public reviewer access | This scoped package is intended for Eagan-lau/CYP_TRACE; project code MIT confirmed | Persistent archival release is separate; consult the actual GitHub commit history for repository version |
| Persistent archival release | Prepared metadata/checksum structure only | Author metadata, release tag and actual archive deposition/DOI |

The normalized human train/evaluation data are included, rather than just a toy
example. The general Figure 2 fixture is numerically complete for that reported
comparison but is NOT a substitute for a representative redistributable
biological training/test dataset for the other general-CYP analyses. The
synthetic evidence test is a software test only and does not close that gap.

## Authorized upstream reconstruction

The inference command `build-evidence-index --dataset-dir PATH --output FILE`
requires three upstream products: `core_sequences.fasta` (sequence-SHA256
headers), `core_reactions.json`, and `core_edges.json`. Those products must be
rebuilt from an authorized project capsule; no public URL is supplied here
because the mixed-source index has not been cleared for redistribution.

For a public release, settle every outstanding row above explicitly. A passed
checksum test cannot establish source rights, scientific reproducibility outside
the listed scope, or actual reviewer access.
