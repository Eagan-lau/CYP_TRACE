# CYP-TRACE: portable saved-score reproduction

This folder is self-contained for the checks below. Python 3.10+ is sufficient;
there are no third-party dependencies. It contains **saved-score evaluation**,
not the complete raw-source reconstruction or a pretrained prediction model.
No BRENDA records, CLEAN weights, sequence strings, reaction structures or source
article text are distributed in the Figure 2 fixture.

## Run from this folder

```
python -m unittest discover -p test_evaluator.py -v
python evaluate_candidates.py example_candidates.tsv --output example_result.json
python evaluate_candidates.py figure2_candidates.tsv.gz --output figure2_result.json
python evaluate_candidates.py figure2_candidates.tsv.gz --ties average --output figure2_tie_average_result.json
python check_external.py
```

The real fixture has 254,790 rows, 95 query panels, 27 groups and 1,341
candidates per panel for each of two methods. Expected MRRs are
0.06887589343729694 (MMseqs2 top-1) and 0.08076905226028033 (ESM cosine top-1).
The shared domain has only three eligible panels from three groups, with MRR
0.75 for both. Every applicable score is 1, so equality follows from the scoring
definition and shared tie rule, not equivalent discriminatory ability. Exact
averaging over all within-score-block orders instead gives 0.6591666666666667
for both. Own-domain end-to-end MRR is 0.0525833874810483 (MMseqs2) and
0.06232376353078993 (ESM); coverage is unchanged. This is a post-evaluation
sensitivity, not a replacement of the original estimates.
`figure2_expected.json` contains the original expected results and
`figure2_tie_average_expected.json` the averaged results;
`figure2_fixture_manifest.json` ties this derivative fixture to saved score
archives. Those archives are needed only to regenerate the fixture, not run it.

`strict_external_scores.tsv` contains the six-isoform, 3,035-label evaluation
and three fixed prediction columns. The independent standard-library check
reproduces tie-aware macro AP and within-isoform top-10%/25% precision and recall.
It does not refit logistic regression or kNN. The labels are source-reported
classifications, not newly certified experimental observations.

## Use with another method

Supply a tab-separated table with one row per method/query/panel/candidate:

| Field | Meaning |
|---|---|
| method | Expert or scoring route identifier |
| query, candidate | Stable opaque identities; same candidates across methods |
| group | Aggregation unit, e.g. protein cluster |
| score | Finite number if applicable; empty/NA/null otherwise |
| applicability | 1 for defined score, 0 outside the expert domain |
| positive | 1 for a documented target, 0 for another supplied alternative |
| panel | Optional repeated panel identifier; default `default` |
| exposure | Optional supplied exposure stratum; default `unspecified` |
| tie_key | Optional deterministic key; default SHA256(candidate) |

Larger scores rank first. `positive=0` is not an experimentally established
negative. Each panel must contain a target. Duplicates, unequal candidate or
target sets, inconsistent group/exposure/tie metadata, and numeric scores
outside applicability are rejected. Route applicability must be determined
without using held-out target labels; the evaluator cannot establish that
provenance merely from the submitted table.

Default ties follow the supplied key in ascending order. The real Figure 2
fixture supplies SHA256(original reaction identifier), not a hash of the opaque
replacement ID. `--ties average` integrates uniform permutations within each
equal-score block analytically. If the highest-scoring block containing positives
has n candidates and m positives, preceded by h higher-scoring candidates,
the first-positive position j has probability C(n-j,m-1)/C(n,m). The evaluator
sums this probability divided by h+j for j=1,...,n-m+1. It does not take the
reciprocal of a mean rank or estimate a bootstrap confidence interval. Exhaustive
positive-position checks for block sizes 1–8 test the implementation.

The hierarchy assigns equal weight to groups, then queries within group, then
panels within query. Conditional MRR renormalizes the original weights over
reachable panels. Common-domain comparisons intersect masks and rebuild that
hierarchy over eligible panels, identically for both methods. No shared positive
means undefined, not zero conditional MRR. When every panel is unreachable,
conditional MRR and the decomposition error are null; end-to-end MRR is zero.
The fixture's exposure is `unspecified`, not certified absence of pretraining
exposure. User-supplied exposure strata are evaluated separately.

No score calibration, model training, uncertainty interval, clinical claim or
biological validation is performed by this evaluator. The manuscript's
task-specific bootstrap analyses are separate.

## Provenance, scope and release

External labels: Ni et al., Scientific Data 12, 1427 (2025),
DOI 10.1038/s41597-025-05753-8; Figshare record
10.6084/m9.figshare.26630515.v4 (CC BY 4.0 in the frozen acquisition metadata).
The scored table is a filtered derivative; scores and opaque fixture IDs are
computational outputs of this study. Credit the original resource as well as
the CYP-TRACE study when reusing it. Source filtering and additional post-hoc
audits do not certify all assay-level labels or original-study independence.

The full local revision includes figure-generation code, frozen plotted data,
normalization audit inputs/code/tables, manuscript sources and render receipts.
Some local raw/development materials have source-specific redistribution terms;
they are not included in this portable supplement. The remote historical raw
capsule is not silently claimed to have been synchronized with this revision.

Project-owned code is MIT (user confirmed 2026-09-15); the included LICENSE
preserves the selected Eagan-lau/CYP_TRACE repository's existing notice. The
destination is https://github.com/Eagan-lau/CYP_TRACE; upload and archival DOI
remain pending. This supplied local review package is real and executable,
but is **not yet a public archived release**. No third-party licence is overridden.
