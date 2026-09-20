# Reviewing the computational results

The [reproduction guide](REPRODUCTION_GUIDE.md) lists installation requirements,
commands and inputs. The [result map](RESULT_REPRODUCTION_MAP.tsv) identifies
the source data and executable checks for each main figure.

From the repository root, run `python run_acceptance.py --output acceptance_run`
after installing the pinned dependencies and `./inference`. The output contains
per-command logs and a machine-readable summary. The suite recalculates
saved-prediction metrics, reconstructs and tests the open biological subset,
rebuilds the human kNN asset and repeats the exploratory logistic comparison.

Historical inputs needed for complete mixed-source upstream refitting are
listed in [the data-access statement](RESTRICTED_DATA_AND_REVIEWER_ACCESS.md).
