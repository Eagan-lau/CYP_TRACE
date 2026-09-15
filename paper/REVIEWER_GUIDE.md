# Reviewer entry point

Use [REPRODUCTION_GUIDE.md](REPRODUCTION_GUIDE.md) for current instructions.
The repository now supplies query-level general-CYP results and an openly
reusable biological subset alongside saved-score evaluation and fixed-human
inference. This replaces the previous attachment-only guide.

From the repository root, run `python paper/reviewer_checks.py --software-root .`
for integrity and saved-score checks, `python paper/reproduce_general_metrics.py`
for assertion/query-level recalculation, and
`python paper/rebuild_open_core.py --output rebuilt_open_core` for biological
reconstruction using frozen admission decisions.

Complete mixed-source raw rebuilding and model refitting have additional
input requirements, listed in the full guide. A passing checksum test does
not certify those unexecuted operations.
