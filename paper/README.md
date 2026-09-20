# CYP-TRACE: paper data and analysis

Companion to *CYP-TRACE: evidence-aware evaluation of cytochrome P450 function
prediction*. Start with [REPRODUCTION_GUIDE.md](REPRODUCTION_GUIDE.md).

This package includes 30 supplementary tables; assertion- and query-level
records for recalculating general-CYP results; original human substrate label
inputs; and an open UniProt–Rhea subset with 375 sequences, 639 reactions and
1,232 edges. Software, normalized human data and inference assets are at the
repository root. The biological subset belongs to the published development
data and retains its original admission decisions.

From the repository root:

```sh
python paper/reviewer_checks.py --software-root .
python paper/reproduce_general_metrics.py
python paper/rebuild_open_core.py --output rebuilt_open_core
python paper/test_open_evidence.py --core rebuilt_open_core
```

The first command uses the standard library, the second requires NumPy, and
the last two also require RDKit. Install `./inference` for the scientific
dependencies. The guide records tested environments, commands and boundaries.

`Tables/` contains machine-readable values and `Table_index.tsv`.
`Source_data/` contains original human inputs, the open biological subset,
derived general-CYP metrics, external audits and figure objects.
`Analysis_code/` preserves scientific implementations and protocols.
`Provenance/` holds technical inventories, checksums and verification records.

See [ATTRIBUTION_OPEN_CORE.md](ATTRIBUTION_OPEN_CORE.md) and
[ATTRIBUTION_HUMAN_DATA.md](ATTRIBUTION_HUMAN_DATA.md) for data terms.
Project-owned code is MIT; that licence does not relicense third-party data.
