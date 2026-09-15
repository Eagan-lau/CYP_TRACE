# Data access and reconstruction scope

This note accompanies *CYP-TRACE: evidence-aware evaluation of cytochrome P450
function prediction*. It distinguishes the publicly runnable analyses from
historical inputs that are not distributed. All paths below are repository-relative.

## Publicly accessible material

The code/data companion provides original human label files, normalized
development and external labels, the frozen human model, candidate-level and
query-level evaluation records, supplementary tables, and analysis code.
`paper/RESULT_REPRODUCTION_MAP.tsv` links the six main figures to their inputs
and executable checks. `paper/Tables/Table_index.tsv` locates all 30 supplementary
table files, including S22. Downloads are accessible without cluster credentials
or a VPN. Source attribution remains separate from the MIT code licence.

The general evidence reconstruction example contains 375 sequences, 639
reactions and 1,232 UniProt–Rhea relationships. It reconstructs admitted evidence
and tests each relationship against the lookup implementation. The full
mixed-source development core in the paper contains 600 sequences, 1,341
reactions and 2,304 relationships. The open example is a source-defined subset
of that core, not an independent test set or a substitute for its full upstream
reconstruction. The full mixed-source index is not distributed.

For the human analysis, both original and normalized public labels are supplied.
`paper/reproduce_logistic.py` rebuilds fingerprints from normalized structures,
selects parameters on the five recorded development folds and then performs
the exploratory external comparison. Selection and evaluation are separate
commands; external labels do not enter the selection command. The original
model-class comparison was introduced after external inspection and retains
that exploratory status. This addition does not change any reported result.

## Inputs not distributed

| Input | Role in this study | Reason for omission and access route | What remains public |
|---|---|---|---|
| Historical P450Rdb v2 reaction and protein files | Mixed-source evidence qualification and general-CYP development core | An explicit redistribution licence or provider-approved reviewer route has not been established. Access is controlled by the P450Rdb provider, via its [download page](https://www.cellknowledge.com.cn/p450rdb_v2/download.html). The authors cannot promise permission or historical-file availability. | File sizes and SHA-256 identifiers in `paper/Provenance/p450_current_download_audit.json`; assertion-decision records, saved candidate/query metrics and reconstruction code. Raw exports and the full biological index are excluded. |
| CLEAN training/model assets | Pretrained EC context, exposure stratification and reaction pruning | Separately obtained research-use assets. Readers should follow the terms and asset instructions in the [CLEAN repository](https://github.com/tttianhao/CLEAN). No redistribution authorization is asserted here. | Analysis implementation and saved pruning comparisons. Re-evaluation from these predictions does not rerun CLEAN. |
| Historical BRENDA bulk archive | Biochemical-context and accession-linked EC audit; not an added core sequence–reaction source | This release does not redistribute the bulk archive. The project records CC BY 4.0 with provider acceptance/DSI notices; omission is not a claim that the data are categorically proprietary. Obtain data under the provider's terms at [BRENDA downloads](https://www.brenda-enzymes.org/download.php). Current downloads need not match the recorded archive. | Audit summaries, source identifiers, acquisition references and code. |
| Historical SABIO-RK exports | Annotation-resolution audit; no assertions admitted to the exact core | Redistribution terms for the historical exports were not certified. Direct acquisition is through [SABIO-RK](https://sabiork.h-its.org/), subject to the provider's access rules. | Audit categories and counts; no bulk raw export. |

The P450Rdb audit on 15 September 2026 found different provider files: 3,849
reaction records versus 3,821 historically, and 1,015 protein records versus
1,012 historically. The byte hashes also differ. The current download is
therefore not silently substituted for the frozen input. An ordinary download
link does not establish access to that historical snapshot.

No controlled-access application service has been established by the authors
for these third-party files. Queries about permission or historical copies
must be directed to the respective provider; approval criteria and response
times are not known. The authors hold no permission to promise confidential
reviewer redistribution of P450Rdb. No such raw files are present in this release.

Some remaining dependencies are technical, not licensing restrictions.
General-CYP sequence searches, feature generation and model refitting require
the complete historical inputs, external tool installations and model assets.
PDB coordinates and ESM resources have their own source routes in Table S1a;
they are not represented here as prohibited material. The current acceptance
suite does not repeat those upstream computations. Earlier plotting scripts
are included, but final author-edited figure layouts are supplied separately.

## What the acceptance result establishes

`python run_acceptance.py --output acceptance_run` checks the packaged file
hashes, recalculates saved-prediction metrics, rebuilds and tests the public
biological evidence subset, rebuilds the frozen human kNN asset, recalculates
its external scores, and repeats logistic selection, fitting and its paired
external interval. It records each command, exit status, software versions,
result files and log hashes. A passing result establishes these named
operations, not an unrestricted end-to-end reconstruction of every analysis.

The public biological example and human datasets contain actual biological
records rather than synthetic stand-ins. The UniProt–Rhea subset has not been
shown to reproduce the full mixed-source general-CYP model conclusions, and
the numerical fixtures alone do not establish that it does. The current
[Digital Discovery data policy](https://www.rsc.org/publishing/publish-with-us/publish-a-journal-article/digital-discovery)
requires reviewer access and, where data are restricted, a public representative
dataset supporting similar conclusions. Whether the supplied material meets
that requirement for the full mixed-source analyses remains an editorial
question, not a claim certified by this acceptance suite.

## Outstanding access actions

The authors still need either a provider-approved route to the exact historical
P450Rdb inputs, or an agreed alternative arrangement with the journal. This
note discloses the gap; it does not replace permission. A persistent archive
and DOI must also be established for the frozen code/data release. No DOI,
embargo, access approval or provider response has been invented.
