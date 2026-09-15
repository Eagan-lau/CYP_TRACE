# Paired development comparison and common-domain reranking

This analysis is specified after the first MAIN/BLAST/ESM point estimates were
viewed. It is exploratory, not a confirmatory test or a replacement for unused
public evidence. Model scores, fitted parameters and partitions are unchanged.

For each named comparison use the identical block, query, target-positive set
and full catalogue. Report (1) existing end-to-end reciprocal rank including
uncovered positives, (2) each method's candidate and positive coverage, and
(3) reranking restricted to the intersection of their explicit candidate masks.
For (3), use the same covered positives and candidate intersection for both
methods. If there is no positive in that intersection, reranking risk/performance
is undefined, not zero. Report how many original queries, panels, groups and
positives remain eligible rather than concealing the conditional denominator.

Paired differences are averaged over chemical panels within a protein, proteins
within a protein group, and protein groups equally. Use 5,000 paired bootstrap
resamples of protein groups with seed 20260909 for descriptive 95% intervals.
These intervals condition on the fixed catalogue, fitted models, observed
publication evidence and particular partitions. They do not quantify chemical
or publication dependence, variation from model refitting, or independent
discovery uncertainty. No multiplicity-adjusted hypothesis claim is made.
For fewer than two contributing groups, the interval is undefined.

Chemical/two-axis dependence-aware inference on edge-decomposable secondary
metrics, repeated training/permutation controls and the contribution of large
components remain separate required analyses. A positive interval from this
conditional resampling alone is not called stable superiority.

Comparisons: conditional ESM residual versus its own chemical prior; ESM nearest
chemical transport versus the same chemical prior and MMseqs2 chemical transport;
ESM nearest-label transfer versus MMseqs2 nearest-label and weighted transfer;
the learned residual versus MMseqs2 weighted transfer; ESM chemical transport
versus analytic uniform expectation. The point-estimate ablation lambda=1 is
reported elsewhere and is not selected because its outer MRR happens to be larger.
