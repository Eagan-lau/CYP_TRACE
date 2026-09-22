# Structure-availability-matched development analysis

Recorded after the 36/600 primary association census, before functional structural
scores. This is a development analysis, not a preregistered independent test.
Keep STRUCTURE_EXPERT_PROTOCOL_V1.md unchanged and exclude all 75 unreconciled
Foldseek entries. No relaxed mapping is used to optimize performance.

Use the 35 frozen outer blocks, unchanged publication purge and 1,341-candidate
MAIN catalogue. Both queries and retained training proteins must belong to the
same 36-member primary set. Targets remain all eligible positives of those
queries. The seen/unseen labels refer to the original full training split, so the
task definition is comparable to the global analysis; also record the number of
labels seen in the restricted training pool. No-training or no-query blocks are
explicitly unevaluable, not silently scored as successful abstention.

Foldseek primary similarity is its 3Di+AA bit score, with E <= 0.001 and bilateral
coordinate-chain coverage >= 0.5. MMseqs2 uses its existing E/coverage-qualified
bit scores on the same training/query subset; ESM uses fresh normalized full
sequence vectors and cosine nearest-neighbour selection. All tools' native
significance/coverage criteria and realized domains remain visible; equal numeric
E-values do not mean equal statistical evidence across tools. Structural TM-score
is approximate and may exceed 1 numerically; do not clamp or interpret as a
probability. No raw score mixing between experts occurs.

Input-audit correction before any functional results: very short, non-significant
alignments can have negative bit scores. Retain and count those source rows;
positive bit weight, significance and coverage are required for transfer. Missing
self results outside the primary subset are reported, not hidden by imputation.
Auxiliary TM/lDDT values may be undefined for short alignments. Count them and
preserve their raw values; they are not inputs to the bit-score expert. All
required score, significance, identity and coordinate fields must remain finite.

For Foldseek and MMseqs2 report nearest-label, bit-weighted label and chemistry
transport experts. ESM reports nearest-label and nearest-protein chemistry
transport. A label is outside a label-transfer expert's domain if it has no
positive support from eligible training neighbours. Transport reaches the whole
catalogue only when there is an eligible neighbour. Stable candidate hash breaks
ties. Single protein reference structures were selected without reaction labels.

Query-independent controls use (i) equal-protein chemical prior in the restricted
training set, (ii) that prior weighted by log(1 + primary compatible-chain count),
and (iii) that prior weighted by log(1 + distinct retained-training publication
count), falling back to unit weight for zero-publication proteins. The latter is
training-only; test publications/annotations are never features. These controls
probe simple deposition/annotation intensity, not full elimination of study bias.
Full-cohort availability and eligible query/edge counts are reported separately.

Report panel->query->protein-group macro metrics. Paired Foldseek/MMseqs2/ESM
comparisons must also intersect candidate domains and covered positives, retaining
undefined outcomes when no common positive exists. Conditional protein-group
bootstrap intervals cannot establish external generalization or distant-homology
benefit. No structural superiority claim from a point estimate alone.
