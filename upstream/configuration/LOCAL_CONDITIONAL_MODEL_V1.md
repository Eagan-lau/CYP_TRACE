# Availability-matched local and global conditional chemistry models

Development protocol recorded before local functional scores. This stage uses
the frozen MAIN catalogue and partitions, not a new independent test set.
It does not replace the complete reaction-centre model or final abstaining router.

Run two separate cohorts: sequence_projected availability (289/600) and
structure_projected availability (23/600) from ordered_site_features_01. In each
cohort, all compared methods have the same eligible query AND training proteins.
Keep original outer/inner edge partitions and publication purges; only then
restrict to available proteins. An empty training or query pool is unevaluable,
not a successful abstention. Keep all 1,341 candidates. Seen/unseen refers to the
original full outer training pool, as in the corrected structural comparison;
also report labels seen after availability restriction. Do not compare the two
cohorts' raw MRRs as a representation effect.

Four protein descriptions have independently fitted linear ridge residuals:

1. Fresh full-sequence ESM2 mean vectors (global).
2. The fixed 35 ordered positions, with 20 standard residue channels, unknown X,
   and a missing-position channel at each position (35 x 22, position-major).
3. Unknown/missing indicators only (35 x 2), with no residue identities.
4. Counts of the 21 amino-acid channels across the 35 positions (composition),
   discarding order. This is a simple information control, not a matched-capacity
   neural model or a complete adjustment for missingness pattern.

Rows are L2-normalized, then centered on each training pool. Apply the existing
train-only chemical backbone, GCV ridge alpha and inner-selected nonnegative
scalar lambda to each description. The same target vectors and prior must agree
across representations. Alpha's zero-residual limit and lambda=0 remain legal.
Neither alpha nor lambda uses outer targets or outer MRR.

The joint model is log P_chem + lambda_global * Delta_global + lambda_local *
Delta_local. Both coefficients minimize the SAME group-weighted inner catalogue
choice log loss. This is a two-feature conditional model, not a fusion of raw
tool scores. Use an analytic gradient and nonnegative L-BFGS-B optimization,
starting from the best zero/single-branch inner solution. The numerical safety
cap is 1,024 (as in the scalar search). An upper-bound solution, failed optimizer
or failed projected-gradient check falls back to the best zero/single-branch
inner solution and is explicitly flagged. Never use outer performance for this
choice. Record coefficient collinearity and do not interpret non-identifiable
coefficients as separate biological effects. Hyperparameter uncertainty and
effective-sample limitations remain open.

Comparators: matched chemical prior, MMseqs2 nearest/weighted direct labels and
chemical transport, ESM nearest labels and chemistry transport, ordered-site
nearest chemistry transport, all four single residual models and the joint
global/local model. Native MMseqs2 applicability masks remain explicit and are
not replaced by zero scores. Other conditional scores require a qualifying
local representation and nonempty training pool in this analysis. The final
coverage router may fall back to global methods outside the local domain, but
that is a separate inference and selective-risk evaluation.

Save every model, inner prediction, selection diagnostic, score/mask, eligible
pool and source hash. Independently verify availability/targets/aggregation,
reconstruct fits by linear solves, check scalar/joint optimality, and report
paired common-domain comparisons with conditional intervals. A fitted local
model is not evidence that mapped sites determine CYP substrate specificity.

Optimizer reference: [SciPy L-BFGS-B documentation](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html).
