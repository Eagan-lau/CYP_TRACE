# Repeated protein-information controls

This development protocol is recorded after global model results and before
control results. It is not a confirmatory preregistration. Run exactly 99
replicates (0..98), seed sequence [20260909, replicate, outer block, stage]. Do
not stop or choose methods according to the observed control distribution.

Control A permutes complete retained-training reaction-label profiles among the
retained-training proteins, separately for each inner training set and each
outer training set. It preserves the multiset of per-protein profiles, every
training candidate count and the equal-protein chemical backbone. It disrupts
protein-to-profile correspondence, not the chemical candidate catalogue. Keep
test labels, public-evidence purges and protein/chemical axes unchanged. Refit
ridge alpha and the nonnegative inner-selected lambda, including lambda=0.
Also evaluate ESM nearest-protein chemical transport from the permuted profiles.
The permutations do not cross a train/test boundary or import held-out labels.

Control B permutes query score rows within each outer block of the already
fitted real conditional and ESM-transport models. Keep each destination query's
true targets, denominator and protein-group identity fixed. This is a fixed-model
query-identity ablation, not retraining. Report singleton blocks and fixed points;
do not pretend a singleton permutation destroys information. The two controls
answer different questions and are never pooled.

Save all per-query rank metrics, group-macro summaries, alpha/lambda receipts,
training-profile assignments and deterministic input hashes. For replicate 0,
also retain inner prediction arrays and outer score arrays for independent
reconstruction. Check chemical-backbone invariance for every fit and profile
multiset invariance for every permutation. No held-out performance selects
hyperparameters, seed or stopping time.

These are diagnostic randomization controls, NOT exact permutation significance
tests: arbitrary within-pool exchangeability is not established for related CYPs
and overlapping chemical/publication groups. Report the empirical control range
and observed-control effect size, not a causal protein contribution or calibrated
p-value. Lineage-preserving, chemical/publication dependence and model-refitting
uncertainty need separate treatment. Passing a control is not independent
biological validation, catalytic probability calibration or submission readiness.
