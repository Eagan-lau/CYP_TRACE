# Fixed-human selective-behavior analysis v1

Status: specified after the external transfer endpoint was opened. This is a
descriptive deployment analysis, not a new confirmatory external test and not a
route-selection opportunity.

## Question

For the fixed human CYP substrate endpoint, how do precision, error among
released cases, positive recall and coverage change as lower-ranked predictions
are withheld? Does the isoform-specific model retain an advantage over pooled
chemistry at matched, label-blind coverage?

## Locked inputs and scope

- Development thresholds use only scaffold-held-out outer predictions from
  `human_substrate_models_01/outer_predictions.jsonl`.
- External evaluation uses the already exposed predictions in
  `external_human_cyp_models_01/external_predictions.jsonl`.
- Results are reported separately for `lineage_eligible` and
  `scaffold_and_lineage_eligible`.
- The six shared externally evaluated isoforms are macro averaged. No external
  label may alter a threshold, ranking, route, coverage target or tie break.

## Two complementary curves

1. **Development-frozen score thresholds.** For each method, isoform and target
   development coverage (100%, 75%, 50%, 25% and 10%), sort the development
   out-of-fold scores using SHA256 of compound InChIKey to break exact ties. The
   score of the last retained row becomes the threshold. Apply that value to the
   external scores without adjustment. Because threshold ties and distribution
   shift can change achieved coverage, report actual development and external
   coverage.
2. **Label-blind matched external coverage.** For each method and isoform, rank
   external rows by score with the same deterministic tie break and retain
   exactly the requested top fraction (ceiling to an integer). This uses the
   external score distribution to control coverage but never reads labels for
   selection. It is the same-candidate comparison of ranking quality, not a
   deployable fixed numerical threshold.

For each curve report selected rows, positive rows, precision, selective risk
(`1 - precision`), positive recall and coverage per isoform and as an unweighted
macro average. Precision concerns source-reported binary labels and is not assay-
harmonized biological precision.

## Uncertainty and interpretation

For label-blind matched coverage, resample external scaffold groups 5,000 times
with replacement and recompute the macro precision difference between the
isoform-specific and pooled models using the already fixed selection masks.
Report percentile 95% intervals. The sampling unit is scaffold, not row.

No new pass/fail gate will be created from these already exposed labels. A
positive curve can support descriptive selective behavior in the fixed-human
module; it cannot validate a threshold prospectively, transfer to unseen
proteins, or rescue any general-reaction route. The released CLI continues to
return ranking scores and applicability information without a binary cutoff.

