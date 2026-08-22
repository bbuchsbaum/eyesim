# Known-item/new-participant reinstatement court

This post-court exploratory analysis asks whether some items have reproducibly
greater eye-movement reinstatement across participants. It deliberately changes
the generalization target from the frozen item-disjoint retrieval court:
information from other participants viewing the **same item** is available,
but the held-out participant never contributes to their own predictor.

The five questions are:

1. Can an item's Transport propensity, estimated from other participants,
   predict a new participant's Transport score?
2. Is that cross-participant item prediction stronger than density?
3. Does a Transport-specific item effect remain after density sigma 80 is
   removed using training participants only?
4. Does response-blind study repeated-viewing propensity transfer to retrieval?
5. Do retrieval, Transport-specific, density, or study item propensities improve
   held-out prediction of the 1--4 newness rating?

## Leakage contract

For every held-out participant, nuisance adjustment and empirical-Bayes item
shrinkage use other participants only. Confidence models use a nested version:
training participants' item predictors also exclude their own gaze, and the
outer evaluation participant is absent from the complete training pool.
Saliency, probe type, effective fixation count, and retained duration are
nuisance predictors. Study item propensity uses only intact repeated-viewing
scores and no retrieval response.

The score endpoint is cross-validated variance explained (`cv_r2`) relative to
the training-only nuisance prediction. Confidence uses held-out ordinal
log-loss gain in bits. Crossed participant/item bootstrap intervals quantify
uncertainty. Private trial predictions remain local and Git-ignored.

Run after loading the package:

```r
source("inst/validation/gaze-weave-item-effects.R")
run_gaze_weave_item_effects()
```

## Current result

The held-out-participant retrieval signal is small and uncertain. Other
participants' Transport propensity explains 0.70% of retrieval-score variance
(crossed-bootstrap 95% interval -1.00% to 2.42%). All three density estimates
have negative point estimates. Transport exceeds density by about 1.0--1.2
percentage points, but every Transport-minus-density interval includes zero.
After removing density sigma 80 using training participants only, the remaining
Transport item propensity explains 0.22% (-0.82% to 1.31%).

Repeated-viewing item propensity is reproducible within the study phase: 7.10%
of score variance (1.50% to 13.15%). Once its scale is calibrated using retrieval
training participants, it transfers only weakly to retrieval Transport: 0.34%
(-1.20% to 1.89%). Descriptively, full-sample item ranks correlate 0.21 between
study and retrieval Transport and 0.29 between study and Transport-specific
retrieval residuals.

No item predictor improves held-out 1--4 newness prediction. For old probes,
the largest point estimate is the Transport-specific item predictor at 0.0037
bits/trial (-0.0139 to 0.0206). For lures, the largest is study item propensity
at 0.0016 bits/trial (-0.0140 to 0.0171). All ordinary and familywise intervals
include zero. These results support modest stimulus-specific gaze organization,
especially during repeated viewing, but not a reliable retrieval or confidence
effect in this sample.
