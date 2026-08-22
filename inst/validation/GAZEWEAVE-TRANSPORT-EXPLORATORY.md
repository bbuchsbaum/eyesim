# Canonical Transport: exploratory recognition heterogeneity

This is a post-court exploratory analysis of the new confidence-rating study. It
cannot revise the frozen GW-18.10 result, and it does not create a confirmatory
recognition claim. The analysis uses the unchanged 1,295 out-of-fold Transport
scores from the exact 0--3000 ms court. Trial-linked outputs remain local and
Git-ignored.

## Questions

1. Does Transport item information vary with saliency, recognition
   correctness, and old-versus-lure probe type?
2. Do those heterogeneous terms improve response prediction in held-out
   participant and item cells?
3. Is the conclusion sensitive to response-blind effective-fixation shrinkage?

The explanatory model is
`gaze_info_bits ~ saliency * correctness * probe_type + gaze_quality + duration`.
Uncertainty uses the crossed participant/item bootstrap with all seven
hierarchical contrasts reported and a family-wise interval.

The predictive court compares a base accuracy model, the base plus an additive
Transport score, and the base plus all hierarchical Transport-by-saliency-by-
probe terms. Models are trained only on rows sharing neither participants nor
items with the evaluation fold. The decision quantity is held-out log-loss
gain in bits, not an in-sample coefficient.

The cheap parameter sensitivity recomputes the final posterior for every
declared reliability value `kappa = 0, 0.25, 0.5, 1, 2, 4, 8, 16` from the
stored pre-shrink candidate posteriors. It reports every value; it does not
select a winner. Chronology, spatial scale, and coverage-prior changes require
new alignments and are not conflated with this exact posterior sensitivity.

## Result

All 20 saliency-by-correctness-by-probe cells are represented in the full
cohort (minimum cell size 5). No crossed interval for the seven explanatory
contrasts excludes zero. In particular, the saliency-by-correctness-by-probe
coefficient is `-0.02603` bits per 20 saliency units (crossed 95% interval
`-0.11604` to `0.07298`; family-wise interval `-0.15372` to `0.12014`).

An additive Transport score changes held-out accuracy information by
`-0.00116` bits (95% interval `-0.00519` to `0.00265`). The complete
Transport-by-saliency-by-probe model changes it by `-0.00651` bits relative to
the no-Transport base (`-0.02384` to `0.00802`) and by `-0.00535` bits relative
to the additive Transport model (`-0.02290` to `0.00942`). The heterogeneous
model therefore adds variance, not predictive evidence.

Across all eight declared reliability values, the three-way interval includes
zero and the heterogeneous held-out gain remains negative (`-0.00596` to
`-0.00679` bits across point estimates). Removing shrinkage does not expose a
recognition effect; increasing shrinkage steadily attenuates both information
and interaction magnitude.

These results do not support a saliency-, correctness-, or probe-specific
recognition rescue for Transport.

## Deeper alignment sensitivity

A ten-cell, one-factor-at-a-time grid is locked in `alignment-grid.csv` before
any further rescoring: the default; temporal weights 0, 1, and 4; chronology
neighbourhoods 1 and 3; spatial-scale multipliers 0.67 and 1.5; and beta
coverage priors (1, 1) and (4, 4). All other settings, the four equal-prior
episodes, candidates, folds, and response-blind calibration remain fixed.

The default cell can reuse the frozen checkpoint; the nine new cells require
fresh alignments. At the observed default runtime this is approximately 5.4
serial hours. If run, every cell must be reported. No behavioral winner may be
selected or promoted to a confirmatory result.

Run after loading the package:

```r
source("inst/validation/gaze-weave-transport-exploratory.R")
run_gaze_weave_transport_exploratory()
```

The hashed fold checkpoints are verified before outcomes are merged. The
canonical source has since received API-only renaming, so this analysis anchors
the frozen measurement to those byte-verified score checkpoints rather than
pretending the old source manifest describes the renamed checkout.
