# GazeWeave recognition sensitivity court

Status: frozen before conditional model fitting
Protocol version: 1
Primary seed: 20260824
Mote: `bd-01M0FQ85Z8T3FW3Y5F1P9FE25V`

## Purpose and evidential boundary

This secondary court asks whether GazeWeave's already cross-fitted item
information is sensitive to retrieval-probe saliency or response correctness
during the complete 0--3000 ms recognition interval. It uses the immutable
trial scores from the completed GW-11 court; it does not refit an alignment,
warp, calibration model, candidate pool, or fixation path.

This analysis was proposed after the GW-11 primary persistence verdict was
known to be `not_supported`. Before this protocol was frozen, the analyst had
also seen task-level results and saliency-stratified delay descriptives, but had
not inspected correctness-conditioned information or the full-trial
saliency/correctness coefficients. The court is therefore explicitly
secondary and partially post hoc. It cannot rescue the persistence endpoint,
select a default engine, establish superiority, or close the independent Wang
gate.

The source fixation exports, GW-11 trial scores, and all outputs from this
court remain local-only and Git-ignored.

## Frozen input and support

Read `inst/validation/gaze-weave-probe-delay-results/trial-scores.csv` and
require its MD5 manifest to pass. Retain rows satisfying:

- `task == "combined"` (probe plus delay, 0--3000 ms);
- `condition == "signal"`;
- `calibrated == TRUE` and `status == "scored"`;
- method in Transport v2, Replay, registered-density ridge, registered
  MultiMatch ridge, or registered elastic ridge.

Every method must have the same 80 participant-item trials, eight participants,
five held-out candidates per trial, finite `gaze_info_bits`, and no duplicate
participant-item key. Frozen support is 59 correct and 21 incorrect trials;
saliency-by-correctness counts are:

| Saliency | Incorrect | Correct |
|---:|---:|---:|
| 20 | 4 | 9 |
| 40 | 6 | 15 |
| 60 | 3 | 11 |
| 80 | 4 | 8 |
| 100 | 4 | 16 |

There are 32 lure trials (16 correct) and 48 old trials (43 correct). Because
probe type and correctness are associated, probe type is a required covariate;
neither condition may be dropped after seeing results.

## Model and estimands

Code predictors as:

```text
saliency_z    = (saliency - 60) / 20
correct_ec    = accuracy - 0.5
probe_type_ec = -0.5 for lure, +0.5 for old
```

For each method, fit the common fixed-effects model

```text
gaze_info_bits ~ saliency_z * correct_ec + probe_type_ec
```

and audit it with the crossed random-intercept model

```text
gaze_info_bits ~ saliency_z * correct_ec + probe_type_ec +
                 (1 | participant) + (1 | item)
```

using `lme4::lmer(..., REML = FALSE)`. Singular random-effects fits are
reported, not silently simplified. Because 46 of 61 items occur once, primary
uncertainty comes from crossed participant-by-item bootstrap weighting of the
fixed-effects model, matching the GW-11 uncertainty policy. The lmer fit is a
hierarchical robustness check, not the source of default p-values.

Report three adjusted contrasts:

1. `saliency_20_to_100`: predicted change from saliency 20 to 100, averaged
   over correctness, equal to four times the `saliency_z` coefficient;
2. `correct_at_60`: correct minus incorrect at saliency 60, equal to the
   `correct_ec` coefficient;
3. `interaction_per_20`: change in the correctness contrast for each 20-point
   saliency increase, equal to the interaction coefficient.

These are associations with held-out item information, not causal effects on
memory accuracy.

## Uncertainty, multiplicity, and stability

- Generate 2,000 deterministic crossed bootstrap draws. For each draw, sample
  eight participants and all observed item identities with replacement, then
  weight each observed trial by the product of its participant and item
  frequencies.
- Reuse the exact draw weights for every method and phase.
- Require at least 95% finite, full-rank fits for each method.
- Report ordinary 95% percentile intervals for every calibrated method.
- Treat the six Transport/Replay method-by-contrast combinations as the
  inferential family. Also report Bonferroni familywise 99.1667% percentile
  intervals. Baseline intervals are comparator diagnostics and do not expand
  the engine family.
- Calculate the fixed-effect design-matrix condition number. If it exceeds 30,
  label all conditional effects unstable.
- Refit after leaving out each participant. A potentially detected contrast
  must retain its full-sample sign in at least seven of eight omissions and
  agree in sign with the corresponding lmer fixed-effect contrast.

The frozen sensitivity label is:

- `supported`: at least one GazeWeave engine contrast has a familywise interval
  excluding zero and passes design, bootstrap, lmer-sign, and leave-one-person
  stability checks;
- `suggestive`: an engine's ordinary 95% interval excludes zero but its
  familywise or stability gate fails;
- `not_supported`: neither engine has an ordinary 95% interval excluding zero,
  or the design/bootstrap integrity gate fails.

Any method may be sensitive; a baseline matching or exceeding GazeWeave must be
reported and precludes a uniqueness claim.

## Phase localization

Apply the identical model and bootstrap draws to `probe` and `delay` signal
scores as secondary localization analyses. These phase results do not create
additional detection claims. If the full-trial court is `not_supported`, they
remain descriptive even if an ordinary interval excludes zero.

Report fitted cell predictions for saliency 20--100, correct and incorrect,
separately by method and phase, with probe type averaged at its effect-coded
zero. Do not filter by observed score, response, saliency, convergence, or
participant after fitting.

## Computational tests

The implementation must test:

- input manifest verification and exact trial/method support;
- invariance to row permutation;
- identical bootstrap weights across methods and phases;
- analytic recovery of known linear coefficients on a deterministic fixture;
- finite handling and explicit rejection of rank-deficient bootstrap draws;
- exact contrast transformations from model coefficients;
- participant and item resampling multiplicities;
- singular lmer reporting without automatic model changes;
- local-only result-manifest integrity and the frozen verdict.
