# Frozen behavioral court for multi-presentation Replay

Status: frozen before behavioral score merge
Protocol version: 1
Seed: 20260829
Measurement dependency: `GAZEWEAVE-PCMRI-MULTI-REPLAY-VERDICT.md`

## Question

Among old probes, does the now-frozen item-specific Replay score improve
out-of-fold prediction of an old response beyond saliency, generic gaze
quality, and candidate-pool size?

The measurement court was completed without response, correctness, probe type,
or saliency columns in its scoring table. It authorized exactly one primary
behavioral score: `replay_all4_shrunk_exhaustive` from the full `[0,3000)`-ms
recognition interval.

## Support and outcome

- Merge behavior only after verifying the frozen measurement manifest.
- Retain old probes with a non-missing response on the 1--4 confidence scale.
- Code responses 1--2 as `said_old = 1` and 3--4 as `said_old = 0`; response 0
  remains missing.
- Preserve the four crossed participant-by-item outer folds used by the
  measurement court.
- Do not inspect lure false alarms, confidence as an ordinal outcome, saliency
  subgroups, alternate windows, or nonlinear effects in this court.

The expected usable old-item support is 989 trials: 892 old responses and 97
new responses.

## Predictive endpoint

For each outer fold, train on trials sharing neither evaluation participants
nor evaluation items. Standardize continuous predictors using training rows
only.

The base logistic model is:

```text
said_old ~ sal10 + z_log_effective_fixations +
  z_total_duration + z_log_candidate_count
```

The full model adds `z_gaze_info_bits`. The primary endpoint is held-out log
loss reduction in bits per trial:

```text
(loss_base - loss_full) / log(2)
```

Use an independent participant-by-item crossed bootstrap with 2,000 draws for
the mean predictive gain. Positive evidence requires the primary Replay 95%
interval to exclude zero. Report the result regardless of sign.

## Prespecified comparators

Run the same base/full prediction for:

1. all-four shrunk Replay (primary);
2. all-four unshrunk Replay;
3. presentation-4 Replay;
4. all-four density sigma 80.

These are descriptive comparator rows, not a multiplicity-free family of new
primary tests. A secondary nested prediction asks whether frozen Replay adds
information beyond density sigma 80 by using `base + density` as the reduced
model and adding Replay as the full model.

## Association audit

Fit the full old-item logistic association on all usable trials and obtain
coefficient intervals by crossed participant-by-item weighted bootstrap. For
the primary Replay model, retain coefficients for:

- frozen Replay information;
- saliency;
- log effective fixation count;
- retained gaze duration;
- log candidate count.

A crossed-random-intercept `glmer` with participant and item intercepts is a
descriptive specification check only. Singularity or convergence failure does
not replace the predictive endpoint and must be reported.

## Interpretation

A positive Replay coefficient or predictive gain is item-specific only to the
extent enforced by the exhaustive contrastive score. Generic fixation support
is reported separately and included in the base model. Density sigma 80 is the
strongest measurement comparator and must remain visible even if Replay is
nominally significant.

No change to the frozen measurement algorithm, reliability parameter,
candidate set, time window, or presentation mixture is permitted after this
behavioral merge.
