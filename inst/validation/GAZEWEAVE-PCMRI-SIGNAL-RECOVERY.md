# GazeWeave pcmri signal-recovery gate

Status: frozen before expanded-support scoring
Protocol version: 1
Seed: 20260827
Mote: `bd-01M0GK8JRQNKXFC6P68WWGGZV7`

## Question

The common-support comparison required at least three retained fixations in
every study presentation and in recognition so that Replay, density,
Transport, and MultiMatch were evaluated on identical trials. That is the
right court for comparing methods, but it discards many trials on which Replay
and density are defined. This gate asks a narrower question:

> Does Replay show item sensitivity or held-out behavioral information when it
> is evaluated on its natural short-path support rather than MultiMatch's
> support?

This is a prespecified measurement-recovery analysis, not permission to search
over windows, outcomes, saliency codings, or model forms. The complete
recognition interval remains the half-open interval `[0, 3000)` ms (500-ms
probe plus 2500-ms delay).

## Frozen support

- Use presentation 4 as the encoding reference.
- Include an old or lure trial when both the presentation-4 study path and the
  complete recognition path contain at least one retained fixation.
- Retain a participant only when each of the two fixed item folds contains at
  least five eligible candidate images.
- Keep five candidates per target, including the target exactly once, using
  the deterministic within-participant candidate construction from the
  full-cohort court.
- Keep the same crossed two-participant-fold by two-item-fold design. Training
  and evaluation participants and items must be disjoint.
- Use the existing spatial crop (the 800 by 600 image rectangle) and retained
  fixation-duration rule. Do not widen either after scoring.

The frozen expected cohort is 2,222 target trials from 45 participants and 120
items (1,102 old and 1,120 lure trials). The pre-retention pool is 2,225 trials
from 46 participants; participant 1023 lacks five candidates in one item fold.

## Methods and endpoints

Score only:

1. GazeWeave Replay with its existing cross-fitted warp and calibration;
2. duration-weighted density similarity at sigma 80 px, calibrated on the same
   folds and candidate sets.

MultiMatch remains in the 1,295-trial common-support court. It is not assigned
synthetic values on one- or two-fixation paths. Transport v2 is omitted because
this gate tests support loss, not the slower symmetric engine.

The primary measurement endpoint is held-out `gaze_info_bits`. Report its
crossed participant-by-item interval, mean rank, and top-one credit. The
primary behavioral endpoint, restricted to old trials with a response, is
held-out log-loss improvement over

```text
said_old ~ sal10
```

after adding training-standardized `gaze_info_bits`. Behavioral models are
trained and evaluated with both participants and items disjoint. A crossed
bootstrap association is secondary and cannot override negative held-out
prediction.

## Interpretation gates

- Expanded support recovers item sensitivity only if Replay's crossed interval
  for mean `gaze_info_bits` excludes zero.
- It recovers behavioral information only if Replay's crossed interval for
  held-out behavioral information excludes zero.
- Replay is better than density only if the paired crossed interval for the
  relevant Replay-minus-density endpoint excludes zero.
- An increase relative to the 1,295-trial court is descriptive because the
  cohorts differ; it is not a paired performance comparison.
- A nominal mixed-model coefficient, saliency subgroup, phase window, ordinal
  response, or nonlinear term is not a recovery unless it passes a separately
  frozen out-of-fold gate.

Participant-linked scores, candidate scores, checkpoints, predictions, and
bootstrap draws stay under the Git-ignored pcmri result directory. Only this
protocol and an aggregate findings document may be published.
