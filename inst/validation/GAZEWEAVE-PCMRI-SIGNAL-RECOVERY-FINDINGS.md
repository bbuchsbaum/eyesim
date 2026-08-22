# GazeWeave pcmri signal-recovery findings

Date: 2026-08-20
Protocol: `GAZEWEAVE-PCMRI-SIGNAL-RECOVERY.md`
Mote: `bd-01M0GK8JRQNKXFC6P68WWGGZV7`
Status: expanded-support gate complete

## Bottom line

The missing behavioral effect is not explained by the 1,295-trial common
MultiMatch support. Expanding Replay and density to all 2,222 trials on which
both paths contain at least one fixation did not produce reliable retrieval
item information or positive held-out recognition prediction.

GazeWeave itself is sensitive: Replay identifies the correct image strongly
across repeated study presentations. The failure is specific to retrieval,
where the item-identifying signal is tiny, single-trial encoding templates and
four-candidate denominators are noisy, and old responses are close to ceiling.

The most defensible next attempt is therefore a reliability improvement:
integrate the four encoding presentations as a mixture of templates, stabilize
the nonmatch denominator, and allow low-information paths to shrink to the
candidate prior. Searching saliency subgroups or alternative response models
is not supported by the diagnostics already run.

## Expanded-support result

The frozen gate retained 2,222 old/lure trials from 45 participants and 120
items. Each trial used five candidates; participant and item folds were both
disjoint between score training and evaluation.

| Endpoint | Replay | Density sigma 80 | Replay minus density |
|---|---:|---:|---:|
| Item information, bits | +0.00481 [-0.00594, +0.01681] | +0.00305 [-0.00488, +0.01098] | +0.00177 [-0.00598, +0.00984] |
| Held-out recognition information, bits/trial | -0.00498 [-0.01989, +0.00867] | -0.02782 [-0.05904, -0.00537] | +0.02284 [+0.00139, +0.05174] |

Replay's positive behavioral contrast with density does not mean that Replay
predicts recognition. Replay itself is slightly worse than the saliency-only
model and its interval includes zero; density generalizes still worse. The
secondary crossed association for Replay was +0.286 standardized log odds,
with interval [-0.056, +0.648]. It also includes zero.

All short paths were genuinely scored. There were no convergence failures, and
Replay took about 27--31 seconds per outer fold on this cohort. The negative
result is not an optimizer or runtime failure.

## What the added trials reveal

The short paths dilute rather than recover item signal.

| Minimum of study/test fixation counts | Trials | Replay bits | Density bits |
|---:|---:|---:|---:|
| 1 | 189 | -0.0486 | -0.0346 |
| 2 | 266 | -0.0276 | -0.0124 |
| 3 or more | 1,767 | +0.0154 | +0.0094 |

Replay coverage rises from 0.45 with one fixation to 0.69 with at least three,
while its mean spatial residual falls from 168 px to 112 px. This is sensible
diagnostic behavior, but the candidate posterior is still too decisive on
some low-information trials. A method that knows a path is weak should return
approximately zero information rather than confidently negative information.

Trials overlapping the previous 1,295-trial court have +0.0172 Replay bits
under the expanded model, while newly admitted trials average -0.0125 bits.
Because calibration and candidate pools change with the cohort, this is a
descriptive localization rather than a randomized support comparison.

## The positive control says the algorithm can see a strong signal

With the same 1,295-trial support and crossed candidate design, Replay gives:

| Study comparison | Replay bits | Crossed 95% interval | Top-one credit |
|---|---:|---:|---:|
| Presentation 1 to 2 | +0.488 | [+0.331, +0.638] | 0.498 |
| Presentation 1 to 3 | +0.366 | [+0.247, +0.484] | 0.452 |
| Presentation 1 to 4 | +0.257 | [+0.150, +0.365] | 0.422 |

Chance top-one credit is 0.20. For presentation 1 to 4, Transport (+0.289
bits), density sigma 80 (+0.271), and Replay (+0.257) are all strongly
positive and close to one another. Thus this dataset validates sensitivity to
repeated-viewing structure, but it does not demonstrate that Replay has more
power than density.

The decline from presentation 2 to presentation 4 is scientifically useful:
later viewing is not a literal replay of initial exploration. It also warns
against treating one study presentation as a noise-free encoding template.
Trial-level Replay scores across the three repeat contrasts correlate only
about 0.31--0.35, although participant means correlate about 0.64--0.69.

## Why the recognition association is hard to recover

### The retrieval item signal is extremely small

Across the common court, every method is near zero: Transport +0.0073 bits,
density +0.0059, MultiMatch position +0.0060, and Replay +0.0032. Expanding
support leaves Replay at +0.0048 bits. This is orders of magnitude below the
study-repeat Replay signal.

### Recognition is near ceiling and clusters, not rows, limit precision

Among 1,077 expanded old trials with a response, 973 are called old and only
104 are called new. Adding trials raises row count, but there are still only 45
participants and 120 crossed items. Fold-specific Replay coefficients range
from negative in one fold to positive in three folds. The crossed interval
implies that only a fairly large standardized association could be detected
reliably.

Using all four confidence categories does not solve this. On common support,
an ordinal mixed model gave Replay `p = .212`; probe-only, delay-only, early-
delay, late-delay, nonlinear, and saliency-interaction checks were also null.
Those checks localize the weak signal and argue against further window or
subgroup search.

### Generic gaze state predicts behavior more readily than item identity

Hits contain somewhat more usable gaze than misses: about 5.36 versus 4.89
retained fixations and 2,136 versus 1,910 ms retained dwell. Retained dwell is a
strong in-sample recognition predictor. Replay information is also correlated
with replay coverage.

The density decomposition has the same pattern. Raw matched similarity predicts
recognition more strongly than matched-minus-nonmatch similarity, and much of
the advantage is carried by similarity to other images. Replay likewise shows
in-sample associations with nonmatch compatibility and coverage. None of these
components improves strict crossed behavioral prediction.

This suggests two distinct constructs:

1. item-specific replay, measured by candidate information;
2. generic retrieval engagement or gaze quality, reflected by dwell,
   coverage, central tendency, and broad compatibility.

The second may be behaviorally meaningful, but it is not episodic
reinstatement and must not be folded into the primary replay score merely to
obtain a stronger coefficient.

### The candidate denominator and encoding template are noisy

Five-candidate scoring uses only four nonmatches. Their within-trial score
spread is large relative to the true-template advantage, so which four images
enter the panel contributes appreciable measurement error. In addition, using
presentation 4 alone discards three observed encoding paths despite modest
trial-level agreement across repetitions. Both choices attenuate any relation
between latent replay strength and memory behavior.

## Recommended algorithmic gate

Implement and test the following as one frozen reliability court before any
new behavioral model is fitted.

### 1. Multi-presentation Replay likelihood

For candidate item `k`, preserve its four study paths separately and integrate
their predictive likelihoods:

```text
p(Y | item k) = sum_r w_r p(Y | study path k,r)
```

Use equal weights first, or learn shared repetition weights only inside the
training fold. Do not concatenate paths, because that invents saccades between
presentations. A later replay graph may pool spatial states and within-trial
transitions while retaining this no-invented-transition contract.

### 2. Stable nonmatch evidence

Replace one K=5 panel with either all eligible same-participant candidates in
the held-out item fold or the average of several deterministic candidate
panels. Candidate identities must remain identical across methods, and
calibration must use the same candidate-pool policy as evaluation. Report the
score's repeated-panel reliability before testing behavior.

### 3. Evidence-aware shrinkage to the prior

Learn, inside training folds, a candidate-invariant reliability gate based on
retrieval path support such as retained fixation count and duration. Combine
the model posterior with the candidate prior:

```text
p_final(k | Y) = q(Y) p_replay(k | Y) + (1 - q(Y)) prior(k)
```

This is an abstention mechanism, not a way to create signal. It should move
uninformative one- and two-fixation trials toward zero bits and must not use
response correctness, saliency outcome contrasts, or knowledge of the true
candidate.

### 4. Freeze the score, then test behavior

Only after the measurement gate is passed should the behavioral court be run.
The null model should include saliency and generic gaze-quality terms such as
retained dwell and fixation count. The specific Replay score must improve
held-out log loss beyond that baseline and survive crossed participant-by-item
uncertainty. Report generic gaze state separately.

### Lower-priority extension

A phase-aware Replay model could retain the full 3000-ms endpoint while using
the known 500-ms boundary to reset transitions and permit different probe and
delay background rates. Existing phase-specific null results make this less
promising than improving template and denominator reliability, but it is more
principled than selecting a favorable subwindow.

## What not to do

- Do not widen the spatial crop after seeing the result; it already matches the
  image rectangle and loses relatively few on-screen fixations.
- Do not freely optimize pair-specific translation or affine warps.
- Do not select the 20% saliency condition, one delay subwindow, or a nonlinear
  response model post hoc.
- Do not call raw or nonmatch similarity "reinstatement."
- Do not interpret Replay's advantage over a negatively generalizing density
  model as positive behavioral information.

## Parsimony and interpretation

Replay remains parsimonious at the inferential level: one item-information
endpoint plus coverage, spatial residual, restart, and warp diagnostics. Its
posterior gaze braid has a clearer meaning than an entropy-regularized
transport coupling, and the repeated-viewing court confirms that the endpoint
responds to real correspondence.

The proposed mixture and abstention add machinery, so they earn their place
only if they improve held-out proper scores and reliability. They should not
create additional primary outcomes. The honest present claim remains:

> GazeWeave is sensitive to repeated viewing and competitive with density and
> MultiMatch, but this recognition dataset contains too little reliable
> item-specific retrieval signal to establish a behavioral replay effect or a
> power advantage for GazeWeave.
