# GazeWeave maximal-cohort recognition court

Status: completed; frozen before full-cohort similarity scoring
Protocol version: 1
Primary seed: 20260825
Mote: `bd-01M0FRA5QSTQZJ4GMNVH6MBN6R`

Execution note, added after scoring: folds 1--2 used the original log-domain
masked Sinkhorn projection. Profiling identified that projection as the dominant
Transport cost. Folds 3--4 used an equivalence-tested standard-domain matrix
scaling backend with automatic log-domain fallback. This changed neither the
objective nor any protocol choice. Randomized, synthetic-alignment, and five
real-path differential checks agreed within floating-point precision; the full
receipt is in `GAZEWEAVE-RECOGNITION-FULL-COHORT-VERDICT.md`.

## Purpose and evidential boundary

This court supersedes the eight-participant GW-12 pilot for claims about
whether held-out gaze reinstatement information varies with probe saliency or
recognition correctness. It uses the complete 0--3000 ms retrieval interval:
the 500-ms visual probe plus the subsequent 2500-ms delay.

The eight-person pilot result was inspected before this protocol. No
full-cohort similarity score or saliency/correctness coefficient was inspected
before freezing this version. This remains an internal, same-study analysis;
it cannot establish external replication, close the Wang-data gate, or by
itself justify a universal default engine.

Raw fixation exports, candidate-level scores, trial-level scores, checkpoints,
and fitted model objects remain local-only and Git-ignored. Aggregate protocol,
test, and verdict documents may be published.

## Data-only maximal cohort

Use the two checksum-verified local Wynn exports declared in
`GAZEWEAVE-PROBE-DELAY-CONTROLS.md`. Eligibility is evaluated without looking
at a similarity score, response correctness, saliency coefficient, or fitted
model.

A participant-item target is eligible when:

- retrieval probe type is `old` or `lure`;
- all four study presentations are present and each retains at least three
  valid fixations in [0, 2500) ms;
- the combined retrieval interval [0, 3000) ms retains at least three valid
  fixations; and
- study and retrieval mappings are unique for participant plus base image.

Assign the 120 base image identities globally to two item folds using the
previously declared data-only item-fold seed 20260820. A participant enters the
primary cohort when both item folds contain at least five eligible targets.
No eligible target from a retained participant is sampled away.

The frozen data-only audit found:

- 46 participants in the exports;
- 1,348 eligible participant-item targets from 45 participants;
- 36 participants with at least five eligible targets in both item folds; and
- 1,295 retained targets, 96.1% of all combined-window eligible targets.

This is the maximal cohort under the fixed five-candidate estimand. Participants
with insufficient candidate support are excluded because changing candidate
pool size changes task difficulty and the interpretation of information gain.

## Folds and candidate sets

Assign the 36 retained participants to two folds using only their eligible
trial counts: process participants from largest to smallest count and allocate
each to the fold with the smaller accumulated count, resolving ties from seed
20260825. Cross participant fold with the two global item folds to make four
outer folds.

For each outer fold:

- evaluation targets have a held-out participant and held-out item identity;
- training excludes every evaluation participant and every evaluation item;
- warp, engine parameters, calibration, and learned comparator weights use
  training rows only; and
- every evaluation target is compared with exactly five study templates from
  the same participant and held-out item fold: its true template and four
  deterministic circular neighbours in a participant-fold-specific seeded
  item permutation.

The same candidate set and ordering are used by every method. Every candidate
set must contain the true item exactly once. Candidate construction may not use
saliency, correctness, probe type, gaze coordinates, duration, or any model
score.

Training calibration may use all eligible candidates in each training
participant's non-held-out item fold. Evaluation always uses the frozen
five-candidate pool. Candidate count, candidate identities, and fold overlap
are written to the local audit artifacts.

## Methods

### GazeWeave engines

- directional GazeWeave Replay;
- symmetric GazeWeave Transport v2.

Use the already frozen GW-11 engine specifications, candidate-invariant
screen-centred contraction plus translation, and training-only calibration.

### Standard density similarities

Compute duration-weighted Gaussian density maps on the 800-by-600 pixel screen
and compare them by cosine similarity without registration at:

- sigma = 80 pixels;
- sigma = 160 pixels.

For each bandwidth, retain both:

- the native raw similarity/rank result; and
- a one-feature ridge scale calibrated strictly inside training folds, yielding
  candidate probabilities and `gaze_info_bits` on the common five-candidate
  task.

The existing registered multiscale density ridge comparator remains in the
court and may use sigma values 30, 60, 80, 120, and 160 pixels. Its weights and
ridge penalty are nested-training estimates.

### MultiMatch landscape

Report all six current `eyesim` MultiMatch metrics separately:

- vector/shape;
- direction;
- length;
- position;
- duration; and
- position EMD.

For each raw metric, retain its native similarity/rank and a separately
training-calibrated one-feature information score. Also retain the existing
registered six-feature MultiMatch ridge composite. Do not average the metrics.
The optional MultiMatch dependencies must be available; a missing comparator
invalidates a non-smoke court rather than disappearing silently.

### Additional comparator

Retain the registered elastic-consensus ridge comparator from the previous
court.

## Outcomes and conditional model

The sole GazeWeave primary endpoint remains held-out item information:

```text
gaze_info_bits = log2(p_true / 0.2)
```

For each calibrated method, fit:

```text
gaze_info_bits ~ saliency_z * correct_ec + probe_type_ec
```

where:

```text
saliency_z    = (saliency - 60) / 20
correct_ec    = accuracy - 0.5
probe_type_ec = -0.5 for lure and +0.5 for old
```

Report:

1. predicted saliency 20-to-100 change, averaged over correctness;
2. correct-minus-incorrect at saliency 60; and
3. their interaction per 20 saliency points.

These are adjusted associations, not causal effects on memory.

## Uncertainty, multiplicity, and decision rule

- Use 2,000 shared crossed participant-by-item bootstrap draws.
- Require at least 95% converged held-out alignments for each GazeWeave engine;
  retain and report every finite nonconverged score rather than filtering it.
- Require at least 95% finite full-rank draws and design condition number at
  most 30.
- Audit the common crossed random-intercept model with `lme4::lmer`, ML,
  `bobyqa`, and participant/item random intercepts. Report singularity without
  changing the model.
- Report leave-one-participant-out signs; a GazeWeave detection must retain its
  sign in at least 80% of omissions and agree with the lmer estimate.
- The six Replay/Transport method-by-contrast tests form the primary family and
  receive Bonferroni 99.1667% intervals.
- Comparator conditional intervals, ranks, top-one credit, information, and
  log loss form a descriptive landscape. Also report Benjamini-Hochberg
  adjusted bootstrap sign probabilities across the calibrated comparator
  method-by-contrast table, but do not use them to rescue the GazeWeave gate.

The GazeWeave sensitivity label is:

- `supported`: an engine familywise interval excludes zero and all integrity
  and stability gates pass;
- `suggestive`: an engine ordinary 95% interval excludes zero but the
  familywise or stability rule fails;
- `not_supported`: neither engine ordinary interval excludes zero, or an
  integrity gate fails.

A comparator matching or exceeding GazeWeave precludes a uniqueness or
superiority claim. Phase-localization analysis is deferred until after the
combined-window verdict and cannot redefine this court.

## Computational and audit gates

Tests must establish:

- checksum verification and exact eligibility support;
- deterministic cohort and candidate construction under row permutation;
- five candidates, one truth, one participant, and one item fold per set;
- zero participant and item overlap in every outer fold;
- every retained target evaluated exactly once;
- analytic and metamorphic correctness of sigma-specific density output;
- availability of every requested native and calibrated MultiMatch metric;
- shared candidate identities across all engines and comparators;
- finite calibrated scores, explicit skipped rows, convergence receipts, and
  runtime/memory receipts;
- shared crossed-bootstrap weights, rank-deficient-draw rejection, and
  singular mixed-model reporting; and
- local result-manifest integrity and Git-ignore protection.
