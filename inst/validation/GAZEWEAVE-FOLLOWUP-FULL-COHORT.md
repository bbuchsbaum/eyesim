# GazeWeave full-cohort repetition and phase follow-up

Status: frozen before follow-up similarity scoring
Protocol version: 1
Seed: 20260826
Mote: `bd-01M0G9SX99EGTK21XBMT8073RT`

## Evidential boundary

This is a secondary follow-up to the completed maximal-cohort recognition
court. It does not alter or rescue that court's `not_supported` verdict. The
follow-up asks two questions motivated by the experimental design:

1. Can held-out gaze identify an image across its four study presentations?
2. Is retrieval information localized to the probe or delay and conditional
   on probe type, saliency, and the participant's response?

The raw fixation exports, checkpoints, candidate scores, and participant-linked
trial scores remain local-only and Git-ignored. Aggregate protocols and
verdicts may be published.

## Frozen support and candidate task

Reuse the score-blind GW-13 cohort, item folds, participant folds, and five-item
candidate construction:

- 36 participants;
- 1,295 participant-item trials;
- 120 base images;
- four crossed participant-by-item outer folds;
- five candidates per evaluated target, including the truth exactly once;
- no evaluated participant or item in training;
- identical candidate identities across methods within a task.

Study presentations 1--4 must each contain at least three retained fixations,
as already required by the GW-13 cohort. Retrieval phase eligibility is based
only on whether the declared interval contains at least one retained fixation;
it may not depend on gaze score, saliency, probe type, response, or accuracy.
Support is common across calibrated GazeWeave and density methods within each
phase. MultiMatch is reported only where its fixation-count contract is met.

## Gate A: study-repeat sensitivity ceiling

Use presentation 1 as the reference and separately score presentations 2, 3,
and 4 as sources:

```text
P1 -> P2
P1 -> P3
P1 -> P4
```

Each contrast receives its own training-only warp, model parameters, and
calibration. Do not concatenate presentations into a scanpath: that would
invent transitions between trials. A multi-presentation latent template is a
later model extension, conditional on this ceiling.

Primary positive-control requirement:

- for Replay or Transport v2, mean held-out information for P1 -> P4 has a
  crossed participant-by-item 95% interval above zero; and
- top-one credit is descriptively above the five-candidate chance rate of 0.20.

Report information, log loss, rank, top-one credit, spatial residual,
chronology, coverage, and warp diagnostics for every repetition contrast.
Estimate the P1-reference linear change from P2 through P4, but do not require
it to be positive: repeated viewing may become more selective rather than more
similar. Also report fixation count and dwell-time changes independently of
similarity.

## Gate B: phase-conditional recognition

Score the following retrieval intervals independently:

```text
probe        [0, 500) ms
delay        [500, 3000) ms
early_delay  [500, 1500) ms
late_delay   [1500, 3000) ms
```

Presentation 4 remains the reference for this gate so phase localization is
not confounded with a new template definition. The already completed combined
window remains an external descriptive reference and is not rescored here.

Represent probe type and response jointly:

```text
old  + correct   = hit
old  + incorrect = miss
lure + correct   = correct_rejection
lure + incorrect = false_alarm
```

Treat saliency as the five observed levels (20, 40, 60, 80, 100) for cell
descriptions. Do not assume a linear dose response. The primary engine family
contains these prespecified contrasts:

1. old hit minus miss at 20% saliency during the full delay;
2. old hit minus miss at 20% saliency during late delay;
3. the delay-minus-probe change in that low-saliency hit-minus-miss contrast;
4. lure false alarm minus correct rejection at 20% saliency during delay;
5. the difference between old hit-minus-miss at 20% and 100% saliency during
   delay.

The first three test partial-probe-triggered reinstatement and persistence. The
lure contrast tests whether an old response, rather than objective accuracy,
tracks apparent reinstatement. The 100% contrast distinguishes memory-dependent
delay information from visibility-related probe matching.

Use shared crossed participant-by-item bootstrap weights. Apply one familywise
interval across the ten Replay/Transport method-by-contrast tests. Report all
cell counts and reject any contrast whose observed cell is empty or whose
bootstrap valid-draw fraction is below 0.95. Density sigma 80 and 160 and the
registered density composite are mandatory comparators. MultiMatch is
secondary and may not define the support of short probe paths.

## Computational and correctness gates

Before non-smoke scoring, tests must establish:

- exact repetition assignment and row-order invariance;
- no invented cross-presentation transitions;
- half-open phase boundaries and duration conservation;
- deterministic phase eligibility and candidate construction;
- one truth and five candidates per evaluated set;
- zero participant and item overlap in every outer fold;
- correct response-class mapping;
- score support shared by mandatory methods within each phase;
- finite calibrated evidence and explicit handling of unavailable MultiMatch;
- reproducible checkpoints, runtime receipts, convergence, and manifests.

Execution is staged for cost control without outcome-dependent method
selection: Replay, mandatory density comparators, and Transport v2 are all
declared before scoring. Transport may run after the faster methods, but the
order does not change the estimand or inclusion rule.
