# GazeWeave implementation decision

Decision date: 2026-08-16
Evidence update: 2026-08-20
Development package version: `0.1.0.9000`
Repository base commit: `2d87fb87aea9f33cd84ca643032d70c6f9c4980e`
Status: no universal engine default; explicit engine choice required pending
external, task-matched imagery evidence.

## Decision

GazeWeave Replay is the directional encoding-to-recall engine. GazeWeave
Transport v2 is the symmetric engine. The real-data court did not justify
silently choosing either one, so `gaze_weave_cv()` requires `engine` unless it
can infer the engine from an explicit specification. The original fused
partial-transport implementation remains available when a legacy
`gaze_weave_spec()` is passed, and is named `transport_v1` in the common
dispatcher.

The sole primary endpoint is held-out item information gain:

```
gaze_info_bits = log2(p_true / prior_true)
```

Candidate probabilities, model parameters, warps, and baseline weights are fit
without the evaluated participant-item units. Coverage, spatial error,
chronology, contraction, template rank, and alignment stability are retained as
diagnostics, not additional primary similarity measures.

This decision is deliberately narrower than a superiority claim. Replay remains
fast, directional, and interpretable. A resolution-invariant calibration
correction repaired its repeated-viewing proper-score failure on a disjoint
participant-item sample, but neither GazeWeave engine achieved positive mean
information in that sample's blank-screen imagery control. The evidence does
not establish that GazeWeave has better power than registered density or a
cross-validated MultiMatch composite. A later maximal-cohort recognition court
likewise found no reliable saliency, correctness, or interaction effect for any
GazeWeave engine and no corrected comparator effect.

## Frozen validation designs

### Synthetic perturbation court

- Protocol: `GAZEWEAVE-V2-COURT.md`
- Seed: `20260816`
- Engines: Replay and Transport v2
- Fair baselines: learned MultiMatch ridge composite; raw and equivalently
  registered density ridge; elastic spatial matching where available
- Manipulations: fixation refinement, contraction, translation, time dilation,
  spatial noise, local swaps, block reorder, reversal, deletion, insertion,
  central bias, fixation count, duration concentration, and candidate-pool
  difficulty
- Outputs: `gaze-weave-v2-results/`
- Human verdict: `GAZEWEAVE-V2-VERDICT.md`

### Real positive and negative controls

- Protocol: `GAZEWEAVE-REAL-CONTROLS.md`
- Seed: `20260817`
- Source data: `test_data/study_fixations_all.csv` and
  `test_data/testdelay_fixations.csv`
- Participants: 104, 121, 124, 128, 315, 317, 319, 327
- Images: 5, 18, 69, 73, 81, 82, 88, 105
- Outer design: two participant folds crossed with two item folds; training for
  each of the four cells excludes every evaluated participant and item
- Candidate set: four within-participant encoding items
- Positive controls: presentation 1 to presentation 4 repeated viewing; fourth
  presentation to the 3000-ms post-test blank-screen imagery interval
- Negative controls: shuffled fixation order, wrong item, and generic
  participant gaze
- Engines: Replay and Transport v2
- Fair baselines: nested-ridge MultiMatch composite, registered-density ridge,
  and elastic-matching ridge
- Uncertainty: 1000-draw crossed participant-by-item bootstrap
- Outputs: `gaze-weave-real-results/`
- Human verdict: `GAZEWEAVE-REAL-VERDICT.md`

The real court produced no skipped comparators and no failed candidate fits.
Six of seven prespecified gates passed. The repeated-viewing gate failed because
Replay mean information was below zero despite strong top-one rank performance.

### Resolution-invariant disjoint replication

- Protocol: `GAZEWEAVE-REPLICATION-CONTROLS.md`
- Seed: `20260818`
- Participants: 9, 18, 111, 130, 300, 302, 316, and 325
- Images: 4, 11, 24, 62, 63, 80, 92, and 101
- Disjointness: no participant or item overlaps the original real court
- Algorithm change: Replay candidate evidence uses mean log likelihood per
  normalized duration bin; the raw total likelihood remains diagnostic
- Evaluation contract: unchanged tasks, folds, candidate pools, controls,
  engines, baselines, endpoints, thresholds, and 1000-draw bootstrap
- Outputs: `gaze-weave-replication-results/`
- Human verdict: `GAZEWEAVE-REPLICATION-VERDICT.md`

The correction was specified before replication scoring. The replication again
passed six of seven gates. Repeated viewing now passed for both engines, while
the blank-screen imagery gate failed because neither engine had positive mean
information. This is an internal replication within the same deposited dataset,
not independent external validation.

### Phase-resolved probe-delay persistence court

- Protocol: `GAZEWEAVE-PROBE-DELAY-CONTROLS.md`, version 2
- Seed: `20260820` with data-only derived seeds through `20260823`
- Source: local, Git-ignored phase-resolved Wynn fixation exports
- Frozen sample: eight participants, 80 participant-item trials, five
  candidates per held-out participant and item fold
- Tasks: repeated viewing; 0--500 ms probe; 500--3000 ms delay; combined;
  early delay; and late delay
- Outputs: local, Git-ignored `gaze-weave-probe-delay-results/`
- Human verdict: `GAZEWEAVE-PROBE-DELAY-VERDICT.md`

The court recovered substantial repeated-viewing point estimates and top-one
credit (0.513 for Transport v2 and registered density; 0.475 for Replay), but
all crossed intervals remained wide. No method carried reliable item
information through the delay. Transport v2 was +0.003 bits over the full delay
and -0.004 bits late; Replay was -0.068 and -0.033 bits. The predeclared verdict
was `not_supported`. This strengthens the no-default, no-superiority decision
and cannot close the independent Wang-data gate.

### Secondary recognition sensitivity court

- Protocol: `GAZEWEAVE-RECOGNITION-SENSITIVITY.md`, version 1
- Seed: `20260824`
- Source: immutable full-trial scores from the phase-resolved court
- Sample: eight participants, 80 signal trials, 59 correct and 21 incorrect
- Estimands: saliency 20-to-100 change, correct-minus-incorrect at saliency 60,
  and their interaction, adjusted for old-versus-lure probe type
- Uncertainty: 2,000 shared crossed participant-by-item bootstrap draws;
  familywise interval over six GazeWeave engine contrasts
- Outputs: local, Git-ignored
  `gaze-weave-recognition-sensitivity-results/`
- Human verdict: `GAZEWEAVE-RECOGNITION-SENSITIVITY-VERDICT.md`

The secondary verdict was `not_supported`. Transport v2 estimated +0.144 bits
from saliency 20 to 100 and -0.067 bits for correct versus incorrect responses
at saliency 60; Replay estimated +0.252 and -0.171 bits. All ordinary and
familywise intervals included zero, as did every probe- and delay-specific
GazeWeave interval. Registered elastic matching had the only ordinary
full-trial comparator interval excluding zero (+0.151 correctness bits), but
its family-comparable interval included zero and its phase localizations were
null. This partially post-hoc court does not alter the no-default,
no-superiority decision or close the independent Wang-data gate.

### Maximal-cohort recognition sensitivity court

- Protocol: `GAZEWEAVE-RECOGNITION-FULL-COHORT.md`, version 1
- Seed: `20260825`; item-fold seed `20260820`
- Source: local, Git-ignored Wynn fixation exports
- Sample: 36 participants and 1,295 combined-window trials, retaining 96.1% of
  score-blind eligible participant-item trials
- Evaluation: four disjoint participant-by-item folds and five identical
  candidates per held-out trial for every method
- Methods: Transport v2, Replay, raw density at sigma 80 and 160 pixels, all
  six raw MultiMatch dimensions, registered learned density and MultiMatch
  composites, and registered elastic matching
- Outputs: local, Git-ignored
  `gaze-weave-recognition-full-cohort-results/`
- Human verdict: `GAZEWEAVE-RECOGNITION-FULL-COHORT-VERDICT.md`

The full-cohort verdict was `not_supported`. Transport estimated +0.023 bits
from saliency 20 to 100, -0.011 bits for correct versus incorrect at saliency
60, and -0.007 interaction bits per 20 saliency points; every ordinary and
familywise interval included zero. Replay was also null, and no comparator
survived its descriptive multiplicity adjustment. Transport had the largest
mean calibrated information (+0.0073 bits), but density sigma 80 (+0.0059),
MultiMatch position (+0.0060), and Replay (+0.0032) were similarly close to
zero. This court supersedes the eight-participant sensitivity pilot but remains
same-study evidence and cannot close the independent Wang-data gate.

During this court, an equivalence-tested standard-domain masked Sinkhorn backend
reduced like-sized Transport fold runtime by 3.4- to 5.8-fold. Five real-path
comparisons agreed with the reference backend to at most `4.44e-16` score
difference. Replay remains far faster, so computational parsimony continues to
favor Replay for directional encoding-to-recall applications.

## Original real-court evidence

| Task | Method | Information bits | Log loss | Top-one credit |
|---|---|---:|---:|---:|
| Repeated viewing | Registered density ridge | 0.467 | 1.062 | 0.547 |
| Repeated viewing | Transport v2 | 0.271 | 1.199 | 0.531 |
| Repeated viewing | MultiMatch ridge | 0.100 | 1.317 | 0.469 |
| Repeated viewing | Elastic ridge | -0.024 | 1.403 | 0.391 |
| Repeated viewing | Replay | -0.326 | 1.612 | 0.609 |
| Blank-screen imagery | Transport v2 | 0.038 | 1.360 | 0.391 |
| Blank-screen imagery | Registered density ridge | 0.032 | 1.364 | 0.313 |
| Blank-screen imagery | Elastic ridge | 0.002 | 1.385 | 0.297 |
| Blank-screen imagery | MultiMatch ridge | -0.005 | 1.390 | 0.297 |
| Blank-screen imagery | Replay | -0.197 | 1.523 | 0.359 |

Every crossed-bootstrap information interval for blank-screen imagery included
zero. Replay minus registered density was -0.511 information bits, with a
crossed-bootstrap interval of -1.478 to 0.123. No superiority rule passed.
MultiMatch showed the clearest sensitivity to shuffled order; registered
density and elastic matching were order-invariant by construction.

Replay completed the primary real-data fits in 8.3 seconds, versus 1970.1
seconds for the current R Transport v2 solver and 111.9 seconds for all baseline
fits combined. Runtime supports Replay's computational parsimony but does not
override its failed proper-score gate.

## Disjoint replication evidence

| Task | Method | Information bits | Log loss | Top-one credit |
|---|---|---:|---:|---:|
| Repeated viewing | Registered density ridge | 0.462 | 1.066 | 0.516 |
| Repeated viewing | Transport v2 | 0.384 | 1.120 | 0.531 |
| Repeated viewing | Replay | 0.229 | 1.227 | 0.453 |
| Repeated viewing | Elastic ridge | 0.015 | 1.376 | 0.203 |
| Repeated viewing | MultiMatch ridge | -0.672 | 1.852 | 0.422 |
| Blank-screen imagery | Elastic ridge | 0.022 | 1.371 | 0.422 |
| Blank-screen imagery | Transport v2 | -0.002 | 1.387 | 0.344 |
| Blank-screen imagery | Registered density ridge | -0.011 | 1.394 | 0.266 |
| Blank-screen imagery | MultiMatch ridge | -0.045 | 1.417 | 0.156 |
| Blank-screen imagery | Replay | -0.072 | 1.436 | 0.328 |

The resolution correction transferred: Replay repeated-viewing information
changed from negative in the original court to +0.229 bits in a disjoint sample,
and its calibration temperatures were finite, non-boundary residual scales.
The imagery bootstrap intervals nevertheless included zero for Replay,
Transport, and registered density. Replay remained 0.147 pooled information
bits behind registered density, with a 95% interval from -0.340 to 0.012. No
superiority rule passed.

Replay completed the replication fits in 12.0 seconds, versus 2516.6 seconds
for Transport v2 and 141.0 seconds for all baselines. This strengthens the
computational-parsimony claim, but not a predictive-superiority claim.

## Interpretation contract

- Replay is directional: recall is modeled as locally ordered replay,
  occasional chunk restarts, and candidate-independent background gaze.
- Replay candidate evidence is calibrated from mean log predictive density per
  duration bin. The raw total grid likelihood remains diagnostic and is not
  exponentiated with a grid-length-dependent unit prior.
- Replay braid weights are posterior correspondence probabilities. Replay
  coverage and background mass therefore have posterior meanings conditional
  on the fitted model.
- Transport v2 is symmetric. Its coupling is a regularized optimized
  correspondence, not posterior psychological uncertainty.
- Registration is candidate-invariant and learned out of fold. The recommended
  nuisance model is screen-centered isotropic contraction plus participant or
  session calibration translation.
- Pair-specific free translation and affine registration are not defaults.
- The candidate set is part of the estimand. `gaze_info_bits` values from
  different candidate pools are not automatically interchangeable.

## Machine-readable evidence

Each result directory contains a `manifest-md5.csv` with file digests. The real
court additionally includes fold assignments, candidate
audits, convergence summaries, perturbation summaries, timing, method
availability, and the frozen decision object. Validation tests verify manifests,
fold disjointness, finite scores, null behaviour, and the declared verdict.

Reproduction commands from the package root are:

```r
source("inst/validation/gaze-weave-v2-court.R")
run_gaze_weave_v2_court("inst/validation/gaze-weave-v2-results")

source("inst/validation/gaze-weave-real-controls.R")
run_gaze_weave_real_court("inst/validation/gaze-weave-real-results")

source("inst/validation/gaze-weave-real-replication.R")
run_gaze_weave_real_replication(
  "inst/validation/gaze-weave-replication-results"
)

source("inst/validation/gaze-weave-probe-delay.R")
run_gaze_weave_probe_delay_court(
  "inst/validation/gaze-weave-probe-delay-results"
)

source("inst/validation/gaze-weave-recognition-sensitivity.R")
run_gaze_weave_recognition_sensitivity(
  "inst/validation/gaze-weave-recognition-sensitivity-results"
)

source("inst/validation/gaze-weave-recognition-full-cohort.R")
run_gaze_weave_recognition_full_cohort(
  output_dir = "inst/validation/gaze-weave-recognition-full-cohort-results",
  workers = 4L
)
```

## Unresolved limitations

1. Both real courts use disjoint subsets of one dataset distributed with
   `eyesim`; neither is independent external replication.
2. Blank-screen imagery evidence is weak for every evaluated method.
3. Participant/session contraction estimates are diagnostically useful but can
   remain confounded with which image regions were fixated on short trials.
4. Replay parameters are selected from a frozen grid rather than a fully
   hierarchical Bayesian model.
5. Transport v2 is substantially faster after the masked-Sinkhorn rewrite but
   remains computationally expensive relative to Replay; its coupling can also
   be non-unique even when its scientific score is stable.
6. Calibration must be rechecked when candidate-pool size, task, population,
   fixation detector, or coordinate system changes.
7. No current evidence supports replacing MultiMatch's diagnostic dimensions
   when the scientific question concerns a particular component such as
   direction or duration rather than item identification.

The next scientific gate is preregistered, external validation with larger
participant and item samples, a task-matched candidate pool, and the same
nested calibration and baseline policy.
