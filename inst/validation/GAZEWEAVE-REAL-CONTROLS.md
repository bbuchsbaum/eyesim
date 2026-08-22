# GazeWeave real-data control court

Status: frozen before outcome scoring
Protocol version: 1
Random seed: 20260817

## Purpose

This court asks whether the provisional synthetic-court decision transfers to
real gaze data. It is not a tuning exercise and it cannot establish universal
superiority. The primary comparison is among GazeWeave Replay, GazeWeave
Transport v2, a nested calibrated MultiMatch composite, a nested registered
multiscale-density composite, and the order-free elastic-consensus baseline.

The court has two directional tasks:

1. **Repeated viewing:** first study presentation predicts fourth study
   presentation. This is the strong local-order positive control.
2. **Blank-screen imagery:** fourth study presentation predicts gaze during
   the subsequent stimulus-free visualization interval. This is the partial
   replay, contraction, and background-gaze test.

Only held-out candidate performance is evidence. Alignment components and
warp parameters are diagnostics.

## Data and provenance

The source files are the repository-deposited
`test_data/study_fixations_all.csv` and
`test_data/testdelay_fixations.csv`. They accompany Wynn, Ryan, and
Buchsbaum (2020), *Eye movements support behavioral pattern completion*,
PNAS 117(11), 6246-6254, DOI 10.1073/pnas.1917586117, and the young/older
adult extension Wynn, Buchsbaum, and Ryan (2021), *Encoding and retrieval eye
movements mediate age differences in pattern completion*, Cognition 214,
104746, DOI 10.1016/j.cognition.2021.104746.

The experiment used an EyeLink II at 500 Hz with nine-point calibration.
EyeLink classified saccades above 0.5 degrees and missing signal lasting at
least three samples; remaining samples were fixations. Images and the blank
post-test field occupied 800 by 600 pixels within the 1024 by 768 display.
Study images appeared for 3000 ms and were repeated four times. A degraded or
intact test image was followed by a 50-ms mask and a 3000-ms blank field during
which participants visualized the image.

The PNAS article deposits the analysis code and data in this repository. The
repository is MIT licensed; no separate data-only licence is recorded, which
is retained as a provenance limitation rather than silently inferred.

## Preprocessing

- Keep subjects represented in both source files and rows whose image is not
  the missing sentinel `.`.
- Retain positive-duration fixations longer than 80 ms, following the
  repository's earlier temporal analysis.
- Retain coordinates in the declared image rectangle: x in [112, 912] and y
  in [84, 684], then translate its upper-left corner to (0, 0).
- Derive study repetition by ranking the four study trial numbers within each
  participant and image. Keep participant-item pairs having all four
  presentations.
- For the primary blank-screen path, keep fixation onset `FixOffset` in
  [0, 3000), truncate duration at 3000 ms, and reapply the 80-ms threshold.
- Use ordinal onset solely to preserve event order. Duration remains the gaze
  mass. Coordinates and all spatial parameters are in pixels because viewing
  distance is not recorded in the deposited tables.
- Require at least three retained fixations in study repetitions 1 and 4 and
  in the primary blank-screen path so every comparator receives the same rows.

The frozen validation sample contains eight participants, balanced four young
and four older adults, and eight images. Eligibility depends only on phase
presence and the preprocessing requirements above. Participants are sampled
with the frozen seed within age group from those with at least 100 eligible
items; images are then sampled from items eligible for every selected
participant. The script writes the resulting identifiers before any model is
fit.

The data-only eligibility pass selected participants 104, 121, 124, 128, 315,
317, 319, and 327, and images 5, 18, 69, 73, 81, 82, 88, and 105. There were
43 eligible young participants, 25 eligible older participants, and 82 images
complete across the selected participants. These identifiers are now frozen.

## Two-way outer cross-fitting

Participants are assigned to two folds within age group and images to two
folds, using the frozen seed. The four Cartesian participant-fold by item-fold
cells are the evaluation folds. For each cell:

- evaluation rows are the intersection of held-out participants and held-out
  items;
- training rows exclude every held-out participant and every held-out item;
- each evaluation path is contrasted only with the four held-out item
  templates from the same participant;
- only unaltered signal rows enter training;
- the same folds and candidate sets are used by every method.

Thus the scored participant and item identities are absent from warp,
background, emission, transition, coverage, temperature, ridge-weight, and
ridge-penalty fitting. The reference encoding paths are observed templates,
not learned nuisance parameters.

## Frozen model families

- Screen: 800 by 600 px, centre (400, 300).
- Registration: global isotropic contraction plus translation, learned only
  in the two-way outer training set and applied identically to all candidates.
- Transport v2: Gaussian mixture bandwidths 30, 60, and 120 px; two ordinal
  successors; coverage 0.5, 0.75, or 1; omission-penalty candidates 0.25,
  0.75, and 1.5; one deterministic start; entropy continuation 0.05 then
  0.015; 80 iterations; four candidate workers where available.
- Replay: 48 duration-mass bins; maximum local skip two; Student-t degrees of
  freedom four; spatial scale floor 6 px; the same frozen transition grid as
  the final synthetic court.
- Replay likelihood temperature: two inner item folds within each
  outer-training set, with a fixed unit-scale Gaussian prior on log temperature;
  score normalization and information are evaluated in log space.
- Density: fixed 30, 60, and 120 px maps on a 24 by 24 grid.
- MultiMatch: all six eyesim dimensions, including position EMD.
- Supervised baselines: conditional ridge composites with lambda in 0.01,
  0.1, 1, or 10, selected using training-only two-fold item splits.
- Elastic consensus: 80-px consensus, 400-px rigidity, and 40-px matching
  radii; at most 25 iterations. These are pixel-scale sensitivity parameters,
  not claimed degree conversions.

## Negative controls

Every task is evaluated with four source conditions:

- `signal`: the observed source path;
- `shuffled_order`: a frozen random permutation of fixation events, retaining
  the exact position-duration multiset;
- `wrong_item`: a different selected item path from the same participant,
  presented under the nominal item's label;
- `generic_gaze`: a compact path reproducing the participant's pooled spatial
  mean and covariance but containing no item label information.

In addition, every candidate denominator consists of within-participant wrong
items. The generic path is constructed only for negative-control evaluation
and never enters parameter fitting.

## Endpoints and uncertainty

Primary method summaries are held-out item information gain, candidate log
loss, top-one identification credit, mean rank, Brier score, and five-bin ECE.
Uncertainty is a frozen 1000-draw crossed participant-by-item bootstrap. The
court also reports:

- power at the empirical 95th percentile of the wrong-item and generic-gaze
  null information, with fractional tie credit;
- false-positive rate at that same threshold;
- signal-minus-shuffled information for the local-order check;
- contraction, translation, replay coverage, background coverage, spatial
  residuals, convergence, elapsed time, and object size;
- results by task, age group, old/lure probe, degradation, and cue duration
  where cell sizes permit;
- a Replay-only sensitivity analysis using every deposited post-test fixation
  rather than the primary 0-3000 ms window.

## Decision rules

The real-data gate passes only if all of the following hold:

1. Replay and Transport each have positive mean signal information and
   top-one credit above 1/4 in repeated viewing.
2. At least one GazeWeave engine has positive mean signal information and
   top-one credit above 1/4 during blank-screen imagery.
3. At least one GazeWeave engine loses information after order shuffling in
   repeated viewing, while registered density is invariant up to numerical
   tolerance.
4. Each engine's pooled empirical false-positive rate is at most 0.075 and
   registration does not increase null top-one credit by more than 0.05 over
   its raw counterpart.
5. Candidate sets, fold audits, and fitting provenance show no participant or
   item overlap between training and evaluation; at least 99% of engine
   candidate sets converge.

Replay retains the directional default only if the gate passes, its pooled
signal log loss is no more than 0.05 nats worse than Transport, its top-one
credit is no more than 0.05 worse, and it is not materially worse than the
best fair baseline (no more than 0.10 nats log-loss loss and no more than 0.10
top-one loss). Otherwise the result is `no_real_data_default`. A superiority
claim requires a crossed-bootstrap 95% interval excluding zero in GazeWeave's
favour for both log loss and information against the strongest fair baseline;
passing the engine-selection gate alone is not superiority evidence.

The blank-screen task is expected to be substantially harder than repeated
viewing. Failure is publishable evidence against the current formulation, not
permission to alter the frozen court. Any post-result change must be labelled
as a new protocol version or an exploratory addendum.

## Pre-final pilot addendum A: distinguish a sequence ablation from a null

A one-fold engineering run was used only to verify calibration and output
construction. Before the final court, the operating-characteristic wording
was narrowed so that `shuffled_order` is not included in the false-positive
null. A shuffled path retains the correct item's complete spatial occupancy
and therefore is deliberately positive for density; treating it as an
item-identity null would make the density comparator fail by construction.
Matched-FPR thresholds use only `wrong_item` and `generic_gaze`.
`shuffled_order` remains the predeclared sequence ablation and is evaluated
only through the paired signal-minus-shuffled check. No final outcome had been
computed when this clarification was made.

## Implementation addendum B: inner-cross-fitted Replay calibration

After the first real court, regeneration of the unchanged synthetic court
exposed a null-calibration failure in the newly added Replay temperature fit.
The temperature had used in-sample outer-training candidate scores. Synthetic
addendum C replaced these with two inner item folds, added a fixed unit-scale
Gaussian prior on log temperature, and moved information/log-loss evaluation
to normalized log probabilities.

The complete real court was then rerun from raw data with seed `20260817`, the
same cohort, folds, candidate pools, controls, baselines, and gate thresholds.
The new calibration caused the repeated-viewing gate to fail and yielded
`no_real_data_default`. No further temperature adjustment is permitted against
this court; such work requires a new protocol or independent data.
