# GazeWeave probe-delay persistence court

Status: frozen before outcome scoring
Protocol version: 2
Random seed: 20260820
Mote: `bd-01M0FFZNNG9PGWQ50KPZEGP4XA`

## Purpose and evidential boundary

This provisional court asks whether item-specific gaze structure evoked by a
brief, possibly degraded retrieval probe remains detectable after the probe
has disappeared. It is a same-parent-study extension of the existing Wynn
real-data court, not an independent replication. It cannot close the external
Wang-data gate, select a universal default engine, or support a superiority
claim by itself.

The directional estimand is:

> How much held-out information does retrieval gaze provide about the identity
> of the participant's corresponding studied image, during the visible probe
> and after the probe has disappeared?

The primary endpoint remains held-out `gaze_info_bits = log2(p_true / prior)`.
Warp, replay coverage, spatial residual, chronology, rank, and calibration are
diagnostics.

## Data and provenance

The frozen inputs are:

- `test_data/wynn_probe_delay/study_fix_input_new.csv`, SHA-256
  `ae4f889074a856376085e5020a4ce20444bd56ae6d471e9606bc404b1f9240a9`;
- `test_data/wynn_probe_delay/testdelay_fix_input_matched.csv`, SHA-256
  `32caac0a991972175f929aeb9f8b474442c60e9bfb04177d55c99239e8e0f244`.

The fixation-level CSVs are local-only, explicitly ignored by Git, and must
not be committed or published. The manifest, protocol, and aggregate court
outputs contain no fixation-level rows.

They contain 46 pseudonymous participants, 120 base image identities, two
image versions, five saliency levels (20, 40, 60, 80, and 100), four intended
study presentations, and retrieval conditions `old`, `lure`, and `newtest`.
For `old`, the retrieval image version matches study. For `lure`, retrieval
uses the alternate version of the same base image. `newtest` has no studied
template and is reserved for negative-control analysis.

The experimental image rectangle is x = 112--912 and y = 84--684, translated
to an 800-by-600 analysis screen. The retrieval sequence contains a 500-ms
probe followed by a 2,500-ms delay. No separate data-only licence is inferred.

## Import and phase construction

- Verify filenames, row/column counts, schemas, and SHA-256 checksums before
  every non-smoke run.
- Keep rows with finite coordinates and timing, then intersect fixation
  intervals with the declared spatial and temporal windows. Translate retained
  coordinates by (-112, -84).
- Study paths use the half-open interval [0, 2500) ms. The 2500-ms boundary is
  supported by the supplied relative timing field and is used only as a
  deterministic truncation rule.
- Order the four study presentations within participant and base image by
  `Run`, then `TrialTotal`, then `Trial`. Require four distinct presentations.
- Retrieval fixation intervals are defined by `FixProbeOnset` and
  `FixDuration`. Split boundary-straddling fixations rather than assigning all
  their mass to one phase.
- Define retrieval tasks:
  - `probe`: [0, 500) ms;
  - `delay`: [500, 3000) ms, with onset reset to zero;
  - `combined`: [0, 3000) ms;
  - `early_delay`: [500, 1500) ms, diagnostic only;
  - `late_delay`: [1500, 3000) ms, diagnostic only.
- After clipping, retain positive-duration fixations longer than 80 ms and
  preserve physical within-window onset. Duration is gaze mass; default
  chronology is ordinal local order.
- The primary reference path is study presentation four. Presentation one to
  presentation four is the strong repeated-viewing positive control.
- Match on participant plus base `ImageNumber`, not `ImageVersion`: this makes
  lure trials a declared pattern-completion condition rather than a label
  error.
- Never filter the primary court by response accuracy, saliency, probe type,
  fitted score, or model convergence.

## Eligibility and frozen cohort rule

Eligibility is determined without examining any similarity score or response
accuracy. A signal participant-item pair must:

- be `old` or `lure` at retrieval;
- have four study presentations, each with at least three retained fixations;
- have at least one retained probe fixation and at least three retained delay
  and combined fixations;
- have at least one retained fixation in each early- and late-delay diagnostic;
- have a unique study and retrieval mapping for participant plus base image.

The initially declared common-item rectangle was infeasible: among all sets of
eight complete, fixation-eligible participants, the largest intersection was
five items. This was discovered by a data-only eligibility audit before any
similarity model was fit. Protocol version 2 therefore preserves eight
participants and ten trials per participant without requiring identical item
sets.

Assign all eligible base-image identities globally to two item folds by
shuffling with seed 20260820. Participants are grouped by their complete set of
old/lure images (four counterbalancing groups). Within each group, sample two
complete participants who have at least five eligible items in both item
folds, using seed 20260821. Within each selected participant and item fold,
sample five eligible items using seed 20260822. This yields 80 paths, five
candidate items per held-out participant and item fold, and preserves global
item-level separation even though participants need not contribute the same
items.

The frozen participant and item-fold assignments are:

| Participant | Fold 1 items | Fold 2 items |
|---|---|---|
| 1015 | 6, 48, 62, 95, 98 | 17, 25, 68, 71, 116 |
| 1018 | 5, 62, 66, 69, 98 | 27, 51, 68, 106, 117 |
| 1019 | 53, 72, 77, 119, 120 | 9, 39, 43, 78, 97 |
| 1026 | 13, 24, 40, 67, 100 | 17, 18, 27, 59, 108 |
| 1038 | 29, 79, 94, 96, 105 | 35, 83, 88, 89, 112 |
| 1039 | 13, 20, 69, 90, 99 | 27, 35, 75, 84, 112 |
| 1042 | 6, 29, 62, 73, 100 | 22, 41, 93, 106, 111 |
| 1044 | 7, 29, 95, 101, 114 | 8, 65, 85, 93, 112 |

The data-only balance is 48 old and 32 lure trials; saliency counts at 20, 40,
60, 80, and 100 are 13, 21, 14, 12, and 20. There are 22 complete participants,
18 meeting the two-fold eligibility rule, and two selected participants from
each counterbalancing group.

## Two-way cross-fitting and candidates

Within each counterbalancing group, assign its two selected participants to
opposite participant folds using seed 20260823. Cross these participant folds
with the frozen global item folds to form four outer evaluation folds. Each
training set excludes every held-out participant and every identity assigned
to the held-out item fold. Each evaluation path is contrasted against five
held-out candidate items from the same participant. The same folds, candidate
sets, and unaltered signal rows are used by every applicable method and phase.

## Frozen engines and baselines

Reuse the model families and numerical settings from protocol version 1 of
`GAZEWEAVE-REAL-CONTROLS.md` without outcome-driven tuning:

- candidate-invariant isotropic contraction plus translation;
- Transport v2 with 30, 60, and 120 pixel spatial scales, two ordinal
  successors, and coverage 0.5, 0.75, or 1;
- Replay with 48 duration-mass bins, maximum local skip two, Student-t
  emissions, restartable replay, and background gaze;
- registered multiscale density and elastic consensus;
- nested registered MultiMatch ridge where both paths contain at least three
  fixations.

MultiMatch is structurally undefined for most 500-ms probe paths and for some
early- and late-delay paths. It is therefore excluded from those phase courts
rather than imputed. Its path coverage is reported. Delay, combined, and
repeated-viewing courts require at least three fixations and include the full
baseline set.

## Controls and outcomes

For repeated viewing, probe, delay, and combined tasks, score observed signal,
frozen order shuffle, a wrong item from the same participant and item fold, and
participant-generic gaze. Early- and late-delay tasks are predeclared
signal-only persistence diagnostics; their null calibration comes from the
full delay task. Order shuffle is a chronology ablation, not an item-identity
null. Wrong item and generic gaze define the operating null.

Report by method and task:

- gaze information bits, log loss, top-one credit, mean rank, Brier score,
  and five-bin ECE;
- 1000-draw crossed participant-by-item bootstrap intervals;
- empirical false-positive rate at the pooled 95th-percentile null threshold;
- contraction, translation, replay/background coverage, convergence, runtime,
  and fold-overlap receipts;
- old versus lure and saliency-level summaries;
- early- versus late-delay information as a persistence diagnostic.

## Interpretation gate

- `supported`: at least one predeclared GazeWeave engine has a crossed-bootstrap
  interval above zero for both delay and late delay, while its pooled null
  false-positive rate is at most 0.075 and fold overlap is zero.
- `suggestive`: mean information is positive for delay and late delay but at
  least one interval includes zero.
- `not supported`: delay or late-delay mean information is nonpositive, or
  required leakage/null checks fail.

Probe and combined results explain when evidence appears but cannot rescue a
failed delay-persistence gate. Density or another baseline may be strongest;
that result must be reported. No post-result parameter changes are permitted
under protocol version 2.
