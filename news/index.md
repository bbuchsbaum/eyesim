# Changelog

## eyesim 0.2.0.9000

- Plot titles, subtitles and captions wrap to the width they are drawn
  at, in lines of even length (new
  [`element_text_wrap()`](https://bbuchsbaum.github.io/eyesim/reference/element_text_wrap.md),
  used by
  [`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)),
  so long GazeWeave titles no longer run off narrow figures. Where a
  layout (such as side-by-side patchwork) reserves less height than the
  wrapped text needs, the text is set smaller to fit rather than
  overprinting. Newlines in a label are kept as line breaks. A user’s
  `theme(plot.title = element_text(...))` merges into it and keeps
  wrapping. The package’s long GazeWeave labels also carry line breaks,
  so they stay readable under other themes.

- [`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md)
  numbers and joins fixations in onset order, even when rows are not
  sorted by onset.

- Fixations whose markers overlap, and so cannot be numbered
  individually, are outlined together (a ring for a single fixation) and
  labelled once (“4, 9, 17”), with a leader that stays clear of every
  other fixation. Numbers still without a position are listed in the
  in-panel note (or counted, when more than six).

- With a finite upper `limit`, colour bars label the limit “≥ limit” (to
  three significant digits), and density at or above it takes the top
  colour. The scale is the same whether or not a plot’s data exceed the
  limit, so plots sharing `limits` share colours and one collected
  colour bar.

## eyesim 0.2.0

- All plots share one visual system. New exports:
  [`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md),
  [`theme_eyesim_spatial()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md),
  [`eyesim_colours()`](https://bbuchsbaum.github.io/eyesim/reference/eyesim_colours.md),
  [`scale_fill_eyesim_density()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md)
  and
  [`scale_colour_eyesim_time()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md).
  Reference/encoding and source/recall gaze, correspondence,
  diagnostics, time and density each keep one colour or ramp across
  [`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md),
  [`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md),
  [`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md)
  and the GazeWeave Replay and Transport
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  methods. Visible changes:

  - `bg_image` is drawn with its stored colours. It was
    contrast-stretched.
  - Stimulus-space plots keep a 1:1 aspect ratio, and `xlim`/`ylim` no
    longer drop fixations outside them, so density contours close.
  - Density is transparent at zero and has a colour bar.
    [`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md)
    and
    [`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md)
    gain `limits` (put several maps on one scale) and `legend`;
    [`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md)
    gains `colours`.
  - `plot.fixation_group(transform = )` now transforms the density
    colour scale, as documented. The fill was always cube-root before,
    so the default (`"identity"`) now fills linearly and low densities
    look lighter. `colours` is now used; it was ignored. `alpha_range`
    defaults to `c(0, 0.9)`. Point area is proportional to duration.
  - `plot.eye_density(alpha = )` is the opacity at maximum density;
    lower densities are more transparent.
  - [`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md)
    titles frames by what it animates over (“Fixation N” or “Onset N”)
    and keeps earlier fixations as faded marks.
  - `type = "density"` and `type = "filled_contour"` draw
    non-overlapping density bands between evenly spaced edges from zero,
    coloured on the shared density scale (so they honour `transform`,
    `colours` and `limits`). The band containing zero is transparent,
    and a stepped colour bar shows the band edges. Densities above
    `limits` take the top colour. Over a `bg_image`, thin white contour
    lines keep density edges visible.
  - Fixation numbers are placed when the plot is drawn, at its actual
    size: next to their own fixation, or further out with a leader line
    that stays clear of every other fixation and leader. Numbers are
    never drawn over another number or marker. The first and last
    fixations are placed first, and a ring marks the first fixation (its
    leader starts at the ring). Numbers that cannot be placed
    unambiguously are omitted, and a note inside the panel (unaffected
    by `labs(caption = )`) says how many. Rank colour bars show no tick
    labels.
  - [`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md)
    gains `aspect = c("equal", "free")`: `"free"` fills the panel at the
    cost of true angles and shapes. With shared `limits` below the peak,
    the stepped colour bar gains a “\> limit” step, and a lower limit is
    ignored (with a warning) for the banded types.
  - GazeWeave overlays have a legend for the reference/source paths.
    Braid widths, arrow opacity and point areas are scaled from zero
    (widths and opacity keep a small visible floor). In both braids,
    point area is fixation mass and source points are filled by the
    share of their mass with a correspondence (faint = unmatched, or
    background in Replay).
  - Requires ggplot2 \>= 4.0.0 and imports viridisLite.

- [`fixation_entropy()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.md)
  on a fixation group accepts `aggregate`. With a vector `sigma`,
  `aggregate = "none"` returns one entropy per scale. It was passed on
  to [`ks::kde()`](https://mvstat.net/ks/reference/kde.html), which
  rejected it.

- `sample_fixations(fast = FALSE)` now returns the same object as the
  default path: a `data.frame` with `duration = 1`, and all four columns
  for an empty `time`. It returned a tibble with `duration` NA before
  the first onset and 0 after, and a column-less tibble for an empty
  `time`.

- [`gaze_transport_spec()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_spec.md)
  gains `revision = c("2026.10", "2026.08")`. The default, `"2026.10"`,
  changes the Transport solver in both the native and the reference
  backend. `"2026.08"` reproduces the frozen August 2026 solver bit for
  bit, and the frozen validation courts now pin it. Specifications saved
  before revisions existed are solved as `"2026.08"`; an unknown
  revision is an error. Under `"2026.10"`:

  - **Estimand change (unreleased).** The chronology residual is
    `1 - 2A / (R + S + kappa)` with `kappa = 1e-3`, one constant shared
    by both backends. It was `1 - 2A / (R + S)`, set to 0 when
    `R + S <= 1e-15`. That definition jumped by up to 0.1 in the
    objective and was 0/0 for one-fixation sources. The new residual is
    continuous, has a bounded gradient, and tends to 1 as selected edge
    mass vanishes. The limit 1 is the neutral value: a source without
    edges gets the same chronology term from every candidate.
  - **Stopping rule (first-order).** `tolerance` is now in objective
    units (nats), with default `1e-6`. A stage stops when three
    conditions hold: the plan is feasible, no negligible-mass cell is
    trapped, and a first-order model predicts at most `tolerance` of
    remaining decrease. The model sums entropic single-cell relaxation
    gains on real cells and a quadratic model on slack cells, both from
    least-squares dual prices. This is a stopping heuristic, not an
    optimality certificate. On a 60-pair review set, a tight
    continuation of stopped nodes lowered the objective by at most
    9.7e-7 at the 95th percentile. On 9 of 692 nodes it still found 1e-5
    to 1.4e-3, along slow directions or away from saddles. Tolerances of
    1e-7 to 1e-9 did not remove those nodes and multiplied stalls.
  - **Traps and steps.** A reviewer showed that zeroing a support cell
    could produce plans that a weighted measure alone certified, 2e-6 to
    25 above the node optimum. Such cells, including cells at exactly
    zero, are now detected from their uncapped relaxed mass and reseeded
    (a 1e-3 mixture with the independent start). All 40 review traps
    then re-converged within 4e-7 or stayed uncertified. Of 332
    exact-zero variants, none certify in place. After reseeding, 24
    settle in a different local optimum 1e-4 to 0.14 above the original.
    Every mirror step is capped at `min(step_size, 50 / max|g|)`, so the
    exponent clamp never binds.
  - **Line searches, stalls and `maxit`.** Both revisions accept a step
    that raises the objective by up to 1e-12 relative. A line search
    therefore fails only when every trial rose by more than that, which
    at 1e-8 projection accuracy is noise-dominated. The native backend
    reported this as a numerical failure, and `backend = "auto"` then
    silently re-solved the pair with the 30-200x slower reference
    backend. A noise-limited line search now ends the stage as
    `"stalled_projection_limited"`. Reaching `maxit` ends the node as
    `"not_converged"`, including under `backend = "optimized"`; it
    previously raised an error on 14 of 60 review pairs at
    `step_size = 0.05`. Both statuses are scored, recorded, and counted
    by
    [`gaze_transport_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_cv.md).
  - **Projection.** When standard Sinkhorn exhausts its iterations, a
    damped log-domain dual Newton projection takes over, then log-domain
    Sinkhorn, in both backends.
  - **Starts and entropy schedule.** Every coverage node is solved from
    the independent start and the adjacent-coverage continuation,
    replacing the reference backend’s silent cold restart.
    `multistart = 2` adds a spatial start and is now native. The default
    entropy schedule is `c(0.15, 0.05, 0.015)`.
  - **Known residual error: local optima.** The objective is non-convex,
    and the default score can fall short of a multistart oracle. The
    error is one-sided (default scores are lower than the oracle’s) and
    depends on the pair mix:
    - On our 60-pair review set, the gap was at most 0.12 nats (0.39
      null-score SD), with 95th percentile 0.0056.
    - On a second reviewer’s 88-pair set (half null), the gap was at
      most 0.117 nats (0.24 SD), with 95th percentile 0.023 (0.047 SD).
      The 95th percentile was 0.033 for null pairs and 0.006 for matched
      pairs. The error is concentrated in null and sparse-source pairs.
      `multistart = 2` never moved a score by more than 2e-5 on the
      reviewer’s set; on ours it lowered the 95th percentile only to
      0.0053 at 1.48 times the runtime. It is therefore opt-in. Scores
      also still depend on `step_size` at this level.
  - The mutual-information term is evaluated in the log domain when the
    product of two marginals underflows.
  - **Backends and fallbacks.** Under `backend = "auto"`, specifications
    the native backend cannot solve (`projection_method = "log"`) go to
    the reference backend for every pair. Any per-pair numerical
    fallback raises a warning of class
    `gaze_transport_backend_fallback`. Every alignment records its
    backend, fallback, status and solver revision in `convergence`.
    [`gaze_transport_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_cv.md)
    reports fallback, stall and not-converged counts per row and in
    `solver`. Validation checkpoints record and verify the solver
    revision.
  - **Polish gap.** The polish Frank-Wolfe gap is `NA`, with a reason,
    where a selected mass with a positive target vanishes. The
    Jensen-Shannon gradient is unbounded there, so the reported gap came
    from the gradient floor. The objective is non-convex, so the gap
    measures stationarity only.

- GazeWeave Replay gains
  `gaze_replay_spec(revision = c("2026.10", "2026.08"))`. The new
  default, `"2026.10"`, changes Replay as follows:

  - The observation unit is one recall fixation. The HMM makes exactly
    one hidden-state visit per recall fixation (background, or one state
    per encoding fixation, with the same transition structure as
    before), and each fixation emits its position and its duration.
    Durations are log-normal: replay states have log-duration mean
    `a + b * log(encoding duration)`, and the background has its own
    log-normal. Both are normalised densities whose parameters are
    shared by every candidate. The candidate score is the total trial
    log likelihood, not a mean per duration bin, so evidence grows with
    the number of recall fixations. `grid_size` is ignored under
    `"2026.10"` (a one-time message says so when it is supplied).
  - When a `screen` is declared, each replay emission is a Student
    density truncated to the screen and normalised by its on-screen
    mass. Under `"2026.08"` the on-screen mass was 0.99 at the centre
    and 0.35 near a corner. Encoding fixations recorded off the screen
    are projected onto it, both in the emissions and in the EM scale
    update. Without a screen, the support is the whole plane, and a
    one-time message says so.
  - The background state’s position density is a screen-truncated kernel
    density of recall fixations (equal weight per fixation within a
    trial), mixed with a uniform screen floor whose weight EM fits.
    `"2026.08"` used a single Student density at the training recall
    mean. The evaluation protocol is declared by `background_support`:
    - Under `"training"` (the default), backgrounds come from training
      recalls only. The `background_by` level is used when it has at
      least three trials, otherwise the population density.
    - Under `"held_out"` (opt-in), a held-out participant’s own recalls
      of items outside the scored row’s candidate pool are used. It
      engages only with an item-level split and a `contrast_on` that
      blocks items; it never engages under participant holdout with
      shared items. It is transductive: it reads other evaluation rows’
      item labels, so it cannot be used on unlabelled recalls.
    - At scoring time the background drops training rows whose exclusion
      key (`setdiff(match_on, background_by)`) matches any candidate in
      the pool, never only the true one, so it is independent of the
      label. When `background_by` is not in `match_on`, that key
      includes the participant, so other participants’ recalls of the
      candidate items remain. If the exclusion would leave fewer than
      three trials, nothing is excluded. An earlier draft excluded only
      the true item, which leaked the label: on an ambiguous null, top-1
      accuracy with two candidates rose from 0.48 to 0.90.
    - Pooled fallbacks are reported by a message and counted in
      [`gaze_replay_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_cv.md)
      provenance (`background_fallback`), as are background levels
      unseen in training (e.g. under participant holdout).
  - The replay scale, the duration parameters, and all four transition
    probabilities are fitted by Baum–Welch EM on training rows only.
    `transition_grid` now supplies only the starting values. The
    background entry, initial and persistence probabilities may all
    approach one, so a fully null recall can be represented. The replay
    scale is capped at one sixth of the shorter screen side: a replay
    state is local, and a near-uniform replay state cannot be
    distinguished from the background. EM reports non-convergence if the
    likelihood ever decreases.
  - On data simulated from the fixation-level model, EM recovers the
    replay scale to within 1%, the background share to within 0.02, and
    the duration slope `b` to within 0.01. `"2026.08"` inflated the
    scale by 60%. With the validation-court transition grid, it also
    could not exceed a stationary background share of 0.67 under a full
    null.
  - An earlier, unreleased draft of `"2026.10"` kept the duration-grid
    observation model. It treated repeated bins of one fixation as
    independent emissions, so replay and background were not identified:
    on a uniform null with 20-30 fixations at 64 bins the background
    share was 0.44-0.62. The fixation-level model fits the same null
    with a background share of 0.98, and the candidate score spread is
    about 1% of a typical true-candidate margin. The draft’s
    repeated-grid message and `provenance$repeated_grid_rows` are
    removed, because `"2026.10"` is unreleased; its results are not
    comparable with this version.
  - Fitted models record `revision`, EM diagnostics (`training$em`), and
    the background model. Their `version` is 5.
  - [`gaze_replay_align()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_align.md)
    and
    [`gaze_replay_align_episode()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_align_episode.md)
    gain `background_key`. Revision `"2026.08"` reproduces the earlier
    fits and scores bit for bit. The frozen validation court scripts pin
    it.
  - Specs and models saved without a revision are treated as
    `"2026.08"`. An unknown revision is an error. Multi-presentation
    court checkpoints record the Replay revision.

- Shared candidate calibration for Transport and Replay under revision
  `"2026.10"` (unreleased), configured by the new
  [`gaze_calibration_control()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_calibration_control.md)
  and the `calibration_control` argument of
  [`gaze_transport_spec()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_spec.md)
  and
  [`gaze_replay_spec()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_spec.md).
  Revision `"2026.08"` is unchanged: its specifications carry no
  control, and its cross-fitted results are
  [`identical()`](https://rdrr.io/r/base/identical.html) to those before
  this change (pinned by
  `tests/testthat/fixtures/calibration_base_0d9fbcf.rds`, fitted at
  master 0d9fbcf). Replay fitted without `contrast_on` still scores at
  temperature one, as before; with too few inner candidates it now
  returns the declared prior instead of temperature one.

  - **Evidence-scaled temperature.** Each held-out row is scored at
    `T_i = T * (n_ref / n_i)^gamma`, where `n_i` counts the row’s
    evidence (Transport: duration-effective recall fixations; Replay:
    recall fixations) and `n_ref` is their geometric mean over the
    calibration rows. `T` and `gamma` in `[0, 1]` are fitted jointly on
    inner out-of-fold rows only. `gamma = 0` is one global temperature.
    The engine’s lower temperature bound caps every row’s confidence.
    Results gain `temperature` (per row), `inverse_temperature` and
    `evidence_count`; `fit$calibration$fold_calibration` records each
    fold’s fit.
  - **Prior.** The prior on `log T` centred at `T = 1` is replaced by a
    zero-centred Gaussian prior (sd 3) on the standardised inverse
    temperature `sigma / T`, where `sigma` is the median within-row
    standard deviation of the inner scores. It pulls only toward the
    declared candidate prior. The old prior pulled Replay’s total trial
    log likelihood, whose natural scale is several nats, toward
    overconfidence.
  - **Stein shrinkage.** The fitted inverse temperature is multiplied by
    `max(0, 1 - p / LR)`, with `LR` the inner likelihood-ratio statistic
    against the declared prior and `p` the number of free parameters.
    With a few dozen inner rows a null fit otherwise returns a positive
    inverse temperature about half the time.
  - **Ranking is calibration-independent.** Calibration changes
    probabilities and bits only. Under `"2026.10"`, `template_rank` and
    `top1_credit` (candidate and per-episode) come from one fixed,
    label-free reference ranking for both engines and both methods:
    `log prior_k + r_k`, with `r_k` the native candidate score after any
    typicality offset (Transport: log mean of the episode scores at
    inverse temperature one; Replay: trial log likelihood), whatever the
    fitted `T`, `gamma` or Stein factor. Previously they followed the
    calibrated posterior. That is not calibration-invariant for
    multi-episode Transport, whose inverse temperature acts inside the
    log mean over episodes: episodes a = (0, -10) and b = (-3, -3) rank
    b first as `1/T -> 0` (arithmetic mean) but a first at `T = 1`. The
    earlier rule for a calibration that returns the prior (prior first,
    then the ranking score) was therefore not the `1/T -> 0` limit, and
    the claim that calibration never changes the ranking was false for
    multi-episode Transport. A test pins identical top-1 and rank AUC
    across `"evidence_scaled"`, `"global"`, an unshrunk fit and a Stein
    factor of zero on a two-episode Transport CV. AUC should be computed
    from ranks or `ranking_score`, not from the calibrated
    probabilities.
  - **Kappa.** The effective-fixation reliability shrink is dropped from
    the default: it could only temper sparse rows, and with `gamma` it
    was a second, weakly identified parameter on the same axis.
    Transport’s `reliability` now defaults to `"none"` under
    `"2026.10"`, and `"effective_fixations"` requires
    `gaze_calibration_control(method = "global")`. That method
    reproduces the pre-change `"2026.10"` probabilities and bits exactly
    only with `reliability = "effective_fixations"` passed explicitly
    for Transport (the new default `"none"` fits a different model); a
    saved `"2026.10"` specification without a control is treated the
    same way. Its rank and top-1 follow the reference ranking above and
    so differ from the pre-change values only where the fitted
    temperature or kappa mixture had reordered candidates (multi-episode
    Transport or a non-uniform prior); the pinned fixture,
    single-episode with a uniform prior, is
    [`identical()`](https://rdrr.io/r/base/identical.html).
  - **Typicality offset** (default `"standardized"` for Transport,
    `"none"` for Replay). `typicality = "mean"` subtracts each
    candidate’s mean score against a seeded, item-balanced subsample of
    the fold’s training recalls of other items (inner rows use their
    inner training rows). Scores are centred within each source and the
    means are shrunk by an empirical-Bayes factor. `"standardized"` also
    divides by each candidate’s shrunk score standard deviation over the
    sources. The offsets never read a held-out recall or label, so
    relabelling a held-out row leaves every candidate score unchanged.
    Adjusted scores are ranking scores, not normalised likelihoods. Mean
    subtraction does not equalise score variances; `fit$calibration`
    reports per-candidate dispersion before and after.
  - Simulated evidence (all synthetic; numbers from
    `inst/validation/gaze-weave-calibration-a3.R`, gates asserted by the
    slow checks in
    `tests/testthat/test_gaze_weave_calibration_revision.R` with
    `EYESIM_SLOW_TESTS=true`). “Before” is the pre-change `"2026.10"`
    calibration.
    - Centre-biased nulls, pooled over four seeds: Replay mean bits
      -0.107 (SE 0.034) before, 0.000 after; Transport -0.015 before,
      0.000 after. Top-1 stays within MC error of chance. Every fold’s
      Stein factor was zero, so the declared prior was returned. With
      the new prior but without the Stein factor, Replay stayed at
      -0.10.
    - Transport with 2- and 8-fixation recalls: held-out log loss 0.885
      before, 0.939 with `gamma = 0`, 0.760 with fitted `gamma` (0.72 to
      1 across folds). Mean max-p minus accuracy moved from +0.13 /
      -0.11 (sparse / rich) to -0.03 / -0.02 on one seed and from +0.04
      / -0.27 to -0.06 / -0.05 on the other. Top-1 and rank-based AUC
      are unchanged.
    - Replay with 6- and 30-fixation recalls: fitted `gamma` is 0 to
      0.3, because the total log likelihood already grows with the
      fixation count. Held-out log loss 0.626 before, 0.631 with
      `gamma = 0`, 0.653 fitted: no calibration gain for Replay on this
      fixture.
    - Typicality, centre-biased null, six candidates per pool (two
      central encodings), argmax share per central candidate (chance
      0.167): Transport 0.396 without an offset, 0.128 with `"mean"`,
      0.160 with `"standardized"` (MC error 0.020); Replay 0.296, 0.271
      and 0.266 (MC error 0.015). The offsets equalise the candidates’
      held-out means (Transport between-candidate SD 0.32 to 0.08) but
      not their variances: the per-candidate SD ratio was 9.3 before and
      6.0 after `"standardized"` (11.1 after `"mean"`) for Transport.
      For Replay it exceeded 1900 before and after, but this is not a
      heavy-tail effect: under these nulls the HMM background absorbs
      most recall fixations, and in 57% of rows every candidate lay
      within 0.01 nats of the others. The SD ratio comes from these
      near-degenerate pools, and Replay’s residual central bias from the
      argmax among near-ties (central candidates won 0.54 of the tie
      rows against a 0.37 share of candidates; few rows). Replay top-1
      and AUC should therefore be computed from `ranking_score` with a
      tie tolerance declared in advance (for example 0.01 nats),
      splitting credit among tied candidates.
    - Discrimination on signal data: Transport top-1 0.781 / AUC 0.942
      without an offset, 0.781 / 0.944 with `"standardized"`, 0.771 /
      0.940 with `"mean"`; Replay 0.938 / 0.985 without, 0.951 / 0.988
      with either offset.

- Canonicalized edge-normalized Transport as the sole public Transport
  method. The API is now
  [`gaze_transport_spec()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_spec.md),
  [`gaze_transport_align()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_align.md),
  [`gaze_transport_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_cv.md),
  and related unversioned helpers;
  [`gaze_weave_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_weave_cv.md)
  accepts `engine = "transport"`. The unreleased original and v2
  implementations, exports, S3 classes, tests, and documentation were
  removed. Their labels remain only in internal and immutable
  development-court provenance.

- Added edge-normalized GazeWeave Transport as the symmetric explanatory
  alignment engine. It integrates bounded coverage, scores separate
  study presentations as fixed equal-prior likelihoods, uses nested
  response-blind calibration, and retains a pure-R oracle alongside the
  optimized native backend. Default tidy output contains only
  identifiers and `gaze_info_bits`; coverage, coordinate-unit spatial
  RMSE, local order, contraction, rank, episode weighting, and solver
  stability are nested diagnostics.

- Added Transport registered overlays, optimized-correspondence gaze
  braids, raw-versus-registered paths, and a one-score evidence ledger.
  Labels state that Transport correspondence is an optimized alignment,
  not posterior replay probability.

- Froze and ran the Transport simulation, full-cohort repeated-viewing,
  and exact 0–3000 ms retrieval courts. Invariants, runtime, and the
  study positive control passed. Retrieval information and held-out
  old-item response prediction did not pass their uncertainty gates, so
  Transport is retained as an explanatory engine without a retrieval or
  superiority claim.

- Added GazeWeave Replay, a directional probability model of locally
  ordered, restartable encoding-to-recall replay mixed with
  candidate-independent background gaze. Candidate probability is
  calibrated on inner-held-out items within each outer-training set from
  mean log predictive density per normalized duration bin, with
  log-domain evidence and a weak fixed residual-temperature prior. Raw
  total HMM likelihood remains available as a diagnostic. The primary
  endpoint is `gaze_info_bits = log2(p_true / prior_true)`.

- Added engine-specific registered overlays, gaze braids, and
  diagnostics. Replay braids show posterior correspondence
  probabilities. Transport braids are labelled as regularized optimized
  correspondences, not posterior uncertainty.

- Added frozen synthetic and real-data validation courts under
  `inst/validation`. The real court compared Replay and the then-current
  Transport development estimator with nested-ridge MultiMatch,
  registered density, and elastic-matching baselines under shared
  participant-by-item folds. The synthetic court passed all gates. The
  real court failed its repeated-viewing gate after calibration exposed
  unstable Replay probabilities, so
  [`gaze_weave_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_weave_cv.md)
  requires an explicit engine choice and no GazeWeave superiority is
  claimed. Registered density was strongest overall and blank-screen
  imagery evidence was weak for every method.

- Added a prospectively frozen participant- and item-disjoint
  replication of the resolution-invariant Replay calibration. The
  correction restored positive held-out Replay information in repeated
  viewing, but neither GazeWeave engine had positive mean information in
  the blank-screen control. The package therefore retains explicit
  engine selection and makes no power or superiority claim.

- [`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md),
  [`fixation_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_similarity.md),
  and
  [`scanpath_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/scanpath_similarity.md)
  now return an `n_perm` column when `permutations > 0`. It records the
  number of permuted non-matching comparisons that actually contributed
  to `perm_sim` for each row, which varies with small `permute_on`
  strata or when fewer candidates than requested are available, and is
  `0` when no baseline could be computed (`perm_sim = NA`). Use it to
  exclude rows with too thin a permutation baseline,
  e.g. `dplyr::filter(res, n_perm >= k)`. The column is additive: output
  for `permutations = 0` is unchanged.

### Correctness fixes from the eyes4s parity audit

- [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  now honours an explicit `weights` vector. It was previously ignored,
  so the map was always unweighted unless `duration_weighted = TRUE`.
  Weights must be finite, non-negative, not all zero, and supplied one
  per fixation; they are subset along with `window`, and they take
  precedence over `duration_weighted`. If the weights of the fixations
  left after `window` sum to zero, or `duration_weighted = TRUE` and
  every duration is zero,
  [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  now warns and returns `NULL` under both backends. The ks path
  previously failed with “missing value where TRUE/FALSE needed”, or,
  with `normalize = FALSE`, returned a map of `NaN`.

- [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  now forwards extra named arguments in `...` to
  [`ks::kde()`](https://mvstat.net/ks/reference/kde.html), as
  documented. Previously any extra argument failed with “unused
  argument”. Arguments that
  [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  sets itself (including `eval.points`), and unknown names, are rejected
  with a clear error, as are extra arguments under `kde_pkg = "MASS"`.

- [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  now requires `kde_pkg` to be `"ks"` or `"MASS"`. Any other value, such
  as a mistyped `"ks "`, previously fell back to MASS silently.

- Weighted densities under `kde_pkg = "MASS"` (explicit `weights` or
  `duration_weighted = TRUE`) now work. The internal weighted kernel
  failed with “invalid ‘times’ argument” and
  [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  returned `NULL`.

- [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  and the geometric warp behind
  [`affine_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)
  and
  [`contract_transform()`](https://bbuchsbaum.github.io/eyesim/reference/latent_pca_transform.md)
  now round maps with a fixed `zapsmall(digits = 7)`. They previously
  used `getOption("digits")`, so the same call returned maps differing
  by up to about 5e-9 when a session changed `options(digits)`. Results
  under the default `digits = 7` are unchanged. The documented rounding
  rule was checked against R 4.5.1.

- The permutation baseline in
  [`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md),
  [`template_similarity_cv()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity_cv.md),
  [`fixation_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_similarity.md),
  and
  [`scanpath_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/scanpath_similarity.md)
  now removes the true match before sampling candidates. It previously
  sampled first and then dropped the match, so with three candidates and
  `permutations = 2` most rows received one control instead of the
  documented two. `n_perm` is now
  `min(permutations, available non-matching candidates)`. Whenever
  `permutations` is smaller than the candidate pool, the sampled
  controls, and hence `perm_sim` and `eye_sim_diff`, differ from earlier
  versions for the same seed.

- When several source rows share a `match_on` key, every copy of that
  key is now removed from each such row’s permutation candidates.
  Previously one copy remained, so a row could be compared with its own
  template in the baseline. This matches
  [`sample_density_time()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density_time.md).

- Permutation candidates are now distinct templates. A template matched
  by several source rows previously appeared once per matching row, so
  it could be drawn twice and was counted twice in `perm_sim` and
  `n_perm`. This applies to
  [`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md),
  [`template_similarity_cv()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity_cv.md),
  [`fixation_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_similarity.md),
  [`scanpath_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/scanpath_similarity.md),
  and
  [`sample_density_time()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density_time.md),
  and changes their permutation columns whenever keys repeat within a
  stratum. Reference rows that no source row matches remain outside the
  candidate set, as before; the documentation now says so.

- Documentation:
  [`?template_similarity`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
  no longer claims that permutation sampling uses a “fixed future seed”.
  Sampling uses the session RNG, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  makes the baseline reproducible. No behaviour change.

- [`template_similarity_cv()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity_cv.md)
  no longer resets the caller’s random number stream. Fold assignment
  and permutation draws still run under `seed`, so its results are
  unchanged, but the session RNG state is restored on exit. The
  documentation now states that permutation controls come only from the
  held-out fold.

- [`gaze_weave_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_weave_cv.md),
  [`gaze_transport_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_transport_cv.md),
  [`gaze_replay_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_replay_cv.md),
  and
  [`gaze_baseline_cv()`](https://bbuchsbaum.github.io/eyesim/reference/gaze_baseline_cv.md)
  likewise no longer reset the caller’s random number stream, and no
  longer create `.Random.seed` in a session that had none. Fold
  assignment still runs under `seed`, so folds and results are
  unchanged.

- `similarity(method = "fisherz")` now returns exactly
  `atanh(1 - .Machine$double.eps)` (about 18.37) for every pair of
  identical maps, and for any correlation within
  `64 * .Machine$double.eps` of 1, such as a rescaled copy of a map;
  correlations that close to -1 give its negative. Identical constant
  maps previously returned 1, and identical non-constant maps returned
  values between about 17.3 and 18.37, depending on rounding error in
  [`cor()`](https://rdrr.io/r/stats/cor.html). The clamp is documented
  in
  [`?template_similarity`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md).

- [`similarity()`](https://bbuchsbaum.github.io/eyesim/reference/similarity.md)
  on two density maps now refuses maps whose lattices (x and y grid
  coordinates) differ, instead of comparing the `z` matrices cell by
  cell as if they were aligned. A numeric `y` must have one value per
  grid cell.
  [`template_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md)
  and related wrappers raise the same error on both the fast cosine and
  the general path, and for multiscale maps, where the per-scale error
  handling previously would have turned it into `NA`.
  [`repetitive_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/repetitive_similarity.md)
  also raises it rather than returning `NaN` with a warning per pair.
  Because
  [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  defaults `xbounds` and `ybounds` to each group’s own data range, maps
  built with default bounds generally lie on different lattices and are
  now refused; pass common `xbounds`, `ybounds`, and `outdim` to every
  map that will be compared.

- [`similarity()`](https://bbuchsbaum.github.io/eyesim/reference/similarity.md)
  for fixation groups with `method = "overlap"`, and therefore
  `fixation_similarity(method = "overlap")`, now defaults to
  `dthresh = 60`, the value documented for
  [`fixation_overlap()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_overlap.md).
  It previously defaulted to 40, so the same pair scored differently
  through the two entry points.

- [`fixation_overlap()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_overlap.md)
  now builds its default `time_samples` grid from the onsets of both
  fixation groups, `seq(0, max(c(x$onset, y$onset)), by = 20)`. It
  previously used `x` only, so swapping the arguments could change the
  result substantially (0.098 versus 0.833 in one case).

- [`sample_fixations()`](https://bbuchsbaum.github.io/eyesim/reference/sample_fixations.md)
  with the default `fast = TRUE` now holds the last fixation for time
  points after its onset, as `fast = FALSE` always did. It previously
  returned `NA` there. It also no longer fails for a group with a single
  fixation. Downstream,
  [`sample_density()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density.md)
  with `times`,
  [`sample_density_time()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density_time.md),
  [`template_sample()`](https://bbuchsbaum.github.io/eyesim/reference/template_sample.md)
  with `time`, and
  [`fixation_overlap()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_overlap.md)
  now score time points after the last onset instead of treating them as
  missing. Two smaller changes on the fast path: fixations with tied
  onsets previously had their coordinates averaged
  (`approx(ties = mean)`) and now yield the last tied fixation, as
  `fast = FALSE` does; and a fixation with `NA` coordinates was
  previously skipped, carrying the preceding fixation forward, and now
  yields `NA`. The `fast = FALSE` path now orders fixations by onset
  first, as the fast path does; for groups whose onsets were not in
  increasing order it previously returned coordinates by row position
  rather than by onset.

- [`rep_fixations()`](https://bbuchsbaum.github.io/eyesim/reference/rep_fixations.md)
  now counts `floor(duration * resolution)` replicates with a
  floating-point tolerance. A duration of 0.29 at resolution 100
  previously gave 28 copies because `0.29 / 0.01` evaluates to
  28.999999999999996; it now gives 29.

- [`sample_density_time()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density_time.md)
  now treats every `time_bins` interval as half-open, `[lower, upper)`,
  including the last, as documented. The final bin was previously
  closed, so a time point equal to the last break (for example t = 3000
  with the default `times` and bins ending at 3000) was averaged into
  the last bin. It now falls outside all bins, which changes the last
  `bin_*`, `perm_bin_*`, and `diff_bin_*` values whenever a sample lies
  on the final break.

- [`multi_match()`](https://bbuchsbaum.github.io/eyesim/reference/multi_match.md)
  now computes `mm_position_emd` over all fixations. It previously used
  the saccade table, which omits the final fixation, so two
  three-fixation paths differing only in their last fixation scored 1.
  Every `mm_position_emd` value changes, including the `_perm` and
  `_diff` columns from
  [`scanpath_similarity()`](https://bbuchsbaum.github.io/eyesim/reference/scanpath_similarity.md).

- [`template_multireg()`](https://bbuchsbaum.github.io/eyesim/reference/template_multireg.md)
  now defaults to `method = "lm"`, as documented. Omitting `method`
  previously failed with “the condition has length \> 1”. Unknown
  methods are rejected by
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

- [`template_regression()`](https://bbuchsbaum.github.io/eyesim/reference/template_regression.md)
  now stops with a clear message when a `baseline_key` value used by
  `source_tab` appears in more than one row of `baseline_tab`. It
  previously failed with “\$ operator is invalid for atomic vectors”.

- [`template_sample()`](https://bbuchsbaum.github.io/eyesim/reference/template_sample.md)
  now drops rows whose template or fixation group is `NULL`, with a
  warning, as its code comment intended. The
  [`is.null()`](https://rdrr.io/r/base/NULL.html) filter on list columns
  was a no-op, so a single `NULL` template made the whole call fail.

- [`fixation_entropy()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.md)
  on a fixation group with a single fixation now returns `NA` for
  `method = "density"` when `sigma` is not supplied, because
  [`suggest_sigma()`](https://bbuchsbaum.github.io/eyesim/reference/suggest_sigma.md)
  cannot estimate a bandwidth from one point. It previously failed with
  “missing values present in assertion”.
  [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  now reports an `NA` or non-finite `sigma` with a clear message.

- [`fixation_entropy()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.md)
  now rejects maps with negative values, such as the difference of two
  densities, with a clear error. Previously an exact difference map
  (total zero) gave `NA` while a signed map with a positive total gave a
  meaningless number. Maps with zero total mass still give `NA`.

- Documentation only, no behaviour change:

  - [`?eye_density.fixation_group`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.fixation_group.md)
    now states that `sigma` is the kernel standard deviation under
    `kde_pkg = "ks"` but the
    [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html) bandwidth
    under `"MASS"`, whose kernel standard deviation is `sigma / 4`.
  - [`?sample_density`](https://bbuchsbaum.github.io/eyesim/reference/sample_density.md)
    now states that points are looked up at the nearest lattice point
    with [`round()`](https://rdrr.io/r/base/Round.html), so exact
    midpoints resolve half to even on the index scale, and that
    off-lattice points are clamped to the edge.
  - [`?fixation_entropy.default`](https://bbuchsbaum.github.io/eyesim/reference/fixation_entropy.default.md)
    now states that the grid method counts fixations outside
    `xbounds`/`ybounds` in the nearest edge cell.
