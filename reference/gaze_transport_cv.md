# Nested-calibrated exhaustive-candidate Transport

Outer folds hold out participant-item cells. Within each outer-training
set, warp fitting and episode-scale temperature/reliability calibration
are repeated on item-grouped inner folds. Every held-out row is scored
against every permitted candidate in its contrast stratum with the
actual declared candidate prior and one common solver, warp, and
equal-episode policy.

## Usage

``` r
gaze_transport_cv(
  ref_tab,
  source_tab,
  match_on,
  contrast_on = NULL,
  refvar = "fixgroup",
  sourcevar = "fixgroup",
  spec = gaze_transport_spec(),
  split_on = match_on,
  n_folds = NULL,
  seed = 20260822L,
  fit_source_filter = NULL,
  eval_source_filter = NULL,
  episode_on = NULL,
  priorvar = NULL,
  keep_alignments = c("true", "all", "none")
)
```

## Arguments

- ref_tab, source_tab:

  Reference and source tables.

- match_on:

  One or more columns defining the true template match.

- contrast_on:

  Optional columns restricting nonmatching candidates.

- refvar, sourcevar:

  List columns containing fixation groups.

- spec:

  A \[gaze_transport_spec()\].

- split_on:

  Columns defining the held-out unit. Defaults to \`match_on\`.

- n_folds:

  Number of cross-fitting folds. \`NULL\` uses up to five; for Replay,
  which scores each held-out row against the candidates in its own fold,
  the default is lowered when needed so that every fold keeps a true
  candidate and a nonmatch in each contrast stratum. Calibration needs
  more: each training fold must hold at least \`2 \* calibration folds\`
  items (4 by default) per contrast stratum, or that fold falls back to
  an uncalibrated score and a warning of class
  \`gaze_calibration_fallback\` is raised.

- seed:

  Fold-assignment seed. \`NULL\` uses the engine-specific default. Folds
  depend on \`seed\` alone: the draw uses R's default generator and
  locale-independent ordering whatever \`RNGkind()\` and locale are in
  effect, and the caller's generator and random number stream are left
  unchanged.

- fit_source_filter, eval_source_filter:

  Optional logical vectors or functions selecting fitting and evaluation
  source rows.

- episode_on:

  Optional columns identifying the separate study presentations within a
  candidate. Their likelihoods receive equal prior weight and are never
  concatenated.

- priorvar:

  Optional reference-table column containing one positive design-prior
  weight per candidate. \`NULL\` declares a uniform design.

- keep_alignments:

  Which candidate alignments to keep in \`results\$alignments\`:
  \`"true"\` (the default) keeps only each row's true candidate, which
  is all the alignment plots and diagnostics read; \`"all"\` keeps every
  candidate's alignment, for auditing; \`"none"\` keeps none. A full
  alignment holds every coverage node's plan and iteration record, so
  \`"all"\` grows with rows times candidates. Scores, the candidate
  tables in \`results\$candidates\` and the solver summary do not depend
  on this choice.

## Value

A \`gaze_transport_fit\` with held-out candidate probabilities, fold
receipts, response-blind quality diagnostics, and calibration checks.
\`results\$backend_fallbacks\`, \`results\$solver_stalled\` and
\`results\$solver_not_converged\` count, per held-out row, the episode
alignments solved by the reference backend after a native failure, those
recorded as \`"stalled_projection_limited"\`, and those in which a
coverage node reached \`maxit\` (\`"not_converged"\`, still scored). The
\`solver\` element summarises the solver revision, the backends used,
and these counts (including inner calibration alignments in
\`backend_fallback_count\`). Per-alignment fallback warnings are
replaced by one summary warning.

## Parallel evaluation

Held-out rows, inner-calibration rows and typicality candidates are
scored in parallel under a \[future\]\[future::plan\] plan with more
than one worker, for example \`future::plan(future::multisession,
workers = 4)\`. The solver draws no random numbers, so the results are
identical to a sequential run. The default sequential plan runs as
before.

## Experimental

GazeWeave is experimental. Its interface, defaults and result columns
may change in future releases; record the package version with any
results.

## Result columns

Every GazeWeave fit (\`gaze_transport_cv()\`, \`gaze_replay_cv()\`,
\`gaze_baseline_cv()\`) has \`results\` with one row per held-out source
row (per source row and method for the baseline court) holding the
source identifiers and these columns:

- \`.cv_fold\`:

  Outer fold that held the row out.

- \`gaze_info_bits\`:

  \`log2(posterior_true / prior_true)\`, the primary endpoint. \`NA\`
  for uncalibrated baseline methods.

- \`posterior_true\`, \`prior_true\`:

  Held-out probability of the true candidate, and its declared prior.

- \`odds_bits\`, \`log_loss\`, \`brier_score\`:

  Other proper scores of the same posterior.

- \`template_rank\`, \`top1_credit\`:

  Rank of the true candidate (1 is best) and credit for ranking it
  first, split over ties.

- \`candidate_count\`:

  Number of candidates scored.

- \`candidates\`:

  Per-row table with at least \`candidate_key\`, \`log_score\`,
  \`prior\`, \`posterior\` and \`is_true\`.

Each engine adds its own diagnostics (for example Transport's solver
counts or Replay's background level), documented on its own page. Every
fit also has \`spec\`, \`keys\`, \`folds\` and \`provenance\`;
\`provenance\$engine\` names the engine, and Transport and Replay record
their estimator \`provenance\$revision\`.

## Examples

``` r
# Example data: two participants study and later recall six images; each
# recall roughly repeats the study path of the same image.
set.seed(1)
make_path <- function(anchor, jitter) {
  fixation_group(
    x = anchor[1] + c(0, 2, 4, 3) + rnorm(4, sd = jitter),
    y = anchor[2] + c(0, 1, -1, 2) + rnorm(4, sd = jitter),
    onset = c(0, 250, 500, 750), duration = c(200, 220, 180, 240)
  )
}
study <- expand.grid(participant = c("p1", "p2"), image = 1:6,
                     stringsAsFactors = FALSE)
anchors <- lapply(study$image, function(i) c(3 * i, 2 * (i %% 3)))
study$fixgroup <- lapply(anchors, make_path, jitter = 0.2)
recall <- study[c("participant", "image")]
recall$fixgroup <- lapply(anchors, make_path, jitter = 0.6)
# \donttest{
fit <- gaze_transport_cv(study, recall,
                         match_on = c("participant", "image"),
                         contrast_on = "participant",
                         # 2 coverage nodes keep the example fast; the
                         # default is 12.
                         spec = gaze_transport_spec(coverage_nodes = 2),
                         n_folds = 2)
fit
#> Cross-fitted GazeWeave Transport
#>   held-out rows: 12 
#>   primary endpoint: gaze_info_bits
#>   finite scores: 12 
#>   mean gaze_info_bits: 2.3153 
#>   episode normalization: fixed equal-prior likelihood mixture; one episode is the identity case
mean(fit$results$gaze_info_bits)
#> [1] 2.315255
# The alignment behind the first held-out row.
gaze_transport_result(fit, row = 1)$diagnostics
#> $matched_coverage
#> [1] 0.4877017
#> 
#> $spatial_rmse
#> [1] 0.6258791
#> 
#> $spatial_rmse_unit
#> [1] "native_coordinate_units"
#> 
#> $local_order_preservation
#> [1] 0.99835
#> 
#> $contraction_scale
#> [1] 1
#> 
#> $template_rank
#> [1] 1
#> 
#> $candidate_count
#> [1] 6
#> 
#> $episode_count
#> [1] 1
#> 
#> $episode_normalization
#> [1] "fixed_equal_prior_likelihood_mixture"
#> 
#> $episodes
#>             episode equal_episode_weight   log_score matched_coverage
#> episode_1 episode_1                    1 -0.07275334        0.4877017
#>           spatial_rmse       spatial_rmse_unit local_order_preservation
#> episode_1    0.6258791 native_coordinate_units                  0.99835
#>           contraction_scale solver_converged       solver_backend stationarity
#> episode_1                 1             TRUE native_rcpparmadillo   3.1132e-05
#>           relative_objective_change
#> episode_1              5.860762e-08
#> 
#> $solver_stability
#> $solver_stability$all_converged
#> [1] TRUE
#> 
#> $solver_stability$backends
#> [1] "native_rcpparmadillo"
#> 
#> $solver_stability$maximum_stationarity
#> [1] 3.1132e-05
#> 
#> $solver_stability$maximum_relative_objective_change
#> [1] 5.860762e-08
#> 
#> 
# }
```
