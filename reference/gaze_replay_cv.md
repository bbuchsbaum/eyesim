# Cross-fitted directional GazeWeave Replay analysis

Fits warp, background, emission, and transition parameters on training
items, then scores every held-out source path against its permitted
encoding candidates.

## Usage

``` r
gaze_replay_cv(
  ref_tab,
  source_tab,
  match_on,
  contrast_on = NULL,
  refvar = "fixgroup",
  sourcevar = "fixgroup",
  spec = gaze_replay_spec(),
  split_on = match_on,
  n_folds = NULL,
  seed = 1,
  fit_source_filter = NULL,
  eval_source_filter = NULL,
  template_on = NULL
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

  A \[gaze_replay_spec()\].

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

- template_on:

  Optional columns identifying separate encoding presentations within
  each matched episode. See \[fit_gaze_replay_model()\].

## Value

A \`gaze_replay_fit\` with \`gaze_info_bits\` as its primary endpoint.

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
fit <- gaze_replay_cv(study, recall, match_on = c("participant", "image"),
                      contrast_on = "participant", n_folds = 2)
#> Warning: Replay: the temperature calibration could not be fitted in 2 of 2 outer folds, so their rows are scored with the declared prior (gaze_info_bits = 0). Reason: fewer than 4 candidates in at least one training contrast stratum. Use fewer folds or more items per contrast stratum.
fit
#> Cross-fitted GazeWeave Replay analysis
#>   rows: 12 
#>   folds: 2 
#>   mean gaze_info_bits: 0 
#>   calibration fallback: 2 of 2 folds (see the gaze_calibration_fallback warning)
head(broom::tidy(fit))
#> # A tibble: 6 × 3
#>   participant image gaze_info_bits
#>   <chr>       <int>          <dbl>
#> 1 p1              1              0
#> 2 p2              1              0
#> 3 p1              2              0
#> 4 p2              2              0
#> 5 p1              3              0
#> 6 p2              3              0
# }
```
