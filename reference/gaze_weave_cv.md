# Cross-fitted GazeWeave analysis

\`gaze_weave_cv()\` is the common entry point for the two GazeWeave
engines. Use Transport for symmetric explanatory alignment and Replay
for a directional encoding-to-recall model. There is no universal engine
default; the supplied engine-specific specification selects the
scientific model.

## Usage

``` r
gaze_weave_cv(
  ref_tab,
  source_tab,
  match_on,
  contrast_on = NULL,
  refvar = "fixgroup",
  sourcevar = "fixgroup",
  spec = NULL,
  engine = NULL,
  split_on = match_on,
  n_folds = NULL,
  seed = NULL,
  fit_source_filter = NULL,
  eval_source_filter = NULL,
  episode_on = NULL,
  priorvar = NULL,
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

  A \[gaze_transport_spec()\] or \[gaze_replay_spec()\].

- engine:

  Either \`"transport"\` or \`"replay"\`. When omitted, it is inferred
  from \`spec\`.

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

  Optional reference columns identifying separate study presentations
  for Transport. Their likelihoods receive equal prior weight.

- priorvar:

  Optional Transport reference column containing one positive
  design-prior weight per candidate.

- template_on:

  Optional reference columns identifying separate encoding templates of
  one candidate for Replay (the counterpart of \`episode_on\`); see
  \[gaze_replay_cv()\].

## Value

An engine-specific fitted object. Use \`broom::tidy()\` for the sole
primary endpoint, \`gaze_info_bits\`.

## Details

Both engines return held-out candidate probabilities and
\`gaze_info_bits = log2(p_true / prior_true)\`. Replay correspondence is
a posterior probability. Transport correspondence is an optimized
alignment and must not be interpreted as statistical uncertainty.

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

## Experimental

GazeWeave is experimental. Its interface, defaults and result columns
may change in future releases; record the package version with any
results.

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
# Which studied image does each recall come from? Candidates are the
# images the same participant studied.
fit <- gaze_weave_cv(study, recall, match_on = c("participant", "image"),
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
head(broom::tidy(fit))
#> # A tibble: 6 × 3
#>   participant image gaze_info_bits
#>   <chr>       <int>          <dbl>
#> 1 p1              1           2.47
#> 2 p2              1           2.35
#> 3 p1              2           2.35
#> 4 p2              2           2.42
#> 5 p1              3           2.38
#> 6 p2              3           2.42
# }
```
