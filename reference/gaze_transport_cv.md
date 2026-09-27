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
  priorvar = NULL
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

  Number of cross-fitting folds.

- seed:

  Fold-assignment seed. \`NULL\` uses the engine-specific default. Folds
  depend on \`seed\` alone; the caller's random number stream is left
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
