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
  candidate and a nonmatch in each contrast stratum.

- seed:

  Fold-assignment seed. \`NULL\` uses the engine-specific default. Folds
  depend on \`seed\` alone; the caller's random number stream is left
  unchanged.

- fit_source_filter, eval_source_filter:

  Optional logical vectors or functions selecting fitting and evaluation
  source rows.

- template_on:

  Optional columns identifying separate encoding presentations within
  each matched episode. See \[fit_gaze_replay_model()\].

## Value

A \`gaze_replay_fit\` with \`gaze_info_bits\` as its primary endpoint.
