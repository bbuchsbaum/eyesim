# Nested cross-validated comparator court for GazeWeave

Every method receives the same outer folds, held-out candidate sets,
candidate priors, and cross-fitted warp. Registered ridge composites
select their penalty using inner folds only. Frozen scores remain
explicitly uncalibrated compatibility scores.

## Usage

``` r
gaze_baseline_cv(
  ref_tab,
  source_tab,
  match_on,
  contrast_on = NULL,
  refvar = "fixgroup",
  sourcevar = "fixgroup",
  spec,
  split_on = match_on,
  n_folds = NULL,
  seed = 1,
  fit_source_filter = NULL,
  eval_source_filter = NULL,
  inner_split_on = split_on
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

  A \[gaze_baseline_spec()\].

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

  Optional logical vectors or functions evaluated on \`source_tab\`.
  They permit a model to be trained on positive-control rows and
  evaluated on both signal and independently generated null rows without
  admitting the null rows to preprocessing or hyperparameter selection.

- inner_split_on:

  Columns defining inner-fold groups. Defaults to \`split_on\`; specify
  a finer training-only unit when an outer fold holds out an entire
  session or replication.

## Value

A \`gaze_baseline_fit\` with one long-format row per held-out trial and
comparator.
