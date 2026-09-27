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

  A \[gaze_transport_spec()\] or \[gaze_replay_spec()\].

- engine:

  Either \`"transport"\` or \`"replay"\`. When omitted, it is inferred
  from \`spec\`.

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

  Optional reference columns identifying separate study presentations
  for Transport. Their likelihoods receive equal prior weight.

- priorvar:

  Optional Transport reference column containing one positive
  design-prior weight per candidate.

## Value

An engine-specific fitted object. Use \`broom::tidy()\` for the sole
primary endpoint, \`gaze_info_bits\`.

## Details

Both engines return held-out candidate probabilities and
\`gaze_info_bits = log2(p_true / prior_true)\`. Replay correspondence is
a posterior probability. Transport correspondence is an optimized
alignment and must not be interpreted as statistical uncertainty.
