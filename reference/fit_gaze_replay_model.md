# Fit a directional GazeWeave replay model on training trials

The fitted warp, background distribution, emission scale, and transition
parameters are learned from the supplied matched training trials. When
\`contrast_on\` is supplied, the likelihood temperature is learned from
inner-held-out item scores with a fixed weak log-temperature prior.
Scored trials must be disjoint from all of these rows. Training strata
too small for two candidate-bearing inner folds retain temperature one
and record the fallback in \`model\$calibration\`.

## Usage

``` r
fit_gaze_replay_model(
  ref_tab,
  source_tab,
  match_on,
  contrast_on = NULL,
  refvar = "fixgroup",
  sourcevar = "fixgroup",
  spec = gaze_replay_spec(),
  template_on = NULL
)
```

## Arguments

- ref_tab, source_tab:

  Matched training reference and source tables.

- match_on:

  Columns defining matched paths.

- contrast_on:

  Optional columns defining training candidate sets for
  likelihood-temperature calibration.

- refvar, sourcevar:

  Fixation-group list columns.

- spec:

  A \[gaze_replay_spec()\].

- template_on:

  Optional columns that uniquely identify separate encoding
  presentations within each \`match_on\` episode. Their likelihoods are
  combined with equal prior weight; paths are never concatenated.

## Value

A fitted \`gaze_replay_model\`.
