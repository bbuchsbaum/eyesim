# Plot a GazeWeave Replay alignment

Plot a GazeWeave Replay alignment

## Usage

``` r
# S3 method for class 'gaze_replay_alignment'
autoplot(
  object,
  type = c("overlay", "braid", "diagnostics"),
  coupling_threshold = 0.005,
  ...
)
```

## Arguments

- object:

  A \`gaze_replay_alignment\`.

- type:

  One of \`"overlay"\`, \`"braid"\`, or \`"diagnostics"\`.

- coupling_threshold:

  Minimum proportion of posterior correspondence drawn in the braid.

- ...:

  Unused.

## Value

A \`ggplot\` object.
