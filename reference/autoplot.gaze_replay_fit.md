# Plot a fitted GazeWeave Replay analysis

Plot a fitted GazeWeave Replay analysis

## Usage

``` r
# S3 method for class 'gaze_replay_fit'
autoplot(
  object,
  row = NULL,
  key = NULL,
  type = c("combined", "overlay", "braid", "diagnostics"),
  coupling_threshold = 0.005,
  ...
)
```

## Arguments

- object:

  A \`gaze_replay_fit\`.

- row:

  Optional result row number.

- key:

  Optional named list selecting one result by identifiers.

- type:

  One of \`"combined"\`, \`"overlay"\`, \`"braid"\`, or
  \`"diagnostics"\`.

- coupling_threshold:

  Minimum posterior correspondence proportion drawn in the braid.

- ...:

  Unused.

## Value

A combined patchwork when available, or an individual \`ggplot\`.
