# Plot a Transport pair alignment

Plot a Transport pair alignment

## Usage

``` r
# S3 method for class 'gaze_transport_alignment'
autoplot(
  object,
  type = c("overlay", "braid", "raw_registered", "diagnostics"),
  coupling_threshold = 0.005,
  ...
)
```

## Arguments

- object:

  A \`gaze_transport_alignment\`.

- type:

  One of \`"overlay"\`, \`"braid"\`, \`"raw_registered"\`, or
  \`"diagnostics"\`.

- coupling_threshold:

  Minimum proportion of transported mass drawn in the braid.

- ...:

  Unused.

## Value

A \`ggplot\` object.
