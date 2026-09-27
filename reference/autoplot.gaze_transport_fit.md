# Plot one fitted Transport result

\`type = "evidence"\` plots the single prespecified score. Alignment
panels display one named study episode; this selection never changes the
score, which always uses the fit's fixed equal-prior episode likelihood
mixture.

## Usage

``` r
# S3 method for class 'gaze_transport_fit'
autoplot(
  object,
  row = NULL,
  key = NULL,
  type = c("combined", "evidence", "overlay", "braid", "raw_registered", "diagnostics"),
  episode = NULL,
  coupling_threshold = 0.005,
  ...
)
```

## Arguments

- object:

  A fitted object returned by \[gaze_transport_cv()\].

- row:

  Optional result row number.

- key:

  Optional named list selecting one result by identifiers.

- type:

  One of \`"combined"\`, \`"evidence"\`, \`"overlay"\`, \`"braid"\`,
  \`"raw_registered"\`, or \`"diagnostics"\`.

- episode:

  Optional retained episode identifier for an alignment panel.

- coupling_threshold:

  Minimum proportion of transported mass drawn in the braid.

- ...:

  Unused.

## Value

A combined patchwork when available, a list of plots otherwise, or one
\`ggplot\` for an individual panel.
