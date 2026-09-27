# Extract one auditable Transport result

Transport has one primary outcome, \`gaze_info_bits\`. Coverage, spatial
error, local order, contraction, rank, and solver behavior are nested
under \`diagnostics\`; they are explanatory properties of the optimized
alignment, not additional evidence scores. When several study
presentations are available, their calibrated likelihoods receive fixed
equal prior weight.

## Usage

``` r
gaze_transport_result(x, row = NULL, key = NULL)
```

## Arguments

- x:

  A fitted object returned by \[gaze_transport_cv()\].

- row:

  Optional result row number.

- key:

  Optional named list selecting one result by identifiers.

## Value

A \`gaze_transport_result\` containing one \`gaze_info_bits\` value,
nested diagnostics, candidate evidence, and provenance.
