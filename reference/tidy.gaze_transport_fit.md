# Tidy a fitted Transport analysis

The default output contains only identifiers and the prespecified
primary endpoint. Use \`diagnostics = TRUE\` to add calibration fields
and one nested diagnostic record per row.

## Usage

``` r
# S3 method for class 'gaze_transport_fit'
tidy(x, diagnostics = FALSE, ...)
```

## Arguments

- x:

  A fitted object returned by \[gaze_transport_cv()\].

- diagnostics:

  Include candidate/calibration columns and a nested \`diagnostics\`
  list-column.

- ...:

  Unused.

## Value

A tibble. The default has identifiers and \`gaze_info_bits\` only.
