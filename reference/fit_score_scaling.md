# Fit and apply group-specific score standardization

Fit and apply group-specific score standardization

## Usage

``` r
fit_score_scaling(
  scores,
  columns = attr(scores, "score_columns"),
  by = c("subject", "draw")
)

apply_score_scaling(scores, scaling)
```

## Arguments

- scores:

  A table with numeric score columns; missing values are allowed.

- columns:

  Score column names. Defaults to the scoring result's metadata.

- by:

  Grouping columns. Use character(0) for a single global fit.

- scaling:

  A fitted scaling object.

## Value

`fit_score_scaling()` returns an `eyesim_score_scaling` list with group
keys, columns, centers, scales, and counts. `apply_score_scaling()`
returns a score table with the fitted object attached as `scaling`.

## Details

Fitting uses finite observations, the arithmetic mean and sample SD. A
constant column maps to zero; an entirely missing column remains
missing. Application rejects unseen groups and never refits. Map
normalization is a separate operation. Scores that are already
standardized cannot be scaled again. Fitting separately by control draw
prevents additional draws from changing earlier draws' scaling.
