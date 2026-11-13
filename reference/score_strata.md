# Score case-control assignments against reference maps

Score case-control assignments against reference maps

## Usage

``` r
score_strata(
  strata,
  refs,
  standardize = c("within_subject", "none"),
  scaling = NULL
)
```

## Arguments

- strata:

  Explicit case-control locations from
  [`sample_controls()`](https://bbuchsbaum.github.io/eyesim/reference/sample_controls.md).

- refs:

  A reference set from
  [`reference_maps()`](https://bbuchsbaum.github.io/eyesim/reference/reference_maps.md).
  Its preparation is fixed; log transforms are not applied again during
  scoring.

- standardize:

  Standardize each reference within subject and control draw, or retain
  raw reference values.

- scaling:

  Optional fitted object from
  [`fit_score_scaling()`](https://bbuchsbaum.github.io/eyesim/reference/fit_score_scaling.md).
  Supplying it applies those parameters instead of estimating them.
  Explicitly combining `standardize = "none"` with a scaling object is
  an error.

## Value

The strata table with one column per reference and attributes
`score_columns`, `reference_spec`, and `scaling` (when applied).
