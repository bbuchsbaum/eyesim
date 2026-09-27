# Declare cross-fitted contraction registration for GazeWeave

Declare cross-fitted contraction registration for GazeWeave

## Usage

``` r
gaze_warp_contraction(
  center = "screen",
  translation = TRUE,
  fit_by = NULL,
  shrink = 1e-06
)
```

## Arguments

- center:

  Either \`"screen"\` or a finite two-element numeric vector.

- translation:

  Whether to estimate translation in addition to isotropic scale.

- fit_by:

  Optional columns defining separate warp strata, such as
  \`"participant"\`.

- shrink:

  Positive stabilization term for moment ratios.

## Value

A \`gaze_warp_spec\` object.
