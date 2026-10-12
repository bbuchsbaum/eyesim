# Declare a multiscale Gaussian spatial model for GazeWeave

The resulting spatial cost is the negative log of a normalized mixture
of Gaussian overlap kernels. Mixture weights are normalized to sum to
one.

## Usage

``` r
gaze_gaussian_mixture(
  sigmas,
  weights = NULL,
  unit = c("native", "px", "deg", "normalized")
)
```

## Arguments

- sigmas:

  Positive Gaussian bandwidths.

- weights:

  Optional non-negative mixture weights.

- unit:

  Coordinate unit for \`sigmas\`.

## Value

A \`gaze_spatial_spec\` object.

## Examples

``` r
# Two spatial scales, in the same unit as the fixation coordinates.
gaze_gaussian_mixture(c(30, 60), unit = "px")
#> $sigmas
#> [1] 30 60
#> 
#> $weights
#> [1] 0.5 0.5
#> 
#> $unit
#> [1] "px"
#> 
#> attr(,"class")
#> [1] "gaze_spatial_spec" "list"             
```
