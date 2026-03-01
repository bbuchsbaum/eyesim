# Suggest Kernel Bandwidth for Density Estimation

Estimates a reasonable kernel bandwidth (sigma) for fixation density
maps using a 2D variant of Silverman's rule of thumb, adapted for
eye-tracking data. The estimate accounts for spatial spread of fixations
and sample size.

## Usage

``` r
suggest_sigma(x, y = NULL, xbounds = NULL, ybounds = NULL)
```

## Arguments

- x:

  A `fixation_group` object, or a numeric vector of x-coordinates.

- y:

  A numeric vector of y-coordinates (only used when `x` is not a
  `fixation_group`).

- xbounds:

  Optional numeric vector of length 2 for the x-axis display bounds.
  Used to scale the estimate relative to display size.

- ybounds:

  Optional numeric vector of length 2 for the y-axis display bounds.
  Used to scale the estimate relative to display size.

## Value

A single numeric value representing the suggested sigma (kernel standard
deviation in coordinate units).

## Details

The function uses a 2D Silverman rule: \\\sigma = n^{-1/6} \cdot
\sqrt{(IQR_x^2 + IQR_y^2)/2} / 1.349\\. When display bounds are
provided, the result is clamped to between 1% and 15% of the mean
display dimension to avoid extreme values.

## Examples

``` r
fg <- fixation_group(x = runif(30, 0, 1280), y = runif(30, 0, 1024),
                     onset = cumsum(rep(200, 30)), duration = rep(200, 30))
suggest_sigma(fg, xbounds = c(0, 1280), ybounds = c(0, 1024))
#> [1] 172.8
```
