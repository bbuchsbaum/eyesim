# Add Scanpath to Dataset

This function adds a scanpath to the input dataset.

## Usage

``` r
add_scanpath(x, ...)
```

## Arguments

- x:

  The input dataset to which the scanpath will be added.

- ...:

  Additional arguments to be passed to the method that adds the
  scanpath.

## Value

A dataset with the added scanpath.

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
tab <- tibble::tibble(trial = 1, fixgroup = list(fg))
sp <- add_scanpath(tab)$scanpath[[1]]
# Saccade amplitude (rho) and direction (theta) from each fixation to the next
as.data.frame(sp)[, c("x", "y", "rho", "theta")]
#>     x   y      rho     theta
#> 1 100 100 111.8034 0.4636476
#> 2 200 150 111.8034 0.4636476
#> 3 300 200   0.0000 0.0000000
```
