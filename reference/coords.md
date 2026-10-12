# extract coordinates

extract coordinates

## Usage

``` r
coords(x)
```

## Arguments

- x:

  the object

## Value

A two-column numeric matrix with columns `x` and `y`, one row per
fixation.

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
coords(fg)
#>        x   y
#> [1,] 100 100
#> [2,] 200 150
#> [3,] 300 200
```
