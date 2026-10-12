# get_density

Get the density matrix from an object.

## Usage

``` r
get_density(x, ...)
```

## Arguments

- x:

  the object

- ...:

  additional arguments

## Value

The density values as a numeric matrix (for an `eye_density`, its `z`
component, with rows indexing `x` and columns `y`).

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
ed <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300),
                  outdim = c(40, 30))
dim(get_density(ed))
#> [1] 40 30
```
