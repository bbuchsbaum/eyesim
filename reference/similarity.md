# Compute Similarity Between Two Objects

This function computes the similarity between two objects using a
specified similarity metric.

## Usage

``` r
similarity(x, y, method, ...)
```

## Arguments

- x:

  The first object to compare.

- y:

  The second object to compare.

- method:

  A character string specifying the similarity metric to be used.

- ...:

  Additional arguments passed to the similarity computation method.

## Value

A numeric value representing the similarity between the two input
objects.

## Details

Density maps are compared cell by cell, so two maps must share a
lattice: the same `x` and `y` grid coordinates, as produced by the same
bounds and `outdim`. Comparing maps on different lattices is an error. A
plain numeric `y` must have one value per grid cell of `x`.

## See also

Other similarity:
[`similarity.scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/similarity.scanpath.md)

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
fg2 <- fixation_group(x = c(110, 190, 320), y = c(90, 160, 210),
                      onset = c(0, 200, 400), duration = c(200, 200, 200))
d1 <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300),
                  outdim = c(40, 30))
d2 <- eye_density(fg2, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300),
                  outdim = c(40, 30))
similarity(d1, d2, method = "pearson")
#> [1] 0.9574041
similarity(d1, d2, method = "cosine")
#> [1] 0.9793084
```
