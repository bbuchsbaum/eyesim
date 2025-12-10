# sample_density

Sample a smooth fixation density map with a set of discrete fixations.

## Usage

``` r
sample_density(x, fix, ...)
```

## Arguments

- x:

  An object representing the smooth fixation density map to be sampled.

- fix:

  A data frame containing the discrete fixations used for sampling the
  density map.

- ...:

  Additional arguments passed on to the method.

## Value

A data frame with columns 'z' (density estimates at fixation locations)
and 'time' (onset time of fixations).

## Details

This function samples a given smooth fixation density map with a set of
discrete fixations to estimate the density at the locations of those
fixations.

## See also

[`sample_density.density`](https://bbuchsbaum.github.io/eyesim/reference/sample_density.density.md)
for an example of a specific method implementation for this generic.
