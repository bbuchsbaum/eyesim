# eye_density

Compute a density map for a set of eye fixations.

## Usage

``` r
eye_density(x, sigma, xbounds, ybounds, outdim, weights, normalize, ...)
```

## Arguments

- x:

  A data frame containing the eye fixations.

- sigma:

  Numeric, the standard deviation of the kernel used for density
  estimation.

- xbounds:

  Numeric vector of length 2, defining the minimum and maximum x-axis
  bounds of the density map.

- ybounds:

  Numeric vector of length 2, defining the minimum and maximum y-axis
  bounds of the density map.

- outdim:

  Numeric vector of length 2, specifying the number of rows and columns
  in the output density map.

- weights:

  Optional, a numeric vector of fixation weights to be used in the
  density estimation.

- normalize:

  Logical, whether to normalize the output density map such that its
  values sum to 1.

- ...:

  Additional arguments passed on to the method.

## Value

An object of class "eye_density", "density", and "list" containing the
computed density map and other relevant information.

## Details

This function uses kernel density estimation to compute a density map of
eye fixations. It takes various parameters to control the computation
and output of the density map.

## See also

[`eye_density.fixation_group`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.fixation_group.md)
for an example of a specific method implementation for this generic.
