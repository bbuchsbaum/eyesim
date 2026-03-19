# sample_density

Sample a smooth fixation density map with a set of discrete fixations.

This function samples a smooth fixation density map represented by the
object `x` with a set of discrete fixations provided in `fix`.

## Usage

``` r
sample_density(x, fix, ...)

# S3 method for class 'density'
sample_density(
  x,
  fix,
  times = NULL,
  normalize = c("none", "max", "sum", "zscore"),
  ...
)
```

## Arguments

- x:

  An object of class "density" representing the smooth fixation density
  map.

- fix:

  A data frame or tibble containing discrete fixations with columns "x",
  "y", and "onset".

- ...:

  Additional arguments passed on to the method, including `normalize`
  (one of `"none"`, `"max"`, `"sum"`, or `"zscore"`) to control density
  map normalization before sampling.

- times:

  A vector of numeric values representing the time points at which the
  density map should be sampled (default is NULL).

- normalize:

  A character string specifying how to normalize the density map before
  sampling. One of:

  `"none"`

  :   No normalization (default). Returns raw density values.

  `"max"`

  :   Divide by the maximum density value, yielding values in \[0, 1\].

  `"sum"`

  :   Divide by the sum of all density values, yielding a probability
      distribution.

  `"zscore"`

  :   Z-score the density map (subtract mean, divide by SD). Useful for
      comparing across maps with different scales.

## Value

A data frame with columns 'z' (density estimates at fixation locations)
and 'time' (onset time of fixations).

A data frame with columns "z" and "time", where "z" contains the sampled
density values and "time" contains the corresponding time points.

## Details

This function samples a given smooth fixation density map with a set of
discrete fixations to estimate the density at the locations of those
fixations.

The function first checks if the `times` parameter is NULL. If so, it
directly samples the density map using the coordinates of the fixations
in the `fix` argument. If the `times` parameter is provided, the
function first calls the `sample_fixations` function to generate a new
fixation sequence with the specified time points, and then samples the
density map using the coordinates of the new fixation sequence. The
result is a data frame containing the sampled density values and the
corresponding time points.

## See also

`sample_density.density` for an example of a specific method
implementation for this generic.
