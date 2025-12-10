# Sample a smooth fixation density map with a set of discrete fixations.

This function samples a smooth fixation density map represented by the
object `x` with a set of discrete fixations provided in `fix`.

## Usage

``` r
# S3 method for class 'density'
sample_density(x, fix, times = NULL)
```

## Arguments

- x:

  An object of class "density" representing the smooth fixation density
  map.

- fix:

  A data frame or tibble containing discrete fixations with columns "x",
  "y", and "onset".

- times:

  A vector of numeric values representing the time points at which the
  density map should be sampled (default is NULL).

## Value

A data frame with columns "z" and "time", where "z" contains the sampled
density values and "time" contains the corresponding time points.

## Details

The function first checks if the `times` parameter is NULL. If so, it
directly samples the density map using the coordinates of the fixations
in the `fix` argument. If the `times` parameter is provided, the
function first calls the `sample_fixations` function to generate a new
fixation sequence with the specified time points, and then samples the
density map using the coordinates of the new fixation sequence. The
result is a data frame containing the sampled density values and the
corresponding time points.
