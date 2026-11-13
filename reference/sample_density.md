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
  ...,
  transform = c("none", "log"),
  floor = NULL
)

# S3 method for class 'list'
sample_density(x, fix, times = NULL, ..., format = c("long", "wide"))

# S3 method for class 'data.frame'
sample_density(
  x,
  fix,
  times = NULL,
  ...,
  map_var = "density",
  name_var = "reference"
)

# S3 method for class 'eyesim_reference_map'
sample_density(
  x,
  fix,
  times = NULL,
  ...,
  normalize = "none",
  transform = "none",
  floor = NULL
)

# S3 method for class 'eyesim_references'
sample_density(x, fix, ..., format = c("wide", "long"))
```

## Arguments

- x:

  A density map, a named list or table of maps, or a reference set.

- fix:

  A table of locations with x/y coordinates. Single maps and map lists
  also use onset (unless times are supplied); reference sets require
  their subject/item key columns.

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

  :   Divide by the maximum density value, yielding values between zero
      and one.

  `"sum"`

  :   Divide by the sum of all density values, yielding a probability
      distribution.

  `"zscore"`

  :   Z-score the density map (subtract mean, divide by SD). Useful for
      comparing across maps with different scales.

- transform:

  Map value transformation, applied after normalization. `log` computes
  `log(z + floor)`; `none` leaves values unchanged.

- floor:

  Explicit positive additive offset required for `log`. For maps
  normalized to unit sum, `1 / (10 * length(x$z))` is a possible choice.
  This is an offset, not a hard lower bound.

- format:

  For a named list of maps, return long rows or one column per map.

- map_var, name_var:

  For a map table, columns containing maps and unique names.

## Value

Single-map sampling returns z and time columns. Named lists and map
tables return fixation/reference/z/time rows, or fixation/time and one
column per map with format = "wide". Reference sets return one column
per reference, aligned with the input locations, or fixation/reference/z
rows with format = "long". Reference values are already prepared and do
not accept further normalization or transformation during lookup.

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

Each point is looked up at the nearest lattice point, separately in x
and y; there is no interpolation. The nearest index is found with
[`round()`](https://rdrr.io/r/base/Round.html), so a coordinate exactly
midway between two grid points is resolved half to even on the index
scale (on the grid 0, 25, 50, 75, 100, the value 12.5 maps to 25 and
37.5 also maps to 25). Coordinates outside the lattice are clamped to
the nearest edge point, so an off-grid fixation takes the edge value
rather than `NA`.

## See also

`sample_density.density` for an example of a specific method
implementation for this generic.
