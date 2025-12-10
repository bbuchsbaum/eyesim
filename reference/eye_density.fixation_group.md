# Compute a density map for a fixation group.

This function computes a density map for a given fixation group using
kernel density estimation.

## Usage

``` r
# S3 method for class 'fixation_group'
eye_density(
  x,
  sigma = 50,
  xbounds = c(min(x$x), max(x$x)),
  ybounds = c(min(x$y), max(x$y)),
  outdim = c(100, 100),
  normalize = TRUE,
  duration_weighted = FALSE,
  window = NULL,
  min_fixations = 2,
  origin = c(0, 0),
  kde_pkg = "ks",
  ...
)
```

## Arguments

- x:

  A fixation_group object.

- sigma:

  The standard deviation(s) of the kernel. Can be a single numeric value
  or a numeric vector. If a vector is provided, a multiscale density
  object (\`eye_density_multiscale\`) will be created. Default is 50.

- xbounds:

  The x-axis bounds. Default is the range of x values in the fixation
  group.

- ybounds:

  The y-axis bounds. Default is the range of y values in the fixation
  group.

- outdim:

  The output dimensions of the density map. Default is c(100, 100).

- normalize:

  Whether to normalize the output map. Default is TRUE.

- duration_weighted:

  Whether to weight the fixations by their duration. Default is FALSE.

- window:

  The temporal window over which to compute the density map. Default is
  NULL.

- min_fixations:

  Minimum number of fixations required to compute a density map. If
  fewer fixations are present after optional filtering, the function
  returns NULL. Default is 2.

- origin:

  The origin of the coordinate system. Default is c(0,0).

## Value

An object of class \`eye_density\` (inheriting from \`density\` and
\`list\`) if \`sigma\` is a single value, or an object of class
\`eye_density_multiscale\` (a list of \`eye_density\` objects) if
\`sigma\` is a vector. Returns \`NULL\` if filtering by \`window\`
leaves fewer than \`min_fixations\` fixations, or if density computation
fails (e.g., due to zero weights).

## Details

The function computes a density map for a given fixation group using
kernel density estimation. If \`sigma\` is a single value, it computes a
standard density map. If \`sigma\` is a vector, it computes a density
map for each value in \`sigma\` and returns them packaged as an
\`eye_density_multiscale\` object, which is a list of individual
\`eye_density\` objects.
