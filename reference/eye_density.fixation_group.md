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
  weights = NULL,
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
  With `kde_pkg = "ks"`, `sigma` is the standard deviation of the
  isotropic Gaussian kernel. With `kde_pkg = "MASS"`, `sigma` is passed
  as the bandwidth `h` of
  [`kde2d`](https://rdrr.io/pkg/MASS/man/kde2d.html), whose kernel
  standard deviation is `h / 4`; the same `sigma` therefore gives a
  kernel four times narrower. Use `4 * sigma` under MASS to approximate
  the ks map.

- xbounds:

  The x-axis bounds. Default is the range of x values in the fixation
  group.

- ybounds:

  The y-axis bounds. Default is the range of y values in the fixation
  group.

- outdim:

  The output dimensions of the density map. Default is c(100, 100).

- weights:

  Optional numeric vector of non-negative fixation weights, one per row
  of `x` (before any `window` filtering). Explicit weights take
  precedence over `duration_weighted`. If NULL and duration_weighted is
  TRUE, uses fixation durations as weights. Default is NULL.

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

- kde_pkg:

  A character string specifying which package to use for kernel density
  estimation. Options are "ks" (default) or "MASS"; any other value is
  an error. Both support weighted estimation; under "MASS" weights use
  an internal weighted version of
  [`MASS::kde2d`](https://rdrr.io/pkg/MASS/man/kde2d.html). Note the
  different meaning of `sigma` under "MASS".

- ...:

  Additional named arguments passed to
  [`kde`](https://mvstat.net/ks/reference/kde.html), for example
  `binned = FALSE`.
  [`eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/eye_density.md)
  sets `x`, `H`, `gridsize`, `xmin`, `xmax`, `w`, and the evaluation
  grid (`eval.points`) itself, so these cannot be supplied. Extra
  arguments are an error when `kde_pkg = "MASS"`.

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

After optional normalization, each map is passed through
`zapsmall(z, digits = 7)`. The precision is fixed and does not depend on
`options(digits)`. In R 4.5.1, where this was verified,
[`zapsmall()`](https://rdrr.io/r/base/zapsmall.html) computes
`round(z, max(0, 7 - log10(max(abs(z)))))`: values are rounded to 7
significant digits relative to the map maximum, so values below roughly
`max(z) * 5e-8` become zero.
