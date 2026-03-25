# Entropy of fixation patterns

These methods quantify how concentrated or diffuse a fixation pattern is
using Shannon entropy. For density maps, entropy is computed from the
normalized density surface. For fixation groups, entropy can be computed
from either a derived density map or a discrete occupancy grid.

## Usage

``` r
# Default S3 method
fixation_entropy(x, ...)

# S3 method for class 'eye_density'
fixation_entropy(x, normalize = TRUE, base = exp(1), ...)

# S3 method for class 'density'
fixation_entropy(x, normalize = TRUE, base = exp(1), ...)

# S3 method for class 'eye_density_multiscale'
fixation_entropy(
  x,
  normalize = TRUE,
  base = exp(1),
  aggregate = c("mean", "none"),
  ...
)

# S3 method for class 'fixation_group'
fixation_entropy(
  x,
  normalize = TRUE,
  base = exp(1),
  method = c("density", "grid"),
  sigma = NULL,
  xbounds = NULL,
  ybounds = NULL,
  outdim = c(50, 50),
  grid = c(10, 10),
  duration_weighted = FALSE,
  ...
)
```

## Arguments

- x:

  The input object.

- ...:

  Additional arguments passed to \`eye_density()\` when \`method =
  "density"\`.

- normalize:

  Logical; if \`TRUE\` (default), divide entropy by the maximum possible
  entropy for the number of valid bins so results lie in \`\[0, 1\]\`.

- base:

  Logarithm base used in the Shannon entropy. Defaults to \`exp(1)\`;
  use \`2\` for bits.

- aggregate:

  For multiscale density objects, one of \`"mean"\` (default) or
  \`"none"\`.

- method:

  For \`fixation_group\` objects, one of \`"density"\` (default) or
  \`"grid"\`.

- sigma:

  Optional bandwidth for density-based entropy on fixation groups. If
  \`NULL\`, \`suggest_sigma()\` is used.

- xbounds, ybounds:

  Optional display bounds for fixation groups. If not supplied, the
  observed fixation ranges are used with a small padding.

- outdim:

  Grid dimensions for density-based entropy from fixation groups.

- grid:

  Grid dimensions for occupancy-grid entropy from fixation groups.

- duration_weighted:

  Logical; if \`TRUE\`, duration-weighted KDE is used for density-based
  entropy.

## Value

A numeric entropy value. Multiscale density objects return either the
mean entropy across scales or a named numeric vector, depending on
\`aggregate\`.
