# Calculate Eye Density by Groups

This function calculates the eye density for fixations, grouped by
specified variables.

## Usage

``` r
density_by(
  x,
  groups,
  sigma = 50,
  xbounds = c(0, 1000),
  ybounds = c(0, 1000),
  outdim = c(100, 100),
  duration_weighted = TRUE,
  window = NULL,
  min_fixations = 2,
  keep_vars = NULL,
  fixvar = "fixgroup",
  result_name = "density",
  ...
)
```

## Arguments

- x:

  A data frame containing fixations and additional grouping variables.

- groups:

  A character vector specifying the grouping variables to use.

- sigma:

  A numeric value or a numeric vector specifying the bandwidth(s) for
  the kernel density estimation. If a vector is provided, multiscale
  densities are computed. Default is 50.

- xbounds:

  A numeric vector of length 2 specifying the x-axis bounds for the
  density calculation (default is c(0, 1000)).

- ybounds:

  A numeric vector of length 2 specifying the y-axis bounds for the
  density calculation (default is c(0, 1000)).

- outdim:

  A numeric vector of length 2 specifying the dimensions of the output
  density matrix (default is c(100, 100)).

- duration_weighted:

  A logical value indicating whether the density should be weighted by
  fixation duration (default is TRUE).

- window:

  A numeric vector of length 2 specifying the time window for selecting
  fixations (default is NULL).

- min_fixations:

  Minimum number of fixations required for computing a density map. Rows
  with fewer fixations after optional filtering will receive \`NULL\` in
  the result column. Default is 2.

- keep_vars:

  A character vector specifying additional variables to keep in the
  output (default is NULL).

- result_name:

  A character string specifying the name for the density result variable
  (default is "density").

- ...:

  Additional arguments passed to the \`eye_density.fixation_group\`
  function.

## Value

A data frame containing the original grouping variables, the fixations,
and the density result. If \`sigma\` is a vector, the column specified
by \`result_name\` will contain \`eye_density_multiscale\` objects.

## Examples

``` r
# Create a data frame with fixations and a grouping variable
fixations <- data.frame(
  subject = rep(c("A", "B"), each = 25),
  x = runif(50, 0, 1000),
  y = runif(50, 0, 1000),
  duration = runif(50, 1, 5),
  onset = seq(1, 50)
)
eyetab <- eye_table("x", "y", "duration", "onset", groupvar=c("subject"), data=fixations)
#> Error in select_at(., c("x", "y", "duration", "onset", groupvar)): could not find function "select_at"

# Calculate eye density by subject
result <- density_by(eyetab, groups = "subject")
#> Error: object 'eyetab' not found
```
