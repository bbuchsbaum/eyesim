# Generate a Simulated Eye-Movement Data Frame

The \`simulate_eye_table\` function generates a simulated \`eye_table\`
object containing eye-movement data with columns for x and y
coordinates, fixation duration, onset time, and an optional grouping
variable.

## Usage

``` r
simulate_eye_table(
  n_fixations,
  n_groups,
  clip_bounds = c(0, 1280, 0, 1280),
  relative_coords = TRUE
)
```

## Arguments

- n_fixations:

  The number of fixations to simulate.

- n_groups:

  The number of groups to simulate.

- clip_bounds:

  A numeric vector of length 4 representing the clip bounds for the
  field of view in the form c(xmin, xmax, ymin, ymax). Default is c(0,
  1280, 0, 1280).

- relative_coords:

  A logical value indicating whether to compute relative coordinates
  (TRUE by default). If TRUE, x and y coordinates will be transformed
  based on the clip_bounds.

## Value

A \`data.frame\` of class "eye_table" with simulated data.

## Examples

``` r
sim_eye_table <- simulate_eye_table(n_fixations = 100, n_groups = 10)
#> Error in select_at(., c("x", "y", "duration", "onset", groupvar)): could not find function "select_at"
```
