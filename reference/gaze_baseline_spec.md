# Specify the fair GazeWeave comparator court

The comparator court preserves the native MultiMatch dimensions, uses a
fixed multiscale density score, and includes order-free
elastic-consensus matching. Ridge composites for registered MultiMatch
and density features are tuned in inner folds of each outer training
split.

## Usage

``` r
gaze_baseline_spec(
  screen,
  density_sigmas,
  density_grid = 32L,
  methods = c("multimatch", "density", "elastic"),
  warp = gaze_warp_none(),
  lambda_grid = c(0, 0.01, 0.1, 1, 10),
  inner_folds = 2L,
  elastic_radii = c(consensus = 2, rigidity = 10, matching = 1),
  elastic_maxit = 25L,
  elastic_tolerance = 1e-05
)
```

## Arguments

- screen:

  Screen geometry created by \[gaze_screen()\].

- density_sigmas:

  Positive density bandwidths in \`screen\$unit\`.

- density_grid:

  Number of grid points per spatial axis.

- methods:

  Any of \`"multimatch"\`, \`"density"\`, and \`"elastic"\`.

- warp:

  A candidate-invariant warp specification.

- lambda_grid:

  Non-negative ridge penalties.

- inner_folds:

  Number of inner folds used to select the ridge penalty.

- elastic_radii:

  Named consensus, rigidity, and matching radii in \`screen\$unit\`.
  Defaults follow the 2, 10, and 1 degree values proposed by Wang,
  Holmqvist, and Alexa when the declared unit is degrees.

- elastic_maxit:

  Maximum elastic-consensus iterations.

- elastic_tolerance:

  Convergence tolerance in screen coordinates.

## Value

A frozen \`gaze_baseline_spec\`.
