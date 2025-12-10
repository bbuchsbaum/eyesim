# Plot a fixation_group object

This function creates a plot of the fixation_group object using ggplot2.
Different plot types and options can be specified to customize the
output.

## Usage

``` r
# S3 method for class 'fixation_group'
plot(
  x,
  type = c("points", "contour", "filled_contour", "density", "raster"),
  bandwidth = 60,
  xlim = range(x$x),
  ylim = range(x$y),
  size_points = TRUE,
  show_points = TRUE,
  show_path = TRUE,
  bins = max(as.integer(length(x$x)/10), 4),
  bg_image = NULL,
  colours = rev(RColorBrewer.brewer.pal(n = 10, "Spectral")),
  alpha_range = c(0.5, 1),
  alpha = 0.8,
  window = NULL,
  transform = c("identity", "sqroot", "curoot", "rank"),
  ...
)
```

## Arguments

- x:

  A fixation_group object.

- type:

  The type of plot to display (default: c("points", "contour",
  "filled_contour", "density", "raster")).

- bandwidth:

  The bandwidth for the kernel density estimator (default: 60).

- xlim:

  The x-axis limits (default: range of x values in the fixation_group
  object).

- ylim:

  The y-axis limits (default: range of y values in the fixation_group
  object).

- size_points:

  Whether to size points according to fixation duration (default: TRUE).

- show_points:

  Whether to show the fixations as points (default: TRUE).

- show_path:

  Whether to show the fixation path (default: TRUE).

- bins:

  Number of bins for density calculations (default:
  max(as.integer(length(x\$x)/10), 4)).

- bg_image:

  An optional background image file name.

- colours:

  Color palette to use for the plot (default:
  rev(RColorBrewer.brewer.pal(n=10, "Spectral"))).

- alpha_range:

  A vector specifying the minimum and maximum alpha values for density
  plots (default: c(0.5, 1)).

- alpha:

  The opacity level for the points (default: 0.8).

- window:

  A vector specifying the time window for selecting fixations (default:
  NULL).

- transform:

  The transformation to apply to the density values (default:
  c("identity", "sqroot", "curoot", "rank")).

- ...:

  Additional arguments (currently unused).

## Value

A ggplot object representing the fixation group plot.

## See also

Other visualization:
[`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md),
[`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md)

## Examples

``` r
# Create a fixation_group object
fg <- fixation_group(x=runif(50, 0, 100), y=runif(50, 0, 100), duration=rep(1,50), onset=seq(1,50))
# Plot the fixation group
plot_fixation_group(fg)
#> Error in plot_fixation_group(fg): could not find function "plot_fixation_group"
```
