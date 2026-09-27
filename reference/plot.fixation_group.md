# Plot a fixation_group object

Plots a fixation group in stimulus coordinates with the shared eyesim
style. \`type = "points"\` draws the scanpath: fixations coloured by
onset and sized by duration, joined in order. The density types draw a
kernel density of the fixations on the eyesim density scale, with the
fixations overlaid in grey. The plot keeps a 1:1 aspect ratio;
\`xlim\`/\`ylim\` set the visible extent (and the extent of
\`bg_image\`) without dropping data. A ring marks the first fixation.

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
  colours = NULL,
  alpha_range = c(0, 0.9),
  alpha = 0.8,
  window = NULL,
  transform = c("identity", "sqroot", "curoot", "rank"),
  legend = TRUE,
  limits = NULL,
  aspect = c("equal", "free"),
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

  Whether to size points according to fixation duration; point area is
  proportional to duration (default: TRUE).

- show_points:

  Whether to show the fixations as points (default: TRUE).

- show_path:

  Whether to show the fixation path (default: TRUE).

- bins:

  Number of density bands for \`type = "density"\` (default:
  max(as.integer(length(x\$x)/10), 4)); \`type = "filled_contour"\` uses
  10 bands.

- bg_image:

  An optional background image file name (or \`cimg\`).

- colours:

  Optional colour ramp. For \`type = "points"\` it replaces the onset
  (time) ramp; for the density types it replaces the density ramp.
  Default \`NULL\` uses the eyesim ramps.

- alpha_range:

  Opacity at zero and at maximum density for the density layers
  (default: c(0, 0.9), so zero density is transparent).

- alpha:

  The opacity level for the points (default: 0.8).

- window:

  A vector specifying the time window for selecting fixations (default:
  NULL).

- transform:

  The transformation applied to the density colour scale (default:
  c("identity", "sqroot", "curoot", "rank")).

- legend:

  Whether to show legends: onset (points) or density (density types),
  plus duration when \`size_points = TRUE\` (default: \`TRUE\`).

- limits:

  Optional density limits for the colour scale of the density types
  (default: from zero to the maximum). Use the same \`limits\` across
  plots to compare them on one scale; for the banded types they also fix
  the band edges. Densities above the upper limit take the top colour.

- aspect:

  \`"equal"\` (default) keeps one data unit the same length on both
  axes, so saccade angles and cluster shapes are true; wide scanpaths
  then leave blank space. \`"free"\` fills the panel and distorts the
  geometry.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot object representing the fixation group plot.

## Details

Fixations are numbered and joined in onset order, whatever the row
order. With fewer than 50 fixations, the numbers are placed when the
plot is drawn, using its actual size: next to their own fixation, or
further out with a leader line that stays clear of every other fixation,
never covering another number or marker. Fixations whose markers overlap
are outlined as a group with one label listing their numbers. Any number
still without a position is listed in a note inside the panel (or
counted, if more than six).

\`type = "density"\` and \`type = "filled_contour"\` draw
non-overlapping density bands between evenly spaced edges from zero,
each coloured at its midpoint; the band containing zero is transparent,
and the stepped colour bar shows the same edges. Over a \`bg_image\`,
thin white contour lines separate the density from the stimulus colours.

## See also

Other visualization:
[`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md),
[`element_text_wrap()`](https://bbuchsbaum.github.io/eyesim/reference/element_text_wrap.md),
[`eyesim_colours()`](https://bbuchsbaum.github.io/eyesim/reference/eyesim_colours.md),
[`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md),
[`scale_fill_eyesim_density()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md),
[`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)

## Examples

``` r
# Create a fixation_group object
fg <- fixation_group(x=runif(50, 0, 100), y=runif(50, 0, 100), duration=rep(1,50), onset=seq(1,50))
# Plot the fixation group using the S3 method
plot(fg)
```
