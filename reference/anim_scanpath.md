# Animate a Fixation Scanpath with gganimate

Creates an animated scanpath: fixations appear one at a time (by order
or by onset), earlier fixations stay visible as faded marks, and colour
encodes onset on the shared eyesim time scale.

## Usage

``` r
anim_scanpath(
  x,
  bg_image = NULL,
  xlim = range(x$x),
  ylim = range(x$y),
  alpha = 1,
  anim_over = c("index", "onset"),
  type = c("points", "raster"),
  time_bin = 1
)
```

## Arguments

- x:

  A \`fixation_group\` object.

- bg_image:

  An optional image file name (or \`cimg\`) to use as the background.

- xlim:

  The range in x coordinates (default: range of x values in the fixation
  group).

- ylim:

  The range in y coordinates (default: range of y values in the fixation
  group).

- alpha:

  The opacity of each dot (default: 1).

- anim_over:

  Animate over index (ordered) or onset (real time) (default: c("index",
  "onset")).

- type:

  Display as points or a raster (default: c("points", "raster")).

- time_bin:

  The size of the time bins (default: 1).

## Value

A gganimate object representing the animated scanpath.

## See also

Other visualization:
[`element_text_wrap()`](https://bbuchsbaum.github.io/eyesim/reference/element_text_wrap.md),
[`eyesim_colours()`](https://bbuchsbaum.github.io/eyesim/reference/eyesim_colours.md),
[`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md),
[`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md),
[`scale_fill_eyesim_density()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md),
[`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)

## Examples

``` r
# Create a fixation group
fg <- fixation_group(x=c(.1,.5,1), y=c(1,.5,1), onset=1:3, duration=rep(1,3))
# Animate the scanpath for the fixation group
if (requireNamespace("gganimate", quietly = TRUE)) {
  anim_sp <- anim_scanpath(fg)
}
```
