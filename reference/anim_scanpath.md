# Animate a Fixation Scanpath with gganimate

This function creates an animated visualization of a fixation scanpath
using gganimate.

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

  An optional image file name to use as the background.

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
[`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md),
[`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md)

## Examples

``` r
# Create a fixation group
fg <- fixation_group(x=c(.1,.5,1), y=c(1,.5,1), onset=1:3, duration=rep(1,3))
# Animate the scanpath for the fixation group
anim_sp <- anim_scanpath(fg)
```
