# Plot Eye Density

Draws a fixation density map on the shared eyesim density scale. Zero
density is transparent, so a background stimulus shows through, and the
map keeps the aspect ratio of its coordinate grid. Over a \`bg_image\`,
thin white contour lines at five evenly spaced levels (from zero to the
upper limit or the map's maximum) keep density edges visible.

## Usage

``` r
# S3 method for class 'eye_density'
plot(
  x,
  alpha = 0.8,
  bg_image = NULL,
  transform = c("identity", "sqroot", "curoot", "rank"),
  colours = NULL,
  legend = TRUE,
  limits = NULL,
  ...
)
```

## Arguments

- x:

  An "eye_density" object.

- alpha:

  Maximum opacity of the density layer (default: 0.8); lower densities
  are progressively more transparent.

- bg_image:

  An optional image file name (or \`cimg\`) to use as the background.

- transform:

  The transformation to apply to the density values (default:
  c("identity", "sqroot", "curoot", "rank")).

- colours:

  Optional vector of colours for the density ramp (default: the eyesim
  density ramp).

- legend:

  Whether to show the density colour bar (default: \`TRUE\`).

- limits:

  Optional density limits for the colour scale (default: from zero to
  the map's maximum). Give several maps the same \`limits\` to compare
  them on one scale.

- ...:

  Additional args

## Value

A ggplot object representing the eye density plot.

## See also

Other visualization:
[`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md),
[`element_text_wrap()`](https://bbuchsbaum.github.io/eyesim/reference/element_text_wrap.md),
[`eyesim_colours()`](https://bbuchsbaum.github.io/eyesim/reference/eyesim_colours.md),
[`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md),
[`scale_fill_eyesim_density()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md),
[`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)

## Examples

``` r
# Create a fixation group and compute eye density
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
ed <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300))
# Plot the eye density
plot(ed)
```
