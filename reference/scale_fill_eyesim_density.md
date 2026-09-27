# eyesim colour scales

\`scale_fill_eyesim_density()\` is the sequential density scale used by
the density plots. Opacity is part of the ramp and is anchored at zero
density, so zero is transparent (a background stimulus shows through)
and the colour bar shows exactly what is drawn. Set shared \`limits\` to
compare conditions on one scale. \`scale_colour_eyesim_time()\` is the
fixation-order/onset scale.

## Usage

``` r
scale_fill_eyesim_density(
  ...,
  name = "Fixation\ndensity",
  limits = c(0, NA),
  transparent = TRUE,
  aesthetics = "fill"
)

scale_colour_eyesim_time(..., aesthetics = "colour")
```

## Arguments

- ...:

  Passed to \[ggplot2::scale_fill_gradientn()\] or
  \[ggplot2::scale_colour_gradientn()\].

- name:

  Legend title.

- limits:

  Scale limits. The default starts the density scale at zero.

- transparent:

  Whether low densities fade to transparent.

- aesthetics:

  Aesthetics the scale applies to.

## Value

A ggplot2 scale.

## See also

Other visualization:
[`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md),
[`element_text_wrap()`](https://bbuchsbaum.github.io/eyesim/reference/element_text_wrap.md),
[`eyesim_colours()`](https://bbuchsbaum.github.io/eyesim/reference/eyesim_colours.md),
[`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md),
[`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md),
[`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
ed <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300))
# Put several maps on one scale by giving them the same limits
plot(ed, limits = c(0, 1e-4))
```
