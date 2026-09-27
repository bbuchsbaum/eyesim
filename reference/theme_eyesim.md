# eyesim ggplot2 theme

\`theme_eyesim()\` is the house theme shared by all eyesim plots: plain
background, light major grid, left-aligned title block, and a small grey
caption for provenance notes. \`theme_eyesim_spatial()\` is the variant
used for stimulus-space plots (scanpaths, density maps): no axes or
grid, and a thin frame marking the stimulus extent.

## Usage

``` r
theme_eyesim(base_size = 10, base_family = "")

theme_eyesim_spatial(base_size = 10, base_family = "")
```

## Arguments

- base_size:

  Base font size in points.

- base_family:

  Base font family.

## Value

A \[ggplot2::theme()\] object.

## Details

Both return ordinary ggplot2 themes, so they can be added to any plot
and further modified with \[ggplot2::theme()\].

## See also

Other visualization:
[`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md),
[`element_text_wrap()`](https://bbuchsbaum.github.io/eyesim/reference/element_text_wrap.md),
[`eyesim_colours()`](https://bbuchsbaum.github.io/eyesim/reference/eyesim_colours.md),
[`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md),
[`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md),
[`scale_fill_eyesim_density()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md)

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_eyesim()
```
