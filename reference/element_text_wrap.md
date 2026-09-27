# Wrapping text element

Like \[ggplot2::element_text()\], but the text is wrapped to the width
available when the plot is drawn. \[theme_eyesim()\] uses it for the
plot title, subtitle and caption.

## Usage

``` r
element_text_wrap(...)
```

## Arguments

- ...:

  Arguments passed to \[ggplot2::element_text()\].

## Value

A theme element.

## See also

Other visualization:
[`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md),
[`eyesim_colours()`](https://bbuchsbaum.github.io/eyesim/reference/eyesim_colours.md),
[`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md),
[`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md),
[`scale_fill_eyesim_density()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md),
[`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() +
  labs(caption = paste(rep("A long caption that wraps.", 8), collapse = " ")) +
  theme(plot.caption = element_text_wrap(hjust = 0))
```
