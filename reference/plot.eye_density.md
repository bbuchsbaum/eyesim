# Plot Eye Density

This function creates a plot of the eye density using ggplot2.

## Usage

``` r
# S3 method for class 'eye_density'
plot(
  x,
  alpha = 0.8,
  bg_image = NULL,
  transform = c("identity", "sqroot", "curoot", "rank"),
  ...
)
```

## Arguments

- x:

  An "eye_density" object.

- alpha:

  The opacity level for the raster layer (default: 0.8).

- bg_image:

  An optional image file name to use as the background.

- transform:

  The transformation to apply to the density values (default:
  c("identity", "sqroot", "curoot", "rank")).

- ...:

  Additional args

## Value

A ggplot object representing the eye density plot.

## See also

Other visualization:
[`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md),
[`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md)

## Examples

``` r
# Assume `ed` is an "eye_density" object
# Plot the eye density
plot_eye_density <- plot.eye_density(ed)
#> Error in plot.eye_density(ed): could not find function "plot.eye_density"
```
