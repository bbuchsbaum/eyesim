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
# Create a fixation group and compute eye density
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
ed <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300))
# Plot the eye density
plot(ed)
```
