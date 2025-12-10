# Construct a Scanpath of a Fixation Group of Related Objects

This function creates a scanpath containing polar coordinates (rho,
theta) along with absolute x and y spatial coordinates for a given
fixation group.

## Usage

``` r
scanpath(x, ...)
```

## Arguments

- x:

  The fixations.

- ...:

  Extra arguments.

## Value

A scanpath object.

## Examples

``` r
# Create a fixation group
fg <- fixation_group(x=c(.1,.5,1), y=c(1,.5,1), onset=1:3, duration=rep(1,3))
# Create a scanpath for the fixation group
sp <- scanpath(fg)
```
