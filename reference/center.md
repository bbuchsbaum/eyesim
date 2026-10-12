# Center Eye-Movements in a New Coordinate System

This function centers the eye-movements in a new coordinate system with
the specified origin.

## Usage

``` r
center(x, origin, ...)

# S3 method for class 'fixation_group'
center(x, origin = NULL, ...)
```

## Arguments

- x:

  The input object containing the eye-movements to be centered.

- origin:

  A vector containing the x and y coordinates of the new coordinate
  system's origin.

- ...:

  Additional arguments to be passed to the centering method.

## Value

An object containing the centered eye-movements in the new coordinate
system.

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
# Centre on the mean fixation position (the default) ...
center(fg)
#> Fixation group: 3 fixations
#>   x range: [-100, 100]
#>   y range: [-50, 50]
#>   duration range: [200, 200]
#>   onset range: [0, 400]
# ... or on a chosen origin, such as the centre of an 800 x 600 screen.
center(fg, origin = c(400, 300))
#> Fixation group: 3 fixations
#>   x range: [-300, -100]
#>   y range: [-200, -100]
#>   duration range: [200, 200]
#>   onset range: [0, 400]
```
