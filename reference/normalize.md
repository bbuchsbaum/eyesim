# Normalize Eye-Movements to Unit Range

This function normalizes the eye-movements to a unit range based on the
specified x and y bounds.

## Usage

``` r
normalize(x, xbounds, ybounds, ...)

# S3 method for class 'fixation_group'
normalize(x, xbounds, ybounds, ...)
```

## Arguments

- x:

  The input object containing the eye-movements to be normalized.

- xbounds:

  A vector containing the minimum and maximum x bounds for
  normalization.

- ybounds:

  A vector containing the minimum and maximum y bounds for
  normalization.

- ...:

  Additional arguments to be passed to the normalization method.

## Value

An object containing the normalized eye-movements in the unit range.

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
normalize(fg, xbounds = c(0, 400), ybounds = c(0, 300))
#> Fixation group: 3 fixations
#>   x range: [0.2, 0.8]
#>   y range: [0.3, 0.7]
#>   duration range: [200, 200]
#>   onset range: [0, 400]
```
