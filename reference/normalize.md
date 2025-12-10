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
if (FALSE) { # \dontrun{
# Example usage of the normalize function
# input_object <- eye-movement data object
x_bounds <- c(0, 1000) # X bounds for normalization
y_bounds <- c(0, 1000) # Y bounds for normalization
normalized_object <- normalize(input_object, xbounds = x_bounds, ybounds = y_bounds)
} # }
```
