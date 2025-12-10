# Center Eye-Movements in a New Coordinate System

This function centers the eye-movements in a new coordinate system with
the specified origin.

## Usage

``` r
center(x, origin, ...)
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
if (FALSE) { # \dontrun{
# Example usage of the center function
new_origin <- c(500, 500) # New coordinate system origin
centered_object <- center(input_object, origin = new_origin)
} # }
```
