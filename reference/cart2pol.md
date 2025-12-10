# Convert Cartesian Coordinates to Polar Coordinates

This function converts Cartesian coordinates (x, y) to polar coordinates
(rho, theta).

## Usage

``` r
cart2pol(x, y)
```

## Arguments

- x:

  A numeric vector representing the x-coordinates.

- y:

  A numeric vector representing the y-coordinates.

## Value

A matrix with two columns, where the first column is rho (the radial
coordinate) and the second column is theta (the angular coordinate).

## Examples

``` r
cart2pol(c(1, 2), c(2, 2))
#> Error in cart2pol(c(1, 2), c(2, 2)): could not find function "cart2pol"
```
