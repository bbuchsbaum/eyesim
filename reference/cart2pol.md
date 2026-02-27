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
#>           rho     theta
#> [1,] 2.236068 1.1071487
#> [2,] 2.828427 0.7853982
```
