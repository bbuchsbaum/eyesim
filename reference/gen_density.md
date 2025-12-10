# This function creates a density object from the provided x, y, and z matrices. The density object is a list containing the x, y, and z values with a class attribute set to "density" and "list".

This function creates a density object from the provided x, y, and z
matrices. The density object is a list containing the x, y, and z values
with a class attribute set to "density" and "list".

## Usage

``` r
gen_density(x, y, z)
```

## Arguments

- x:

  A numeric vector representing the x-axis values of the density map.

- y:

  A numeric vector representing the y-axis values of the density map.

- z:

  A matrix representing the density values at each (x, y) coordinate.

## Value

A density object which is a list containing the x, y, and z values with
a class attribute set to "density" and "list".

## Details

The function first checks if the dimensions of the z matrix are equal to
the length of the x and y vectors. If not, it throws an error. Then, it
creates a list containing the x, y, and z values and sets the class
attribute of the list to "density" and "list".
