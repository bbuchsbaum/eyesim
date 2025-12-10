# Convert an eye_density object to a data.frame.

This function converts an eye_density object into a data.frame with x,
y, and z values.

## Usage

``` r
# S3 method for class 'eye_density'
as.data.frame(x, ...)
```

## Arguments

- x:

  An eye_density object to be converted into a data.frame.

- ...:

  Additional arguments passed to the method (currently not used).

## Value

A data.frame with columns x, y, and z representing the x-axis, y-axis,
and density values, respectively.

## Details

The function extracts the x and y values from the eye_density object,
then creates a data.frame with all possible combinations of x and y
using purrr::cross_df(). It then adds a new column 'z' to the data.frame
with the density values from the eye_density object.
