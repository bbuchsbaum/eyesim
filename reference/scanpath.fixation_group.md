# Create a Scanpath for a Fixation Group

This function creates a scanpath for a fixation group.

## Usage

``` r
# S3 method for class 'fixation_group'
scanpath(x, ...)
```

## Arguments

- x:

  A fixation group object.

- ...:

  Additional arguments (currently unused).

## Value

A scanpath object.

## Examples

``` r
# Create a fixation group
fg <- fixation_group(x = 1:5, y = 6:10, duration = rep(0.2, 5), onset = 1:5)
# Create a scanpath for the fixation group using the S3 generic
scanpath_obj <- scanpath(fg)
```
