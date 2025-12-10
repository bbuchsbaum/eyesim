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

## Value

A scanpath object.

## Examples

``` r
# Create a fixation group
fixgroup <- data.frame(x = 1:5, y = 6:10)
# Create a scanpath for the fixation group
scanpath_obj <- scanpath.fixation_group(fixgroup)
#> Error in scanpath.fixation_group(fixgroup): could not find function "scanpath.fixation_group"
```
