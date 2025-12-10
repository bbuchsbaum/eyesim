# Add Scanpath to a Data Frame

This function adds a scanpath to a data frame.

## Usage

``` r
# S3 method for class 'data.frame'
add_scanpath(x, outvar = "scanpath", fixvar = "fixgroup")
```

## Arguments

- x:

  A data frame.

- outvar:

  The output variable name for the scanpath. Defaults to "scanpath".

- fixvar:

  The fixation group variable name. Defaults to "fixgroup".

## Value

A data frame with the added scanpath.

## Examples

``` r
# Create a data frame with a fixation group
df <- data.frame(x = 1:5, y = 6:10, fixgroup = rep(1, 5))
# Add a scanpath to the data frame
df <- add_scanpath.data.frame(df)
#> Error in add_scanpath.data.frame(df): could not find function "add_scanpath.data.frame"
```
