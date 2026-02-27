# Add Scanpath to a Data Frame

This function adds a scanpath to a data frame.

## Usage

``` r
# S3 method for class 'data.frame'
add_scanpath(x, outvar = "scanpath", fixvar = "fixgroup", ...)
```

## Arguments

- x:

  A data frame.

- outvar:

  The output variable name for the scanpath. Defaults to "scanpath".

- fixvar:

  The fixation group variable name. Defaults to "fixgroup".

- ...:

  Additional arguments (currently unused).

## Value

A data frame with the added scanpath.

## Examples

``` r
fg <- fixation_group(x = 1:5, y = 6:10, duration = rep(0.2, 5), onset = 1:5)
df <- tibble::tibble(fixgroup = list(fg))
df <- add_scanpath(df)
```
