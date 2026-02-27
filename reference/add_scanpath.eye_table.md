# Add Scanpath to an Eye Table

This function adds a scanpath to an eye table.

## Usage

``` r
# S3 method for class 'eye_table'
add_scanpath(x, outvar = "scanpath", fixvar = "fixgroup", ...)
```

## Arguments

- x:

  An eye table object.

- outvar:

  The output variable name for the scanpath. Defaults to "scanpath".

- fixvar:

  The fixation group variable name. Defaults to "fixgroup".

- ...:

  Additional arguments (currently unused).

## Value

An eye table object with the added scanpath.

## Examples

``` r
# Create an eye table with a fixation group
fg <- fixation_group(x = 1:5, y = 6:10, duration = rep(0.2, 5), onset = 1:5)
et <- as_eye_table(tibble::tibble(id = 1, fixgroup = list(fg)))
et <- add_scanpath(et)
```
