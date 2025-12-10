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
df <- data.frame(x = 1:5, y = 6:10, fixgroup = rep(1, 5))
eye_table_df <- as_eye_table(df)
# Add a scanpath to the eye table
eye_table_df <- add_scanpath.eye_table(eye_table_df)
#> Error in add_scanpath.eye_table(eye_table_df): could not find function "add_scanpath.eye_table"
```
