# Reapply the 'eye_table' Class to an Object

The \`as_eye_table\` function is a simple utility function that
reapplies the 'eye_table' class to a given object if it is not already
part of its class attribute. This function is not intended for serious
use.

## Usage

``` r
as_eye_table(x)
```

## Arguments

- x:

  The input object to which the 'eye_table' class should be (re)applied.

## Value

The input object with the 'eye_table' class (re)applied to its class
attribute.

## Examples

``` r
# Create a simple data.frame
df <- data.frame(x = 1:5, y = 6:10)

# Apply the 'eye_table' class to df
eye_table_df <- as_eye_table(df)
class(eye_table_df)
#> [1] "eye_table"  "data.frame"
```
