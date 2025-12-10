# Subset an 'eye_table' Object

The \`\[.eye_table\` function is an S3 method for subsetting 'eye_table'
objects. It uses \`NextMethod()\` to perform the subsetting operation
and then reapplies the 'eye_table' class to the result using the
\`as_eye_table\` function. This ensures that the returned object still
has the 'eye_table' class after subsetting.

## Usage

``` r
# S3 method for class 'eye_table'
x[i, j, drop = FALSE]
```

## Arguments

- x:

  The 'eye_table' object to be subsetted.

- i:

  Row indices for subsetting.

- j:

  Column indices for subsetting.

- drop:

  A logical value indicating whether to drop dimensions that have only
  one level after subsetting. Defaults to \`FALSE\`.

## Value

An 'eye_table' object after subsetting.

## Examples

``` r
# Create an 'eye_table' object
df <- data.frame(x = 1:5, y = 6:10)
eye_table_df <- as_eye_table(df)

# Subset the 'eye_table' object
subset_eye_table <- eye_table_df[1:3,]
class(subset_eye_table)
#> [1] "eye_table"  "data.frame"
```
