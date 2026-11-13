# Arithmetic on eye density maps

Addition returns the mean of two maps. Subtraction returns their signed
difference; division returns their log ratio (zeros can yield infinite
or undefined cells). Difference and log-ratio surfaces have class
\`eye_contrast\`, not \`density\`, and cannot be passed to density-only
methods such as entropy. Maps must share grid coordinates and smoothing
bandwidth. Smoothing metadata is retained; the estimator of each operand
is recorded in \`source_estimators\`.

## Usage

``` r
# S3 method for class 'eye_density'
Ops(e1, e2)
```

## Arguments

- e1, e2:

  Eye density maps on the same grid and at the same bandwidth.

## Value

A density for addition, or an \`eye_contrast\` for subtraction/division.
