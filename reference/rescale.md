# rescale

Rescale spatial coordinates.

## Usage

``` r
rescale(x, sx, sy)

# S3 method for class 'fixation_group'
rescale(x, sx, sy)
```

## Arguments

- x:

  An object containing spatial coordinates to be rescaled.

- sx:

  A numeric value representing the x-axis scale factor.

- sy:

  A numeric value representing the y-axis scale factor.

## Value

An object with rescaled spatial coordinates.

## Details

This function rescales the spatial coordinates of an object by the given
scale factors.

## See also

`rescale.fixation_group` for an example of a specific method
implementation for this generic.
