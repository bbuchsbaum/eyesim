# rep_fixations

Replicate a fixation sequence.

## Usage

``` r
rep_fixations(x, resolution)

# S3 method for class 'fixation_group'
rep_fixations(x, resolution = 100)
```

## Arguments

- x:

  An object representing a fixation sequence.

- resolution:

  A numeric value representing the temporal resolution of the replicated
  fixations.

## Value

An object containing the replicated fixation sequence with the specified
temporal resolution.

## Details

This function replicates a fixation sequence with a specified temporal
resolution. It can be useful when working with fixation data that needs
to be resampled or when creating fixation sequences with consistent
temporal spacing.

## See also

`rep_fixations.fixation_group` for an example of a specific method
implementation for this generic.
