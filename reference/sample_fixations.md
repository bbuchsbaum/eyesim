# sample_fixations

Sample fixation coordinates at arbitrary time points.

## Usage

``` r
sample_fixations(x, time, ...)

# S3 method for class 'fixation_group'
sample_fixations(x, time, fast = TRUE, ...)
```

## Arguments

- x:

  the fixation group

- time:

  the continuous time points to sample

- ...:

  Additional arguments passed to methods.

- fast:

  Logical. If TRUE (default), uses a vectorized lookup; if FALSE,
  evaluates each time point in turn. Both paths first order the
  fixations by onset (keeping the input order among tied onsets) and
  drop fixations with a missing onset, so they return the same
  coordinates for any input order; of tied onsets, the last one is used.

## Value

A data frame of class `sampled_fixation_group` with one row per time
point and columns `x`, `y` (`NA` before the first onset), `onset` (the
time point) and `duration` (1).

## Details

Each time point takes the coordinates of the most recent fixation whose
onset is at or before it. Fixation durations are not used: gaze during a
saccade is attributed to the preceding fixation, and the last fixation
is held for all later time points. Time points before the first onset
give `NA`.

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
sample_fixations(fg, time = c(-50, 0, 250, 500))
#> Fixation group: 4 fixations
#>   x range: [100, 300]
#>   y range: [100, 200]
#>   duration range: [1, 1]
#>   onset range: [-50, 500]
```
