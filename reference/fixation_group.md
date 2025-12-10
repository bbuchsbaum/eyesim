# Create a Fixation Group Object

This function creates a fixation group object containing fixation data
with x/y coordinates, duration, onset times, and an optional group
index.

## Usage

``` r
fixation_group(x, y, duration, onset, group = 0)
```

## Arguments

- x:

  A numeric vector of x-coordinates for each fixation.

- y:

  A numeric vector of y-coordinates for each fixation.

- duration:

  A numeric vector of fixation durations. If missing, computed from
  onset differences.

- onset:

  A numeric vector of fixation onset times.

- group:

  An optional group identifier (default is 0).

## Value

A tibble of class "fixation_group" with columns: index, x, y, duration,
onset, group_index.

## Examples

``` r
fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                     onset = c(0, 200, 400), duration = c(200, 200, 200))
```
