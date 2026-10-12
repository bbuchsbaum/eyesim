# Report comparator availability

Optional MultiMatch dependencies are checked explicitly. Density and the
built-in elastic-consensus reference implementation have no optional
dependency.

## Usage

``` r
gaze_baseline_availability(spec)
```

## Arguments

- spec:

  A \[gaze_baseline_spec()\].

## Value

A data frame with one row per requested comparator.

## Examples

``` r
spec <- gaze_baseline_spec(screen = gaze_screen(20, 10, unit = "deg"),
                           density_sigmas = c(1, 2))
gaze_baseline_availability(spec)
#>       method available reason
#> 1 multimatch      TRUE   <NA>
#> 2    density      TRUE   <NA>
#> 3    elastic      TRUE   <NA>
```
