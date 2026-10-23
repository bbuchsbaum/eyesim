# Revisitation of previously fixated locations

Scores, fixation by fixation, whether gaze returns to a location already
fixated earlier in the same viewing episode ("revisitation"; Voss et
al., 2011; Kragel et al., 2021), and summarises each episode as counts,
rates and, optionally, rates expected under a spatial null model.

## Usage

``` r
fixation_revisits(x, ...)

# S3 method for class 'fixation_group'
fixation_revisits(
  x,
  radius,
  min_lag = 2,
  max_lag = Inf,
  require_departure = TRUE,
  reference = NULL,
  dist_method = c("euclidean", "manhattan"),
  ...
)
```

## Arguments

- x:

  A `fixation_group` or an `eye_table`.

- ...:

  Further arguments passed to methods.

- radius:

  Distance within which two fixations count as the same location.

- min_lag:

  Minimum number of fixations back for an eligible earlier fixation
  (default 2).

- max_lag:

  Maximum number of fixations back (default `Inf`).

- require_departure:

  Logical; require the preceding fixation to lie at least `radius` away
  (default `TRUE`).

- reference:

  Optional two-column matrix or data frame of `x`, `y` locations used
  for the spatial null model.

- dist_method:

  `"euclidean"` (default) or `"manhattan"`.

## Value

For a `fixation_group`, a tibble with one row per fixation in onset
order: `index`, `onset`, `opportunity`, `revisit`, `lag` (fixations back
to the most recent matched fixation), `match_index`, `distance` (to the
nearest eligible earlier fixation) and `expected` (null probability of a
revisit, or `NA` without `reference`). Use
[`revisit_summary()`](https://bbuchsbaum.github.io/eyesim/reference/revisit_summary.md)
for one row per episode. For an `eye_table`, see
`fixation_revisits.eye_table`.

## Details

Fixation \\i\\ is a *revisit* when it lies within `radius` of at least
one earlier fixation \\j\\ that is at least `min_lag` fixations back
(and at most `max_lag`). The default `min_lag = 2` means a fixation
never revisits the fixation immediately before it. With
`require_departure = TRUE` (default) gaze must also have left: the
preceding fixation must lie at least `radius` away from fixation \\i\\,
so a run of small refixations on one spot is not counted as a return. A
fixation is an *opportunity* when it could have been a revisit: it has
an eligible earlier fixation and, under `require_departure`, arrives
from at least `radius` away.

Revisit counts rise mechanically when an observer makes many fixations
within a small area, so raw counts partly restate fixation number and
spatial dispersion. Supplying `reference` locations gives a spatial
null: for each opportunity, the expected probability of a revisit is the
proportion of reference locations that would count as a revisit of the
same eligible earlier fixations, among those reference locations that
also satisfy the departure condition. A natural reference is the
observer's fixations from other trials, which preserves the observer's
spatial bias but not the order within the episode (see the `eye_table`
method, `null = "subject"`).

`radius` is in the units of the fixation coordinates. Kragel et al.
(2021) used about 2 degrees of visual angle; convert degrees to pixels
for your display geometry before calling.

## References

Voss, J. L., Warren, D. E., Gonsalves, B. D., Federmeier, K. D., Tranel,
D., & Cohen, N. J. (2011). Spontaneous revisitation during visual
exploration as a link among strategic behavior, learning, and the
hippocampus. *PNAS*, 108(31), E402-E409.

Kragel, J. E., Schuele, S., VanHaerents, S., Rosenow, J. M., & Voss, J.
L. (2021). Rapid coordination of effective learning by the human
hippocampus. *Science Advances*, 7(25), eabf7144.

## Examples

``` r
# A-B-C-A: the fourth fixation returns to the first location
fg <- fixation_group(x = c(100, 400, 400, 105), y = c(100, 100, 400, 98),
                     onset = c(0, 250, 500, 750), duration = rep(200, 4))
rv <- fixation_revisits(fg, radius = 50)
rv
#> Revisits: 1 of 2 opportunities (4 fixations; radius 50, lag 2-Inf, departure required)
#> # A tibble: 4 × 8
#>   index onset opportunity revisit   lag match_index distance expected
#>   <int> <dbl> <lgl>       <lgl>   <int>       <int>    <dbl>    <dbl>
#> 1     1     0 FALSE       FALSE      NA          NA    NA          NA
#> 2     2   250 FALSE       FALSE      NA          NA    NA          NA
#> 3     3   500 TRUE        FALSE      NA          NA   424.         NA
#> 4     4   750 TRUE        TRUE        3           1     5.39       NA
revisit_summary(rv)
#> # A tibble: 1 × 8
#>   n_fix opportunities revisits revisit_rate expected_revisits expected_rate
#>   <int>         <int>    <int>        <dbl>             <dbl>         <dbl>
#> 1     4             2        1          0.5                NA            NA
#> # ℹ 2 more variables: excess_revisits <dbl>, excess_rate <dbl>
```
