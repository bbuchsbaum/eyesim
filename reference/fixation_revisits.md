# Revisitation of previously fixated locations

Scores, fixation by fixation, whether gaze returns to a location already
fixated earlier in the same viewing episode, and summarises each episode
as counts, rates and, optionally, rates expected under a spatial
reference distribution. This is a spatial-return measure related to
Kragel et al. (2021), not the backtracking-sequence classification of
Voss et al. (2011), which classified A-B-C-A as "other" rather than
revisitation.

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

  Maximum number of fixations back, a whole number at least `min_lag`,
  or `Inf` (default) for no limit.

- require_departure:

  Logical; require the preceding fixation to lie at least `radius` from
  the current fixation (default `TRUE`).

- reference:

  Optional two-column matrix or data frame of `x`, `y` finite numeric
  locations used for the spatial reference distribution.

- dist_method:

  `"euclidean"` (default) or `"manhattan"`.

## Value

For a `fixation_group`, a tibble with one row per fixation in onset
order: `index`, `onset`, `opportunity`, `revisit`, `lag` (fixations back
to the most recent matched fixation), `match_index` (that fixation's
original index), `nearest_distance` (to the nearest eligible earlier
fixation), `match_distance` (to `match_index`, which need not be
nearest), `distance` (a compatibility alias for `nearest_distance`) and
`expected` (null probability of a revisit, or `NA` without reference
support). Distances are `NA` outside opportunities; `match_distance` is
also `NA` without a match. Attributes retain `radius`, `min_lag`,
`max_lag`, `require_departure`, `dist_method` and `null` (whether a
reference was supplied). Use
[`revisit_summary()`](https://bbuchsbaum.github.io/eyesim/reference/revisit_summary.md)
for one row per episode. For an `eye_table`, see
`fixation_revisits.eye_table`.

## Details

Fixation \\i\\ is a *revisit* when it lies within `radius` of at least
one earlier fixation \\j\\ that is at least `min_lag` fixations back
(and at most `max_lag`). The default `min_lag = 2` means a fixation
never revisits the fixation immediately before it. With
`require_departure = TRUE` (default), the incoming step must be at least
`radius`: the preceding fixation must lie at least `radius` away from
fixation \\i\\. This does not require leaving and re-entering a region
centred on the earlier fixation. For example, with radius 10, positions
0, -9, 9 on one axis count as a return, although both later positions
are within radius of the original position. A fixation is an
*opportunity* when it could have been a revisit: it has an eligible
earlier fixation and, under `require_departure`, arrives from at least
`radius` away.

Revisit counts rise mechanically when an observer makes many fixations
within a small area, so raw counts partly restate fixation number and
spatial dispersion. Supplying `reference` locations gives a spatial
null: for each opportunity, the expected probability of a revisit is the
proportion of reference locations that would count as a revisit of the
same eligible earlier fixations, among those reference locations that
also satisfy the departure condition. A natural reference is the
observer's fixations from other trials, which preserves the observer's
spatial bias but not the order within the episode (see the `eye_table`
method, `null = "subject"`). Observed minus expected is a descriptive
correction under this reference distribution, not an inferential test or
a generative model of scanpaths.

Coordinates and onsets must be finite numeric values. Invalid values
cause an error; dropping fixations would change adjacency, lags and
denominators. Fixations are scored in onset order, preserving input
order for tied onsets. Lag bounds are whole numbers; large finite bounds
are supported without conversion to 32-bit integers.

Reference coverage is updated as fixations enter and leave the eligible
lag window. For \\n\\ fixations, \\m\\ reference points and maximum
eligible history size \\L\\, scoring uses \\O(n + nL + nm)\\ time and
\\O(n + m)\\ storage, including result rows. No pairwise distance matrix
is allocated.

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
#> # A tibble: 4 × 10
#>   index onset opportunity revisit   lag match_index distance nearest_distance
#>   <int> <dbl> <lgl>       <lgl>   <int>       <int>    <dbl>            <dbl>
#> 1     1     0 FALSE       FALSE      NA          NA    NA               NA   
#> 2     2   250 FALSE       FALSE      NA          NA    NA               NA   
#> 3     3   500 TRUE        FALSE      NA          NA   424.             424.  
#> 4     4   750 TRUE        TRUE        3           1     5.39             5.39
#> # ℹ 2 more variables: match_distance <dbl>, expected <dbl>
revisit_summary(rv)
#> # A tibble: 1 × 8
#>   n_fix opportunities revisits revisit_rate expected_revisits expected_rate
#>   <int>         <int>    <int>        <dbl>             <dbl>         <dbl>
#> 1     4             2        1          0.5                NA            NA
#> # ℹ 2 more variables: excess_revisits <dbl>, excess_rate <dbl>
```
