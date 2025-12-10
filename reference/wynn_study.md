# Eye-tracking study data from Wynn et al.

Eye-tracking fixation data from a study phase experiment.

## Usage

``` r
wynn_study
```

## Format

A data frame or tibble containing eye-tracking fixation data with
columns for x/y coordinates, fixation duration, onset times, and
subject/trial identifiers.

## Source

Wynn et al. eye-tracking study.

## Examples

``` r
data(wynn_study)
head(wynn_study)
#> # A tibble: 6 × 8
#> # Groups:   ImageNumber, Subject [6]
#>   Subject fixgroup            ImageVersion Version Trial ImageNumber blockstudy
#>     <dbl> <list>              <chr>        <chr>   <int> <fct>            <int>
#> 1       5 <fxtn_grp [39 × 6]> 1_A          A          18 1                    1
#> 2       6 <fxtn_grp [42 × 6]> 1_A          A          13 1                    1
#> 3       7 <fxtn_grp [26 × 6]> 1_A          A          25 1                    1
#> 4       8 <fxtn_grp [30 × 6]> 1_A          A          11 1                    1
#> 5       9 <fxtn_grp [32 × 6]> 1_A          A          15 1                    1
#> 6      10 <fxtn_grp [41 × 6]> 1_A          A          12 1                    1
#> # ℹ 1 more variable: ImageNumberS <chr>
```
