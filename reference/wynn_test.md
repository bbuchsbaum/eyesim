# Eye-tracking test data from Wynn et al.

Eye-tracking fixation data from a test phase experiment.

## Usage

``` r
wynn_test
```

## Format

A data frame or tibble containing eye-tracking fixation data with
columns for x/y coordinates, fixation duration, onset times, and
subject/trial identifiers.

## Source

Wynn et al. eye-tracking study.

## Examples

``` r
data(wynn_test)
head(wynn_test)
#> # A tibble: 6 × 11
#> # Groups:   ImageNumber, Subject [6]
#>   Subject fixgroup   ImageVersion Saliency Accuracy Version Trial Duration
#>     <dbl> <list>     <chr>           <int>    <dbl> <chr>   <int>    <int>
#> 1       5 <fxtn_grp> 1_A                20        1 A         129      250
#> 2       6 <fxtn_grp> 1_A                20        1 A         140      500
#> 3       7 <fxtn_grp> 1_A                20        1 A         135      750
#> 4       8 <fxtn_grp> 1_A                40        0 A         143      250
#> 5       9 <fxtn_grp> 1_A                40        0 A         136      500
#> 6      10 <fxtn_grp> 1_A                60        0 A         124      750
#> # ℹ 3 more variables: ImageNumber <int>, ImageNumberS <chr>, probe_type <chr>
```
