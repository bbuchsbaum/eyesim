# Revisits for every episode of an eye_table

Applies
[`fixation_revisits`](https://bbuchsbaum.github.io/eyesim/reference/fixation_revisits.md)
to each row (episode) of an `eye_table` and appends per-episode summary
columns.

## Usage

``` r
# S3 method for class 'eye_table'
fixation_revisits(
  x,
  radius,
  min_lag = 2,
  max_lag = Inf,
  require_departure = TRUE,
  null = c("none", "subject", "item"),
  subject_var = NULL,
  item_var = NULL,
  var = "fixgroup",
  keep_fixations = FALSE,
  dist_method = c("euclidean", "manhattan"),
  ...
)
```

## Arguments

- x:

  An `eye_table`.

- radius, min_lag, max_lag, require_departure, dist_method:

  See
  [`fixation_revisits`](https://bbuchsbaum.github.io/eyesim/reference/fixation_revisits.md).

- null:

  The spatial null model:

  `"none"`

  :   (default) no null; expected columns are `NA`.

  `"subject"`

  :   every fixation of the same subject from the *other* rows of the
      table (leave-episode-out). Preserves the observer's spatial bias
      but not the content of the viewed image.

  `"item"`

  :   every fixation from other rows showing the same item (`item_var`).
      Preserves where the image draws gaze, but not the observer's own
      bias. If `subject_var` is given, all of the same subject's rows
      are excluded, which matters when an observer sees an item more
      than once.

- subject_var:

  Column identifying the subject; required when `null = "subject"`,
  optional for `null = "item"`.

- item_var:

  Column identifying the viewed item (image); required when
  `null = "item"`.

- var:

  Name of the fixation-group column (default `"fixgroup"`).

- keep_fixations:

  Logical; keep the per-fixation tables in a list column
  `revisit_fixations` (default `FALSE`).

- ...:

  Ignored.

## Value

`x` with columns from
[`revisit_summary`](https://bbuchsbaum.github.io/eyesim/reference/revisit_summary.md)
(`n_fix`, `opportunities`, `revisits`, `revisit_rate`,
`expected_revisits`, `expected_rate`, `excess_revisits`, `excess_rate`).

## Examples

``` r
set.seed(1)
et <- simulate_eye_table(n_fixations = 120, n_groups = 6)
et$subject <- rep(c("s1", "s2"), each = 3)
et$image <- rep(c("a", "b", "c"), times = 2)
res <- fixation_revisits(et, radius = 150, null = "subject", subject_var = "subject")
tibble::as_tibble(res)[, c("subject", "n_fix", "revisits", "expected_revisits", "excess_revisits")]
#> # A tibble: 6 × 5
#>   subject n_fix revisits expected_revisits excess_revisits
#>   <chr>   <int>    <int>             <dbl>           <dbl>
#> 1 s1         20        5              5.41          -0.406
#> 2 s1         20        6              7.97          -1.97 
#> 3 s1         20        9              6.00           3.00 
#> 4 s2         20        5              6.52          -1.52 
#> 5 s2         20        6              4.78           1.22 
#> 6 s2         20        2              6.06          -4.06 
# same-image null from the other subject's viewing
res_item <- fixation_revisits(et, radius = 150, null = "item", item_var = "image", subject_var = "subject")
tibble::as_tibble(res_item)[, c("subject", "image", "revisits", "expected_revisits")]
#> # A tibble: 6 × 4
#>   subject image revisits expected_revisits
#>   <chr>   <chr>    <int>             <dbl>
#> 1 s1      a            5              5.32
#> 2 s1      b            6              6.12
#> 3 s1      c            9              5.67
#> 4 s2      a            5              5.69
#> 5 s2      b            6              4.82
#> 6 s2      c            2              7.08
```
