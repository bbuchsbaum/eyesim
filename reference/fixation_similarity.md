# Fixation Similarity

Compute the similarity between each fixation group in a `source_tab` and
a matching fixation group in `ref_tab`.

## Usage

``` r
fixation_similarity(
  ref_tab,
  source_tab,
  match_on,
  permutations = 0,
  permute_on = NULL,
  method = c("sinkhorn", "overlap"),
  refvar = "fixgroup",
  sourcevar = "fixgroup",
  window = NULL,
  ...
)
```

## Arguments

- ref_tab:

  The reference table containing the fixation groups to compare.

- source_tab:

  The source table containing the fixation groups to compare.

- match_on:

  The column name in both tables used to match fixation groups.

- permutations:

  The number of permutations to perform for permutation tests (default
  is 0, no permutations).

- permute_on:

  The column name on which to permute for permutation tests (default is
  NULL).

- method:

  The similarity metric to use; options are "sinkhorn" and "overlap"
  (default is "sinkhorn").

- refvar:

  The name of the column containing fixation groups in the reference
  table (default is "fixgroup").

- sourcevar:

  The name of the column containing fixation groups in the source table
  (default is "fixgroup").

- window:

  The temporal window over which to compute similarity (default is
  NULL).

- ...:

  Additional arguments to pass to the similarity metric function.

## Value

A table containing the computed similarities between fixation groups.

## Details

Permutation handling and units follow `template_similarity`:

- Candidate sets are defined by `permute_on`. Every copy of the true
  match is removed first; sampling is then without replacement when
  `permutations` is smaller than the number of remaining candidates.

- When `permutations` is greater than or equal to the available
  non-matching candidates, all candidates are used (exhaustive
  baseline).

- When permutations are requested, the result includes `eye_sim`,
  `perm_sim` (mean permuted similarity),
  `eye_sim_diff = eye_sim - perm_sim` (all on the scale of `method`),
  and `n_perm` (the number of permuted comparisons that contributed to
  `perm_sim` for that row; `0` when no baseline could be computed). If
  `method = "fisherz"`, convert to correlations via
  [`tanh()`](https://rdrr.io/r/base/Hyperbolic.html) if desired.

For `method = "overlap"`, pass `time_samples` through `...`; the
distance threshold `dthresh` defaults to 60, as in
[`fixation_overlap`](https://bbuchsbaum.github.io/eyesim/reference/fixation_overlap.md).

## Examples

``` r
set.seed(1)
make_tab <- function(jitter) {
  tibble::tibble(trial = 1:3, fixgroup = lapply(1:3, function(i) {
    fixation_group(x = c(100, 200, 300) * i + rnorm(3, 0, jitter),
                   y = c(100, 150, 200) + rnorm(3, 0, jitter),
                   onset = c(0, 200, 400), duration = c(200, 200, 200))
  }))
}
ref <- make_tab(0)
src <- make_tab(20)
fixation_similarity(ref, src, match_on = "trial", method = "overlap",
                    time_samples = seq(0, 600, by = 20), permutations = 2)
#> fixation_similarity: similarity metric is overlap
#> # A tibble: 3 × 7
#>   trial fixgroup           matchind eye_sim perm_sim eye_sim_diff n_perm
#>   <int> <list>                <int>   <dbl>    <dbl>        <dbl>  <int>
#> 1     1 <fxtn_grp [3 × 6]>        1       1        0            1      2
#> 2     2 <fxtn_grp [3 × 6]>        2       1        0            1      2
#> 3     3 <fxtn_grp [3 × 6]>        3       1        0            1      2
```
