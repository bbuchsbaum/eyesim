# Run similarity analysis for fixation data

This function compares the similarity between each fixation group in a
`source_tab` with a matching fixation group in `ref_tab` using the
specified similarity metric. Optionally, permutation tests can be
performed for assessing the significance of similarity values.

## Usage

``` r
run_similarity_analysis(
  ref_tab,
  source_tab,
  match_on,
  permutations,
  permute_on = NULL,
  method,
  refvar,
  sourcevar,
  window = NULL,
  multiscale_aggregation = "mean",
  ...
)
```

## Arguments

- ref_tab:

  A data frame containing the reference fixation groups.

- source_tab:

  A data frame containing the source fixation groups to be compared with
  the reference fixation groups.

- match_on:

  A column name in both `ref_tab` and `source_tab` used for matching the
  fixation groups.

- permutations:

  The number of permutations to perform for assessing the significance
  of similarity values (default: 0, no permutation tests).

- permute_on:

  An optional column name for limiting the matching indices in
  permutation tests (default: NULL).

- method:

  The similarity metric to use for comparing fixation groups (e.g.,
  "sinkhorn", "overlap").

- refvar:

  A column name in `ref_tab` containing the reference fixation groups.

- sourcevar:

  A column name in `source_tab` containing the source fixation groups.

- window:

  An optional numeric vector specifying the temporal window for
  computing similarity (default: NULL).

- ...:

  Extra arguments passed to the similarity function.
