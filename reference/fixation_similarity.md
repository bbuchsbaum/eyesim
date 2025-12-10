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

- Candidate sets are defined by `permute_on`; sampling is without
  replacement when `permutations` is smaller than the number of
  candidates.

- When `permutations` is greater than or equal to the available
  non-matching candidates, all candidates are used (exhaustive
  baseline).

- When permutations are requested, the result includes `eye_sim`,
  `perm_sim` (mean permuted similarity), and
  `eye_sim_diff = eye_sim - perm_sim`, all on the scale of `method`. If
  `method = "fisherz"`, convert to correlations via
  [`tanh()`](https://rdrr.io/r/base/Hyperbolic.html) if desired.

## Examples

``` r
# Example usage of the fixation_similarity function
ref_table <- # reference table data
source_table <- # source table data
match_column <- # column name to match fixation groups
similarity_results <- fixation_similarity(ref_table, source_table, match_column)
#> fixation_similarity: similarity metric is sinkhornoverlap
#> Error: object 'source_table' not found
```
