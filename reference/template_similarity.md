# template_similarity

Compute similarity between each density map in a `source_tab` with a
matching ("template") density map in `ref_tab`.

## Usage

``` r
template_similarity(
  ref_tab,
  source_tab,
  match_on,
  permute_on = NULL,
  refvar = "density",
  sourcevar = "density",
  method = c("spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov", "emd"),
  permutations = 10,
  multiscale_aggregation = "mean",
  ...
)
```

## Arguments

- ref_tab:

  A data frame or tibble containing reference density maps.

- source_tab:

  A data frame or tibble containing source density maps.

- match_on:

  A character string representing the variable used to match density
  maps between `ref_tab` and `source_tab`.

- permute_on:

  A character string representing the variable used to stratify
  permutations (default is NULL).

- refvar:

  A character string representing the name of the variable containing
  density maps in the reference table (default is "density").

- sourcevar:

  A character string representing the name of the variable containing
  density maps in the source table (default is "density").

- method:

  A character string specifying the similarity method to use. Possible
  values are "spearman", "pearson", "fisherz", "cosine", "l1",
  "jaccard", and "dcov" (default is "spearman").

- permutations:

  A numeric value specifying the number of permutations for the baseline
  map (default is 10).

- multiscale_aggregation:

  If the density maps are multiscale (i.e., \`eye_density_multiscale\`
  objects), this specifies how to aggregate similarities from different
  scales. Options: "mean" (default, returns the average similarity
  across scales), "none" (returns a list or vector of similarities, one
  per scale, within the result columns). See
  \`similarity.eye_density_multiscale\`.

- ...:

  Extra arguments to pass to the \`similarity\` function.

## Value

A data frame or tibble containing the source table and additional
columns with the similarity scores and permutation results.

## Details

Permutation baseline and exhaustive behavior:

- The set of permutation candidates is determined by `permute_on`. If
  `permute_on` is provided, candidates are restricted within that
  stratum (e.g., within-participant); otherwise all reference items are
  candidates.

- If `permutations` is less than the number of available non-matching
  candidates, a random subset of that size is drawn (without
  replacement) for each trial. Internally, sampling is performed with a
  fixed future seed to aid reproducibility.

- If `permutations` is greater than or equal to the number of available
  non-matching candidates, the procedure uses all candidates (excluding
  the true match). In other words, the permutation baseline is
  exhaustive when possible.

- For small-N designs, you can set `permutations` to a large number to
  trigger exhaustive behavior. For example, with 3 images per
  participant and `permute_on = participant`, there are only 2
  non-matching candidates per trial; any `permutations >= 2` will result
  in using both.

Returned columns and units:

- `eye_sim`: the observed similarity for the matched pair, on the scale
  of `method`.

- `perm_sim`: the mean similarity across permuted non-matching pairs
  (same scale as `eye_sim`).

- `eye_sim_diff`: `eye_sim - perm_sim`. Units match `method`.

Notes on `method` and interpretation:

- If `method = "fisherz"`, values are Fisher z (atanh of Pearson *r*).
  Convert back to *r* via `tanh(z)` for reporting on the correlation
  scale.

- If `method = "pearson"` or `"spearman"`, values are correlations
  (roughly in `[-1, 1]`).

- Other methods (e.g., `"emd"`, `"cosine"`) produce scores on their
  respective scales.
