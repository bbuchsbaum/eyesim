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
  similarity_transform = NULL,
  similarity_transform_args = list(),
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

- similarity_transform:

  Optional preprocessing hook applied before similarity is computed.
  Should be a function that accepts (`ref_tab`, `source_tab`,
  `match_on`, `refvar`, `sourcevar`) and returns a list with updated
  tables/column names. See \`latent_pca_transform\`,
  \`coral_transform\`, and \`cca_transform\`.

- similarity_transform_args:

  Named list of extra arguments passed to \`similarity_transform\`.

- ...:

  Extra arguments to pass to the \`similarity\` function.

## Value

A data frame or tibble containing the source table and additional
columns with the similarity scores and permutation results.

## Details

Permutation baseline and exhaustive behavior:

- The set of permutation candidates is determined by `permute_on`.
  Candidates are the distinct reference rows matched by at least one
  source row, within the same `permute_on` stratum if given (e.g.,
  within-participant). Each candidate counts once, however many source
  rows match it, so `n_perm` counts distinct templates. A reference row
  that no source row matches is never a candidate.

- The true match is removed from the candidate set before any sampling.
  If several source rows share the same `match_on` key, every copy of
  that key is removed, so a row is never compared with its own template
  in the baseline.

- If `permutations` is less than the number of available non-matching
  candidates, a random subset of that size is drawn (without
  replacement) for each trial.

- Sampling uses the session random number generator; there is no
  internal fixed seed. Call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) immediately before
  the call to make the baseline reproducible. The call advances the
  session RNG. With `method = "cosine"`, the default
  `multiscale_aggregation = "mean"`, no `window` or extra arguments, and
  every reference and source map on one lattice, a vectorized path
  samples with [`sample()`](https://rdrr.io/r/base/sample.html)
  directly; other methods draw per-row streams through
  `furrr::furrr_options(seed = TRUE)`, which are derived from the
  session RNG. The same seed can therefore select different controls for
  different methods.

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

- `n_perm`: the number of permuted non-matching comparisons that
  contributed to `perm_sim` for that row. This varies across rows (e.g.,
  small `permute_on` strata, or fewer candidates than requested) and is
  `0` when no baseline could be computed (`perm_sim = NA`). Use it to
  drop rows with too few permutations, e.g.
  `dplyr::filter(res, n_perm >= k)`.

Notes on `method` and interpretation:

- If `method = "fisherz"`, values are Fisher z (atanh of Pearson *r*).
  Convert back to *r* via `tanh(z)` for reporting on the correlation
  scale. To keep z finite, *r* is clamped to
  `[-1 + .Machine$double.eps, 1 - .Machine$double.eps]`. Identical maps,
  whether constant or not, and any *r* within `64 * .Machine$double.eps`
  of 1 (for example a rescaled copy of a map) are treated as *r* = 1 and
  give `atanh(1 - .Machine$double.eps)` (about 18.37) rather than `Inf`;
  *r* near -1 is handled symmetrically.

- If `method = "pearson"` or `"spearman"`, values are correlations
  (roughly in `[-1, 1]`).

- Other methods (e.g., `"emd"`, `"cosine"`) produce scores on their
  respective scales.
