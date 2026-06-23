# eyesim (development version)

* `template_similarity()`, `fixation_similarity()`, and `scanpath_similarity()`
  now return an `n_perm` column when `permutations > 0`. It records the number of
  permuted non-matching comparisons that actually contributed to `perm_sim` for
  each row, which varies with small `permute_on` strata or when fewer candidates
  than requested are available, and is `0` when no baseline could be computed
  (`perm_sim = NA`). Use it to exclude rows with too thin a permutation baseline,
  e.g. `dplyr::filter(res, n_perm >= k)`. The column is additive: output for
  `permutations = 0` is unchanged.
