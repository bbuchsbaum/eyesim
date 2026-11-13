# Draw replayable mismatched-template assignments

Draw replayable mismatched-template assignments

## Usage

``` r
mismatch_assignments(strata, refs, n = 100, trial = "trial", seed = 1)
```

## Arguments

- strata:

  A case-control table from
  [`sample_controls()`](https://bbuchsbaum.github.io/eyesim/reference/sample_controls.md).

- refs:

  A reference set.

- n:

  Number of null draws.

- trial:

  Trial column, nested within the reference subject column.

- seed:

  Seed, or NULL to use the session RNG.

## Value

A data frame containing original subject/item/trial keys, `null_draw`,
and `replacement_item`. Each trial gets one other valid item from the
same subject per draw. Assignments can be saved and passed to
[`mismatch_draws()`](https://bbuchsbaum.github.io/eyesim/reference/mismatch_draws.md).
