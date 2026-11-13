# Summarize case minus mean-control reference scores

Summarize case minus mean-control reference scores

## Usage

``` r
case_control_contrast(
  scored,
  by = c("fixation", "trial"),
  subject = NULL,
  trial = "trial",
  columns = attr(scored, "score_columns")
)
```

## Arguments

- scored:

  A table returned by
  [`score_strata()`](https://bbuchsbaum.github.io/eyesim/reference/score_strata.md).

- by:

  Return one contrast per fixation/control draw or an equally weighted
  mean of fixation contrasts per subject/trial/control draw.

- subject, trial:

  Subject and trial column names. Subject defaults to the reference
  specification when available.

- columns:

  Score columns, defaulting to the score metadata.

## Value

A data frame with one contrast per reference. Control draws remain
separate. Missing scores propagate rather than silently changing the
pool.
