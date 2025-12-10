# Template Regression

This function performs a template regression on the provided reference
and source tables to estimate the beta weights for the baseline and
source maps.

## Usage

``` r
template_regression(
  ref_tab,
  source_tab,
  match_on,
  baseline_tab,
  baseline_key,
  method = c("lm", "rlm", "rank")
)
```

## Arguments

- ref_tab:

  A data frame containing the reference maps.

- source_tab:

  A data frame containing the source maps.

- match_on:

  A character string specifying the column name to be used for matching
  between the reference and source tables.

- baseline_tab:

  A data frame containing the baseline maps.

- baseline_key:

  A character string specifying the column name to be used for matching
  between the baseline table and the source table.

- method:

  A character vector of available regression methods. Default is c("lm",
  "rlm", "rank"). The selected method will be used for the regression
  analysis. - "lm": Linear regression (default). - "rlm": Robust linear
  regression. - "rank": Rank-based correlation.

## Value

A data frame with the source table augmented with two new columns,
beta_baseline and beta_source, representing the estimated beta weights
for the baseline and source maps, respectively.
