# Template Multiple Regression

This function performs a multiple regression on the provided source
table with the specified response and covariates using the chosen
regression method.

## Usage

``` r
template_multireg(
  source_tab,
  response,
  covars,
  method = c("lm", "rlm", "nnls", "logistic"),
  intercept = TRUE
)
```

## Arguments

- source_tab:

  A data frame containing the source maps.

- response:

  A character string specifying the column name of the response
  variable.

- covars:

  A character vector specifying the column names of the covariates.

- method:

  A character vector of available regression methods. Default is c("lm",
  "rlm", "nnls", "logistic"). The selected method will be used for the
  regression analysis. - "lm": Linear regression (default). - "rlm":
  Robust linear regression. - "nnls": Non-negative least squares. -
  "logistic": Logistic regression.

- intercept:

  A logical value indicating whether to include an intercept term in the
  model. Default is TRUE.

## Value

A data frame with the source table augmented with a new column,
multireg, containing the results of the multiple regression analysis.
