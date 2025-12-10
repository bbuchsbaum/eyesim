# Match Scaling Parameters for Fixation Data

This function matches the scaling parameters of fixation data between a
reference table and a source table based on a common matching variable.

## Usage

``` r
match_scale(
  ref_tab,
  source_tab,
  match_on,
  refvar = "fixgroup",
  sourcevar = "fixgroup",
  window,
  ...
)
```

## Arguments

- ref_tab:

  A data frame containing the reference fixation data.

- source_tab:

  A data frame containing the source fixation data.

- match_on:

  A string specifying the variable name in both \`ref_tab\` and
  \`source_tab\` to match on.

- refvar:

  A string specifying the variable name in \`ref_tab\` containing the
  reference fixations (default is "fixgroup").

- sourcevar:

  A string specifying the variable name in \`source_tab\` containing the
  source fixations (default is "fixgroup").

- window:

  A numeric vector of length 2 specifying the time window to restrict
  the fixations in the source fixation data (default is NULL, which
  considers all fixations).

- ...:

  Additional arguments passed to the \`estimate_scale\` function.

## Value

A data frame containing the original source fixation data with
additional columns for the matched scaling parameters.

## Examples

``` r
# Example usage of the match_scale function
ref_tab <- # reference fixation data
source_tab <- # source fixation data
matched_data <- match_scale(ref_tab, source_tab, match_on = "subject_id")
#> Error in match_scale(ref_tab, source_tab, match_on = "subject_id"): could not find function "match_scale"
```
