# Scanpath Similarity

Compute similarity between scanpaths in a source table and matching
scanpaths in a reference table.

## Usage

``` r
scanpath_similarity(
  ref_tab,
  source_tab,
  match_on,
  permutations = 0,
  permute_on = NULL,
  method = c("multimatch"),
  refvar = "scanpath",
  sourcevar = "scanpath",
  window = NULL,
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

- permutations:

  A numeric value specifying the number of permutations for the baseline
  map (default is 10).

- permute_on:

  A character string representing the variable used to stratify
  permutations (default is NULL).

- method:

  The similarity method to use. Currently only "multimatch" is
  supported.

- refvar:

  The name of the column containing scanpaths in the reference table
  (default is "scanpath").

- sourcevar:

  The name of the column containing scanpaths in the source table
  (default is "scanpath").

- window:

  An optional temporal window for restricting the scanpath comparison.

- ...:

  Extra arguments to pass to the \`similarity\` function.

## Value

A table containing the computed similarities between scanpaths.
