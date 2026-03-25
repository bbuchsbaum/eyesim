# Cross-Fitted Template Similarity

Compute template similarity on held-out folds while fitting any optional
domain-adaptation transform only on training rows.

## Usage

``` r
template_similarity_cv(
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
  split_on = match_on,
  n_folds = NULL,
  seed = 1,
  fit_source_filter = NULL,
  eval_source_filter = NULL,
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

- split_on:

  Character vector of source-table columns used to assign folds. All
  rows sharing the same \`split_on\` values are held out together.
  Defaults to \`match_on\`.

- n_folds:

  Number of folds. Defaults to \`min(5, n_unique_groups)\`.

- seed:

  Random seed used for fold assignment.

- fit_source_filter:

  Optional logical vector or function selecting which source rows are
  eligible for transform fitting. Functions receive \`source_tab\` and
  must return a logical vector with one value per row.

- eval_source_filter:

  Optional logical vector or function selecting which source rows are
  scored. Functions receive \`source_tab\` and must return a logical
  vector with one value per row.

- ...:

  Extra arguments to pass to the \`similarity\` function.

## Value

A tibble containing only held-out evaluation rows from \`source_tab\`,
augmented with similarity columns and a \`.cv_fold\` column. Fold
metadata is stored in \`attr(x, "similarity_cv")\`.

## Details

This function is the leakage-safe counterpart to
\`template_similarity()\` when using learned transforms such as
\`latent_pca_transform\`, \`coral_transform\`, or \`cca_transform\`. For
each fold, it:

1.  assigns held-out source rows using \`split_on\`,

2.  excludes held-out \`match_on\` keys from transform fitting,

3.  fits the transform on the remaining training rows only,

4.  applies the fitted transform to held-out source rows and their
    matched reference rows, and

5.  computes similarity only on the held-out rows.
