# Latent-space transforms for template-based similarity

These helpers convert densities to low-dimensional vectors before
similarity is computed. They are intended to be passed via the
\`similarity_transform\` argument of \`template_similarity\`.

## Usage

``` r
latent_pca_transform(
  ref_tab,
  source_tab,
  match_on,
  refvar = "density",
  sourcevar = "density",
  comps = 30,
  center = TRUE,
  scale. = FALSE,
  ...
)

contract_transform(
  ref_tab,
  source_tab,
  match_on,
  refvar = "density",
  sourcevar = "density",
  shrink = 1e-06,
  fit_by = NULL,
  unique_match_only = FALSE,
  ...
)

affine_transform(
  ref_tab,
  source_tab,
  match_on,
  refvar = "density",
  sourcevar = "density",
  shrink = 1e-06,
  fit_by = NULL,
  unique_match_only = FALSE,
  ...
)

coral_transform(
  ref_tab,
  source_tab,
  match_on,
  refvar = "density",
  sourcevar = "density",
  comps = 30,
  center = TRUE,
  scale. = FALSE,
  shrink = 0.001,
  fit_by = NULL,
  ...
)

cca_transform(
  ref_tab,
  source_tab,
  match_on,
  refvar = "density",
  sourcevar = "density",
  comps = 10,
  center = TRUE,
  scale. = FALSE,
  shrink = 0.001,
  fit_by = NULL,
  unique_match_only = FALSE,
  ...
)
```

## Arguments

- ref_tab, source_tab:

  Data frames passed from \`template_similarity\`.

- match_on:

  Column used for matching rows between tables.

- refvar, sourcevar:

  Density column names to transform.

- comps:

  Number of latent components to retain (post-PCA).

- center, scale.:

  Logical flags passed to \`stats::prcomp\`.

- ...:

  Additional arguments (currently unused; reserved for extension).

- shrink:

  Ridge term for covariance regularization.

- fit_by:

  Optional character vector of column names used to fit CORAL, CCA, or
  geometric transforms separately within strata shared by \`ref_tab\`
  and \`source_tab\`.

- unique_match_only:

  Logical; if \`TRUE\`, supervised fitting uses at most one matched pair
  per \`match_on\` key within each fitting stratum.

## Value

A list with \`ref_tab\`, \`source_tab\`, \`refvar\`, \`sourcevar\`, and
\`info\`.
