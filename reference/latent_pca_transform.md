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

## Value

A list with \`ref_tab\`, \`source_tab\`, \`refvar\`, \`sourcevar\`, and
\`info\`.
