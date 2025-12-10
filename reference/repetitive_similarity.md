# Repetitive Similarity Analysis for Density Maps

Computes within-condition and between-condition similarity for density
maps. For each density map (trial), this function calculates its average
similarity to all other maps within the same condition (\`repsim\`) and
its average similarity to all maps from different conditions
(\`othersim\`).

## Usage

``` r
repetitive_similarity(
  tab,
  density_var = "density",
  condition_var,
  method = c("spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov", "emd"),
  pairwise = FALSE,
  multiscale_aggregation = "mean",
  ...
)
```

## Arguments

- tab:

  A data frame or tibble containing the density maps and condition
  identifiers.

- density_var:

  A character string specifying the name of the column containing the
  density maps (must be of class "density" or compatible). Default is
  "density".

- condition_var:

  A character string specifying the name of the column identifying the
  condition for each trial.

- method:

  A character string specifying the similarity method to use, passed to
  \`similarity.density\`. Possible values include "spearman", "pearson",
  "fisherz", "cosine", "l1", "jaccard", "dcov". Default is "spearman".

- pairwise:

  A logical value indicating whether to return the raw pairwise
  similarity scores within the same condition for each trial. Default is
  FALSE. If TRUE, a list column named \`pairwise_repsim\` will be added.

- multiscale_aggregation:

  If densities are multiscale (i.e., \`eye_density_multiscale\`
  objects), this specifies how to aggregate similarities from different
  scales. Options include "mean" (default) or "none" (to get a list
  column of similarity vectors in \`pairwise_repsim\`). Note: \`repsim\`
  and \`othersim\` always report the mean similarity across scales.

- ...:

  Additional arguments passed to the \`similarity.density\` function.

## Value

The input tibble \`tab\` augmented with the following columns: -
\`repsim\`: The average similarity of the trial's density map to other
maps within the same condition. - \`othersim\`: The average similarity
of the trial's density map to maps from all other conditions. -
\`pairwise_repsim\` (optional, if \`pairwise = TRUE\`): A list column
containing vectors of similarity scores between the trial and each other
trial within the same condition.

## Examples

``` r
# \donttest{
  # Generate a small synthetic dataset of density maps across two conditions.
  # Each "density_map" is created from normally-distributed random samples.
  set.seed(123)
  n_trials   <- 20
  conditions <- rep(c("A", "B"), each = n_trials / 2)

  my_data <- tibble::tibble(
    subject         = rep(1:4, length.out = n_trials),
    trial_condition = conditions,
    density_map     = purrr::map(seq_len(n_trials), function(i) {
      x <- rnorm(100,
                 mean = ifelse(conditions[i] == "A", 0, 2),
                 sd   = 1)
      stats::density(x)
    })
  )

  # Compute within- and between-condition similarity.
  result <- repetitive_similarity(
    my_data,
    density_var   = "density_map",
    condition_var = "trial_condition",
    method        = "cosine"
  )
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.

  # Optionally, return the raw pairwise similarities.
  result_pairwise <- repetitive_similarity(
    my_data,
    density_var   = "density_map",
    condition_var = "trial_condition",
    method        = "cosine",
    pairwise      = TRUE
  )
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
#> Warning: Less than 2 common valid data points for similarity calculation.
# }
```
