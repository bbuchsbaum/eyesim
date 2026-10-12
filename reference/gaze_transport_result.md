# Extract one auditable Transport result

Transport has one primary outcome, \`gaze_info_bits\`. Coverage, spatial
error, local order, contraction, rank, and solver behavior are nested
under \`diagnostics\`; they are explanatory properties of the optimized
alignment, not additional evidence scores. When several study
presentations are available, their calibrated likelihoods receive fixed
equal prior weight.

## Usage

``` r
gaze_transport_result(x, row = NULL, key = NULL)
```

## Arguments

- x:

  A fitted object returned by \[gaze_transport_cv()\].

- row:

  Optional result row number.

- key:

  Optional named list selecting one result by identifiers.

## Value

A \`gaze_transport_result\` containing one \`gaze_info_bits\` value,
nested diagnostics, candidate evidence, and provenance.

## Examples

``` r
# Example data: two participants study and later recall six images; each
# recall roughly repeats the study path of the same image.
set.seed(1)
make_path <- function(anchor, jitter) {
  fixation_group(
    x = anchor[1] + c(0, 2, 4, 3) + rnorm(4, sd = jitter),
    y = anchor[2] + c(0, 1, -1, 2) + rnorm(4, sd = jitter),
    onset = c(0, 250, 500, 750), duration = c(200, 220, 180, 240)
  )
}
study <- expand.grid(participant = c("p1", "p2"), image = 1:6,
                     stringsAsFactors = FALSE)
anchors <- lapply(study$image, function(i) c(3 * i, 2 * (i %% 3)))
study$fixgroup <- lapply(anchors, make_path, jitter = 0.2)
recall <- study[c("participant", "image")]
recall$fixgroup <- lapply(anchors, make_path, jitter = 0.6)
# \donttest{
fit <- gaze_transport_cv(study, recall,
                         match_on = c("participant", "image"),
                         contrast_on = "participant",
                         # 2 coverage nodes keep the example fast; the
                         # default is 12.
                         spec = gaze_transport_spec(coverage_nodes = 2),
                         n_folds = 2)
gaze_transport_result(fit, row = 1)$diagnostics
#> $matched_coverage
#> [1] 0.4877017
#> 
#> $spatial_rmse
#> [1] 0.6258791
#> 
#> $spatial_rmse_unit
#> [1] "native_coordinate_units"
#> 
#> $local_order_preservation
#> [1] 0.99835
#> 
#> $contraction_scale
#> [1] 1
#> 
#> $template_rank
#> [1] 1
#> 
#> $candidate_count
#> [1] 6
#> 
#> $episode_count
#> [1] 1
#> 
#> $episode_normalization
#> [1] "fixed_equal_prior_likelihood_mixture"
#> 
#> $episodes
#>             episode equal_episode_weight   log_score matched_coverage
#> episode_1 episode_1                    1 -0.07275334        0.4877017
#>           spatial_rmse       spatial_rmse_unit local_order_preservation
#> episode_1    0.6258791 native_coordinate_units                  0.99835
#>           contraction_scale solver_converged       solver_backend stationarity
#> episode_1                 1             TRUE native_rcpparmadillo   3.1132e-05
#>           relative_objective_change
#> episode_1              5.860762e-08
#> 
#> $solver_stability
#> $solver_stability$all_converged
#> [1] TRUE
#> 
#> $solver_stability$backends
#> [1] "native_rcpparmadillo"
#> 
#> $solver_stability$maximum_stationarity
#> [1] 3.1132e-05
#> 
#> $solver_stability$maximum_relative_objective_change
#> [1] 5.860762e-08
#> 
#> 
# }
```
