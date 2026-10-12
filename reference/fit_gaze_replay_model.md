# Fit a directional GazeWeave replay model on training trials

The fitted warp, background distribution, emission scale, and transition
parameters are learned from the supplied matched training trials. When
\`contrast_on\` is supplied, the likelihood temperature is learned from
inner-held-out item scores with a fixed weak log-temperature prior.
Scored trials must be disjoint from all of these rows. Training strata
too small for two candidate-bearing inner folds retain temperature one
and record the fallback in \`model\$calibration\`.

## Usage

``` r
fit_gaze_replay_model(
  ref_tab,
  source_tab,
  match_on,
  contrast_on = NULL,
  refvar = "fixgroup",
  sourcevar = "fixgroup",
  spec = gaze_replay_spec(),
  template_on = NULL
)
```

## Arguments

- ref_tab, source_tab:

  Matched training reference and source tables.

- match_on:

  Columns defining matched paths.

- contrast_on:

  Optional columns defining training candidate sets for
  likelihood-temperature calibration.

- refvar, sourcevar:

  Fixation-group list columns.

- spec:

  A \[gaze_replay_spec()\].

- template_on:

  Optional columns that uniquely identify separate encoding
  presentations within each \`match_on\` episode. Their likelihoods are
  combined with equal prior weight; paths are never concatenated.

## Value

A fitted \`gaze_replay_model\`.

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
model <- fit_gaze_replay_model(study, recall,
                               match_on = c("participant", "image"),
                               contrast_on = "participant")
#> Replay revision 2026.10: no screen declared, so emission densities are normalised over the whole plane. Declare gaze_screen() for screen-truncated emissions. (Shown once per session.)
model
#> $spec
#> $max_skip
#> [1] 2
#> 
#> $student_df
#> [1] 4
#> 
#> $scale_floor
#> [1] 0.05
#> 
#> $transition_grid
#> $transition_grid$background
#> [1] 0.02 0.08
#> 
#> $transition_grid$restart
#> [1] 0.01 0.05 0.12
#> 
#> $transition_grid$advance
#> [1] 0.15 0.35
#> 
#> $transition_grid$background_stay
#> [1] 0.85 0.95
#> 
#> 
#> $temperature_bounds
#> [1] 5e-02 1e+02
#> 
#> $reliability
#> [1] "none"
#> 
#> $reliability_kappa_bounds
#> [1]   0 100
#> 
#> $warp
#> $type
#> [1] "none"
#> 
#> $fit_by
#> NULL
#> 
#> attr(,"class")
#> [1] "gaze_warp_spec" "list"          
#> 
#> $screen
#> NULL
#> 
#> $chronology
#> $type
#> [1] "local_order"
#> 
#> $horizon
#> [1] 0.2
#> 
#> attr(,"class")
#> [1] "gaze_chronology_spec" "list"                
#> 
#> $calibration
#> $calibration$folds
#> [1] 2
#> 
#> $calibration$seed
#> [1] 104729
#> 
#> $calibration$log_temperature_prior_sd
#> [1] 1
#> 
#> $calibration$log1p_kappa_prior_sd
#> [1] 1
#> 
#> $calibration$control
#> GazeWeave calibration control (revision 2026.10)
#>   method: evidence_scaled 
#>   gamma bounds: 0 to 1 
#>   inverse-temperature prior sd: 3 
#>   typicality offset: none 
#> 
#> 
#> $version
#> [1] 4
#> 
#> $revision
#> [1] "2026.10"
#> 
#> $background_by
#> NULL
#> 
#> $background_support
#> [1] "training"
#> 
#> attr(,"class")
#> [1] "gaze_replay_spec" "list"            
#> 
#> $warp
#> $type
#> [1] "none"
#> 
#> $fit_by
#> NULL
#> 
#> $group_models
#> $group_models$all
#> $group_models$all$A
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> $group_models$all$translation
#> [1] 0 0
#> 
#> $group_models$all$center
#> [1] 0 0
#> 
#> $group_models$all$shift
#> [1] 0 0
#> 
#> $group_models$all$scale
#> [1] 1
#> 
#> 
#> 
#> $info
#> $info$type
#> [1] "none"
#> 
#> $info$groups
#> [1] "all"
#> 
#> 
#> attr(,"class")
#> [1] "gaze_warp_model" "list"           
#> 
#> $emission_models
#> $emission_models$all
#> $emission_models$all$revision
#> [1] "2026.10"
#> 
#> $emission_models$all$replay_scale
#> [1] 0.5252298
#> 
#> $emission_models$all$duration
#> $emission_models$all$duration$replay_intercept
#> [1] -1.631652e-08
#> 
#> $emission_models$all$duration$replay_slope
#> [1] 1
#> 
#> $emission_models$all$duration$replay_sd
#> [1] 0.05
#> 
#> $emission_models$all$duration$background_mean
#> [1] 5.35813
#> 
#> $emission_models$all$duration$background_sd
#> [1] 0.05442562
#> 
#> 
#> $emission_models$all$screen
#> $emission_models$all$screen$xlim
#> [1] -Inf  Inf
#> 
#> $emission_models$all$screen$ylim
#> [1] -Inf  Inf
#> 
#> $emission_models$all$screen$area
#> [1] Inf
#> 
#> $emission_models$all$screen$source
#> [1] "plane"
#> 
#> 
#> $emission_models$all$training_fixations
#> [1] 48
#> 
#> 
#> 
#> $parameters
#> $parameters$background
#> [1] 1e-06
#> 
#> $parameters$restart
#> [1] 9.99999e-07
#> 
#> $parameters$advance
#> [1] 0.999999
#> 
#> $parameters$background_stay
#> [1] 1e-06
#> 
#> 
#> $temperature
#> [1] 7.643808
#> 
#> $reliability_kappa
#> [1] 0
#> 
#> $calibration
#> $inverse_temperature
#> [1] 0.1308249
#> 
#> $temperature
#> [1] 7.643808
#> 
#> $gamma
#> [1] 0
#> 
#> $evidence_reference
#> [1] 4
#> 
#> $score_scale
#> [1] 20.78655
#> 
#> $standardized_inverse_temperature
#> [1] 2.719397
#> 
#> $log_loss
#> [1] 0.02698156
#> 
#> $penalized_objective
#> [1] 0.06121822
#> 
#> $global_temperature
#> [1] 7.053444
#> 
#> $global_log_loss
#> [1] 0.01964426
#> 
#> $uniform_log_loss
#> [1] 1.098612
#> 
#> $training_n
#> [1] 12
#> 
#> $gamma_bounds
#> [1] 0 1
#> 
#> $temperature_bounds
#> [1] 5e-02 1e+02
#> 
#> $inverse_temperature_prior_sd
#> [1] 3
#> 
#> $at_confidence_bound
#> [1] FALSE
#> 
#> $at_prior
#> [1] FALSE
#> 
#> $unshrunk_standardized_inverse_temperature
#> [1] 2.947007
#> 
#> $likelihood_ratio
#> [1] 25.89523
#> 
#> $free_parameters
#> [1] 2
#> 
#> $stein_factor
#> [1] 0.9227657
#> 
#> $kappa
#> [1] 0
#> 
#> $scheme
#> [1] "inner_cross_fitted_evidence_scaled_temperature"
#> 
#> $evidence_measure
#> [1] "recall_fixations"
#> 
#> $score_semantics
#> [1] "total_trial_log_likelihood"
#> 
#> $folds
#> $folds[[1]]
#> $folds[[1]]$fold
#> [1] 1
#> 
#> $folds[[1]]$train_rows
#> [1]  2  6  7  9 11 12
#> 
#> $folds[[1]]$eval_rows
#> [1]  1  3  4  5  8 10
#> 
#> $folds[[1]]$train_match_keys
#> [1] "2:p2|1:1" "2:p2|1:3" "2:p1|1:4" "2:p1|1:5" "2:p1|1:6" "2:p2|1:6"
#> 
#> $folds[[1]]$eval_match_keys
#> [1] "2:p1|1:1" "2:p1|1:2" "2:p2|1:2" "2:p1|1:3" "2:p2|1:4" "2:p2|1:5"
#> 
#> $folds[[1]]$overlap_match_n
#> [1] 0
#> 
#> 
#> $folds[[2]]
#> $folds[[2]]$fold
#> [1] 2
#> 
#> $folds[[2]]$train_rows
#> [1]  1  3  4  5  8 10
#> 
#> $folds[[2]]$eval_rows
#> [1]  2  6  7  9 11 12
#> 
#> $folds[[2]]$train_match_keys
#> [1] "2:p1|1:1" "2:p1|1:2" "2:p2|1:2" "2:p1|1:3" "2:p2|1:4" "2:p2|1:5"
#> 
#> $folds[[2]]$eval_match_keys
#> [1] "2:p2|1:1" "2:p2|1:3" "2:p1|1:4" "2:p1|1:5" "2:p1|1:6" "2:p2|1:6"
#> 
#> $folds[[2]]$overlap_match_n
#> [1] 0
#> 
#> 
#> 
#> $seed
#> [1] 104729
#> 
#> attr(,"class")
#> [1] "gaze_evidence_calibration" "list"                     
#> 
#> $template_on
#> NULL
#> 
#> $template_prior
#> [1] "equal"
#> 
#> $training
#> $training$source_n
#> [1] 12
#> 
#> $training$match_keys
#>  [1] "2:p1|1:1" "2:p2|1:1" "2:p1|1:2" "2:p2|1:2" "2:p1|1:3" "2:p2|1:3"
#>  [7] "2:p1|1:4" "2:p2|1:4" "2:p1|1:5" "2:p2|1:5" "2:p1|1:6" "2:p2|1:6"
#> 
#> $training$template_count
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $training$em
#> $training$em$duration
#> $training$em$duration$replay_intercept
#> [1] -1.631652e-08
#> 
#> $training$em$duration$replay_slope
#> [1] 1
#> 
#> $training$em$duration$replay_sd
#> [1] 0.05
#> 
#> $training$em$duration$background_mean
#> [1] 5.35813
#> 
#> $training$em$duration$background_sd
#> [1] 0.05442562
#> 
#> 
#> $training$em$log_likelihood_trace
#>  [1] -384.0999 -290.5843 -276.6010 -276.0123 -275.8176 -275.7500 -275.7267
#>  [8] -275.7187 -275.7160 -275.7151 -275.7147 -275.7146 -275.7146
#> 
#> $training$em$iterations
#> [1] 13
#> 
#> $training$em$converged
#> [1] TRUE
#> 
#> $training$em$monotone
#> [1] TRUE
#> 
#> $training$em$floor_weight
#> [1] 0.001000004
#> 
#> $training$em$background_fraction
#> [1] 1.335245e-08
#> 
#> $training$em$pair_background_coverage
#>  [1] 2.708974e-09 5.269253e-08 1.294903e-08 5.908431e-09 9.316896e-09
#>  [6] 1.267299e-08 2.491733e-08 5.756674e-09 6.305082e-09 3.640834e-09
#> [11] 1.927637e-08 4.084307e-09
#> 
#> 
#> $training$background_level
#>  [1] "participant" "participant" "participant" "participant" "participant"
#>  [6] "participant" "participant" "participant" "participant" "participant"
#> [11] "participant" "participant"
#> 
#> 
#> $revision
#> [1] "2026.10"
#> 
#> $background
#> $background$table
#>            x           y mass trial background_key exclude_key
#> 1   2.234045  1.62777999 0.25     1            all    2:p1|1:1
#> 2   4.656041  3.02526952 0.25     1            all    2:p1|1:1
#> 3   6.265232  0.45344701 0.25     1            all    2:p1|1:1
#> 4   5.715960  4.09481726 0.25     1            all    2:p1|1:1
#> 5   2.607249  2.23051121 0.25     2            all    2:p2|1:1
#> 6   6.060372  4.00930565 0.25     2            all    2:p2|1:1
#> 7   7.430024  0.61855813 0.25     2            all    2:p2|1:1
#> 8   6.546105  3.72301316 0.25     2            all    2:p2|1:1
#> 9   6.859369  3.80800428 0.25     3            all    2:p1|1:2
#> 10  7.609582  4.83253202 0.25     3            all    2:p1|1:2
#> 11  9.875572  3.29651300 0.25     3            all    2:p1|1:2
#> 12  8.764315  5.89360171 0.25     3            all    2:p1|1:2
#> 13  5.696426  3.93988556 0.25     4            all    2:p2|1:2
#> 14  8.805823  5.42759978 0.25     4            all    2:p2|1:2
#> 15  9.871252  2.95586136 0.25     4            all    2:p2|1:2
#> 16  8.892266  5.97741950 0.25     4            all    2:p2|1:2
#> 17  8.591004  0.31889772 0.25     5            all    2:p1|1:3
#> 18 10.805438  0.08896355 0.25     5            all    2:p1|1:3
#> 19 13.036096 -0.81606528 0.25     5            all    2:p1|1:3
#> 20 11.646663  1.07813011 0.25     5            all    2:p1|1:3
#> 21  8.819414 -1.14861566 0.25     6            all    2:p2|1:3
#> 22 10.683032  1.70594999 0.25     6            all    2:p2|1:3
#> 23 12.608743 -1.99898346 0.25     6            all    2:p2|1:3
#> 24 11.965862  1.72188176 0.25     6            all    2:p2|1:3
#> 25 11.330448  1.22821968 0.25     7            all    2:p1|1:4
#> 26 13.549509  2.01563668 0.25     7            all    2:p1|1:4
#> 27 17.252300  1.27011226 0.25     7            all    2:p1|1:4
#> 28 15.010437  3.98886410 0.25     7            all    2:p1|1:4
#> 29 11.809159  2.60001728 0.25     8            all    2:p2|1:4
#> 30 13.442383  2.62723998 0.25     8            all    2:p2|1:4
#> 31 15.107524  0.16934389 0.25     8            all    2:p2|1:4
#> 32 14.354885  5.12157437 0.25     8            all    2:p2|1:4
#> 33 15.255060  3.62845417 0.25     9            all    2:p1|1:5
#> 34 16.856812  6.32366148 0.25     9            all    2:p1|1:5
#> 35 19.635090  2.84698378 0.25     9            all    2:p1|1:5
#> 36 18.531854  5.14530321 0.25     9            all    2:p1|1:5
#> 37 14.913360  4.27419928 0.25    10            all    2:p2|1:5
#> 38 17.124523  4.95370824 0.25    10            all    2:p2|1:5
#> 39 20.384787  2.79959949 0.25    10            all    2:p2|1:5
#> 40 18.063481  5.97916438 0.25    10            all    2:p2|1:5
#> 41 18.472584 -0.73879405 0.25    11            all    2:p1|1:6
#> 42 21.245147  1.59033734 0.25    11            all    2:p1|1:6
#> 43 22.616435 -0.86804512 0.25    11            all    2:p1|1:6
#> 44 21.724745  1.11964998 0.25    11            all    2:p1|1:6
#> 45 18.312614 -0.25812705 0.25    12            all    2:p2|1:6
#> 46 19.904747  0.44433430 0.25    12            all    2:p2|1:6
#> 47 22.878752 -1.10626238 0.25    12            all    2:p2|1:6
#> 48 20.540351  2.24120707 0.25    12            all    2:p2|1:6
#> 
#> $background$screen
#> $background$screen$xlim
#> [1] -Inf  Inf
#> 
#> $background$screen$ylim
#> [1] -Inf  Inf
#> 
#> $background$screen$area
#> [1] Inf
#> 
#> $background$screen$source
#> [1] "plane"
#> 
#> 
#> $background$control
#> $background$control$max_iterations
#> [1] 200
#> 
#> $background$control$tolerance
#> [1] 1e-06
#> 
#> $background$control$probability_bounds
#> [1] 0.000001 0.999999
#> 
#> $background$control$background_floor
#> [1] 0.05
#> 
#> $background$control$background_floor_bounds
#> [1] 0.001000 0.999999
#> 
#> $background$control$background_min_trials
#> [1] 3
#> 
#> $background$control$quadrature_nodes
#> [1] 40
#> 
#> $background$control$max_scale_fraction
#> [1] 0.1666667
#> 
#> $background$control$duration_sd_floor
#> [1] 0.05
#> 
#> 
#> $background$bandwidth_floor
#> [1] 0.05
#> 
#> $background$floor_center
#> [1] 12.79923  2.38043
#> 
#> $background$floor_sd
#> [1] 17.0373
#> 
#> $background$background_by
#> NULL
#> 
#> $background$exclude_on
#> [1] "participant" "image"      
#> 
#> $background$support
#> [1] "training"
#> 
#> $background$floor_weight
#> [1] 0.001000004
#> 
#> 
#> $match_on
#> [1] "participant" "image"      
#> 
#> $version
#> [1] 5
#> 
#> attr(,"class")
#> [1] "gaze_replay_model" "list"             
# }
```
