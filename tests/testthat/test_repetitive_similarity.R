library(testthat)
library(dplyr)
library(tibble)
library(purrr)

context("Repetitive Similarity")

test_that("repetitive_similarity works with standard (non-multiscale) density objects", {
  set.seed(101)
  n_conditions <- 2
  n_trials_per_cond <- 4
  n_total_trials <- n_conditions * n_trials_per_cond
  single_sigma <- 50

  # Create fixation data with conditions
  fix_data <- tibble::tibble(
    trial_id = rep(1:n_total_trials, each = 20),
    condition = rep(c("Cond1", "Cond2"), each = n_trials_per_cond * 20),
    x = runif(n_total_trials * 20, 0, 1024),
    y = runif(n_total_trials * 20, 0, 768),
    duration = rnorm(n_total_trials * 20, 220, 25),
    onset = rep(seq(0, 950, by = 50), n_total_trials) + rep((0:(n_total_trials-1))*1000, each=20)
  )
  # Introduce slight condition difference
  fix_data <- fix_data %>%
    dplyr::mutate(x = ifelse(condition == "Cond2", x - 40, x),
                  y = ifelse(condition == "Cond2", y - 20, y))

  # Create eye table grouped by trial and condition
  et_per_trial <- eye_table(x = "x", y = "y", duration = "duration", onset = "onset",
                          groupvar = c("trial_id", "condition"), data = fix_data,
                          clip_bounds = c(0, 1024, 0, 768))

  # Generate standard (single-scale) densities per trial
  dens_per_trial <- density_by(et_per_trial, groups = c("trial_id", "condition"),
                             sigma = single_sigma, xbounds = c(0, 1024), ybounds = c(0, 768),
                             result_name = "standard_density")

  # Now dens_per_trial has n_total_trials rows, with n_trials_per_cond per condition.
  # Test repetitive_similarity - it returns input table augmented with new columns
  rep_sim_standard <- repetitive_similarity(
    dens_per_trial, # Use dens_per_trial
    density_var = "standard_density",
    condition_var = "condition",
    pairwise = TRUE, # We need pairwise results for detailed checks
    method = "pearson" # Using pearson for this test
  )

  # Check output structure and types
  expect_true(is.data.frame(rep_sim_standard))
  expect_equal(nrow(rep_sim_standard), n_total_trials) # Should return input rows + new columns
  # Check expected columns: original from dens_per_trial + repsim, othersim, pairwise_repsim
  expect_true(all(c(names(dens_per_trial), "repsim", "othersim", "pairwise_repsim") %in% names(rep_sim_standard)))

  # Check calculation result types
  expect_true(is.numeric(rep_sim_standard$repsim))
  expect_true(is.numeric(rep_sim_standard$othersim))
  expect_true(is.list(rep_sim_standard$pairwise_repsim))

  # Check pairwise_repsim structure
  expect_equal(length(rep_sim_standard$pairwise_repsim), n_total_trials) # One list element per input row

  # For each row, the pairwise_repsim list should contain similarities to other trials in the same condition
  expected_pairwise_count_per_row <- n_trials_per_cond - 1

  # Check the first row's pairwise results (assuming first row is Cond1)
  # Need to be careful if data gets reordered, but test data is sequential
  # Access the list of pairwise comparisons for the first row
  first_row_pairwise_list <- rep_sim_standard$pairwise_repsim[[1]]
  expect_true(is.list(first_row_pairwise_list))
  expect_equal(length(first_row_pairwise_list), expected_pairwise_count_per_row)
  # Check that the first pairwise comparison result is numeric
  expect_true(is.numeric(first_row_pairwise_list[[1]]))

  # Check a row from the second condition (e.g., row 5, assuming 4 trials per cond)
  second_cond_row_idx <- n_trials_per_cond + 1
  # Access the list of pairwise comparisons for this row
  second_cond_pairwise_list <- rep_sim_standard$pairwise_repsim[[second_cond_row_idx]]
  expect_true(is.list(second_cond_pairwise_list))
  expect_equal(length(second_cond_pairwise_list), expected_pairwise_count_per_row)
  # Check that the first pairwise comparison result is numeric
  expect_true(is.numeric(second_cond_pairwise_list[[1]]))
}) 