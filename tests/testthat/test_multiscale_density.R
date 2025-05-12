library(testthat)
library(dplyr)
library(tibble)
library(purrr)

context("Multiscale Density Analysis Workflow")

test_that("density_by creates eye_density_multiscale objects with vector sigma", {
  set.seed(123)
  n_trials <- 6
  sigmas_vec <- c(30, 60, 90)

  # Create simple fixation data
  fix_data <- tibble(
    trial = rep(1:n_trials, each = 10),
    condition = rep(c("A", "B"), each = n_trials / 2 * 10),
    x = runif(n_trials * 10, 0, 1000),
    y = runif(n_trials * 10, 0, 1000),
    duration = rnorm(n_trials * 10, 300, 50),
    onset = rep(seq(0, 450, by = 50), n_trials) + rep((0:(n_trials-1))*500, each=10)
  )

  et <- eye_table(x = "x", y = "y", duration = "duration", onset = "onset",
                  groupvar = c("trial", "condition"), data = fix_data, clip_bounds = c(0, 1000, 0, 1000))

  # Generate multiscale densities
  multiscale_densities <- density_by(et, groups = c("trial", "condition"),
                                     sigma = sigmas_vec,
                                     xbounds = c(0, 1000), ybounds = c(0, 1000),
                                     result_name = "multiscale_map")

  # --- Assertions ---
  expect_true("multiscale_map" %in% names(multiscale_densities))
  expect_equal(nrow(multiscale_densities), n_trials)

  # Check first non-null result
  first_map <- multiscale_densities$multiscale_map[[1]]
  expect_s3_class(first_map, "eye_density_multiscale")
  expect_true(is.list(first_map))
  expect_equal(length(first_map), length(sigmas_vec))

  # Check sigmas stored correctly
  stored_sigmas <- sapply(first_map, `[[`, "sigma")
  expect_equal(stored_sigmas, sigmas_vec)

  # Check individual elements are eye_density objects
  expect_s3_class(first_map[[1]], "eye_density")
  expect_s3_class(first_map[[2]], "eye_density")
  expect_s3_class(first_map[[3]], "eye_density")

})

test_that("template_similarity works with multiscale objects and aggregation methods", {
  set.seed(456)
  n_trials <- 4
  sigmas_vec <- c(40, 80)

  # Create two sets of fixation data
  fix_data1 <- tibble(
    trial = rep(1:n_trials, each = 15),
    x = runif(n_trials * 15, 0, 800),
    y = runif(n_trials * 15, 0, 600),
    duration = rnorm(n_trials * 15, 250, 40),
    onset = rep(seq(0, 700, by = 50), n_trials) + rep((0:(n_trials-1))*750, each=15)
  )
  # Make second set slightly different
  fix_data2 <- fix_data1 %>% mutate(x = x + rnorm(n(), 0, 50), y = y + rnorm(n(), 0, 50))

  et1 <- eye_table(x = "x", y = "y", duration = "duration", onset = "onset",
                  groupvar = "trial", data = fix_data1, clip_bounds = c(0, 800, 0, 600))
  et2 <- eye_table(x = "x", y = "y", duration = "duration", onset = "onset",
                   groupvar = "trial", data = fix_data2, clip_bounds = c(0, 800, 0, 600))

  # Generate multiscale densities
  dens1 <- density_by(et1, groups = "trial", sigma = sigmas_vec, xbounds = c(0, 800), ybounds = c(0, 600), result_name = "ms_dens")
  dens2 <- density_by(et2, groups = "trial", sigma = sigmas_vec, xbounds = c(0, 800), ybounds = c(0, 600), result_name = "ms_dens")

  # --- Test with aggregation = 'mean' ---
  tsim_mean <- template_similarity(dens1, dens2, match_on="trial", refvar = "ms_dens", sourcevar = "ms_dens",
                                   method="cosine", permutations = 0, multiscale_aggregation = "mean")

  expect_true("eye_sim" %in% names(tsim_mean))
  expect_equal(nrow(tsim_mean), n_trials)
  expect_true(is.numeric(tsim_mean$eye_sim)) # Should be single value per row
  expect_false(is.list(tsim_mean$eye_sim))  # Ensure it's not a list column

  # --- Test with aggregation = 'none' ---
  tsim_none <- template_similarity(dens1, dens2, match_on="trial", refvar = "ms_dens", sourcevar = "ms_dens",
                                   method="cosine", permutations = 0, multiscale_aggregation = "none")

  expect_true("eye_sim" %in% names(tsim_none))
  expect_equal(nrow(tsim_none), n_trials)
  expect_true(is.list(tsim_none$eye_sim))     # Should be a list column
  expect_equal(length(tsim_none$eye_sim[[1]]), length(sigmas_vec)) # Each list element has length = n_sigmas
  expect_true(is.numeric(tsim_none$eye_sim[[1]])) # Elements within the list are numeric vectors

})

test_that("repetitive_similarity works with multiscale objects and aggregation methods", {
  set.seed(789)
  n_conditions <- 2
  n_trials_per_cond <- 3
  n_total_trials <- n_conditions * n_trials_per_cond
  sigmas_vec <- c(25, 50)

  # Create fixation data with conditions
  fix_data <- tibble(
    trial_id = rep(1:n_total_trials, each = 12),
    condition = rep(LETTERS[1:n_conditions], each = n_trials_per_cond * 12),
    x = runif(n_total_trials * 12, 0, 1200),
    y = runif(n_total_trials * 12, 0, 900),
    duration = rnorm(n_total_trials * 12, 200, 30),
    onset = rep(seq(0, 550, by = 50), n_total_trials) + rep((0:(n_total_trials-1))*600, each=12)
  )
  # Introduce slight condition difference (shift condition B slightly)
  fix_data <- fix_data %>%
    mutate(x = ifelse(condition == "B", x + 50, x),
           y = ifelse(condition == "B", y + 30, y))

  et <- eye_table(x = "x", y = "y", duration = "duration", onset = "onset",
                  groupvar = c("trial_id", "condition"), data = fix_data, clip_bounds = c(0, 1200, 0, 900))

  # Generate multiscale densities
  dens <- density_by(et, groups = c("trial_id", "condition"),
                     sigma = sigmas_vec, xbounds = c(0, 1200), ybounds = c(0, 900),
                     result_name = "ms_map")

  # --- Test with aggregation = 'mean' and pairwise=TRUE ---
  rep_sim_mean <- repetitive_similarity(dens, density_var = "ms_map", condition_var = "condition",
                                        method = "spearman", multiscale_aggregation = "mean", pairwise = TRUE)

  expect_true(all(c("repsim", "othersim", "pairwise_repsim") %in% names(rep_sim_mean)))
  expect_equal(nrow(rep_sim_mean), n_total_trials)
  expect_true(is.numeric(rep_sim_mean$repsim))   # repsim/othersim are always mean aggregated
  expect_true(is.numeric(rep_sim_mean$othersim))
  expect_false(is.list(rep_sim_mean$repsim))
  expect_false(is.list(rep_sim_mean$othersim))

  # pairwise_repsim should be a list, but contain single numeric means
  expect_true(is.list(rep_sim_mean$pairwise_repsim))
  # Each element of the list is a list of length n_trials_per_cond - 1
  expect_equal(length(rep_sim_mean$pairwise_repsim[[1]]), n_trials_per_cond - 1)
  # Each of *those* list elements should be a single numeric value (the mean sim)
  expect_true(is.numeric(rep_sim_mean$pairwise_repsim[[1]][[1]]))
  expect_equal(length(rep_sim_mean$pairwise_repsim[[1]][[1]]), 1)


  # --- Test with aggregation = 'none' and pairwise=TRUE ---
  rep_sim_none <- repetitive_similarity(dens, density_var = "ms_map", condition_var = "condition",
                                        method = "spearman", multiscale_aggregation = "none", pairwise = TRUE)

  expect_true(all(c("repsim", "othersim", "pairwise_repsim") %in% names(rep_sim_none)))
  expect_equal(nrow(rep_sim_none), n_total_trials)
  expect_true(is.numeric(rep_sim_none$repsim))   # repsim/othersim are still mean aggregated
  expect_true(is.numeric(rep_sim_none$othersim))

  # pairwise_repsim should be a list containing vectors/lists of per-scale similarities
  expect_true(is.list(rep_sim_none$pairwise_repsim))
  expect_equal(length(rep_sim_none$pairwise_repsim[[1]]), n_trials_per_cond - 1)
  # Each of *those* elements should be a numeric vector of per-scale similarities
  expect_true(is.numeric(rep_sim_none$pairwise_repsim[[1]][[1]]))
  expect_equal(length(rep_sim_none$pairwise_repsim[[1]][[1]]), length(sigmas_vec))

}) 