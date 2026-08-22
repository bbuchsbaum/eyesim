make_baseline_test_spec <- function(
    methods = c("multimatch", "density", "elastic"),
    warp = gaze_warp_none()) {
  gaze_baseline_spec(
    screen = gaze_screen(100, 100, unit = "deg"),
    density_sigmas = c(3, 8),
    density_grid = 8,
    methods = methods,
    warp = warp,
    lambda_grid = c(0.1, 1),
    inner_folds = 2,
    elastic_radii = c(consensus = 5, rigidity = 20, matching = 3),
    elastic_maxit = 15,
    elastic_tolerance = 1e-3
  )
}

make_baseline_cv_tables <- function() {
  paths <- lapply(seq_len(6), function(item) {
    anchor <- c(item * 12, (item %% 2) * 35 + 15)
    coords <- rbind(
      anchor,
      anchor + c(3, 7),
      anchor + c(9, 1),
      anchor + c(6, 10)
    )
    make_gaze_fixations(
      coords, duration = c(1, 2, 1, 1), onset = c(0, 1, 3, 4)
    )
  })
  list(
    reference = tibble::tibble(image_id = seq_len(6), fixgroup = paths),
    source = tibble::tibble(image_id = seq_len(6), fixgroup = paths)
  )
}

test_that("comparator availability reports missing optional dependencies", {
  spec <- make_baseline_test_spec(methods = c("multimatch", "density"))
  testthat::local_mocked_bindings(
    gaze_namespace_available = function(package) FALSE,
    .package = "eyesim"
  )
  availability <- gaze_baseline_availability(spec)

  expect_false(availability$available[availability$method == "multimatch"])
  expect_match(
    availability$reason[availability$method == "multimatch"],
    "missing optional dependency"
  )
  expect_true(availability$available[availability$method == "density"])
})

test_that("multiscale density retains native per-scale diagnostics", {
  reference <- make_gaze_fixations(rbind(
    c(20, 20), c(30, 40), c(50, 25), c(65, 55)
  ))
  near <- make_gaze_fixations(rbind(
    c(21, 20), c(31, 39), c(49, 26), c(64, 54)
  ))
  far <- make_gaze_fixations(rbind(
    c(75, 75), c(82, 70), c(70, 85), c(90, 90)
  ))
  spec <- make_baseline_test_spec(methods = "density")
  near_features <- eyesim:::baseline_density_features(reference, near, spec)
  far_features <- eyesim:::baseline_density_features(reference, far, spec)

  expect_named(near_features, c("density_sigma_3", "density_sigma_8"))
  expect_true(all(near_features > far_features))
  expect_true(all(near_features >= 0 & near_features <= 1))
})

test_that("elastic consensus is order-free and penalizes unrelated geometry", {
  reference <- make_gaze_fixations(rbind(
    c(20, 20), c(30, 40), c(50, 25), c(65, 55)
  ))
  translated <- make_gaze_fixations(rbind(
    c(17, 18), c(27, 38), c(47, 23), c(62, 53)
  ))
  reversed <- make_gaze_fixations(
    cbind(rev(translated$x), rev(translated$y))
  )
  unrelated <- make_gaze_fixations(rbind(
    c(75, 75), c(82, 70), c(70, 85), c(90, 90)
  ))
  spec <- make_baseline_test_spec(methods = "elastic")
  matched <- elastic_consensus_align(reference, translated, spec)
  reordered <- elastic_consensus_align(reference, reversed, spec)
  wrong <- elastic_consensus_align(reference, unrelated, spec)

  expect_true(matched$convergence$converged)
  expect_gt(matched$soft_similarity, wrong$soft_similarity)
  expect_equal(matched$soft_similarity, reordered$soft_similarity, tolerance = 1e-10)
  expect_equal(matched$score, reordered$score, tolerance = 1e-10)
  expect_gte(matched$encoding_coverage, 0)
  expect_lte(matched$encoding_coverage, 1)
  expect_gte(matched$local_deformation, 0)
})

test_that("ridge ranking learns feature combinations without an intercept", {
  make_set <- function(first_true = TRUE) {
    features <- if (first_true) {
      matrix(c(2, 0.1, -1, 0.2), 2, byrow = TRUE)
    } else {
      matrix(c(-1, 0.2, 2, 0.1), 2, byrow = TRUE)
    }
    colnames(features) <- c("signal", "noise")
    list(
      features = features,
      true_index = if (first_true) 1L else 2L
    )
  }
  sets <- list(make_set(TRUE), make_set(FALSE), make_set(TRUE), make_set(FALSE))
  model <- eyesim:::fit_baseline_ridge(sets, c("signal", "noise"), lambda = 0.1)
  predictions <- lapply(sets, function(set) {
    eyesim:::predict_baseline_ridge(model, set)
  })

  expect_equal(model$convergence, 0)
  expect_gt(model$coefficients[["signal"]], 0)
  expect_true(all(vapply(seq_along(sets), function(i) {
    which.max(predictions[[i]]) == sets[[i]]$true_index
  }, logical(1))))
})

test_that("non-finite native metrics become explicit skipped rows", {
  set <- list(
    features = matrix(c(NA_real_, NA_real_), ncol = 1,
                      dimnames = list(NULL, "registered_mm_vector")),
    candidate_key = c("a", "b"),
    true_index = 1L,
    contrast_key = "all"
  )
  result <- eyesim:::baseline_scored_evidence(
    set,
    method = "multimatch_mm_vector_registered",
    court = "frozen_registered",
    score = c(NA_real_, NA_real_),
    calibrated = FALSE
  )

  expect_identical(result$status, "skipped")
  expect_match(result$reason, "non-finite")
  expect_true(is.na(result$compatibility_bits))
  expect_true(is.na(result$gaze_info_bits))
  expect_true(all(result$candidates$status == "skipped"))
})

test_that("fair court uses disjoint nested folds and identical candidate sets", {
  skip_if_not_installed("igraph")
  if (!any(vapply(
    c("emdist", "T4transport", "transport"),
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  ))) {
    skip("No MultiMatch EMD backend is installed")
  }
  tables <- make_baseline_cv_tables()
  fit <- gaze_baseline_cv(
    tables$reference,
    tables$source,
    match_on = "image_id",
    split_on = "image_id",
    n_folds = 3,
    seed = 29,
    spec = make_baseline_test_spec()
  )

  expect_s3_class(fit, "gaze_baseline_fit")
  expect_true(all(vapply(fit$folds, `[[`, integer(1), "overlap_match_n") == 0L))
  expect_true(all(unlist(lapply(fit$folds, function(fold) {
    vapply(fold$inner_folds, `[[`, integer(1), "overlap_match_n")
  })) == 0L))

  expected_methods <- c(
    "density_raw", "density_registered",
    "density_sigma_3_raw", "density_sigma_8_raw",
    "density_ridge_registered",
    "density_sigma_3_raw_calibrated", "density_sigma_8_raw_calibrated",
    "elastic_consensus_registered", "elastic_ridge_registered",
    "multimatch_ridge_registered",
    "multimatch_mm_vector_raw_calibrated",
    "multimatch_mm_direction_raw_calibrated",
    "multimatch_mm_length_raw_calibrated",
    "multimatch_mm_position_raw_calibrated",
    "multimatch_mm_duration_raw_calibrated",
    "multimatch_mm_position_emd_raw_calibrated"
  )
  expect_true(all(expected_methods %in% unique(fit$results$method)))
  counts <- table(fit$results$method)
  expect_true(all(counts == 6L))
  expect_true(all(fit$results$candidate_count == 2L))
  expect_true(all(fit$results$status == "scored"))
  expect_true(all(is.finite(fit$results$compatibility_bits)))
  expect_true(all(fit$results$calibrated[
    fit$results$court == "supervised_registered"
  ]))
  expect_false(any(fit$results$calibrated[
    grepl("^frozen", fit$results$court)
  ]))
  expect_true(all(fit$results$calibrated[
    fit$results$court == "supervised_raw"
  ]))
  expect_true(all(is.finite(fit$results$gaze_info_bits[fit$results$calibrated])))
  expect_true(all(is.na(fit$results$gaze_info_bits[!fit$results$calibrated])))

  one_trial <- fit$results[fit$results$image_id == 1L, ]
  candidate_keys <- lapply(one_trial$candidates, function(candidates) {
    candidates$candidate_key
  })
  expect_true(all(vapply(candidate_keys[-1], identical, logical(1), candidate_keys[[1]])))
  diagnostic_row <- one_trial$candidates[[which(
    one_trial$method == "multimatch_ridge_registered"
  )]]
  expect_true(all(c(
    "registered_mm_position", "registered_mm_duration",
    "registered_density_sigma_3", "native_diagnostics"
  ) %in% names(diagnostic_row)))

  composite_folds <- lapply(fit$folds, `[[`, "composites")
  expect_true(all(vapply(composite_folds, function(composites) {
    composites$multimatch_ridge_registered$lambda %in% c(0.1, 1) &&
      composites$density_ridge_registered$lambda %in% c(0.1, 1) &&
      composites$elastic_ridge_registered$lambda %in% c(0.1, 1)
  }, logical(1))))
  expect_true(all(vapply(composite_folds, function(composites) {
    composites$multimatch_ridge_registered$convergence == 0L &&
      composites$density_ridge_registered$convergence == 0L &&
      composites$elastic_ridge_registered$convergence == 0L
  }, logical(1))))
})
