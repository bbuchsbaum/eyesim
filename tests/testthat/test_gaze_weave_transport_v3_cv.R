make_transport_v3_cv_fixture <- function() {
  make_path <- function(anchor, episode = 1, source = FALSE) {
    offset <- if (source) c(0.1, 0.05) else c(0, 0)
    coords <- rbind(
      anchor + offset + c(0, 0),
      anchor + offset + c(1, 0.3 + 0.03 * episode),
      anchor + offset + c(2, -0.2),
      anchor + offset + c(3, 0.5)
    )
    make_gaze_fixations(
      coords, duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
    )
  }
  references <- expand.grid(
    participant = c("p1", "p2"),
    item = 1:4,
    presentation = 1:4,
    stringsAsFactors = FALSE
  )
  references$prior_weight <- references$item
  references$fixgroup <- lapply(seq_len(nrow(references)), function(index) {
    anchor <- c(
      references$item[[index]] * 3,
      if (references$participant[[index]] == "p1") 0 else 10
    )
    make_path(anchor, references$presentation[[index]])
  })
  sources <- expand.grid(
    participant = c("p1", "p2"), item = 1:4,
    stringsAsFactors = FALSE
  )
  sources$fixgroup <- lapply(seq_len(nrow(sources)), function(index) {
    anchor <- c(
      sources$item[[index]] * 3,
      if (sources$participant[[index]] == "p1") 0 else 10
    )
    make_path(anchor, source = TRUE)
  })
  list(references = references, sources = sources)
}

make_transport_v3_cv_spec <- function(reliability = "effective_fixations") {
  gaze_transport_spec(
    coverage_nodes = 2,
    entropy_schedule = 0.03,
    maxit = 40,
    tolerance = 1e-3,
    projection_maxit = 300,
    projection_tolerance = 1e-7,
    backend = "optimized",
    reliability = reliability,
    calibration_folds = 2,
    calibration_seed = 20260822
  )
}

test_that("episode-scale calibration contains the temperature-only boundary", {
  profile_sets <- replicate(
    8,
    list(
      item_a = c(0.8, 0.7, 0.9, 0.75),
      item_b = c(0.1, 0.2, 0.05, 0.15),
      item_c = c(-0.2, -0.1, -0.3, -0.15)
    ),
    simplify = FALSE
  )
  truth <- rep(c(1L, 2L), 4)
  prior <- rep(list(c(0.5, 0.3, 0.2)), 8)
  effective <- rep(c(1, 4), 4)
  calibration <- eyesim:::fit_transport_v3_calibration(
    profile_sets, truth, effective, prior,
    make_transport_v3_cv_spec()
  )

  expect_lte(calibration$reliability_worsening, 1e-10)
  expect_gte(calibration$kappa, 0)
  expect_true(calibration$temperature >= 0.05)
  expect_true(calibration$temperature <= 20)
  expect_equal(
    eyesim:::transport_v3_candidate_score(c(1, 3), 2),
    eyesim:::transport_v3_log_mean_exp(c(0.5, 1.5)),
    tolerance = 1e-14
  )
  expect_lt(
    eyesim:::transport_v3_reliability(1, 2),
    eyesim:::transport_v3_reliability(4, 2)
  )
})

test_that("null scores recover the actual candidate prior", {
  pool_size <- rep(c(2L, 3L, 5L), 3)
  profile_sets <- lapply(pool_size, function(size) {
    stats::setNames(
      rep(list(rep(0, 4)), size), paste0("item_", seq_len(size))
    )
  })
  truth <- vapply(seq_along(pool_size), function(index) {
    (index - 1L) %% pool_size[[index]] + 1L
  }, integer(1))
  prior <- lapply(pool_size, function(size) {
    weights <- seq_len(size)
    weights / sum(weights)
  })
  calibration <- eyesim:::fit_transport_v3_calibration(
    profile_sets, truth, rep(c(1, 2, 4), 3), prior,
    make_transport_v3_cv_spec()
  )

  expect_equal(
    vapply(calibration$evidence, `[[`, numeric(1), "gaze_info_bits"),
    rep(0, length(pool_size)), tolerance = 1e-12
  )
  for (index in seq_along(calibration$evidence)) {
    expect_equal(calibration$evidence[[index]]$candidates$posterior,
                 prior[[index]],
                 tolerance = 1e-12)
  }
})

test_that("nested Transport keeps held-out cells outside every fit", {
  fixture <- make_transport_v3_cv_fixture()
  fit <- gaze_transport_cv(
    fixture$references,
    fixture$sources,
    match_on = c("participant", "item"),
    contrast_on = "participant",
    spec = make_transport_v3_cv_spec(),
    n_folds = 2,
    seed = 20260822,
    episode_on = "presentation",
    priorvar = "prior_weight"
  )

  expect_s3_class(fit, "gaze_transport_fit")
  expect_equal(nrow(fit$results), 8)
  expect_equal(fit$results$candidate_count, rep(4L, 8))
  expect_equal(fit$results$common_episode_count, rep(4L, 8))
  expect_equal(
    fit$results$gaze_info_bits,
    log2(fit$results$posterior_true / fit$results$prior_true),
    tolerance = 1e-12
  )
  expected_prior <- fit$results$item / sum(1:4)
  expect_equal(fit$results$prior_true, expected_prior, tolerance = 1e-14)
  expect_true(all(fit$results$all_converged))

  for (receipt in fit$folds) {
    expect_equal(receipt$overlap_match_n, 0L)
    expect_length(intersect(receipt$train_rows, receipt$eval_rows), 0L)
    expect_false(receipt$heldout_cell_contributed_to_fit)
    expect_lte(receipt$calibration$reliability_worsening, 1e-10)
    for (inner in receipt$inner_receipts) {
      expect_equal(inner$overlap_match_n, 0L)
      expect_length(intersect(inner$train_rows, inner$eval_rows), 0L)
      expect_identical(
        inner$candidate_policy,
        "exhaustive_permitted_with_actual_prior"
      )
    }
  }
})

test_that("all episodes and candidates share warp and solver policy", {
  fixture <- make_transport_v3_cv_fixture()
  fit <- gaze_transport_cv(
    fixture$references,
    fixture$sources,
    match_on = c("participant", "item"),
    contrast_on = "participant",
    spec = make_transport_v3_cv_spec(reliability = "none"),
    n_folds = 2,
    episode_on = "presentation",
    priorvar = "prior_weight"
  )
  candidates <- fit$results$alignments[[1L]]
  episode_results <- unlist(lapply(candidates, function(candidate) {
    candidate$alignment$episodes
  }), recursive = FALSE)
  warp_scale <- vapply(episode_results, function(result) {
    result$diagnostics$warp_scale
  }, numeric(1))
  solver_backend <- vapply(episode_results, function(result) {
    result$convergence$backend
  }, character(1))

  expect_equal(
    unname(warp_scale), rep(warp_scale[[1L]], length(warp_scale))
  )
  expect_equal(
    unname(solver_backend),
    rep(solver_backend[[1L]], length(solver_backend))
  )
  for (candidate in candidates) {
    expect_equal(
      unname(candidate$diagnostics$equal_episode_weights),
      rep(0.25, 4), tolerance = 0
    )
  }
})

test_that("quality and calibration remain explicit separate diagnostics", {
  fixture <- make_transport_v3_cv_fixture()
  fit <- gaze_transport_cv(
    fixture$references,
    fixture$sources,
    match_on = c("participant", "item"),
    contrast_on = "participant",
    spec = make_transport_v3_cv_spec(),
    n_folds = 2,
    episode_on = "presentation",
    priorvar = "prior_weight"
  )

  expect_true(all(c(
    "raw_fixation_count", "effective_fixations", "total_duration",
    "duration_concentration", "duration_entropy_bits",
    "normalized_duration_entropy", "spatial_dispersion",
    "generic_gaze_quality"
  ) %in% names(fit$results)))
  expect_identical(
    fit$provenance$generic_gaze_quality,
    "separate_response_blind_diagnostic_channel"
  )
  expect_true(is.finite(fit$calibration$heldout_log_loss))
  expect_true(is.data.frame(fit$calibration$heldout_reliability_curve))
  expect_true(is.finite(fit$calibration$heldout_ece))
  expect_named(
    fit$calibration$candidate_pool_size_sensitivity,
    c("candidate_count", "gaze_info_bits", "log_loss")
  )
})
