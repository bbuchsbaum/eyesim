make_replay_test_spec <- function(warp = gaze_warp_none(), grid_size = 16L) {
  gaze_replay_spec(
    grid_size = grid_size,
    max_skip = 2L,
    student_df = 4,
    scale_floor = 0.08,
    transition_grid = list(
      background = 0.03,
      restart = 0.02,
      advance = 0.45,
      background_stay = 0.9
    ),
    warp = warp
  )
}

test_that("duration grids are invariant to adjacent identical splitting", {
  original <- make_gaze_fixations(
    rbind(c(0, 0), c(2, 1)),
    duration = c(2, 2),
    onset = c(0, 2)
  )
  split <- make_gaze_fixations(
    rbind(c(0, 0), c(0, 0), c(2, 1)),
    duration = c(1, 1, 2),
    onset = c(0, 1, 2)
  )
  chronology <- gaze_local_order()
  original_grid <- eyesim:::gaze_duration_grid(
    eyesim:::as_gaze_measure(original, chronology), 32
  )
  split_grid <- eyesim:::gaze_duration_grid(
    eyesim:::as_gaze_measure(split, chronology), 32
  )

  expect_equal(original_grid$coords, split_grid$coords, tolerance = 0)
  expect_equal(original_grid$time, split_grid$time, tolerance = 0)
})

test_that("Replay evidence scale is invariant to duration-grid replication", {
  base_scores <- list(
    c(-0.4, -0.9, -1.2),
    c(-1.1, -0.3, -0.8),
    c(-0.7, -1.0, -0.2),
    c(-0.2, -0.6, -1.3)
  )
  truth <- c(1L, 2L, 3L, 1L)
  normalized_16 <- lapply(base_scores, function(score) {
    vapply(score * 16, eyesim:::gaze_replay_mean_log_score,
           numeric(1), grid_size = 16)
  })
  normalized_64 <- lapply(base_scores, function(score) {
    vapply(score * 64, eyesim:::gaze_replay_mean_log_score,
           numeric(1), grid_size = 64)
  })

  expect_equal(normalized_16, base_scores, tolerance = 0)
  expect_equal(normalized_64, base_scores, tolerance = 0)
  fit_16 <- eyesim:::fit_gaze_temperature(normalized_16, truth)
  fit_64 <- eyesim:::fit_gaze_temperature(normalized_64, truth)
  expect_equal(fit_16$temperature, fit_64$temperature, tolerance = 1e-12)
})

test_that("Replay retains raw likelihood while scoring mean log density", {
  reference <- make_gaze_fixations(rbind(c(0, 0), c(1, 1), c(2, 0)))
  training <- tibble::tibble(image_id = 1L, fixgroup = list(reference))
  model <- fit_gaze_replay_model(
    training, training, match_on = "image_id",
    spec = make_replay_test_spec(grid_size = 20)
  )
  result <- gaze_replay_align(reference, reference, model)

  expect_equal(
    result$log_score,
    result$alignment$log_likelihood / model$spec$grid_size,
    tolerance = 0
  )
  expect_equal(
    result$provenance$raw_log_likelihood,
    result$alignment$log_likelihood,
    tolerance = 0
  )
  expect_identical(
    result$provenance$score_semantics,
    "mean_log_likelihood_per_duration_bin"
  )
})

test_that("Replay transitions normalize and preserve declared components", {
  parameters <- list(
    background = 0.05,
    restart = 0.1,
    advance = 0.3,
    background_stay = 0.9
  )
  transition <- eyesim:::gaze_replay_transition(
    c(0.2, 0.3, 0.5), parameters, max_skip = 2
  )

  expect_equal(rowSums(transition$transition), rep(1, 4), tolerance = 1e-14)
  expect_equal(sum(transition$initial), 1, tolerance = 1e-14)
  expect_equal(transition$transition[-1, 1], rep(0.05, 3), tolerance = 1e-14)
  expect_equal(
    rowSums(transition$restart_component)[-1],
    rep(0.1, 3),
    tolerance = 1e-14
  )
})

test_that("log-space forward likelihood matches exhaustive state enumeration", {
  initial <- c(0.35, 0.65)
  transition <- matrix(c(0.8, 0.2, 0.25, 0.75), 2, byrow = TRUE)
  log_emission <- log(matrix(c(
    0.7, 0.2,
    0.4, 0.8,
    0.9, 0.3
  ), 3, 2, byrow = TRUE))
  fit <- eyesim:::gaze_hmm_forward_backward(log_emission, initial, transition)
  paths <- expand.grid(rep(list(1:2), 3))
  path_probability <- apply(paths, 1, function(states) {
    probability <- initial[[states[[1]]]] * exp(log_emission[1, states[[1]]])
    for (time in 2:3) {
      probability <- probability *
        transition[states[[time - 1]], states[[time]]] *
        exp(log_emission[time, states[[time]]])
    }
    probability
  })

  expect_equal(exp(fit$log_likelihood), sum(path_probability), tolerance = 1e-14)
  expect_equal(rowSums(fit$posterior), rep(1, 3), tolerance = 1e-14)
  expect_equal(
    apply(fit$transition_posterior, 1, sum),
    rep(1, 2),
    tolerance = 1e-14
  )
})

test_that("low-rank Replay recursion matches the dense HMM oracle", {
  parameters <- list(
    background = 0.06,
    restart = 0.08,
    advance = 0.3,
    background_stay = 0.88
  )
  transition <- eyesim:::gaze_replay_transition(
    c(0.25, 0.35, 0.4), parameters, max_skip = 2
  )
  set.seed(33)
  log_emission <- matrix(stats::rnorm(28), nrow = 7, ncol = 4)
  dense <- eyesim:::gaze_hmm_forward_backward(
    log_emission, transition$initial, transition$transition
  )
  low_rank <- eyesim:::gaze_replay_forward_backward(log_emission, transition)

  expect_equal(low_rank$log_likelihood, dense$log_likelihood, tolerance = 1e-12)
  expect_equal(low_rank$posterior, dense$posterior, tolerance = 1e-12)
  expect_match(low_rank$complexity, "max_skip")
  expect_lte(
    sum(transition$local_component[-1, -1] > 0),
    3 * (2 + 1)
  )
})

test_that("Replay is directional and penalizes complete reversal", {
  reference <- make_gaze_fixations(rbind(
    c(0, 0), c(0.8, 1.1), c(2, 0.3), c(3.2, 1.7)
  ))
  ref_tab <- tibble::tibble(image_id = 1L, fixgroup = list(reference))
  source_tab <- tibble::tibble(image_id = 1L, fixgroup = list(reference))
  model <- fit_gaze_replay_model(
    ref_tab,
    source_tab,
    match_on = "image_id",
    spec = make_replay_test_spec(grid_size = 24)
  )
  forward <- gaze_replay_align(reference, reference, model, "forward")
  reverse_path <- make_gaze_fixations(
    cbind(rev(reference$x), rev(reference$y))
  )
  reverse <- gaze_replay_align(reference, reverse_path, model, "reverse")

  expect_gt(forward$log_score, reverse$log_score)
  expect_gt(forward$diagnostics$replay_coverage, 0.8)
  expect_true(forward$convergence$converged)
})

test_that("Replay background absorbs gaze far from every encoding state", {
  reference <- make_gaze_fixations(rbind(c(0, 0), c(1, 1), c(2, 0)))
  background_training <- make_gaze_fixations(rbind(c(8, 8), c(9, 9), c(10, 8)))
  model <- fit_gaze_replay_model(
    tibble::tibble(image_id = 1L, fixgroup = list(reference)),
    tibble::tibble(image_id = 1L, fixgroup = list(background_training)),
    match_on = "image_id",
    spec = make_replay_test_spec(grid_size = 20)
  )
  result <- gaze_replay_align(reference, background_training, model)

  expect_gt(result$diagnostics$background_coverage, 0.5)
  expect_lt(result$diagnostics$replay_coverage, 0.5)
})

test_that("Replay candidate evidence is invariant to candidate evaluation order", {
  true_path <- make_gaze_fixations(rbind(c(0, 0), c(1, 1), c(2, 0)))
  wrong_path <- make_gaze_fixations(rbind(c(5, 5), c(6, 4), c(7, 5)))
  model <- fit_gaze_replay_model(
    tibble::tibble(item = "true", fixgroup = list(true_path)),
    tibble::tibble(item = "true", fixgroup = list(true_path)),
    match_on = "item",
    spec = make_replay_test_spec()
  )
  true_result <- gaze_replay_align(true_path, true_path, model, "true")
  wrong_result <- gaze_replay_align(wrong_path, true_path, model, "wrong")
  forward <- eyesim:::score_gaze_engine_results(
    list(true_result, wrong_result), "true"
  )
  reversed <- eyesim:::score_gaze_engine_results(
    list(wrong_result, true_result), "true"
  )

  expect_equal(forward$gaze_info_bits, reversed$gaze_info_bits, tolerance = 1e-14)
  expect_equal(forward$posterior_true, reversed$posterior_true, tolerance = 1e-14)
})

test_that("Replay uses an explicit identity fallback for tiny calibration sets", {
  first <- make_gaze_fixations(rbind(c(0, 0), c(1, 1), c(2, 0)))
  second <- make_gaze_fixations(rbind(c(5, 5), c(6, 4), c(7, 5)))
  reference <- tibble::tibble(
    participant = "p1", item = c("a", "b"),
    fixgroup = list(first, second)
  )
  source <- reference
  model <- fit_gaze_replay_model(
    reference, source,
    match_on = c("participant", "item"),
    contrast_on = "participant",
    spec = make_replay_test_spec()
  )
  results <- list(
    gaze_replay_align(first, first, model, "a"),
    gaze_replay_align(second, first, model, "b")
  )
  evidence <- eyesim:::score_gaze_engine_results(
    results, "a", temperature = model$temperature
  )

  expect_s3_class(model$calibration, "gaze_temperature_fit")
  expect_identical(model$calibration$scheme, "identity_fallback")
  expect_equal(model$temperature, 1)
  expect_equal(evidence$temperature, model$temperature)
  expect_equal(model$version, 3L)
})

test_that("Replay temperature calibration is inner-cross-fitted by item", {
  tabs <- make_gaze_weave_cv_tables()
  keep <- tabs$ref_tab$participant == "p1"
  model <- fit_gaze_replay_model(
    tabs$ref_tab[keep, ],
    tabs$source_tab[tabs$source_tab$participant == "p1", ],
    match_on = c("participant", "image_id"),
    contrast_on = "participant",
    spec = make_replay_test_spec(grid_size = 12)
  )

  expect_identical(
    model$calibration$scheme,
    "inner_cross_fitted_candidate_log_loss"
  )
  expect_true(all(vapply(
    model$calibration$folds, `[[`, integer(1), "overlap_match_n"
  ) == 0L))
  expect_true(is.finite(model$temperature) && model$temperature > 0)
  expect_false(model$calibration$at_boundary)
})

test_that("cross-fitted Replay uses disjoint item keys and normalized candidates", {
  tabs <- make_gaze_weave_cv_tables()
  spec <- make_replay_test_spec(
    warp = gaze_warp_contraction(
      center = c(0, 0),
      translation = TRUE,
      fit_by = "participant"
    ),
    grid_size = 12
  )
  fit <- gaze_replay_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = c("participant", "image_id"),
    contrast_on = "participant",
    split_on = c("participant", "image_id"),
    n_folds = 2,
    seed = 17,
    spec = spec
  )

  expect_s3_class(fit, "gaze_replay_fit")
  expect_true(all(vapply(fit$folds, `[[`, integer(1), "overlap_match_n") == 0L))
  expect_true(all(fit$results$candidate_count == 3L))
  expect_true(all(fit$results$all_converged))
  expect_true(all(is.finite(fit$results$gaze_info_bits)))
  expect_true(all(vapply(fit$results$candidates, function(candidates) {
    abs(sum(candidates$posterior) - 1) < 1e-12
  }, logical(1))))
  expect_true(all(fit$results$template_rank == 1))
  expect_true(all(vapply(fit$folds, function(fold) {
    inherits(fold$calibration, "gaze_temperature_fit") &&
      is.finite(fold$calibration$temperature)
  }, logical(1))))

  scales <- unlist(lapply(fit$folds, function(fold) {
    vapply(fold$warp$groups, `[[`, numeric(1), "scale")
  }))
  expect_equal(scales, rep(tabs$scale, length(scales)), tolerance = 1e-5)
})

test_that("ordered training selects the lower restart candidate", {
  path <- make_gaze_fixations(rbind(
    c(0, 0), c(1, 0.8), c(2, 0.1), c(3, 1), c(4, 0.2)
  ))
  spec <- gaze_replay_spec(
    grid_size = 30,
    max_skip = 2,
    scale_floor = 0.05,
    transition_grid = list(
      background = 0.02,
      restart = c(0.01, 0.4),
      advance = 0.5,
      background_stay = 0.9
    )
  )
  model <- fit_gaze_replay_model(
    tibble::tibble(image_id = 1L, fixgroup = list(path)),
    tibble::tibble(image_id = 1L, fixgroup = list(path)),
    match_on = "image_id",
    spec = spec
  )

  expect_equal(model$parameters$restart, 0.01)
})

test_that("splitting an identical initial state preserves one-step evidence", {
  parameters <- list(
    background = 0.05,
    restart = 0.05,
    advance = 0.3,
    background_stay = 0.9
  )
  compact <- eyesim:::gaze_replay_transition(c(0.5, 0.5), parameters, 2)
  refined <- eyesim:::gaze_replay_transition(c(0.25, 0.25, 0.5), parameters, 2)
  compact_emission <- matrix(log(c(0.1, 0.8, 0.3)), nrow = 1)
  refined_emission <- matrix(log(c(0.1, 0.8, 0.8, 0.3)), nrow = 1)
  compact_fit <- eyesim:::gaze_replay_forward_backward(
    compact_emission, compact
  )
  refined_fit <- eyesim:::gaze_replay_forward_backward(
    refined_emission, refined
  )

  expect_equal(
    compact_fit$log_likelihood,
    refined_fit$log_likelihood,
    tolerance = 1e-14
  )
})
