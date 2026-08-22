make_episode_test_spec <- function(reliability = "none", grid_size = 12L) {
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
    reliability = reliability,
    reliability_kappa_bounds = c(0, 100)
  )
}

make_episode_paths <- function(offset = c(0, 0)) {
  lapply(0:3, function(presentation) {
    make_gaze_fixations(
      rbind(c(0, 0), c(1, 1), c(2, 0)) +
        matrix(offset + c(0.03, -0.02) * presentation,
               nrow = 3L, ncol = 2L, byrow = TRUE),
      duration = c(1, 2, 1),
      onset = c(0, 1, 3)
    )
  })
}

make_episode_tables <- function() {
  a <- make_episode_paths(c(0, 0))
  b <- make_episode_paths(c(8, 7))
  reference <- tibble::tibble(
    item = rep(c("a", "b"), each = 4L),
    presentation = rep(1:4, 2L),
    fixgroup = c(a, b)
  )
  source <- tibble::tibble(
    item = c("a", "b"),
    fixgroup = list(a[[4L]], b[[4L]])
  )
  list(reference = reference, source = source, a = a, b = b)
}

test_that("equal template mixture matches an analytic oracle", {
  observed <- eyesim:::gaze_replay_log_mixture(log(c(0.2, 0.8)))
  shifted <- eyesim:::gaze_replay_log_mixture(log(c(0.2, 0.8)) + 100)

  expect_equal(observed$log_score, log(0.5), tolerance = 1e-14)
  expect_equal(observed$posterior_weight, c(0.2, 0.8), tolerance = 1e-14)
  expect_equal(shifted$log_score, observed$log_score + 100, tolerance = 1e-13)
  expect_equal(
    shifted$posterior_weight, observed$posterior_weight,
    tolerance = 1e-13
  )
})

test_that("path quality is invariant to adjacent fixation refinement", {
  original <- make_gaze_fixations(
    rbind(c(0, 0), c(2, 1)),
    duration = c(2, 2), onset = c(0, 2)
  )
  split <- make_gaze_fixations(
    rbind(c(0, 0), c(0, 0), c(2, 1)),
    duration = c(1, 1, 2), onset = c(0, 1, 2)
  )
  first <- eyesim:::gaze_replay_path_quality(original)
  second <- eyesim:::gaze_replay_path_quality(split)

  expect_identical(first$raw_fixation_count, 2L)
  expect_identical(second$raw_fixation_count, 3L)
  expect_identical(first$coalesced_fixation_count, 2L)
  expect_identical(second$coalesced_fixation_count, 2L)
  expect_equal(first$effective_fixations, 2, tolerance = 0)
  expect_equal(second$effective_fixations, first$effective_fixations,
               tolerance = 0)
  expect_equal(second$total_duration, first$total_duration, tolerance = 0)
})

test_that("reliability shrinkage is monotone and has exact boundaries", {
  support <- c(1, 2, 5, 10)
  weight <- eyesim:::gaze_replay_reliability(support, kappa = 3)

  expect_equal(weight[[1L]], 0.25, tolerance = 0)
  expect_true(all(diff(weight) > 0))
  expect_equal(
    eyesim:::gaze_replay_reliability(support, kappa = 0),
    rep(1, length(support)), tolerance = 0
  )
  expect_error(
    eyesim:::gaze_replay_reliability(c(0, 1), kappa = 1),
    "positive"
  )
})

test_that("zero reliability returns the declared candidate prior", {
  base <- eyesim:::score_gaze_candidates(
    c(2, 0), true_index = 1L, prior = c(0.7, 0.3), reliability = 1
  )
  abstained <- eyesim:::score_gaze_candidates(
    c(2, 0), true_index = 1L, prior = c(0.7, 0.3), reliability = 0
  )
  partial <- eyesim:::score_gaze_candidates(
    c(2, 0), true_index = 1L, prior = c(0.7, 0.3), reliability = 0.4
  )

  expect_equal(abstained$candidates$posterior, c(0.7, 0.3), tolerance = 0)
  expect_equal(abstained$gaze_info_bits, 0, tolerance = 1e-14)
  expect_equal(
    partial$candidates$posterior,
    0.4 * base$candidates$posterior + 0.6 * c(0.7, 0.3),
    tolerance = 1e-14
  )
  expect_equal(partial$reliability, 0.4, tolerance = 0)
})

test_that("joint calibration can abstain on selectively noisy paths", {
  score_sets <- c(
    rep(list(c(-4, 0)), 12L),
    rep(list(c(4, 0)), 12L)
  )
  support <- c(rep(1, 12L), rep(8, 12L))
  fit <- eyesim:::fit_gaze_replay_reliability(
    score_sets,
    true_index = rep(1L, length(score_sets)),
    effective_fixations = support,
    temperature_bounds = c(0.05, 20),
    kappa_bounds = c(0, 100)
  )

  expect_s3_class(fit, "gaze_reliability_fit")
  expect_gt(fit$kappa, 0)
  expect_lte(fit$log_loss, fit$temperature_only_log_loss + 1e-8)
  expect_lt(
    eyesim:::gaze_replay_reliability(1, fit$kappa),
    eyesim:::gaze_replay_reliability(8, fit$kappa)
  )
})

test_that("calibration parameter guard clamps boundary probes", {
  lower <- c(log(0.05), 0)
  upper <- c(log(100), log1p(100))

  expect_null(eyesim:::gaze_replay_calibration_parameters(
    c(NA_real_, 0), lower, upper
  ))
  expect_null(eyesim:::gaze_replay_calibration_parameters(
    c(0, Inf), lower, upper
  ))
  observed <- eyesim:::gaze_replay_calibration_parameters(
    c(0, -.Machine$double.eps), lower, upper
  )
  expect_equal(observed$temperature, 1, tolerance = 0)
  expect_equal(observed$kappa, 0, tolerance = 0)
})

test_that("one-template episode is exactly the ordinary Replay score", {
  tabs <- make_episode_tables()
  model <- fit_gaze_replay_model(
    tabs$reference[tabs$reference$presentation == 4L, ],
    tabs$source,
    match_on = "item",
    spec = make_episode_test_spec()
  )
  ordinary <- gaze_replay_align(
    tabs$a[[4L]], tabs$a[[4L]], model, candidate_key = "a"
  )
  episode <- gaze_replay_align_episode(
    list(tabs$a[[4L]]), tabs$a[[4L]], model,
    candidate_key = "a", template_key = "p4"
  )

  expect_equal(episode$log_score, ordinary$log_score, tolerance = 0)
  expect_equal(episode$diagnostics$replay_coverage,
               ordinary$diagnostics$replay_coverage, tolerance = 0)
  expect_identical(episode$alignment$template_count, 1L)
  expect_equal(episode$alignment$template_posterior, 1, tolerance = 0)
})

test_that("four-presentation episodes preserve separate paths", {
  tabs <- make_episode_tables()
  model <- fit_gaze_replay_model(
    tabs$reference, tabs$source,
    match_on = "item", template_on = "presentation",
    spec = make_episode_test_spec()
  )
  episode <- gaze_replay_align_episode(
    tabs$a, tabs$a[[4L]], model,
    candidate_key = "a", template_key = paste0("p", 1:4)
  )

  expect_identical(model$template_on, "presentation")
  expect_identical(model$version, 4L)
  expect_identical(episode$alignment$template_count, 4L)
  expect_length(episode$alignment$components, 4L)
  expect_equal(sum(episode$alignment$template_posterior), 1,
               tolerance = 1e-14)
  expect_identical(
    vapply(
      episode$alignment$components,
      function(value) value$diagnostics$template_state_count,
      integer(1)
    ),
    rep(3L, 4L)
  )
  expect_equal(episode$diagnostics$template_effective_count,
               exp(-sum(episode$alignment$template_posterior *
                         log(episode$alignment$template_posterior))),
               tolerance = 1e-14)
})

test_that("episode score is invariant to presentation row order", {
  tabs <- make_episode_tables()
  model <- fit_gaze_replay_model(
    tabs$reference, tabs$source,
    match_on = "item", template_on = "presentation",
    spec = make_episode_test_spec()
  )
  forward <- gaze_replay_align_episode(
    tabs$a, tabs$a[[4L]], model,
    candidate_key = "a", template_key = paste0("p", 1:4)
  )
  order <- c(4, 2, 1, 3)
  permuted <- gaze_replay_align_episode(
    tabs$a[order], tabs$a[[4L]], model,
    candidate_key = "a", template_key = paste0("p", 1:4)[order]
  )

  expect_equal(forward$log_score, permuted$log_score, tolerance = 1e-14)
  expect_equal(
    forward$alignment$template_posterior[
      match(permuted$alignment$template_key,
            forward$alignment$template_key)
    ],
    permuted$alignment$template_posterior,
    tolerance = 1e-14
  )
})

test_that("row scoring groups four templates into one item candidate", {
  tabs <- make_episode_tables()
  model <- fit_gaze_replay_model(
    tabs$reference, tabs$source,
    match_on = "item", template_on = "presentation",
    spec = make_episode_test_spec()
  )
  scored <- eyesim:::score_gaze_replay_row(
    tabs$source[1, ], tabs$reference,
    match_on = "item", contrast_on = NULL,
    refvar = "fixgroup", sourcevar = "fixgroup", model = model
  )

  expect_identical(scored$evidence$candidate_count, 2L)
  expect_identical(nrow(scored$candidates), 2L)
  expect_true(scored$candidates$is_true[scored$candidates$item == "a"])
  expect_identical(scored$alignment$template_count, 4L)
  expect_identical(scored$alignment$template_key, as.character(1:4))
  expect_identical(nrow(scored$candidate_components), 8L)
  expect_true(all(c(
    "template_raw_fixation_count", "template_effective_fixations",
    "template_total_duration"
  ) %in% names(scored$candidate_components)))
  expect_true(all(scored$candidate_components$template_effective_fixations > 0))
  expect_true(scored$all_converged)
})

test_that("episode template keys must be complete and unique", {
  tabs <- make_episode_tables()
  model <- fit_gaze_replay_model(
    tabs$reference, tabs$source,
    match_on = "item", template_on = "presentation",
    spec = make_episode_test_spec()
  )

  expect_error(
    gaze_replay_align_episode(
      tabs$a, tabs$a[[4L]], model,
      candidate_key = "a", template_key = rep("same", 4L)
    ),
    "unique"
  )
})

test_that("cross-fitted episode Replay keeps items disjoint and templates grouped", {
  tabs <- make_gaze_weave_cv_tables()
  references <- dplyr::bind_rows(lapply(1:4, function(presentation) {
    part <- tabs$ref_tab
    part$presentation <- presentation
    part$fixgroup <- lapply(part$fixgroup, function(path) {
      shifted <- path
      shifted$x <- shifted$x + 0.01 * presentation
      shifted$y <- shifted$y - 0.005 * presentation
      shifted
    })
    part
  }))
  fit <- gaze_replay_cv(
    references,
    tabs$source_tab,
    match_on = c("participant", "image_id"),
    contrast_on = "participant",
    split_on = c("participant", "image_id"),
    template_on = "presentation",
    n_folds = 2,
    seed = 41,
    spec = make_episode_test_spec(
      reliability = "effective_fixations", grid_size = 8L
    )
  )

  expect_s3_class(fit, "gaze_replay_fit")
  expect_true(all(vapply(
    fit$folds, `[[`, integer(1), "overlap_match_n"
  ) == 0L))
  expect_true(all(fit$results$candidate_count == 3L))
  expect_true(all(fit$results$reliability >= 0 &
                  fit$results$reliability <= 1))
  expect_true(all(is.finite(fit$results$base_gaze_info_bits)))
  expect_true(all(is.finite(fit$results$gaze_info_bits)))
  expect_true(all(fit$results$effective_fixations >= 1))
  expect_true(all(vapply(fit$results$alignment, function(alignment) {
    inherits(alignment, "gaze_replay_episode_alignment") &&
      alignment$template_count == 4L
  }, logical(1))))
  expect_true(all(vapply(fit$folds, function(fold) {
    inherits(fold$calibration, "gaze_reliability_fit") &&
      fold$calibration$kappa >= 0
  }, logical(1))))
})
