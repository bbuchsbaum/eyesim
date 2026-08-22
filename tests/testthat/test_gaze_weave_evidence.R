test_that("candidate probabilities are stable and respect declared priors", {
  prior <- c(0.6, 0.3, 0.1)
  unchanged <- eyesim:::score_gaze_candidates(
    log_score = c(1000, 1000, 1000),
    true_index = 1,
    candidate_key = letters[1:3],
    prior = prior
  )
  concentrated <- eyesim:::gaze_candidate_probabilities(
    c(1000, 0, -1000), prior = prior
  )

  expect_equal(unchanged$candidates$posterior, prior, tolerance = 1e-12)
  expect_equal(unchanged$gaze_info_bits, 0, tolerance = 1e-12)
  expect_equal(unchanged$odds_bits, 0, tolerance = 1e-12)
  expect_equal(sum(concentrated), 1, tolerance = 1e-14)
  expect_equal(concentrated, c(1, 0, 0), tolerance = 1e-14)
})

test_that("information and odds bits obey their probability identities", {
  evidence <- eyesim:::score_gaze_candidates(
    log_score = log(c(4, 1, 1, 1)),
    true_index = 1,
    prior = rep(0.25, 4)
  )
  p_true <- evidence$posterior_true

  expect_equal(evidence$gaze_info_bits, log2(p_true / 0.25), tolerance = 1e-14)
  expect_equal(
    evidence$odds_bits,
    log2((p_true / (1 - p_true)) / (0.25 / 0.75)),
    tolerance = 1e-14
  )
  expect_equal(evidence$odds_bits, log2(4), tolerance = 1e-14)
  expect_equal(evidence$log_loss, -log(p_true), tolerance = 1e-14)
  expect_equal(evidence$candidate_count, 4L)
})

test_that("rank and top-one credit use tolerance-aware ties", {
  evidence <- eyesim:::score_gaze_candidates(
    log_score = c(0, 1e-14, -1e-14, -2),
    true_index = 1,
    tie_tolerance = 1e-10
  )

  expect_equal(evidence$template_rank, 2)
  expect_equal(evidence$top1_credit, 1 / 3)
  expect_equal(evidence$tied_candidates, 3L)
})

test_that("temperature fitting recovers held-out score scale", {
  set.seed(20260815)
  n <- 2000L
  true_temperature <- 2
  margin <- stats::rnorm(n, sd = 2)
  probability_one <- stats::plogis(margin / true_temperature)
  is_one <- stats::rbinom(n, 1, probability_one) == 1
  scores <- cbind(margin / 2, -margin / 2)
  truth <- ifelse(is_one, 1L, 2L)

  fit <- eyesim:::fit_gaze_temperature(scores, truth)

  expect_s3_class(fit, "gaze_temperature_fit")
  expect_equal(fit$temperature, true_temperature, tolerance = 0.25)
  expect_lt(fit$log_loss, fit$uncalibrated_log_loss)
  expect_false(fit$at_boundary)
})

test_that("temperature regularization avoids a separable small-sample boundary", {
  scores <- matrix(
    c(20, 0, 18, 0, 22, 0, 19, 0),
    ncol = 2,
    byrow = TRUE
  )
  fit <- eyesim:::fit_gaze_temperature(scores, rep(1L, 4))

  expect_false(fit$at_boundary)
  expect_gt(fit$temperature, fit$bounds[[1]])
  expect_true(is.finite(fit$penalized_objective))
})

test_that("information remains finite when posterior probabilities underflow", {
  evidence <- eyesim:::score_gaze_candidates(
    log_score = c(0, -2000), true_index = 2
  )

  expect_equal(evidence$posterior_true, 0)
  expect_true(is.finite(evidence$gaze_info_bits))
  expect_true(is.finite(evidence$log_loss))
  expect_lt(evidence$gaze_info_bits, -1000)
})

test_that("cross-fitted calibration has disjoint train and evaluation rows", {
  set.seed(91)
  n <- 400L
  margin <- stats::rnorm(n, sd = 1.5)
  truth <- ifelse(stats::runif(n) < stats::plogis(margin / 1.8), 1L, 2L)
  scores <- cbind(margin / 2, -margin / 2)
  folds <- rep(1:4, length.out = n)
  calibrated <- eyesim:::crossfit_gaze_temperature(scores, truth, folds)

  expect_s3_class(calibrated, "gaze_temperature_cv")
  expect_length(calibrated$evidence, n)
  expect_true(all(is.finite(calibrated$temperature)))
  expect_true(all(vapply(calibrated$folds, `[[`, integer(1), "overlap_n") == 0L))
  expect_true(all(vapply(calibrated$evidence, function(x) {
    is.finite(x$gaze_info_bits) && is.finite(x$log_loss)
  }, logical(1))))
})

test_that("reliability tables audit every candidate probability", {
  evidence <- lapply(seq_len(12), function(i) {
    eyesim:::score_gaze_candidates(
      log_score = c(i / 12, 0, -i / 24),
      true_index = 1
    )
  })
  reliability <- eyesim:::gaze_reliability_table(evidence, bins = 5)

  expect_equal(sum(reliability$n), 36L)
  expect_true(all(reliability$mean_probability >= 0 &
                    reliability$mean_probability <= 1))
  expect_true(all(reliability$observed_frequency >= 0 &
                    reliability$observed_frequency <= 1))
})
