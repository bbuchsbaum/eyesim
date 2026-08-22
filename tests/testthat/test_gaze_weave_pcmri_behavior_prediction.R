source(
  gaze_weave_test_inst_path(
    "validation", "gaze-weave-pcmri-behavior-prediction.R"
  ),
  local = TRUE
)

test_that("behavioral log loss is finite and rewards the correct prediction", {
  good <- pcmri_behavior_log_loss(c(1, 0), c(.9, .1))
  bad <- pcmri_behavior_log_loss(c(1, 0), c(.1, .9))
  expect_true(all(is.finite(good)))
  expect_true(all(good < bad))
})

test_that("behavioral fold fitting forbids participant and item leakage", {
  train <- data.frame(
    participant = rep(c("a", "b"), each = 4),
    item = rep(1:4, 2),
    said_old = rep(c(0, 1), 4),
    sal10 = rep(c(-1, 1), 4),
    gaze_info_bits = seq(-.4, .3, length.out = 8)
  )
  evaluate <- data.frame(
    participant = rep("c", 4), item = 5:8,
    said_old = c(0, 1, 0, 1), sal10 = c(-1, -1, 1, 1),
    gaze_info_bits = c(-.2, .1, -.1, .2)
  )
  result <- suppressWarnings(pcmri_behavior_fit_fold(
    train, evaluate, "replay", 1L
  ))
  expect_equal(nrow(result), 4L)
  expect_true(all(is.finite(result$behavior_info_bits)))
  expect_equal(unique(result$train_n), 8L)
  expect_error(
    pcmri_behavior_fit_fold(
      train, transform(evaluate, participant = "a"), "replay", 1L
    ),
    "overlaps"
  )
})

test_that("crossed weighted means preserve paired method differences", {
  tab <- data.frame(
    participant = rep(c("a", "b"), each = 2),
    item = rep(1:2, 2),
    trial_key = paste(rep(c("a", "b"), each = 2), rep(1:2, 2), sep = ":")
  )
  plan <- recognition_bootstrap_plan(tab, draws = 40L, seed = 91L)
  first <- pcmri_behavior_weighted_mean(rep(1, 4), plan$weights)
  second <- pcmri_behavior_weighted_mean(rep(.75, 4), plan$weights)
  expect_equal(first - second, rep(.25, 40), tolerance = 1e-12)
})
