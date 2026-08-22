item_effect_source <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    gaze_weave_test_inst_path("validation", "gaze-weave-item-effects.R"),
    envir = environment
  )
  environment
}

test_that("empirical-Bayes item shrinkage has analytic boundary behavior", {
  court <- item_effect_source()
  no_signal <- court$item_effect_shrink(rep(0, 8), rep(letters[1:4], each = 2))
  signal <- court$item_effect_shrink(
    rep(c(-1, 1), each = 4), rep(c("a", "b"), each = 4)
  )

  expect_true(all(no_signal$item_propensity == 0))
  expect_true(all(no_signal$shrinkage_weight == 0))
  expect_equal(signal$item_propensity, c(-1, 1), tolerance = 1e-12)
  expect_equal(signal$shrinkage_weight, c(1, 1), tolerance = 1e-12)
})

item_effect_fixture <- function(seed = 11L) {
  set.seed(seed)
  tab <- expand.grid(
    participant = paste0("p", 1:12),
    item = 1:8,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  tab$probe_type <- rep(c("old", "lure"), length.out = nrow(tab))
  tab$degradation <- rep(c(20, 40, 60, 80), length.out = nrow(tab))
  tab$saliency_z <- (tab$degradation - 60) / 20
  tab$effective_fixations <- 5
  tab$total_duration <- 2000
  tab$trial_key <- paste(tab$participant, tab$item, sep = ":")
  item_value <- seq(-0.7, 0.7, length.out = 8)
  tab$score_transport <- item_value[tab$item] + stats::rnorm(nrow(tab), sd = 0.1)
  tab$score_density_sigma_80 <- stats::rnorm(nrow(tab), sd = 0.1)
  tab
}

test_that("cross-participant item prediction excludes self and finds an oracle", {
  court <- item_effect_source()
  fixture <- item_effect_fixture()
  observed <- court$item_effect_crossparticipant(
    fixture, "score_transport"
  )
  shuffled <- fixture
  for (participant in unique(shuffled$participant)) {
    index <- shuffled$participant == participant
    shuffled$item[index] <- sample(shuffled$item[index])
  }
  shuffled$trial_key <- paste(shuffled$participant, shuffled$item, sep = ":")
  permuted <- court$item_effect_crossparticipant(
    shuffled, "score_transport"
  )

  expect_equal(nrow(observed), nrow(fixture))
  expect_true(all(observed$self_participant_overlap == 0L))
  expect_true(all(observed$source_participants == 11L))
  expect_gt(sum(observed$squared_gain), 0)
  expect_gt(sum(observed$squared_gain), sum(permuted$squared_gain))
})

test_that("held-out values cannot change their own item propensity", {
  court <- item_effect_source()
  fixture <- item_effect_fixture()
  observed <- court$item_effect_crossparticipant(fixture, "score_transport")
  mutated <- fixture
  mutated$score_transport[mutated$participant == "p1"] <-
    mutated$score_transport[mutated$participant == "p1"] + 100
  rerun <- court$item_effect_crossparticipant(mutated, "score_transport")

  expect_equal(
    observed$item_propensity[observed$participant == "p1"],
    rerun$item_propensity[rerun$participant == "p1"],
    tolerance = 1e-12
  )
})

test_that("density residualization preserves a Transport-specific item effect", {
  court <- item_effect_source()
  fixture <- item_effect_fixture()
  fixture$score_transport <- 2 * fixture$score_density_sigma_80 +
    rep(seq(-0.7, 0.7, length.out = 8), each = 12) +
    stats::rnorm(nrow(fixture), sd = 0.05)
  observed <- court$item_effect_crossparticipant(
    fixture, "score_transport", density = "score_density_sigma_80"
  )

  expect_gt(sum(observed$squared_gain), 0)
  expect_gt(stats::cor(observed$item_propensity, observed$residual), 0.8)
})

test_that("study transfer learns its scale from training participants", {
  court <- item_effect_source()
  retrieval <- item_effect_fixture()
  item_value <- seq(-0.7, 0.7, length.out = 8)
  retrieval$score_transport <- item_value[retrieval$item] +
    stats::rnorm(nrow(retrieval), sd = 0.05)
  study <- retrieval
  study$score_study <- 10 * item_value[study$item] +
    stats::rnorm(nrow(study), sd = 0.1)
  study$study_effective_fixations <- study$effective_fixations

  observed <- court$item_effect_study_transfer(retrieval, study)

  expect_gt(sum(observed$squared_gain), 0)
  expect_gt(stats::cor(observed$item_propensity, observed$residual), 0.8)
  expect_true(all(observed$self_participant_overlap == 0L))
})

test_that("safe standardization handles zero item variance", {
  court <- item_effect_source()
  train <- data.frame(value = rep(2, 5))
  evaluate <- data.frame(value = c(1, 3))
  result <- court$item_effect_safe_standardize(
    train, evaluate, "value", "z_value"
  )

  expect_identical(result$train$z_value, rep(0, 5))
  expect_identical(result$evaluate$z_value, c(0, 0))
})

test_that("local repeated-viewing and retrieval supports align", {
  retrieval_dir <- gaze_weave_test_inst_path(
    "validation", "gaze-weave-transport-v3-retrieval-results"
  )
  comparator_dir <- gaze_weave_test_inst_path(
    "validation", "gaze-weave-recognition-full-cohort-results"
  )
  study_dir <- gaze_weave_test_inst_path(
    "validation", "gaze-weave-transport-v3-repeated-viewing-results"
  )
  testthat::skip_if_not(
    file.exists(file.path(retrieval_dir, "checkpoint-manifest.csv")) &&
      file.exists(file.path(comparator_dir, "manifest-md5.csv")) &&
      file.exists(file.path(
        study_dir, "checkpoint-study_p3_p4-fold-04.rds"
      )),
    "Local-only item-effect inputs are unavailable."
  )
  court <- item_effect_source()
  panel <- court$transport_density_read_panel(retrieval_dir, comparator_dir)
  retrieval <- court$transport_density_wide(panel)
  study <- court$item_effect_read_study(study_dir, retrieval)

  expect_equal(nrow(retrieval), 1295L)
  expect_equal(nrow(study), 1295L)
  expect_identical(
    sort(retrieval$trial_key), sort(study$trial_key)
  )
  expect_equal(length(unique(study$item)), 120L)
})
