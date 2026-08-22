source(
  gaze_weave_test_inst_path(
    "validation", "gaze-weave-pcmri-multi-replay-behavior.R"
  ),
  local = TRUE
)

test_that("behavior response coding preserves missing and confidence bins", {
  expect_equal(
    multi_replay_behavior_response(c(0, 1, 2, 3, 4, NA)),
    c(NA_integer_, 1L, 1L, 0L, 0L, NA_integer_)
  )
})

test_that("behavioral standardization uses training rows only", {
  train <- data.frame(value = c(2, 4, 6))
  evaluate <- data.frame(value = c(4, 8))
  observed <- multi_replay_behavior_standardize(
    train, evaluate, "value", "z_value"
  )

  expect_equal(observed$center, 4, tolerance = 0)
  expect_equal(observed$scale, 2, tolerance = 0)
  expect_equal(observed$train$z_value, c(-1, 0, 1), tolerance = 0)
  expect_equal(observed$evaluate$z_value, c(0, 2), tolerance = 0)
})

test_that("behavioral fold prediction enforces crossed support", {
  make_rows <- function(participants, items) {
    grid <- expand.grid(
      participant = participants, item = items,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    grid$sal10 <- sample(c(-2, -1, 0, 1, 2), nrow(grid), replace = TRUE)
    grid$gaze_info_bits <- stats::rnorm(nrow(grid))
    grid$effective_fixations <- stats::runif(nrow(grid), 1, 8)
    grid$total_duration <- stats::runif(nrow(grid), 500, 2900)
    grid$candidate_count <- sample(5:12, nrow(grid), replace = TRUE)
    probability <- stats::plogis(
      -0.2 + 0.3 * grid$sal10 + 0.4 * grid$gaze_info_bits
    )
    grid$said_old <- stats::rbinom(nrow(grid), 1, probability)
    grid
  }
  set.seed(91)
  train <- make_rows(paste0("p", 1:5), 1:10)
  evaluate <- make_rows(c("p6", "p7"), 11:15)
  observed <- multi_replay_behavior_fit_fold(
    train, evaluate, "synthetic", 1L
  )

  expect_identical(nrow(observed), nrow(evaluate))
  expect_true(all(is.finite(observed$behavior_info_bits)))
  expect_true(all(observed$probability_base > 0 &
                  observed$probability_base < 1))
  expect_true(all(observed$probability_full > 0 &
                  observed$probability_full < 1))
  expect_identical(observed$train_participants[[1L]], 5L)
  expect_identical(observed$eval_participants[[1L]], 2L)

  expect_error(
    multi_replay_behavior_fit_fold(
      train, train[1:3, ], "synthetic", 1L
    ),
    "overlaps"
  )
})

test_that("frozen local measurement yields the prespecified old-item support", {
  root <- normalizePath(testthat::test_path("..", ".."))
  measurement_dir <- file.path(
    root, "inst", "validation",
    "gaze-weave-pcmri-catalog-replication-results", "multi-replay"
  )
  skip_if_not(
    file.exists(file.path(measurement_dir, "measurement-verdict.csv")),
    "Frozen local measurement court unavailable"
  )
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  scores <- multi_replay_behavior_scores(measurement_dir, verify = TRUE)
  support <- table(scores$method)

  expect_true(all(support[multi_replay_behavior_methods] == 989L))
  primary <- scores[scores$method == multi_replay_behavior_primary, ]
  expect_identical(sum(primary$said_old == 0L), 97L)
  expect_identical(sum(primary$said_old == 1L), 892L)
  expect_true(all(primary$probe_type == "old"))
})

test_that("local behavioral results obey their manifest", {
  result_dir <- file.path(
    testthat::test_path("..", ".."), "inst", "validation",
    "gaze-weave-pcmri-catalog-replication-results", "multi-replay",
    "behavior"
  )
  manifest_path <- file.path(result_dir, "manifest-md5.csv")
  skip_if_not(file.exists(manifest_path), "Behavioral court not finalized locally")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  paths <- file.path(result_dir, manifest$file)

  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), manifest$md5)
})
