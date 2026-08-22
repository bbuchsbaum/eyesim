v2_validation_script <- system.file(
  "validation", "gaze-weave-v2-court.R", package = "eyesim"
)
if (!nzchar(v2_validation_script)) {
  v2_validation_script <- testthat::test_path(
    "..", "..", "inst", "validation", "gaze-weave-v2-court.R"
  )
}
source(v2_validation_script, local = TRUE)

test_that("the final synthetic court is frozen before execution", {
  protocol <- system.file(
    "validation", "GAZEWEAVE-V2-COURT.md", package = "eyesim"
  )
  if (!nzchar(protocol)) {
    protocol <- testthat::test_path(
      "..", "..", "inst", "validation", "GAZEWEAVE-V2-COURT.md"
    )
  }
  text <- paste(readLines(protocol, warn = FALSE), collapse = "\n")

  expect_match(text, "Final seed: `20260816`", fixed = TRUE)
  expect_match(text, "No universal-superiority claim", fixed = TRUE)
  expect_error(
    run_gaze_weave_v2_court(seed = 99L, smoke = FALSE),
    "frozen to seed 20260816"
  )
})

test_that("court generators provide matched signal and independent null rows", {
  court <- v2_court_data(
    n_participants = 1, n_items = 4, n_replications = 2,
    n_fixations = 5, seed = 17
  )

  expect_equal(nrow(court$reference), 3 * 2 * 4)
  expect_equal(nrow(court$source), 2 * nrow(court$reference))
  expect_equal(as.integer(table(court$source$condition)), c(24L, 24L))
  expect_equal(court$truth$scale, 0.78)

  neutral <- court$reference$fixgroup[
    court$reference$family == "neutral_misspecified" &
      court$reference$replicate == 1
  ]
  signatures <- lapply(neutral, function(path) {
    coords <- cbind(path$x, path$y)
    coords[do.call(order, as.data.frame(coords)), , drop = FALSE]
  })
  expect_true(all(vapply(signatures[-1], identical, logical(1), signatures[[1]])))
})

test_that("fractional perturbations balance every declared main effect", {
  design <- v2_fractional_design()

  expect_equal(nrow(design), 16L)
  expect_equal(ncol(design), 13L)
  expect_true(all(vapply(design, function(value) {
    identical(sort(unique(value)), 0:1) && sum(value) == 8L
  }, logical(1))))
})

test_that("dimensionless spatial costs are exactly unit invariant", {
  first <- rbind(c(1.2, 3.4), c(5.6, 7.8), c(9.1, 2.3))
  second <- rbind(c(2.2, 2.4), c(5.1, 8.8))
  spatial <- gaze_gaussian_mixture(c(0.75, 1.5, 3))

  native <- eyesim:::gaze_spatial_cost(first, second, spatial)
  converted <- eyesim:::gaze_spatial_cost(
    first * 50, second * 50,
    gaze_gaussian_mixture(c(0.75, 1.5, 3) * 50)
  )

  expect_identical(native, converted)
})

test_that("baseline filters exclude null rows from all training folds", {
  court <- v2_court_data(
    n_participants = 1, n_items = 4, n_replications = 2,
    n_fixations = 5, seed = 17
  )
  spec <- gaze_baseline_spec(
    gaze_screen(24, 18, "deg"), c(1, 2), density_grid = 8,
    methods = "density", lambda_grid = c(0.1, 1), inner_folds = 2
  )
  match_on <- c("participant", "family", "replicate", "item")
  fit <- gaze_baseline_cv(
    court$reference, court$source, match_on,
    contrast_on = c("participant", "family"),
    split_on = c("participant", "family", "replicate"),
    inner_split_on = match_on,
    n_folds = 2, seed = 17, spec = spec,
    fit_source_filter = function(tab) tab$condition == "signal"
  )

  expect_equal(sort(unique(fit$results$condition)), c("null", "signal"))
  expect_true(all(vapply(fit$folds, function(fold) {
    all(court$source$condition[fold$train_rows] == "signal")
  }, logical(1))))
  expect_true(all(vapply(fit$folds, function(fold) {
    fold$overlap_match_n == 0L &&
      all(vapply(fold$inner_folds, `[[`, integer(1), "overlap_match_n") == 0L)
  }, logical(1))))
  expect_true(all(fit$results$candidate_count == 4L))
})

test_that("resource profiling is bounded and non-instrumenting", {
  profiled <- v2_profile_fit("constant", sum(seq_len(10)))

  expect_equal(profiled$value, 55)
  expect_named(
    profiled$resources,
    c("method", "elapsed_seconds", "gc_used_mb", "gc_high_water_mb", "fit_size_mb")
  )
  expect_true(all(profiled$resources[-1] >= 0))
})

test_that("the persisted final court passes its frozen gates and manifest", {
  result_dir <- system.file(
    "validation", "gaze-weave-v2-results", package = "eyesim"
  )
  if (!nzchar(result_dir)) {
    result_dir <- testthat::test_path(
      "..", "..", "inst", "validation", "gaze-weave-v2-results"
    )
  }
  skip_if_not(dir.exists(result_dir), "final synthetic court has not been run")
  configuration <- utils::read.csv(
    file.path(result_dir, "configuration.csv"), stringsAsFactors = FALSE
  )
  gates <- utils::read.csv(
    file.path(result_dir, "gate-verdict.csv"), stringsAsFactors = FALSE
  )
  manifest <- utils::read.csv(
    file.path(result_dir, "manifest-md5.csv"), stringsAsFactors = FALSE
  )
  perturbations <- utils::read.csv(
    file.path(result_dir, "perturbations.csv"), stringsAsFactors = FALSE
  )

  expect_identical(configuration$final_seed, 20260816L)
  expect_false(configuration$smoke)
  expect_true(configuration$advance)
  expect_identical(configuration$provisional_default, "replay")
  expect_true(all(gates$passed))
  expect_true(all(names(v2_fractional_design()) %in% names(perturbations)))

  paths <- file.path(result_dir, manifest$file)
  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), manifest$md5)
})
