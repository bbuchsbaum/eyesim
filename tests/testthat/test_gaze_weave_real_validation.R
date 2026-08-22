real_result_path <- function(name) {
  installed <- system.file(
    "validation", "gaze-weave-real-results", name, package = "eyesim"
  )
  if (nzchar(installed)) return(installed)
  testthat::test_path(
    "..", "..", "inst", "validation", "gaze-weave-real-results", name
  )
}

test_that("real-control artifacts are complete and checksum-valid", {
  manifest <- utils::read.csv(real_result_path("manifest-md5.csv"))
  paths <- vapply(manifest$file, real_result_path, character(1))
  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), manifest$md5)
  expect_true(all(c(
    "configuration.csv", "gate-verdict.csv", "method-summary.csv",
    "bootstrap-intervals.csv", "trial-scores.csv", "scientific-results.rds"
  ) %in% manifest$file))
})

test_that("real court preserves the frozen cohort and two-way split", {
  configuration <- utils::read.csv(real_result_path("configuration.csv"))
  audit <- utils::read.csv(real_result_path("fold-audit.csv"))

  expect_identical(configuration$seed, 20260817L)
  expect_identical(
    configuration$participants,
    "104;121;124;128;315;317;319;327"
  )
  expect_identical(configuration$items, "5;18;69;73;81;82;88;105")
  expect_true(all(audit$participant_overlap == 0L))
  expect_true(all(audit$item_overlap == 0L))
  expect_true(all(audit$training_conditions == "signal"))
  expect_equal(nrow(audit), 8L)
})

test_that("real court compares every declared calibrated method fairly", {
  scores <- utils::read.csv(real_result_path("trial-scores.csv"))
  calibrated <- c(
    "transport_v2", "replay", "multimatch_ridge_registered",
    "density_ridge_registered", "elastic_ridge_registered"
  )
  expected <- expand.grid(
    task = c("blank_screen_imagery", "repeated_viewing"),
    condition = c("generic_gaze", "shuffled_order", "signal", "wrong_item"),
    method = calibrated,
    stringsAsFactors = FALSE
  )
  counts <- aggregate(
    participant ~ task + condition + method,
    scores[scores$method %in% calibrated, ], length
  )
  names(counts)[[4L]] <- "n"

  expect_equal(nrow(counts), nrow(expected))
  expect_true(all(counts$n == 64L))
  expect_true(all(scores$status == "scored"))
  expect_true(all(scores$converged))
  expect_true(all(scores$candidate_count == 4L))
})

test_that("real court locks no default after the repeated-viewing gate failure", {
  configuration <- utils::read.csv(real_result_path("configuration.csv"))
  gates <- utils::read.csv(real_result_path("gate-verdict.csv"))
  comparison <- utils::read.csv(real_result_path("paired-comparison.csv"))

  expect_false(configuration$advance)
  expect_identical(configuration$provisional_default, "no_real_data_default")
  expect_identical(
    configuration$strongest_baseline, "density_ridge_registered"
  )
  expect_false(configuration$superiority_supported)
  expect_false(gates$passed[gates$gate == "repeated_viewing"])
  expect_true(all(gates$passed[gates$gate != "repeated_viewing"]))
  expect_true(comparison$lower_95[comparison$metric == "delta_info"] < 0)
  expect_true(comparison$upper_95[comparison$metric == "delta_info"] > 0)
  expect_true(comparison$lower_95[comparison$metric == "delta_log_loss"] < 0)
  expect_true(comparison$upper_95[comparison$metric == "delta_log_loss"] > 0)
})

test_that("real controls retain density order invariance and null-safe registration", {
  order <- utils::read.csv(real_result_path("order-controls.csv"))
  density <- order[
    order$task == "repeated_viewing" &
      order$method == "density_registered",
  ]
  registration <- utils::read.csv(
    real_result_path("density-null-registration.csv")
  )

  expect_lte(density$max_absolute_pair_difference, 1e-10)
  expect_lte(
    registration$top1_credit[registration$method == "density_registered"] -
      registration$top1_credit[registration$method == "density_raw"],
    0.05
  )
})

test_that("real results preserve calibration failure, weak imagery, and runtime contrast", {
  summary <- utils::read.csv(real_result_path("method-summary.csv"))
  bootstrap <- utils::read.csv(real_result_path("bootstrap-intervals.csv"))
  resources <- utils::read.csv(real_result_path("resources.csv"))
  imagery <- summary[
    summary$task == "blank_screen_imagery" & summary$condition == "signal" &
      summary$method %in% c("replay", "transport_v2"),
  ]
  imagery_ci <- bootstrap[
    bootstrap$task == "blank_screen_imagery" &
      bootstrap$condition == "signal" &
      bootstrap$metric == "gaze_info_bits" &
      bootstrap$method %in% c("replay", "transport_v2"),
  ]

  expect_gt(imagery$mean_info_bits[imagery$method == "transport_v2"], 0)
  expect_lt(imagery$mean_info_bits[imagery$method == "replay"], 0)
  expect_true(all(imagery_ci$lower_95 < 0 & imagery_ci$upper_95 > 0))
  repeated_replay <- summary[
    summary$task == "repeated_viewing" & summary$condition == "signal" &
      summary$method == "replay",
  ]
  expect_lt(repeated_replay$mean_info_bits, 0)
  expect_gt(repeated_replay$top1_credit, 0.25)
  elapsed <- aggregate(elapsed_seconds ~ method, resources, sum)
  expect_gt(
    elapsed$elapsed_seconds[elapsed$method == "transport_v2"],
    100 * elapsed$elapsed_seconds[elapsed$method == "replay"]
  )
})
