replication_result_path <- function(name) {
  installed <- system.file(
    "validation", "gaze-weave-replication-results", name,
    package = "eyesim"
  )
  if (nzchar(installed)) return(installed)
  testthat::test_path(
    "..", "..", "inst", "validation",
    "gaze-weave-replication-results", name
  )
}

test_that("replication artifacts are complete and checksum-valid", {
  manifest <- utils::read.csv(replication_result_path("manifest-md5.csv"))
  paths <- vapply(manifest$file, replication_result_path, character(1))

  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), manifest$md5)
  expect_true(all(c(
    "configuration.csv", "gate-verdict.csv", "method-summary.csv",
    "bootstrap-intervals.csv", "trial-scores.csv", "scientific-results.rds"
  ) %in% manifest$file))
})

test_that("replication cohort is participant- and item-disjoint", {
  configuration <- utils::read.csv(
    replication_result_path("configuration.csv")
  )
  audit <- utils::read.csv(replication_result_path("fold-audit.csv"))

  expect_identical(configuration$protocol_version, 2L)
  expect_identical(configuration$seed, 20260818L)
  expect_identical(
    configuration$participants, "111;130;18;300;302;316;325;9"
  )
  expect_identical(configuration$items, "4;11;24;62;63;80;92;101")
  expect_length(
    intersect(
      strsplit(configuration$participants, ";", fixed = TRUE)[[1]],
      strsplit(configuration$excluded_participants, ";", fixed = TRUE)[[1]]
    ),
    0L
  )
  expect_length(
    intersect(
      strsplit(configuration$items, ";", fixed = TRUE)[[1]],
      strsplit(configuration$excluded_items, ";", fixed = TRUE)[[1]]
    ),
    0L
  )
  expect_identical(
    configuration$replay_score_scale,
    "mean_log_likelihood_per_duration_bin"
  )
  expect_true(all(audit$participant_overlap == 0L))
  expect_true(all(audit$item_overlap == 0L))
  expect_true(all(audit$training_conditions == "signal"))
  expect_equal(nrow(audit), 8L)
})

test_that("replication preserves the full fair comparator court", {
  scores <- utils::read.csv(replication_result_path("trial-scores.csv"))
  calibrated <- c(
    "transport_v2", "replay", "multimatch_ridge_registered",
    "density_ridge_registered", "elastic_ridge_registered"
  )
  counts <- aggregate(
    participant ~ task + condition + method,
    scores[scores$method %in% calibrated, ],
    length
  )

  expect_equal(nrow(counts), 2L * 4L * length(calibrated))
  expect_true(all(counts$participant == 64L))
  expect_true(all(scores$status == "scored"))
  expect_true(all(scores$converged))
  expect_true(all(scores$candidate_count == 4L))
})

test_that("resolution correction transfers but imagery gate remains negative", {
  configuration <- utils::read.csv(
    replication_result_path("configuration.csv")
  )
  gates <- utils::read.csv(replication_result_path("gate-verdict.csv"))
  summary <- utils::read.csv(replication_result_path("method-summary.csv"))
  bootstrap <- utils::read.csv(
    replication_result_path("bootstrap-intervals.csv")
  )

  replay_repeated <- summary[
    summary$task == "repeated_viewing" & summary$condition == "signal" &
      summary$method == "replay",
  ]
  imagery <- summary[
    summary$task == "blank_screen_imagery" &
      summary$condition == "signal" &
      summary$method %in% c("replay", "transport_v2"),
  ]
  imagery_ci <- bootstrap[
    bootstrap$task == "blank_screen_imagery" &
      bootstrap$condition == "signal" &
      bootstrap$metric == "gaze_info_bits" &
      bootstrap$method %in% c("replay", "transport_v2"),
  ]

  expect_gt(replay_repeated$mean_info_bits, 0)
  expect_gt(replay_repeated$top1_credit, 0.25)
  expect_true(all(imagery$mean_info_bits < 0))
  expect_true(all(imagery_ci$lower_95 < 0 & imagery_ci$upper_95 > 0))
  expect_true(gates$passed[gates$gate == "repeated_viewing"])
  expect_false(gates$passed[gates$gate == "blank_screen_imagery"])
  expect_true(all(gates$passed[gates$gate != "blank_screen_imagery"]))
  expect_false(configuration$advance)
  expect_identical(configuration$provisional_default, "no_real_data_default")
  expect_false(configuration$superiority_supported)
})

test_that("replication retains null safety, density order invariance, and runtime contrast", {
  operating <- utils::read.csv(
    replication_result_path("operating-characteristics.csv")
  )
  order <- utils::read.csv(replication_result_path("order-controls.csv"))
  registration <- utils::read.csv(
    replication_result_path("density-null-registration.csv")
  )
  resources <- utils::read.csv(replication_result_path("resources.csv"))
  density <- order[
    order$task == "repeated_viewing" &
      order$method == "density_registered",
  ]
  elapsed <- aggregate(elapsed_seconds ~ method, resources, sum)

  expect_true(all(operating$null_error <= 0.075))
  expect_lte(density$max_absolute_pair_difference, 1e-10)
  expect_lte(
    registration$top1_credit[registration$method == "density_registered"] -
      registration$top1_credit[registration$method == "density_raw"],
    0.05
  )
  expect_gt(
    elapsed$elapsed_seconds[elapsed$method == "transport_v2"],
    100 * elapsed$elapsed_seconds[elapsed$method == "replay"]
  )
})
