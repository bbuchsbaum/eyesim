probe_delay_validation_path <- function(...) {
  gaze_weave_test_inst_path("validation", ...)
}

probe_delay_local_data_path <- function(...) {
  testthat::test_path("..", "..", "test_data", "wynn_probe_delay", ...)
}

probe_delay_local_result_path <- function(...) {
  testthat::test_path(
    "..", "..", "inst", "validation",
    "gaze-weave-probe-delay-results", ...
  )
}

probe_delay_source_runner <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    probe_delay_validation_path("gaze-weave-probe-delay.R"),
    envir = environment
  )
  environment
}

probe_delay_skip_without_restricted_inputs <- function() {
  paths <- probe_delay_local_data_path(c(
    "study_fix_input_new.csv", "testdelay_fix_input_matched.csv"
  ))
  testthat::skip_if_not(
    all(file.exists(paths)),
    "Restricted Wynn fixation exports are not present in this checkout."
  )
}

probe_delay_skip_without_local_results <- function() {
  testthat::skip_if_not(
    file.exists(probe_delay_local_result_path("manifest-md5.csv")),
    "Local-only probe-delay court results are not present."
  )
}

test_that("probe-delay windows split boundary fixations deterministically", {
  court <- probe_delay_source_runner()
  rows <- data.frame(
    FixX = c(512, 512), FixY = c(384, 384),
    FixEndTime = c(700, 700), stringsAsFactors = FALSE
  )
  probe <- court$probe_delay_clip_rows(
    rows, start = c(400, 490), duration = c(300, 300), window = c(0, 500)
  )
  delay <- court$probe_delay_clip_rows(
    rows, start = c(400, 490), duration = c(300, 300), window = c(500, 3000)
  )

  expect_equal(probe$path_duration, 100)
  expect_equal(probe$path_onset, 400)
  expect_equal(delay$path_duration, c(200, 290))
  expect_equal(delay$path_onset, c(0, 0))
})

test_that("signal-only diagnostics do not manufacture order-control rows", {
  court <- probe_delay_source_runner()
  scored <- data.frame(
    task = "late_delay", method = "replay", condition = "signal",
    participant = "p1", item = 1L, gaze_info_bits = 0,
    compatibility_bits = 0, calibrated = TRUE,
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(court$real_order_checks(scored)), 0L)
})

test_that("restricted fixation exports are ignored and checksum-pinned", {
  probe_delay_skip_without_restricted_inputs()
  root <- testthat::test_path("..", "..")
  ignore <- readLines(file.path(root, ".gitignore"), warn = FALSE)
  expected_ignored <- c(
    "test_data/wynn_probe_delay/study_fix_input_new.csv",
    "test_data/wynn_probe_delay/testdelay_fix_input_matched.csv"
  )
  expect_true(all(expected_ignored %in% ignore))

  court <- probe_delay_source_runner()
  manifest <- court$probe_delay_verify_inputs(
    probe_delay_local_data_path()
  )
  expect_identical(manifest$rows, c(101991L, 29699L))
  expect_identical(manifest$columns, c(16L, 23L))
})

test_that("probe-delay cohort and folds match the prospective contract", {
  probe_delay_skip_without_restricted_inputs()
  court <- probe_delay_source_runner()
  raw <- court$probe_delay_read_inputs(probe_delay_local_data_path())
  cohort <- court$probe_delay_select_cohort(raw)
  observed_items <- lapply(
    split(cohort$pairs$item, cohort$pairs$participant),
    function(value) as.numeric(sort(value))
  )
  observed_fold1 <- lapply(split(
    cohort$pairs$item[cohort$pairs$item_fold == 1L],
    cohort$pairs$participant[cohort$pairs$item_fold == 1L]
  ), function(value) as.numeric(sort(value)))

  expect_identical(cohort$participants, names(court$probe_delay_expected_items))
  expect_identical(observed_items, court$probe_delay_expected_items)
  expect_identical(observed_fold1, court$probe_delay_expected_fold1_items)
  expect_equal(nrow(cohort$pairs), 80L)
  expect_identical(as.integer(cohort$balance$probe_type), c(32L, 48L))
  expect_identical(
    as.integer(cohort$balance$degradation), c(13L, 21L, 14L, 12L, 20L)
  )

  plan <- court$probe_delay_fold_plan(cohort)
  expect_identical(
    as.integer(table(plan$participant_map$participant_fold)), c(4L, 4L)
  )
  expect_true(all(table(
    plan$participant_map$design_signature,
    plan$participant_map$participant_fold
  ) == 1L))
})

test_that("every frozen phase has five candidates per evaluation participant", {
  probe_delay_skip_without_restricted_inputs()
  court <- probe_delay_source_runner()
  raw <- court$probe_delay_read_inputs(probe_delay_local_data_path())
  cohort <- court$probe_delay_select_cohort(raw)
  plan <- court$probe_delay_fold_plan(cohort)
  tasks <- c(
    "repeated_viewing", "probe", "delay", "combined",
    "early_delay", "late_delay"
  )

  for (task in tasks) {
    tables <- court$probe_delay_task_tables(raw, cohort, task)
    expect_equal(nrow(tables$reference), 80L, info = task)
    expect_equal(nrow(tables$signal), 80L, info = task)
    expect_equal(nrow(tables$source), 320L, info = task)
    for (fold in plan$folds) {
      subsets <- court$real_fold_subsets(tables, fold)
      expect_equal(nrow(subsets$source_train), 20L, info = task)
      expect_equal(nrow(subsets$ref_eval), 20L, info = task)
      candidate_n <- table(subsets$source_eval$participant[
        subsets$source_eval$condition == "signal"
      ])
      expect_true(all(candidate_n == 5L), info = task)
      expect_equal(length(intersect(
        unique(subsets$source_train$participant),
        unique(subsets$source_eval$participant)
      )), 0L, info = task)
      expect_equal(length(intersect(
        unique(subsets$source_train$item),
        unique(subsets$source_eval$item)
      )), 0L, info = task)
    }
  }
})

test_that("local-only probe-delay artifacts are complete and checksum-valid", {
  probe_delay_skip_without_local_results()
  manifest <- utils::read.csv(
    probe_delay_local_result_path("manifest-md5.csv"),
    stringsAsFactors = FALSE
  )
  paths <- vapply(
    manifest$file, probe_delay_local_result_path, character(1)
  )

  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), manifest$md5)
  expect_true(all(c(
    "configuration.csv", "persistence-verdict.csv", "path-coverage.csv",
    "method-summary.csv", "bootstrap-intervals.csv", "trial-scores.csv",
    "scientific-results.rds"
  ) %in% manifest$file))
})

test_that("full local court preserves the frozen no-persistence verdict", {
  probe_delay_skip_without_local_results()
  configuration <- utils::read.csv(
    probe_delay_local_result_path("configuration.csv"),
    stringsAsFactors = FALSE
  )
  verdict <- utils::read.csv(
    probe_delay_local_result_path("persistence-verdict.csv"),
    stringsAsFactors = FALSE
  )
  audit <- utils::read.csv(probe_delay_local_result_path("fold-audit.csv"))
  scores <- utils::read.csv(
    probe_delay_local_result_path("trial-scores.csv"),
    stringsAsFactors = FALSE
  )
  coverage <- utils::read.csv(probe_delay_local_result_path("path-coverage.csv"))

  expect_identical(configuration$protocol_version, 2L)
  expect_identical(configuration$seed, 20260820L)
  expect_false(configuration$smoke)
  expect_identical(verdict$verdict, "not_supported")
  expect_true(verdict$leakage_safe)
  expect_true(verdict$engine_convergence)
  expect_false(verdict$independent_replication)
  expect_false(verdict$closes_wang_gate)
  expect_equal(nrow(audit), 24L)
  expect_true(all(audit$participant_overlap == 0L))
  expect_true(all(audit$item_overlap == 0L))
  expect_equal(nrow(scores), 9600L)
  expect_true(all(scores$status == "scored"))
  expect_true(all(scores$candidate_count == 5L))
  expect_true(all(scores$converged[
    scores$method %in% c("transport_v2", "replay")
  ]))
  expect_equal(
    coverage$multimatch_pair_coverage[coverage$task == "probe"], 0.075
  )
})
