source(
  gaze_weave_test_inst_path(
    "validation", "gaze-weave-transport-v3-retrieval.R"
  ),
  local = TRUE
)

v3_retrieval_local_data <- function() {
  all(file.exists(testthat::test_path(
    "..", "..", "test_data", "wynn_probe_delay",
    c("study_fix_input_new.csv", "testdelay_fix_input_matched.csv")
  )))
}

test_that("retrieval measurement tables are behavior-blind and four-episode", {
  skip_if_not(v3_retrieval_local_data(), "Local restricted inputs unavailable")
  old <- setwd(testthat::test_path("..", ".."))
  on.exit(setwd(old), add = TRUE)
  raw <- full_recognition_read_inputs(verify = TRUE)
  cohort <- full_recognition_select_cohort(raw)
  tables <- v3_retrieval_score_blind_tables(raw, cohort)
  expect_identical(nrow(tables$source), 1295L)
  expect_identical(nrow(tables$reference), 4L * 1295L)
  expect_false(any(v3_retrieval_forbidden_score_columns %in%
                     names(tables$source)))
  expect_false(any(v3_retrieval_forbidden_score_columns %in%
                     names(tables$reference)))
  expect_true(all(table(paste(
    tables$reference$participant, tables$reference$item, sep = ":"
  )) == 4L))
})

test_that("retrieval pool retains K=5 candidates and four episodes", {
  skip_if_not(v3_retrieval_local_data(), "Local restricted inputs unavailable")
  old <- setwd(testthat::test_path("..", ".."))
  on.exit(setwd(old), add = TRUE)
  raw <- full_recognition_read_inputs(verify = TRUE)
  cohort <- full_recognition_smoke_cohort(full_recognition_select_cohort(raw))
  tables <- v3_retrieval_score_blind_tables(raw, cohort)
  plan <- full_recognition_candidate_plan(cohort, v3_retrieval_seed)
  value <- v3_retrieval_expand_pool(tables$reference, tables$source, plan)
  expect_true(all(table(value$reference$candidate_set_id) == 20L))
  expect_true(all(table(paste(
    value$reference$candidate_set_id,
    value$reference$candidate_item, sep = "|"
  )) == 4L))
})

test_that("response coding and log loss preserve the frozen meanings", {
  expect_identical(
    v3_retrieval_response(c(0, 1, 2, 3, 4, NA)),
    c(NA_integer_, 1L, 1L, 0L, 0L, NA_integer_)
  )
  expect_equal(v3_retrieval_log_loss(c(1, 0), c(0.8, 0.2)),
               rep(-log(0.8), 2))
})

test_that("behavior fold keeps quality separate and evaluates held-out cells", {
  make <- function(participants, items) {
    tab <- expand.grid(
      participant = participants, item = items,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    tab$said_old <- rep(0:1, length.out = nrow(tab))
    tab$sal10 <- rep(c(-2, 0, 2), length.out = nrow(tab))
    tab$gaze_info_bits <- seq(-1, 1, length.out = nrow(tab)) +
      rep(c(0, 0.1), length.out = nrow(tab))
    tab$effective_fixations <- seq(3, 8, length.out = nrow(tab))
    tab$total_duration <- seq(1000, 2400, length.out = nrow(tab))
    tab$Response <- ifelse(tab$said_old == 1L, 1L, 4L)
    tab
  }
  train <- make(1:4, 1:4)
  evaluate <- make(5:6, 5:6)
  value <- v3_retrieval_fit_behavior_fold(
    train, evaluate, "transport_v3", 1L, "old"
  )
  expect_identical(nrow(value), nrow(evaluate))
  expect_true(all(is.finite(value$behavior_info_bits)))
  expect_length(intersect(train$participant, value$participant), 0L)
  expect_length(intersect(train$item, value$item), 0L)
})

test_that("retrieval output is explicitly Git-ignored", {
  old <- setwd(testthat::test_path("..", ".."))
  on.exit(setwd(old), add = TRUE)
  skip_if_not(
    dir.exists(".git"),
    "Git-ignore policy is checked only in a source checkout"
  )
  status <- system2(
    "git", c("check-ignore", "-q",
             "inst/validation/gaze-weave-transport-v3-retrieval-results/.probe"),
    stdout = FALSE, stderr = FALSE
  )
  expect_identical(status, 0L)
})
