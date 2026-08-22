source(
  gaze_weave_test_inst_path(
    "validation", "gaze-weave-transport-v3-repeated-viewing.R"
  ),
  local = TRUE
)

test_that("repeated-viewing protocol includes initial and adjacent contrasts", {
  expect_identical(
    v3_repeated_contrasts$task,
    c(
      "study_p1_p2", "study_p1_p3", "study_p1_p4",
      "study_p2_p3", "study_p3_p4"
    )
  )
  expect_identical(v3_repeated_seed, 20260823L)
  expect_identical(v3_repeated_shuffle_seed, 20260824L)
})

test_that("reversal and shuffle preserve spatial-duration density", {
  path <- fixation_group(
    x = c(1, 4, 9, 16), y = c(2, 3, 5, 8),
    duration = c(10, 20, 30, 40), onset = c(0, 10, 30, 60)
  )
  reversed <- v3_repeated_perturb_path(path, "reversed", "a")
  shuffled <- v3_repeated_perturb_path(path, "shuffled", "a")
  shuffled_again <- v3_repeated_perturb_path(path, "shuffled", "a")
  signature <- function(x) {
    sort(paste(x$x, x$y, x$duration, sep = ":"))
  }
  expect_identical(signature(reversed), signature(path))
  expect_identical(signature(shuffled), signature(path))
  expect_identical(shuffled, shuffled_again)
  expect_identical(reversed$x, rev(path$x))
})

test_that("candidate expansion retains the frozen K=5 truth contract", {
  reference <- tibble::tibble(
    participant = 1L, item = 1:5,
    fixgroup = lapply(1:5, function(x) {
      fixation_group(
        x = c(x, x + 1), y = c(x, x + 2),
        duration = c(1, 1), onset = c(0, 1)
      )
    })
  )
  source <- reference[1, ]
  source$task <- "study_p1_p4"
  plan <- data.frame(
    participant = 1L, target_item = 1L, item_fold = 1L,
    candidate_set_id = "1:1", candidate_position = 1:5,
    candidate_item = 1:5, is_true = c(TRUE, rep(FALSE, 4))
  )
  value <- v3_repeated_expand_pool(reference, source, plan)
  expect_identical(nrow(value$reference), 5L)
  expect_identical(nrow(value$source), 1L)
  expect_identical(sum(value$reference$is_true), 1L)
  expect_true(all(value$reference$prior_weight == 1))
})

test_that("candidate strata are separate from response-blind calibration folds", {
  source <- tibble::tibble(
    participant = rep(1:2, each = 2), item = 1:4,
    candidate_set_id = paste0("set_", 1:4),
    fixgroup = lapply(1:4, function(x) {
      fixation_group(
        x = c(x, x + 1), y = c(x, x + 2),
        duration = c(1, 1), onset = c(0, 1)
      )
    })
  )
  folds <- eyesim:::make_gaze_weave_folds(
    source, split_on = c("participant", "item"),
    contrast_on = NULL, n_folds = 2L, seed = 1L
  )
  expect_identical(folds$n_folds, 2L)
  expect_identical(sort(unique(folds$fold_id)), 1:2)
})

test_that("advance rule requires both positive estimates and one positive lower bound", {
  make_scores <- function(information = 1, control = 0) {
    expand <- expand.grid(
      participant = 1:4, item = 1:4,
      control = c("intact", "reversed", "shuffled"),
      stringsAsFactors = FALSE
    )
    expand$task <- "study_p1_p4"
    expand$gaze_info_bits <- ifelse(
      expand$control == "intact", information, control
    )
    expand
  }
  cohort <- list(pairs = expand.grid(participant = 1:4, item = 1:4))
  passing <- v3_repeated_advance(make_scores(), cohort, draws = 20L)
  tied <- v3_repeated_advance(make_scores(1, 1), cohort, draws = 20L)
  expect_true(attr(passing, "advance"))
  expect_false(attr(tied, "advance"))
})

test_that("restricted repeated-viewing output is ignored", {
  root <- normalizePath(testthat::test_path("..", ".."))
  old <- setwd(root)
  on.exit(setwd(old), add = TRUE)
  skip_if_not(
    dir.exists(".git"),
    "Git-ignore policy is checked only in a source checkout"
  )
  status <- system2(
    "git", c("check-ignore", "-q",
             paste0(
               "inst/validation/",
               "gaze-weave-transport-v3-repeated-viewing-results/.probe"
             )),
    stdout = FALSE, stderr = FALSE
  )
  expect_identical(status, 0L)
})
