source(
  gaze_weave_test_inst_path(
    "validation", "gaze-weave-pcmri-multi-replay.R"
  ),
  local = TRUE
)

multi_replay_local_data <- function() {
  all(file.exists(testthat::test_path(
    "..", "..", "test_data", "wynn_probe_delay",
    c("study_fix_input_new.csv", "testdelay_fix_input_matched.csv")
  )))
}

multi_replay_cached_design <- local({
  value <- NULL
  function() {
    if (!is.null(value)) return(value)
    if (!multi_replay_local_data()) return(NULL)
    old <- setwd(testthat::test_path("..", ".."))
    on.exit(setwd(old), add = TRUE)
    raw <- full_recognition_read_inputs(verify = TRUE)
    cohort <- multi_replay_select_cohort(raw)
    candidates <- multi_replay_candidate_plan(cohort)
    folds <- full_recognition_fold_plan(cohort, multi_replay_seed)
    value <<- list(
      raw = raw,
      cohort = cohort,
      candidates = candidates,
      folds = folds,
      design = multi_replay_validate_design(
        cohort, candidates, folds, strict = TRUE
      ),
      tables = multi_replay_task_tables(raw, cohort)
    )
    value
  }
})

test_that("multi-presentation cohort matches the frozen private support", {
  skip_if_not(multi_replay_local_data(), "Local restricted inputs unavailable")
  design <- multi_replay_cached_design()

  expect_identical(
    as.integer(design$design$observed),
    as.integer(multi_replay_expected)
  )
  expect_identical(nrow(design$cohort$pairs), 2055L)
  expect_identical(length(design$cohort$participants), 45L)
  expect_false("1023" %in% design$cohort$participants)
  expect_identical(
    as.integer(table(design$cohort$pairs$probe_type)), c(1044L, 1011L)
  )
})

test_that("candidate pools are exhaustive, deterministic, and exact", {
  skip_if_not(multi_replay_local_data(), "Local restricted inputs unavailable")
  design <- multi_replay_cached_design()
  observed <- multi_replay_candidate_plan(design$cohort)

  expect_identical(observed, design$candidates)
  expect_true(all(tapply(
    observed$is_true, observed$candidate_set_id, sum
  ) == 1L))
  pool_sizes <- vapply(split(
    observed$candidate_item,
    interaction(observed$participant, observed$item_fold, drop = TRUE)
  ), function(value) length(unique(value)), integer(1))
  expect_equal(range(pool_sizes), c(5L, 33L))
  expect_equal(stats::median(pool_sizes), 25)
  expect_equal(mean(pool_sizes), 22.83, tolerance = 0.01)
})

test_that("crossed folds preserve four paths without expanding by target", {
  skip_if_not(multi_replay_local_data(), "Local restricted inputs unavailable")
  design <- multi_replay_cached_design()
  audit <- design$design$fold_audit
  subset <- multi_replay_fold_subsets(
    design$tables, design$candidates, design$folds$folds[[1L]]
  )

  expect_identical(audit$eval_trials, c(510L, 513L, 522L, 510L))
  expect_identical(audit$train_trials, c(510L, 522L, 513L, 510L))
  expect_true(all(audit$participant_overlap == 0L))
  expect_true(all(audit$item_overlap == 0L))
  expect_identical(nrow(subset$ref_train), 4L * nrow(subset$source_train))
  expect_identical(nrow(subset$ref_eval), 4L * nrow(subset$source_eval))
  expect_true(all(table(paste(
    subset$ref_eval$participant, subset$ref_eval$item, sep = ":"
  )) == 4L))
  first <- subset$source_eval[1, , drop = FALSE]
  candidates <- multi_replay_candidate_references(subset$ref_eval, first)
  expect_true(all(table(candidates$item) == 4L))
  expect_true(first$item[[1L]] %in% candidates$item)
})

test_that("density episode averages four normalized maps", {
  paths <- lapply(1:4, function(presentation) {
    make_gaze_fixations(
      rbind(c(100, 100), c(180, 160), c(250, 220)) + presentation,
      duration = c(1, 2, 1), onset = c(0, 1, 3)
    )
  })
  reference <- tibble::tibble(
    participant = rep("p1", 4),
    item = rep(1L, 4),
    item_fold = rep(1L, 4),
    presentation = 1:4,
    fixgroup = paths
  )
  spec <- multi_replay_specs(smoke = TRUE)$density
  observed <- multi_replay_density_episode_table(reference, spec)
  component <- lapply(
    reference$fixgroup,
    eyesim:::baseline_density_signature,
    spec = spec
  )
  expected <- lapply(seq_along(spec$density_sigmas), function(i) {
    Reduce(`+`, lapply(component, `[[`, i)) / 4
  })

  expect_identical(nrow(observed), 1L)
  expect_equal(observed$signature[[1L]], expected, tolerance = 0)
  expect_equal(vapply(observed$signature[[1L]], sum, numeric(1)), c(1, 1))
})

test_that("study reversal changes order but preserves occupancy", {
  path <- make_gaze_fixations(
    rbind(c(10, 20), c(30, 40), c(50, 60), c(70, 80)),
    duration = c(1, 2, 3, 4), onset = c(0, 1, 3, 6)
  )
  shuffled <- multi_replay_reverse_path(path)

  expect_identical(shuffled$x, rev(path$x))
  expect_equal(sort(shuffled$x), sort(path$x), tolerance = 0)
  expect_equal(sort(shuffled$y), sort(path$y), tolerance = 0)
  expect_equal(sort(shuffled$duration), sort(path$duration), tolerance = 0)
  expect_equal(sum(shuffled$duration), sum(path$duration), tolerance = 0)
  expect_equal(shuffled$onset, c(0, head(cumsum(shuffled$duration), -1L)))
  expect_identical(multi_replay_reverse_path(path), shuffled)
})

test_that("local smoke checkpoints retain the frozen score contract", {
  result_dir <- file.path(tempdir(), "eyesim-multi-replay-smoke")
  path <- multi_replay_checkpoint(result_dir, "all4", 1L, smoke = TRUE)
  skip_if_not(file.exists(path), "Smoke checkpoint not present in this process")
  checkpoint <- readRDS(path)

  expect_true(all(checkpoint$scored$method %in% c(
    "replay_all4_unshrunk_exhaustive",
    "replay_all4_shrunk_exhaustive"
  )))
  expect_true(all(is.finite(checkpoint$scored$gaze_info_bits)))
  expect_true(all(checkpoint$scored$candidate_count >= 2L))
  expect_true(all(checkpoint$calibration$shrinkage_not_worse))
})
