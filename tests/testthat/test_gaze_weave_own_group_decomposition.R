own_group_source <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    gaze_weave_test_inst_path(
      "validation", "gaze-weave-own-group-decomposition.R"
    ),
    envir = environment
  )
  environment
}

test_that("candidate margin has an analytic truth-minus-nonmatch value", {
  court <- own_group_source()
  candidates <- data.frame(
    log_score = c(2, 1, 0, -1, -2),
    is_true = c(TRUE, FALSE, FALSE, FALSE, FALSE)
  )

  expect_equal(court$own_group_candidate_margin(candidates), 2.5)
  expect_error(
    court$own_group_candidate_margin(transform(candidates, is_true = TRUE)),
    "one finite truth"
  )
})

test_that("legacy frozen specification is rewrapped without changing values", {
  court <- own_group_source()
  legacy <- structure(
    list(version = 3L, estimand = "edge_normalized_episode_transport", x = 2),
    class = c("gaze_transport_v3_spec", "list")
  )
  canonical <- court$own_group_spec_from_frozen(legacy)

  expect_s3_class(canonical, "gaze_transport_spec")
  expect_identical(unclass(canonical), unclass(legacy))
  expect_error(
    court$own_group_spec_from_frozen(structure(
      list(version = 2L, estimand = "old"),
      class = c("gaze_transport_v3_spec", "list")
    )),
    "not the canonical"
  )
})

test_that("donor selection is deterministic, row invariant, and excludes self", {
  court <- own_group_source()
  study <- expand.grid(
    participant = paste0("p", 1:7), presentation = 1:4,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  study$item <- 2L
  study$study_image_version <- "2A"
  study$nfix <- 4L
  study$fixgroup <- lapply(seq_len(nrow(study)), function(index) index)
  first <- court$own_group_select_donors(study, "p1", 2L, "2A", 41L)
  second <- court$own_group_select_donors(
    study[sample(seq_len(nrow(study))), ], "p1", 2L, "2A", 41L
  )

  expect_identical(first$donor_participant, second$donor_participant)
  expect_false("p1" %in% first$donor_participant)
  expect_equal(length(unique(first$donor_participant)), 4L)
  expect_identical(first$episode_id, 1:4)
})

test_that("sparse donor cells retain four presentations without using self", {
  court <- own_group_source()
  study <- expand.grid(
    participant = c("p1", "p2", "p3"), presentation = 1:4,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  study$item <- 2L
  study$study_image_version <- "2A"
  study$nfix <- 4L
  study$fixgroup <- lapply(seq_len(nrow(study)), function(index) index)
  selected <- court$own_group_select_donors(
    study, "p1", 2L, "2A", 41L
  )

  expect_identical(selected$episode_id, 1:4)
  expect_equal(length(unique(selected$donor_participant)), 2L)
  expect_false("p1" %in% selected$donor_participant)
})

test_that("newtest candidate plans are deterministic and item-fold local", {
  court <- own_group_source()
  item_map <- data.frame(item = 1:12, item_fold = rep(1:2, each = 6))
  source <- data.frame(participant = "p1", item = 4L)
  first <- court$own_group_newtest_plan(source, item_map)
  second <- court$own_group_newtest_plan(source, item_map[sample(12), ])

  expect_identical(first, second)
  expect_equal(nrow(first), 5L)
  expect_equal(sum(first$is_true), 1L)
  expect_true(all(first$candidate_item %in% 1:6))
})

test_that("frozen warp reconstruction reproduces canonical scoring", {
  retrieval_dir <- gaze_weave_test_inst_path(
    "validation", "gaze-weave-transport-v3-retrieval-results"
  )
  testthat::skip_if_not(
    file.exists(file.path(
      retrieval_dir, "checkpoint-transport-v3-fold-01.rds"
    )),
    "Local canonical Transport checkpoint is unavailable."
  )
  court <- own_group_source()
  old <- setwd(testthat::test_path("..", ".."))
  on.exit(setwd(old), add = TRUE)
  checkpoint <- readRDS(file.path(
    retrieval_dir, "checkpoint-transport-v3-fold-01.rds"
  ))
  warp <- court$own_group_warp_from_info(checkpoint$warp)
  group <- warp$group_models[[1L]]
  frozen <- checkpoint$warp$groups[[1L]]

  expect_equal(group$scale, frozen$scale, tolerance = 1e-12)
  expect_equal(group$translation, frozen$translation, tolerance = 1e-12)
  expect_equal(group$A, diag(frozen$scale, 2L), tolerance = 1e-12)

  context <- court$own_group_context()
  source <- context$tables$source
  source <- source[
    as.character(source$participant) ==
      as.character(checkpoint$scored$participant[[1L]]) &
      source$item == checkpoint$scored$item[[1L]],
    , drop = FALSE
  ]
  expanded <- court$v3_retrieval_expand_pool(
    context$tables$reference, source, context$candidate_plan
  )
  reproduced <- court$v3_retrieval_score_row(
    expanded$source, expanded$reference, context$spec, warp,
    checkpoint$calibration
  )
  expect_equal(
    reproduced$gaze_info_bits, checkpoint$scored$gaze_info_bits[[1L]],
    tolerance = 2e-5
  )
  expect_equal(
    reproduced$candidates[[1L]]$log_score,
    checkpoint$scored$candidates[[1L]]$log_score,
    tolerance = 2e-5
  )
})

test_that("local inputs support four independent donors and newtest", {
  data_dir <- gaze_weave_test_inst_path("validation")
  testthat::skip_if_not(
    file.exists(file.path(
      data_dir, "gaze-weave-transport-v3-retrieval-results",
      "freeze-private", "score-blind-tables.rds"
    )),
    "Local own-group inputs are unavailable."
  )
  court <- own_group_source()
  old <- setwd(testthat::test_path("..", ".."))
  on.exit(setwd(old), add = TRUE)
  context <- court$own_group_context()
  folds <- court$own_group_fold_plan(context)
  newtest <- dplyr::bind_rows(lapply(folds, function(fold) {
    court$own_group_newtest_source(context, fold)
  }))

  expect_gt(nrow(newtest), 800L)
  expect_true(all(newtest$probe_type == "newtest"))
  expect_equal(anyDuplicated(court$own_group_trial_key(
    newtest$participant, newtest$item
  )), 0L)
  sample_row <- context$tables$source[1, , drop = FALSE]
  plan <- court$own_group_existing_plan(context$candidate_plan, sample_row)
  lookup <- court$own_group_version_lookup(
    context$study, sample_row$participant[[1L]]
  )
  pool <- court$own_group_build_pool(
    context$study, sample_row, plan, lookup
  )
  expect_equal(nrow(pool), 20L)
  expect_false(as.character(sample_row$participant) %in% pool$donor_participant)
})
