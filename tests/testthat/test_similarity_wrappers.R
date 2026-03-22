make_similarity_wrapper_fixgroup <- function(x_offset = 0, y_offset = 0) {
  fixation_group(
    x = c(100, 180, 260, 340) + x_offset,
    y = c(120, 200, 180, 260) + y_offset,
    onset = c(0, 120, 260, 420),
    duration = c(120, 140, 160, 180)
  )
}

make_similarity_wrapper_tables <- function() {
  fixgroups <- list(
    make_similarity_wrapper_fixgroup(0, 0),
    make_similarity_wrapper_fixgroup(80, 40),
    make_similarity_wrapper_fixgroup(-60, 90)
  )

  tibble::tibble(
    trial = c("t1", "t2", "t3"),
    fixgroup = fixgroups,
    scanpath = lapply(fixgroups, scanpath)
  )
}

test_that("fixation_similarity overlap returns numeric similarity values", {
  ref_tab <- make_similarity_wrapper_tables()
  source_tab <- make_similarity_wrapper_tables()

  res <- fixation_similarity(
    ref_tab,
    source_tab,
    match_on = "trial",
    method = "overlap",
    time_samples = seq(0, 420, by = 60)
  )

  expect_true(is.numeric(res$eye_sim))
  expect_false(anyNA(res$eye_sim))
  expect_equal(res$eye_sim, rep(1, nrow(res)))
})

test_that("fixation_similarity sinkhorn accepts fixation_group inputs", {
  if (!requireNamespace("T4transport", quietly = TRUE)) {
    skip("T4transport not available")
  }

  ref_tab <- make_similarity_wrapper_tables()
  source_tab <- make_similarity_wrapper_tables()

  res <- fixation_similarity(
    ref_tab,
    source_tab,
    match_on = "trial",
    method = "sinkhorn"
  )

  expect_true(is.numeric(res$eye_sim))
  expect_false(anyNA(res$eye_sim))
  expect_true(all(res$eye_sim > 0))
  expect_equal(res$eye_sim, rep(res$eye_sim[[1]], nrow(res)))
})

test_that("scanpath_similarity expands multimatch metrics without permutations", {
  ref_tab <- make_similarity_wrapper_tables()
  source_tab <- make_similarity_wrapper_tables()

  res <- scanpath_similarity(
    ref_tab,
    source_tab,
    match_on = "trial",
    method = "multimatch",
    screensize = c(800, 600)
  )

  expected_cols <- c(
    "mm_vector",
    "mm_direction",
    "mm_length",
    "mm_position",
    "mm_duration",
    "mm_position_emd"
  )

  expect_true(all(expected_cols %in% names(res)))
  expect_false("eye_sim" %in% names(res))
  expect_false(any(vapply(res[expected_cols], is.list, logical(1))))
  expect_true(all(vapply(res[expected_cols], is.numeric, logical(1))))
  expect_true(all(res$mm_vector > 0.999))
})

test_that("scanpath_similarity expands multimatch permutation summaries", {
  ref_tab <- make_similarity_wrapper_tables()
  source_tab <- make_similarity_wrapper_tables()

  res <- scanpath_similarity(
    ref_tab,
    source_tab,
    match_on = "trial",
    method = "multimatch",
    permutations = 2,
    screensize = c(800, 600)
  )

  expected_cols <- c(
    "mm_vector",
    "mm_direction",
    "mm_length",
    "mm_position",
    "mm_duration",
    "mm_position_emd"
  )
  perm_cols <- paste0(expected_cols, "_perm")
  diff_cols <- paste0(expected_cols, "_diff")

  expect_true(all(c(expected_cols, perm_cols, diff_cols) %in% names(res)))
  expect_false(any(vapply(res[c(expected_cols, perm_cols, diff_cols)], is.list, logical(1))))
  expect_true(all(vapply(res[c(expected_cols, perm_cols, diff_cols)], is.numeric, logical(1))))
  expect_true(all(res$mm_vector > 0.999))
})

test_that("multi_match returns all six metrics for short scanpaths", {
  fg_short <- fixation_group(
    x = c(100, 180),
    y = c(120, 200),
    onset = c(0, 120),
    duration = c(120, 140)
  )
  fg_long <- make_similarity_wrapper_fixgroup()

  res <- expect_warning(
    multi_match(scanpath(fg_short), scanpath(fg_long), screensize = c(800, 600)),
    "requires 3 or more coordinates"
  )

  expected_cols <- c(
    "mm_vector",
    "mm_direction",
    "mm_length",
    "mm_position",
    "mm_duration",
    "mm_position_emd"
  )

  expect_equal(names(res), expected_cols)
  expect_length(res, length(expected_cols))
  expect_true(all(is.na(res)))
})

test_that("scanpath_similarity keeps six multimatch columns for short-window permutations", {
  ref_tab <- make_similarity_wrapper_tables()
  source_tab <- make_similarity_wrapper_tables()

  res <- expect_warning(
    scanpath_similarity(
      ref_tab,
      source_tab,
      match_on = "trial",
      method = "multimatch",
      permutations = 10,
      window = c(0, 250),
      screensize = c(800, 600)
    ),
    "requires 3 or more coordinates"
  )

  expected_cols <- c(
    "mm_vector",
    "mm_direction",
    "mm_length",
    "mm_position",
    "mm_duration",
    "mm_position_emd"
  )
  perm_cols <- paste0(expected_cols, "_perm")
  diff_cols <- paste0(expected_cols, "_diff")

  expect_true(all(c(expected_cols, perm_cols, diff_cols) %in% names(res)))
  expect_true(all(vapply(res[c(expected_cols, perm_cols, diff_cols)], is.numeric, logical(1))))
  expect_true(all(is.na(unlist(res[c(expected_cols, perm_cols, diff_cols)]))))
})
