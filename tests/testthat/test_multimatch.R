library(testthat)

make_scanpath <- function(x, y, onset, duration) {
  df <- data.frame(x = x, y = y, onset = onset, duration = duration)
  scanpath.fixation_group(df)
}

test_that("multi_match enforces minimum fixation count", {
  scan1 <- make_scanpath(c(0, 1), c(0, 1), c(0, 1), c(100, 100))
  expect_warning(
    res <- multi_match(scan1, scan1, c(100, 100)),
    "requires 3 or more coordinates"
  )
  expect_named(res, c("mm_vector", "mm_direction", "mm_length",
                      "mm_position", "mm_duration", "mm_position_emd"))
  expect_true(all(is.na(res)))
})

test_that("multi_match returns perfect similarity for identical scanpaths", {
  scan1 <- make_scanpath(c(0, 50, 100), c(0, 50, 100), c(0, 1, 2), c(100, 100, 100))
  res <- multi_match(scan1, scan1, c(200, 200))
  expect_equal(unname(res), rep(1, length(res)), tolerance = 1e-8)
})

test_that("EMD position similarity accounts for all fixations", {
  scan1 <- make_scanpath(c(0, 50, 100), c(0, 50, 100), c(0, 1, 2), c(100, 100, 100))
  scan2 <- make_scanpath(c(0, 50, 150), c(0, 50, 150), c(0, 1, 2), c(100, 100, 100))
  res <- multi_match(scan1, scan2, c(200, 200))
  expect_lt(res[["mm_position_emd"]], 1)
  expect_gt(res[["mm_position_emd"]], 0)
})
