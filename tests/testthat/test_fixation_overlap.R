library(testthat)

context("fixation_overlap")

test_that("fixation_overlap returns structure and valid ranges", {
  x <- seq(0, 90, by = 10)
  y <- x
  fg1 <- fixation_group(x = x, y = y, duration = rep(1, length(x)), onset = x)
  fg2 <- fixation_group(x = x, y = y, duration = rep(1, length(x)), onset = x)

  res <- fixation_overlap(fg1, fg2, time_samples = x, dthresh = 1e-6)

  expect_true(is.list(res))
  expect_equal(names(res), c("overlap", "perc"))
  expect_true(is.numeric(res$overlap))
  expect_equal(length(res$overlap), 1)
  expect_true(res$overlap >= 0 && res$overlap <= length(x))
  expect_true(is.numeric(res$perc))
  expect_equal(length(res$perc), 1)
  expect_equal(res$perc, res$overlap/length(x))
  expect_true(res$perc >= 0 && res$perc <= 1)
})
