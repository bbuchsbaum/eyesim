library(testthat)

test_that("EMD similarity equals 1 for identical maps", {
  d <- list(x = 1:2, y = 1:2, z = matrix(c(0.25,0.25,0.25,0.25), nrow=2))
  class(d) <- c("eye_density", "density", "list")
  expect_equal(similarity(d, d, method="emd"), 1)
})

test_that("signed EMD returns 0 for identical residuals", {
  skip_if_not_installed("emdist")
  d <- list(x = 1:2, y = 1:2, z = matrix(c(0.3,0.2,0.2,0.3), nrow=2))
  class(d) <- c("eye_density", "density", "list")
  sal <- list(x = 1:2, y = 1:2, z = matrix(0.25, nrow=2, ncol=2))
  class(sal) <- c("eye_density", "density", "list")
  expect_equal(similarity(d, d, method="emd", saliency_map=sal), 0)
})

test_that("emdw normalises weights, so backends agree on unequal totals", {
  x <- cbind(c(0, 100, 200), 0)
  y <- cbind(c(0, 100, 300), 0)
  wx <- c(1, 1, 1)
  wy <- c(1, 1, 4)
  expected <- 350 / 3
  if (requireNamespace("emdist", quietly = TRUE)) {
    expect_equal(emdw(x, wx, y, wy), expected, tolerance = 1e-6)
  }
  skip_if_not_installed("transport")
  expect_equal(
    transport::wasserstein(transport::wpp(x, wx / 3), transport::wpp(y, wy / 6), p = 1),
    expected, tolerance = 1e-6
  )
})

test_that("emdw ignores zero-mass points and handles empty mass", {
  skip_if_not_installed("emdist")
  x <- cbind(c(0, 50, 100), 0)
  expect_equal(emdw(x, c(1, 0, 1), x, c(1, 0, 1)), 0)
  expect_equal(emdw(x, c(0, 0, 0), x, c(0, 0, 0)), 0)
  expect_true(is.na(emdw(x, c(0, 0, 0), x, c(1, 0, 0))))
})
