library(testthat)

test_that("EMD similarity equals 1 for identical maps", {
  d <- list(x = 1:2, y = 1:2, z = matrix(c(0.25,0.25,0.25,0.25), nrow=2))
  class(d) <- c("eye_density", "density", "list")
  expect_equal(similarity(d, d, method="emd"), 1)
})

test_that("signed EMD returns 0 for identical residuals", {
  d <- list(x = 1:2, y = 1:2, z = matrix(c(0.3,0.2,0.2,0.3), nrow=2))
  class(d) <- c("eye_density", "density", "list")
  sal <- list(x = 1:2, y = 1:2, z = matrix(0.25, nrow=2, ncol=2))
  class(sal) <- c("eye_density", "density", "list")
  expect_equal(similarity(d, d, method="emd", saliency_map=sal), 0)
})
