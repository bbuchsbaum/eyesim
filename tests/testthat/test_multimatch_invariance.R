test_that("translation invariance holds for non-position dimensions and matches Python", {
  testthat::skip_on_cran()
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    testthat::skip("reticulate not available")
  }
  if (!reticulate::py_module_available("multimatch_gaze")) {
    testthat::skip("Python module 'multimatch_gaze' not available")
  }

  mmgaze <- reticulate::import("multimatch_gaze")

  set.seed(7)
  n <- 6
  x <- cumsum(runif(n, -50, 50)) + 300
  y <- cumsum(runif(n, -50, 50)) + 200
  dur <- runif(n, 0.05, 0.3)
  ons <- cumsum(runif(n, 0.05, 0.4))
  fg1 <- fixation_group(x, y, dur, ons)

  # Translate by a constant offset
  dx <- 37; dy <- -21
  fg2 <- fixation_group(x + dx, y + dy, dur, ons)

  sp1 <- scanpath(fg1)
  sp2 <- scanpath(fg2)
  screensize <- c(640, 480)

  r_mm <- multi_match(sp1, sp2, screensize)
  py_mm <- eyesim:::py_multi_match(fg1, fg2, screensize, grouping = FALSE, tdir = 0, tdur = 0, tamp = 0)

  # Non-position dimensions should be ~1 for pure translation
  expect_gt(r_mm[["mm_vector"]], 0.999)
  expect_gt(r_mm[["mm_direction"]], 0.999)
  expect_gt(r_mm[["mm_length"]], 0.999)
  expect_gt(r_mm[["mm_duration"]], 0.999)
  # Position should drop depending on magnitude of (dx,dy)
  expect_lt(r_mm[["mm_position"]], 1)

  # Parity with Python on the five metrics
  r_subset <- r_mm[c("mm_vector", "mm_direction", "mm_length", "mm_position", "mm_duration")]
  testthat::expect_equal(unname(r_subset), unname(py_mm), tolerance = 1e-6)
})

test_that("direction is scale-invariant and matches Python (non-grouping)", {
  testthat::skip_on_cran()
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    testthat::skip("reticulate not available")
  }
  if (!reticulate::py_module_available("multimatch_gaze")) {
    testthat::skip("Python module 'multimatch_gaze' not available")
  }

  set.seed(8)
  n <- 6
  x <- cumsum(runif(n, -40, 40)) + 200
  y <- cumsum(runif(n, -40, 40)) + 300
  dur <- runif(n, 0.05, 0.3)
  ons <- cumsum(runif(n, 0.05, 0.4))
  fg1 <- fixation_group(x, y, dur, ons)

  s <- 1.8
  fg2 <- fixation_group(x[1] + s*(x - x[1]), y[1] + s*(y - y[1]), dur, ons)

  sp1 <- scanpath(fg1)
  sp2 <- scanpath(fg2)
  screensize <- c(800, 600)

  r_mm <- multi_match(sp1, sp2, screensize)
  py_mm <- eyesim:::py_multi_match(fg1, fg2, screensize, grouping = FALSE, tdir = 0, tdur = 0, tamp = 0)

  # Direction should be ~1 despite scaling; vector/length/position may differ
  expect_gt(r_mm[["mm_direction"]], 0.999)

  r_subset <- r_mm[c("mm_vector", "mm_direction", "mm_length", "mm_position", "mm_duration")]
  testthat::expect_equal(unname(r_subset), unname(py_mm), tolerance = 1e-6)
})

