test_that("multi_match matches Python multimatch_gaze on simple paths (no grouping)", {
  testthat::skip_on_cran()
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    testthat::skip("reticulate not available")
  }
  if (!reticulate::py_module_available("multimatch_gaze")) {
    testthat::skip("Python module 'multimatch_gaze' not available")
  }

  # Construct two small scanpaths with clear geometry
  fg1 <- fixation_group(
    x = c(100, 200, 260, 300),
    y = c(100, 140, 120, 180),
    duration = c(0.2, 0.25, 0.2, 0.3),
    onset = c(0.1, 0.3, 0.6, 1.0)
  )
  fg2 <- fixation_group(
    x = c(120, 210, 255, 310),
    y = c(110, 150, 130, 175),
    duration = c(0.22, 0.24, 0.21, 0.28),
    onset = c(0.05, 0.35, 0.7, 1.2)
  )

  sp1 <- scanpath(fg1)
  sp2 <- scanpath(fg2)
  screensize <- c(500, 500)

  # R implementation
  r_mm <- multi_match(sp1, sp2, screensize)

  # Python reference via reticulate (no simplification/grouping)
  py_mm <- eyesim:::py_multi_match(fg1, fg2, screensize,
                                   grouping = FALSE, tdir = 0, tdur = 0, tamp = 0)

  # Compare the five matching dimensions
  r_subset <- r_mm[c("mm_vector", "mm_direction", "mm_length", "mm_position", "mm_duration")]
  testthat::expect_equal(unname(r_subset), unname(py_mm), tolerance = 1e-6)
})

test_that("multi_match matches Python multimatch_gaze on random paths (no grouping)", {
  testthat::skip_on_cran()
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    testthat::skip("reticulate not available")
  }
  if (!reticulate::py_module_available("multimatch_gaze")) {
    testthat::skip("Python module 'multimatch_gaze' not available")
  }

  set.seed(42)
  for (k in 1:3) {
    nfix1 <- sample(4:7, 1)
    nfix2 <- sample(4:7, 1)

    # Make random increasing onsets and positions
    dur1 <- runif(nfix1, 0.05, 0.3)
    dur2 <- runif(nfix2, 0.05, 0.3)
    fg1 <- fixation_group(
      x = cumsum(runif(nfix1, -40, 40)) + 250,
      y = cumsum(runif(nfix1, -40, 40)) + 250,
      duration = dur1,
      onset = cumsum(runif(nfix1, 0.05, 0.4))
    )
    fg2 <- fixation_group(
      x = cumsum(runif(nfix2, -40, 40)) + 250,
      y = cumsum(runif(nfix2, -40, 40)) + 250,
      duration = dur2,
      onset = cumsum(runif(nfix2, 0.05, 0.4))
    )

    sp1 <- scanpath(fg1)
    sp2 <- scanpath(fg2)
    screensize <- c(640, 480)

    r_mm <- multi_match(sp1, sp2, screensize)
    py_mm <- eyesim:::py_multi_match(fg1, fg2, screensize,
                                     grouping = FALSE, tdir = 0, tdur = 0, tamp = 0)

    r_subset <- r_mm[c("mm_vector", "mm_direction", "mm_length", "mm_position", "mm_duration")]
    testthat::expect_equal(unname(r_subset), unname(py_mm), tolerance = 1e-6)
  }
})

