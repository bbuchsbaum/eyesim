# Regression tests for changes made to fix bugs and remove deprecated code.

test_that("eye_table constructs correctly with groupvar and vars", {
  set.seed(42)
  df <- data.frame(
    xpos = runif(20, 0, 500),
    ypos = runif(20, 0, 500),
    dur = abs(rnorm(20, 200, 30)),
    ons = rep(cumsum(abs(rnorm(10, 300, 50))), 2),
    trial = rep(c("A", "B"), each = 10),
    cond = rep(c("old", "new"), each = 10)
  )

  et <- eye_table("xpos", "ypos", "dur", "ons",
                  groupvar = "trial", vars = "cond", data = df,
                  clip_bounds = c(0, 500, 0, 500))

  expect_s3_class(et, "eye_table")
  expect_equal(nrow(et), 2)
  expect_true("fixgroup" %in% names(et))
  expect_true("cond" %in% names(et))

  # Each fixgroup should be a fixation_group

  expect_s3_class(et$fixgroup[[1]], "fixation_group")
  expect_s3_class(et$fixgroup[[2]], "fixation_group")
})

test_that("simulate_eye_table generates per-group onsets", {
  set.seed(1)
  et <- simulate_eye_table(n_fixations = 40, n_groups = 4,
                           clip_bounds = c(0, 500, 0, 500))

  expect_s3_class(et, "eye_table")
  expect_equal(nrow(et), 4)

  # Each group's fixation_group should have onsets starting near 0
  # (not continuing from the previous group)
  for (i in seq_len(nrow(et))) {
    fg <- et$fixgroup[[i]]
    expect_true(fg$onset[1] < 1000,
                info = paste("Group", i, "onset starts at", fg$onset[1],
                             "- should be near 0, not cumulative"))
  }
})

test_that("normalize.fixation_group scales to [0,1] for arbitrary bounds", {
  fg <- fixation_group(
    x = c(100, 500, 900),
    y = c(200, 600, 1000),
    onset = c(0, 200, 400),
    duration = c(200, 200, 200)
  )

  normed <- normalize(fg, xbounds = c(100, 900), ybounds = c(200, 1000))

  # Should be in [0, 1]
  expect_equal(min(normed$x), 0)
  expect_equal(max(normed$x), 1)
  expect_equal(min(normed$y), 0)
  expect_equal(max(normed$y), 1)

  # Middle value should be 0.5
  expect_equal(normed$x[2], 0.5)
  expect_equal(normed$y[2], 0.5)
})

test_that("estimate_scale filters y (not x) when window is provided", {
  # Create two fixation groups with different time ranges
  fg_ref <- fixation_group(
    x = c(100, 200, 300),
    y = c(100, 200, 300),
    onset = c(0, 100, 200),
    duration = c(100, 100, 100)
  )

  fg_source <- fixation_group(
    x = c(100, 200, 300, 400),
    y = c(100, 200, 300, 400),
    onset = c(0, 100, 200, 300),
    duration = c(100, 100, 100, 100)
  )

  # With a window, only source fixations in that window should be used

  result <- eyesim:::estimate_scale(fg_ref, fg_source, window = c(0, 250))
  expect_true(is.list(result))
  expect_true(length(result$par) == 2)
  # Should not error (the old bug caused subsetting of x instead of y)
})

test_that("eye_density rejects non-positive sigma", {
  fg <- fixation_group(
    x = c(100, 200, 300),
    y = c(100, 150, 200),
    onset = c(0, 200, 400),
    duration = c(200, 200, 200)
  )

  expect_error(eye_density(fg, sigma = 0, xbounds = c(0, 400), ybounds = c(0, 300)),
               "sigma must be a positive")
  expect_error(eye_density(fg, sigma = -10, xbounds = c(0, 400), ybounds = c(0, 300)),
               "sigma must be a positive")
  # Valid sigma should work
  expect_s3_class(
    eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300)),
    "eye_density"
  )
})

test_that("as.data.frame.eye_density returns correct grid", {
  fg <- fixation_group(
    x = c(100, 200, 300),
    y = c(100, 150, 200),
    onset = c(0, 200, 400),
    duration = c(200, 200, 200)
  )
  ed <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300),
                    outdim = c(10, 10))
  df <- as.data.frame(ed)

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 100)  # 10 x 10 grid
  expect_true(all(c("x", "y", "z") %in% names(df)))
  expect_equal(length(unique(df$x)), 10)
  expect_equal(length(unique(df$y)), 10)
})
