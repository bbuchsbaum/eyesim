test_that("fixation_overlap matches deterministic overlap counts across metrics", {
  fg_ref <- fixation_group(
    x = c(0, 10, 20),
    y = c(0, 0, 0),
    onset = c(0, 100, 200),
    duration = c(100, 100, 100)
  )

  fg_shift <- fixation_group(
    x = c(2, 15, 40),
    y = c(1, 0, 0),
    onset = c(0, 100, 200),
    duration = c(100, 100, 100)
  )

  times <- c(0, 100, 200)

  euclidean <- fixation_overlap(fg_ref, fg_shift, dthresh = 6,
                                time_samples = times,
                                dist_method = "euclidean")
  manhattan <- fixation_overlap(fg_ref, fg_shift, dthresh = 6,
                                time_samples = times,
                                dist_method = "manhattan")

  expect_equal(euclidean$overlap, 2)
  expect_equal(euclidean$perc, 2 / 3)
  expect_equal(manhattan$overlap, 2)
  expect_equal(manhattan$perc, 2 / 3)
})

test_that("suggest_sigma enforces input requirements and display-based clamp", {
  expect_error(
    suggest_sigma(c(1, 2, 3)),
    "'y' must be provided"
  )

  fg_tight <- fixation_group(
    x = rep(400, 20),
    y = rep(300, 20),
    onset = seq(0, by = 50, length.out = 20),
    duration = rep(50, 20)
  )

  sigma <- suggest_sigma(fg_tight, xbounds = c(0, 1000), ybounds = c(0, 1000))

  expect_equal(sigma, 10)
})

test_that("c.fixation_group concatenates sequentially and rejects invalid inputs", {
  fg1 <- fixation_group(
    x = c(10, 20),
    y = c(5, 10),
    onset = c(0, 100),
    duration = c(100, 100)
  )

  fg2 <- fixation_group(
    x = c(30, 40),
    y = c(15, 20),
    onset = c(0, 100),
    duration = c(100, 100)
  )

  combined <- c(fg1, fg2)

  expect_s3_class(combined, "fixation_group")
  expect_equal(nrow(combined), 4)
  expect_equal(combined$index, 1:4)
  expect_equal(combined$onset, c(0, 100, 200, 300))

  expect_error(
    c(fg1, data.frame(x = 1)),
    "must be fixation_group"
  )
})
