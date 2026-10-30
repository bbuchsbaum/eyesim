test_that("empty fixation groups can infer their durations", {
  inferred <- fixation_group(x = numeric(), y = numeric(), onset = numeric())
  explicit <- fixation_group(x = numeric(), y = numeric(), onset = numeric(),
                             duration = numeric())

  expect_s3_class(inferred, "fixation_group")
  expect_equal(nrow(inferred), 0L)
  expect_identical(inferred, explicit)
})

test_that("inferred durations end at zero for nonempty fixation groups", {
  single <- fixation_group(x = 1, y = 2, onset = 100)
  several <- fixation_group(x = 1:3, y = 4:6, onset = c(100, 110, 135))

  expect_identical(single$duration, 0)
  expect_identical(several$duration, c(10, 25, 0))
})
