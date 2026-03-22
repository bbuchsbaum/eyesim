library(testthat)
library(tibble)

test_that("add_scanpath builds scanpaths row-wise instead of reusing the first row", {
  fg1 <- fixation_group(
    x = c(1, 2, 3),
    y = c(1, 2, 3),
    onset = c(0, 100, 200),
    duration = c(100, 100, 100)
  )
  fg2 <- fixation_group(
    x = c(5, 6, 7),
    y = c(5, 4, 3),
    onset = c(0, 100, 200),
    duration = c(100, 100, 100)
  )

  df <- tibble(id = 1:2, fixgroup = list(fg1, fg2))
  out_df <- add_scanpath(df)

  expect_s3_class(out_df$scanpath[[1]], "scanpath")
  expect_s3_class(out_df$scanpath[[2]], "scanpath")
  expect_equal(out_df$scanpath[[1]]$x, fg1$x)
  expect_equal(out_df$scanpath[[2]]$x, fg2$x)
  expect_false(identical(out_df$scanpath[[1]]$x, out_df$scanpath[[2]]$x))
  expect_false(identical(out_df$scanpath[[1]]$theta, out_df$scanpath[[2]]$theta))

  et <- as_eye_table(df)
  out_et <- add_scanpath(et)

  expect_s3_class(out_et$scanpath[[1]], "scanpath")
  expect_s3_class(out_et$scanpath[[2]], "scanpath")
  expect_equal(out_et$scanpath[[1]]$x, fg1$x)
  expect_equal(out_et$scanpath[[2]]$x, fg2$x)
})
