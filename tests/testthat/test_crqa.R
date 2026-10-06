context("crqa wrapper")

skip_if_not_installed("crqa")

library(testthat)

# Create simple fixation groups
fg1 <- data.frame(x = 1:5, y = 1:5)
fg2 <- data.frame(x = 2:6, y = 3:7)

# Use package function directly for reference
nr <- min(nrow(fg1), nrow(fg2))
ref <- crqa::crqa(as.matrix(fg1[1:nr, 1:2]),
                  as.matrix(fg2[1:nr, 1:2]),
                  method = "mdcrqa", radius = 60)

# Call wrapper
res <- crqa(fg1, fg2, radius = 60)

expect_equal(res, ref)


test_that("crqa uses the x/y columns of a fixation_group, not the index", {
  fga <- fixation_group(x = c(10, 200, 400, 30), y = c(500, 20, 300, 80),
                        onset = c(0, 200, 400, 600), duration = rep(200, 4))
  fgb <- fixation_group(x = c(15, 210, 390, 600), y = c(505, 25, 310, 900),
                        onset = c(0, 200, 400, 600), duration = rep(200, 4))
  ref <- crqa::crqa(as.matrix(fga[, c("x", "y")]), as.matrix(fgb[, c("x", "y")]),
                    method = "mdcrqa", radius = 60)
  expect_equal(crqa(fga, fgb, radius = 60), ref)
  expect_error(crqa(fga[0, ], fgb), "non-empty")
})
