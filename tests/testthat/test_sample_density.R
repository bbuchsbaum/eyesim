context("sample_density vectorization")

library(testthat)

# Simple density map
x_grid <- seq(0, 10, by = 1)
y_grid <- seq(0, 10, by = 1)
z_mat <- outer(x_grid, y_grid, function(a,b) a + b)
dens <- list(x = x_grid, y = y_grid, z = z_mat)
class(dens) <- c("eye_density", "density", "list")

# fixation group
fix <- fixation_group(x = c(1,5,10),
                      y = c(2,5,9),
                      onset = c(0,50,100),
                      duration = rep(1,3))

test_that("times argument matches direct sampling", {
  direct <- sample_density(dens, fix)
  timed  <- sample_density(dens, fix, times = fix$onset)
  expect_equal(direct, timed)
})
