library(testthat)
library(tibble)

make_entropy_density <- function(z) {
  structure(
    list(
      z = z,
      x = seq_len(nrow(z)),
      y = seq_len(ncol(z)),
      sigma = 1
    ),
    class = c("eye_density", "density", "list")
  )
}

test_that("fixation_entropy returns 0 for a point mass density", {
  dens <- make_entropy_density(matrix(c(1, 0, 0, 0), nrow = 2, byrow = TRUE))

  expect_equal(fixation_entropy(dens, normalize = FALSE), 0)
  expect_equal(fixation_entropy(dens, normalize = TRUE), 0)
})

test_that("fixation_entropy returns maximal entropy for a uniform density", {
  dens <- make_entropy_density(matrix(1, nrow = 2, ncol = 2))

  expect_equal(fixation_entropy(dens, normalize = TRUE), 1)
  expect_equal(fixation_entropy(dens, normalize = FALSE, base = 2), 2)
})

test_that("fixation_entropy is invariant to density scaling", {
  dens1 <- make_entropy_density(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  dens2 <- make_entropy_density(matrix(c(10, 20, 30, 40), nrow = 2, byrow = TRUE))

  expect_equal(
    fixation_entropy(dens1, normalize = TRUE),
    fixation_entropy(dens2, normalize = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    fixation_entropy(dens1, normalize = FALSE),
    fixation_entropy(dens2, normalize = FALSE),
    tolerance = 1e-12
  )
})

test_that("grid fixation entropy matches analytic occupancy probabilities", {
  fg <- fixation_group(
    x = c(0.1, 0.2, 0.8, 0.85),
    y = c(0.1, 0.2, 0.8, 0.85),
    onset = c(0, 100, 200, 300),
    duration = rep(100, 4)
  )

  # Two occupied cells with probability 0.5 each
  expect_equal(
    fixation_entropy(fg, method = "grid", grid = c(2, 2), xbounds = c(0, 1), ybounds = c(0, 1),
                     normalize = FALSE, base = 2),
    1
  )
  expect_equal(
    fixation_entropy(fg, method = "grid", grid = c(2, 2), xbounds = c(0, 1), ybounds = c(0, 1),
                     normalize = TRUE, base = 2),
    0.5
  )
})

test_that("grid fixation entropy is invariant to fixation order", {
  fg <- fixation_group(
    x = c(0.1, 0.2, 0.8, 0.85, 0.55),
    y = c(0.1, 0.2, 0.8, 0.85, 0.45),
    onset = c(0, 100, 200, 300, 400),
    duration = rep(100, 5)
  )
  shuffled <- fg[c(5, 3, 1, 4, 2), ]

  expect_equal(
    fixation_entropy(fg, method = "grid", grid = c(3, 3), xbounds = c(0, 1), ybounds = c(0, 1)),
    fixation_entropy(shuffled, method = "grid", grid = c(3, 3), xbounds = c(0, 1), ybounds = c(0, 1))
  )
})

test_that("fixation group density entropy matches derived eye density entropy", {
  fg <- fixation_group(
    x = c(0.2, 0.25, 0.75, 0.8),
    y = c(0.2, 0.25, 0.75, 0.8),
    onset = c(0, 100, 200, 300),
    duration = rep(100, 4)
  )
  dens <- eye_density(fg, sigma = 0.08, xbounds = c(0, 1), ybounds = c(0, 1), outdim = c(30, 30))

  expect_equal(
    fixation_entropy(fg, method = "density", sigma = 0.08, xbounds = c(0, 1), ybounds = c(0, 1), outdim = c(30, 30)),
    fixation_entropy(dens),
    tolerance = 1e-12
  )
})

test_that("multiscale fixation entropy supports mean and none aggregation", {
  fg <- fixation_group(
    x = c(0.2, 0.25, 0.75, 0.8),
    y = c(0.2, 0.25, 0.75, 0.8),
    onset = c(0, 100, 200, 300),
    duration = rep(100, 4)
  )
  dens_ms <- eye_density(fg, sigma = c(0.05, 0.15), xbounds = c(0, 1), ybounds = c(0, 1), outdim = c(25, 25))

  ent_none <- fixation_entropy(dens_ms, aggregate = "none")
  ent_mean <- fixation_entropy(dens_ms, aggregate = "mean")

  expect_length(ent_none, 2)
  expect_true(all(names(ent_none) %in% c("sigma_0.05", "sigma_0.15")))
  expect_equal(ent_mean, mean(ent_none), tolerance = 1e-12)
})

test_that("fixation_entropy returns NA for zero-mass densities", {
  dens <- make_entropy_density(matrix(0, nrow = 2, ncol = 2))
  expect_true(is.na(fixation_entropy(dens)))
})

test_that("fixation_entropy stays within [0, 1] when normalized", {
  set.seed(123)
  fg <- fixation_group(
    x = runif(25),
    y = runif(25),
    onset = seq(0, by = 100, length.out = 25),
    duration = rep(100, 25)
  )

  ent_grid <- fixation_entropy(fg, method = "grid", grid = c(5, 5), xbounds = c(0, 1), ybounds = c(0, 1))
  ent_dens <- fixation_entropy(fg, method = "density", sigma = 0.1, xbounds = c(0, 1), ybounds = c(0, 1),
                               outdim = c(25, 25))

  expect_gte(ent_grid, 0)
  expect_lte(ent_grid, 1)
  expect_gte(ent_dens, 0)
  expect_lte(ent_dens, 1)
})
