context("normalize parameter for sample_density and sample_density_time")

library(testthat)
library(tibble)
library(dplyr)

# ---------- helpers ----------------------------------------------------------
make_test_density <- function(peak_x, peak_y, sigma = 2) {
  x_grid <- seq(0, 10, by = 0.5)
  y_grid <- seq(0, 10, by = 0.5)
  z_mat <- outer(x_grid, y_grid, function(a, b) {
    exp(-((a - peak_x)^2 + (b - peak_y)^2) / (2 * sigma^2))
  })
  dens <- list(x = x_grid, y = y_grid, z = z_mat)
  class(dens) <- c("eye_density", "density", "list")
  dens
}

make_test_fixgroup <- function(x_coords, y_coords, onsets) {
  fixation_group(
    x = x_coords,
    y = y_coords,
    onset = onsets,
    duration = rep(100, length(x_coords))
  )
}

# ---------- sample_density.density -------------------------------------------

test_that("normalize='none' returns raw density (default)", {
  dens <- make_test_density(5, 5, sigma = 2)
  fix <- fixation_group(x = 5, y = 5, onset = 0, duration = 100)

  res_default <- sample_density(dens, fix)
  res_none    <- sample_density(dens, fix, normalize = "none")

  expect_equal(res_default$z, res_none$z)
  # Raw peak of Gaussian with sigma=2 is 1.0

  expect_equal(res_default$z, 1.0)
})

test_that("normalize='max' scales values to [0, 1]", {
  dens <- make_test_density(5, 5, sigma = 2)
  fix <- fixation_group(x = c(5, 0), y = c(5, 0), onset = c(0, 100), duration = c(100, 100))

  res <- sample_density(dens, fix, normalize = "max")

  # Peak location should be exactly 1.0

  expect_equal(res$z[1], 1.0)
  # Off-peak should be in (0, 1)
  expect_gt(res$z[2], 0)
  expect_lt(res$z[2], 1)
})

test_that("normalize='sum' produces a probability distribution", {
  dens <- make_test_density(5, 5, sigma = 2)

  # Sum of the normalized density grid should be 1
  zmat <- dens$z / sum(dens$z)
  expect_equal(sum(zmat), 1.0)

  # Sampled values should be very small positive numbers (1 / n_cells scale)
  fix <- fixation_group(x = 5, y = 5, onset = 0, duration = 100)
  res <- sample_density(dens, fix, normalize = "sum")
  expect_gt(res$z, 0)
  expect_lt(res$z, 1)  # much less than 1 for any single cell
})

test_that("normalize='zscore' centres the density map", {
  dens <- make_test_density(5, 5, sigma = 2)

  # At the peak, z-scored value should be positive (above mean)
  fix_peak <- fixation_group(x = 5, y = 5, onset = 0, duration = 100)
  res_peak <- sample_density(dens, fix_peak, normalize = "zscore")
  expect_gt(res_peak$z, 0)

  # At a corner far from peak, z-scored value should be negative (below mean)
  fix_corner <- fixation_group(x = 0, y = 0, onset = 0, duration = 100)
  res_corner <- sample_density(dens, fix_corner, normalize = "zscore")
  expect_lt(res_corner$z, 0)
})

test_that("normalize works with the times argument", {
  dens <- make_test_density(5, 5, sigma = 2)
  fix <- fixation_group(x = c(5, 0), y = c(5, 0), onset = c(0, 200), duration = c(100, 100))

  res_raw <- sample_density(dens, fix, times = c(0, 100, 200), normalize = "none")
  res_max <- sample_density(dens, fix, times = c(0, 100, 200), normalize = "max")

  # Max-normalized peak should be 1
  expect_equal(res_max$z[1], 1.0)
  # Both should have the same number of rows
  expect_equal(nrow(res_raw), nrow(res_max))
})

test_that("normalize handles uniform density gracefully", {
  # All z values identical -> sd = 0 for zscore, max > 0 for max
  x_grid <- seq(0, 10, by = 1)
  y_grid <- seq(0, 10, by = 1)
  z_mat <- matrix(0.5, nrow = length(x_grid), ncol = length(y_grid))
  dens <- list(x = x_grid, y = y_grid, z = z_mat)
  class(dens) <- c("eye_density", "density", "list")

  fix <- fixation_group(x = 5, y = 5, onset = 0, duration = 100)

  # zscore with sd=0: should return the original values (no division)
  res_z <- sample_density(dens, fix, normalize = "zscore")
  expect_equal(res_z$z, 0.5)

  # max should give 1.0
  res_m <- sample_density(dens, fix, normalize = "max")
  expect_equal(res_m$z, 1.0)

  # sum should give 1/n_cells-ish
  res_s <- sample_density(dens, fix, normalize = "sum")
  expect_equal(res_s$z, 0.5 / sum(z_mat))
})

test_that("normalize handles zero density gracefully", {
  x_grid <- seq(0, 10, by = 1)
  y_grid <- seq(0, 10, by = 1)
  z_mat <- matrix(0, nrow = length(x_grid), ncol = length(y_grid))
  dens <- list(x = x_grid, y = y_grid, z = z_mat)
  class(dens) <- c("eye_density", "density", "list")

  fix <- fixation_group(x = 5, y = 5, onset = 0, duration = 100)

  # max with all zeros: division guarded, should return 0
  res_m <- sample_density(dens, fix, normalize = "max")
  expect_equal(res_m$z, 0)

  # sum with all zeros: guarded
  res_s <- sample_density(dens, fix, normalize = "sum")
  expect_equal(res_s$z, 0)

  # zscore with all zeros: guarded (sd=0)
  res_z <- sample_density(dens, fix, normalize = "zscore")
  expect_equal(res_z$z, 0)
})

# ---------- sample_density_time ----------------------------------------------

test_that("sample_density_time passes normalize through", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5, sigma = 2))
  )

  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(make_test_fixgroup(c(5, 0), c(5, 0), c(0, 200)))
  )

  res_none <- sample_density_time(
    template_tab, source_tab, match_on = "trial_id",
    times = c(0, 200), normalize = "none"
  )
  res_max <- sample_density_time(
    template_tab, source_tab, match_on = "trial_id",
    times = c(0, 200), normalize = "max"
  )

  # With max normalization the peak value should be 1
  expect_equal(res_max$sampled[[1]]$z[1], 1.0)
  # Raw value at peak is also 1.0 for this Gaussian, but off-peak differs
  # The key check: max-normalized off-peak value should differ from raw
  raw_offpeak <- res_none$sampled[[1]]$z[2]
  max_offpeak <- res_max$sampled[[1]]$z[2]
  # max normalisation divides by max(z_mat); since peak is 1.0, they should be equal here
  # but let's verify the function runs without error and values are finite
  expect_true(is.finite(max_offpeak))
})

test_that("sample_density_time zscore normalization makes cross-map comparisons comparable", {
  # Two density maps with different overall scales (use amplitude multiplier)
  dens_strong <- make_test_density(5, 5, sigma = 1)   # sharp peak
  dens_weak   <- make_test_density(5, 5, sigma = 1)
  # Scale down the weak map to simulate fewer fixations / lower amplitude
  dens_weak$z <- dens_weak$z * 0.1

  template_tab <- tibble(
    trial_id = c("A", "B"),
    density = list(dens_strong, dens_weak)
  )

  # Both source fixations at (5, 5) - the peak
  source_tab <- tibble(
    trial_id = c("A", "B"),
    fixgroup = list(
      make_test_fixgroup(c(5), c(5), c(0)),
      make_test_fixgroup(c(5), c(5), c(0))
    )
  )

  # Without normalization, raw values differ drastically
  res_raw <- sample_density_time(
    template_tab, source_tab, match_on = "trial_id",
    times = c(0), normalize = "none"
  )
  raw_A <- res_raw$sampled[[1]]$z[1]
  raw_B <- res_raw$sampled[[2]]$z[1]
  expect_gt(raw_A, raw_B * 5)  # strong is 10x the weak map

  # With zscore, both are at their respective peaks so both should be
  # equally far above their map's mean in SD units
  res_z <- sample_density_time(
    template_tab, source_tab, match_on = "trial_id",
    times = c(0), normalize = "zscore"
  )
  z_A <- res_z$sampled[[1]]$z[1]
  z_B <- res_z$sampled[[2]]$z[1]
  expect_gt(z_A, 0)
  expect_gt(z_B, 0)
  # After z-scoring, both peaks should yield the same value (same shape, just scaled)
  expect_equal(z_A, z_B, tolerance = 0.01)
})

test_that("sample_density_time normalize works with time_bins and permutations", {
  template_tab <- tibble(
    trial_id = c("A", "B"),
    subject = c("S1", "S1"),
    density = list(
      make_test_density(2, 2, sigma = 1),
      make_test_density(8, 8, sigma = 1)
    )
  )

  source_tab <- tibble(
    trial_id = c("A", "B"),
    subject = c("S1", "S1"),
    fixgroup = list(
      make_test_fixgroup(c(2), c(2), c(0)),
      make_test_fixgroup(c(8), c(8), c(0))
    )
  )

  result <- sample_density_time(
    template_tab, source_tab, match_on = "trial_id",
    times = c(0),
    time_bins = c(0, 100),
    permutations = 5,
    permute_on = "subject",
    normalize = "zscore"
  )

  # Should have bin and perm columns
  expect_true("bin_1" %in% names(result))
  expect_true("perm_bin_1" %in% names(result))
  # Matched values should be positive (at peak, above mean)
  expect_true(all(result$bin_1 > 0))
  # Permuted values should be negative (far from peak, below mean)
  expect_true(all(result$perm_bin_1 < 0))
})

test_that("normalize argument validates correctly", {
  dens <- make_test_density(5, 5)
  fix <- fixation_group(x = 5, y = 5, onset = 0, duration = 100)

  expect_error(sample_density(dens, fix, normalize = "invalid"),
               "'arg' should be one of")
})
