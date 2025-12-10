context("sample_density_time")

library(testthat)
library(tibble)
library(dplyr)

# Create simple test density maps
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

# Create test fixation group
make_test_fixgroup <- function(x_coords, y_coords, onsets) {
  fixation_group(
    x = x_coords,
    y = y_coords,
    onset = onsets,
    duration = rep(100, length(x_coords))
  )
}

test_that("sample_density_time basic functionality works", {
  # Create template table with density maps
  template_tab <- tibble(
    trial_id = c("A", "B"),
    density = list(
      make_test_density(2, 2),
      make_test_density(8, 8)
    )
  )

  # Create source table with fixation groups
  source_tab <- tibble(
    trial_id = c("A", "B"),
    fixgroup = list(
      make_test_fixgroup(c(2, 3, 4), c(2, 3, 4), c(0, 100, 200)),
      make_test_fixgroup(c(8, 7, 6), c(8, 7, 6), c(0, 100, 200))
    )
  )

  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0, 100, 200)
  )

  expect_s3_class(result, "tbl_df")
  expect_true("sampled" %in% names(result))
  expect_equal(nrow(result), 2)


  # Check sampled is a list column with correct structure
  expect_true(is.list(result$sampled))
  expect_true(all(sapply(result$sampled, function(x) "z" %in% names(x))))
  expect_true(all(sapply(result$sampled, function(x) "time" %in% names(x))))
})

test_that("sample_density_time with time_bins aggregates correctly", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5))
  )

  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(
      make_test_fixgroup(
        c(5, 5, 5, 5, 5, 5),
        c(5, 5, 5, 5, 5, 5),
        c(0, 50, 100, 150, 200, 250)
      )
    )
  )

  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = seq(0, 300, by = 50),
    time_bins = c(0, 150, 300)
  )

  expect_true("bin_1" %in% names(result))
  expect_true("bin_2" %in% names(result))
  expect_equal(nrow(result), 1)

  # Values should be numeric
  expect_true(is.numeric(result$bin_1))
  expect_true(is.numeric(result$bin_2))
})

test_that("sample_density_time with permutations adds baseline columns", {
  template_tab <- tibble(
    trial_id = c("A", "B", "C"),
    subject = c("S1", "S1", "S1"),
    density = list(
      make_test_density(2, 2),
      make_test_density(5, 5),
      make_test_density(8, 8)
    )
  )

  source_tab <- tibble(
    trial_id = c("A", "B", "C"),
    subject = c("S1", "S1", "S1"),
    fixgroup = list(
      make_test_fixgroup(c(2, 3), c(2, 3), c(0, 100)),
      make_test_fixgroup(c(5, 6), c(5, 6), c(0, 100)),
      make_test_fixgroup(c(8, 7), c(8, 7), c(0, 100))
    )
  )

  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0, 100),
    time_bins = c(0, 100, 200),
    permutations = 10,
    permute_on = "subject"
  )

  expect_true("perm_sampled" %in% names(result))
  expect_true("perm_bin_1" %in% names(result))
  expect_true("diff_bin_1" %in% names(result))
})

test_that("sample_density_time handles missing matches", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5))
  )

  source_tab <- tibble(
    trial_id = c("A", "B"),  # B has no match
    fixgroup = list(
      make_test_fixgroup(c(5), c(5), c(0)),
      make_test_fixgroup(c(5), c(5), c(0))
    )
  )

  expect_warning(
    result <- sample_density_time(
      template_tab = template_tab,
      source_tab = source_tab,
      match_on = "trial_id",
      times = c(0, 100)
    ),
    "Did not find matching template"
  )

  # Should only have 1 row (the matched one)
  expect_equal(nrow(result), 1)
})

test_that("sample_density_time validates inputs", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5))
  )

  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(make_test_fixgroup(c(5), c(5), c(0)))
  )

  # Wrong match_on column

  expect_error(
    sample_density_time(
      template_tab = template_tab,
      source_tab = source_tab,
      match_on = "nonexistent",
      times = c(0, 100)
    ),
    "not found"
  )

  # Invalid time_bins
  expect_error(
    sample_density_time(
      template_tab = template_tab,
      source_tab = source_tab,
      match_on = "trial_id",
      times = c(0, 100),
      time_bins = c(100, 50)  # Not monotonic
    ),
    "monotonically increasing"
  )
})

test_that("sample_density_time preserves source_tab columns", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5))
  )

  source_tab <- tibble(
    trial_id = "A",
    subject = "S1",
    condition = "test",
    extra_col = 42,
    fixgroup = list(make_test_fixgroup(c(5), c(5), c(0)))
  )

  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0, 100)
  )

  expect_true("subject" %in% names(result))
  expect_true("condition" %in% names(result))
  expect_true("extra_col" %in% names(result))
  expect_equal(result$extra_col, 42)
})

test_that("sample_density_time with custom aggregate_fun works", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5))
  )

  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(
      make_test_fixgroup(
        c(5, 5, 5, 5),
        c(5, 5, 5, 5),
        c(0, 50, 100, 150)
      )
    )
  )

  result_mean <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = seq(0, 200, by = 50),
    time_bins = c(0, 200),
    aggregate_fun = mean
  )

  result_max <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = seq(0, 200, by = 50),
    time_bins = c(0, 200),
    aggregate_fun = max
  )

  # Max should be >= mean for the same data
  expect_true(result_max$bin_1 >= result_mean$bin_1)
})

# ============================================================================
# VALUE CORRECTNESS TESTS
# ============================================================================

test_that("sampled values are correct for known density peaks", {
  # Create a Gaussian density peaked at (5, 5)
  # At the peak, value should be ~1.0; away from peak, lower
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5, sigma = 2))
  )

  # Fixations: one at peak (5,5), one away (0,0), one in between (3,3)
  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(
      make_test_fixgroup(
        x_coords = c(5, 0, 3),
        y_coords = c(5, 0, 3),
        onsets = c(0, 100, 200)
      )
    )
  )

  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0, 100, 200)
  )

  sampled <- result$sampled[[1]]

  # Value at peak (5,5) should be highest (close to 1)
  val_at_peak <- sampled$z[sampled$time == 0]
  val_at_origin <- sampled$z[sampled$time == 100]
  val_at_middle <- sampled$z[sampled$time == 200]

  expect_gt(val_at_peak, 0.9)  # Should be close to 1 at peak
  expect_lt(val_at_origin, 0.1)  # Should be near 0 far from peak
  expect_gt(val_at_middle, val_at_origin)  # Middle should be higher than origin
  expect_lt(val_at_middle, val_at_peak)  # Middle should be lower than peak
})

test_that("temporal interpolation samples correct coordinates between fixations", {
  # Create density with distinct values at different locations
  # Peak at (2, 2) - fixation starts here
  # We'll sample at times between fixations to test interpolation
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(2, 2, sigma = 2))
  )

  # Three fixations: (2,2) at t=0, (8,8) at t=200, (8,8) at t=400
  # Interpolation at t=100 should still be at (2,2) due to forward-fill
  # t=300 should be at (8,8) due to forward-fill from t=200
  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(
      make_test_fixgroup(
        x_coords = c(2, 8, 8),
        y_coords = c(2, 8, 8),
        onsets = c(0, 200, 400)
      )
    )
  )

  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0, 100, 200, 300)  # 100 is between fix 1&2, 300 is between fix 2&3
  )

  sampled <- result$sampled[[1]]

  # At t=0, we're at (2,2) - peak of density
  val_t0 <- sampled$z[sampled$time == 0]
  # At t=100, forward-fill should keep us at (2,2)
  val_t100 <- sampled$z[sampled$time == 100]
  # At t=200, we move to (8,8) - far from peak
  val_t200 <- sampled$z[sampled$time == 200]
  # At t=300, forward-fill keeps us at (8,8)
  val_t300 <- sampled$z[sampled$time == 300]

  # t=0 and t=100 should both be at peak (high value)
  expect_equal(val_t0, val_t100, tolerance = 0.01)
  expect_gt(val_t0, 0.9)

  # t=200 and t=300 should both be far from peak (low value)
  expect_equal(val_t200, val_t300, tolerance = 0.01)
  expect_lt(val_t200, 0.1)
})

test_that("permutation baseline uses non-matching templates", {
  # Create templates with very different peaks
  # A peaks at (2,2), B peaks at (8,8)
  template_tab <- tibble(
    trial_id = c("A", "B"),
    subject = c("S1", "S1"),
    density = list(
      make_test_density(2, 2, sigma = 1),  # Sharp peak at (2,2)
      make_test_density(8, 8, sigma = 1)   # Sharp peak at (8,8)
    )
  )

  # Source A has fixations at (2,2) - matches template A's peak
  # Source B has fixations at (8,8) - matches template B's peak
  source_tab <- tibble(
    trial_id = c("A", "B"),
    subject = c("S1", "S1"),
    fixgroup = list(
      make_test_fixgroup(c(2), c(2), c(0)),
      make_test_fixgroup(c(8), c(8), c(0))
    )
  )

  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0),
    time_bins = c(0, 100),
    permutations = 10,
    permute_on = "subject"
  )

  # For trial A: fixations at (2,2), matched template peaks at (2,2)
  # So bin_1 should be high (~1), but perm uses template B which peaks at (8,8)
  # So perm_bin_1 should be low (~0)
  expect_gt(result$bin_1[1], 0.9)
  expect_lt(result$perm_bin_1[1], 0.1)
  expect_gt(result$diff_bin_1[1], 0.8)  # Large positive difference

  # For trial B: fixations at (8,8), matched template peaks at (8,8)
  # So bin_1 should be high, perm uses template A which peaks at (2,2)
  expect_gt(result$bin_1[2], 0.9)
  expect_lt(result$perm_bin_1[2], 0.1)
  expect_gt(result$diff_bin_1[2], 0.8)
})

test_that("NULL density in template is handled gracefully",
{
  template_tab <- tibble(
    trial_id = c("A", "B"),
    density = list(
      make_test_density(5, 5),
      NULL  # NULL density for trial B
    )
  )

  source_tab <- tibble(
    trial_id = c("A", "B"),
    fixgroup = list(
      make_test_fixgroup(c(5), c(5), c(0)),
      make_test_fixgroup(c(5), c(5), c(0))
    )
  )

  # Should not error, but B should have NA values
  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0, 100),
    time_bins = c(0, 200)
  )

  # Trial A should have valid values
  expect_false(is.na(result$bin_1[result$trial_id == "A"]))

  # Trial B should have NA (NULL template)
  expect_true(is.na(result$bin_1[result$trial_id == "B"]))
})

test_that("multi-subject permutation stays within subject", {
  # Two subjects, each with 2 trials
  # Subject 1: templates peak at (2,2) and (3,3)

  # Subject 2: templates peak at (7,7) and (8,8)
  template_tab <- tibble(
    trial_id = c("S1_A", "S1_B", "S2_A", "S2_B"),
    subject = c("S1", "S1", "S2", "S2"),
    density = list(
      make_test_density(2, 2, sigma = 1),
      make_test_density(3, 3, sigma = 1),
      make_test_density(7, 7, sigma = 1),
      make_test_density(8, 8, sigma = 1)
    )
  )

  # Each source fixates at its matched template's peak
  source_tab <- tibble(
    trial_id = c("S1_A", "S1_B", "S2_A", "S2_B"),
    subject = c("S1", "S1", "S2", "S2"),
    fixgroup = list(
      make_test_fixgroup(c(2), c(2), c(0)),
      make_test_fixgroup(c(3), c(3), c(0)),
      make_test_fixgroup(c(7), c(7), c(0)),
      make_test_fixgroup(c(8), c(8), c(0))
    )
  )

  set.seed(123)  # For reproducibility
  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0),
    time_bins = c(0, 100),
    permutations = 10,
    permute_on = "subject"
  )

  # All matched values should be high (at peak)
  expect_true(all(result$bin_1 > 0.9))

  # For S1 trials: permuted should use other S1 template (nearby peak)
  # So perm values should still be moderate (not near 0)
  s1_perm <- result$perm_bin_1[result$subject == "S1"]
  expect_true(all(s1_perm > 0.3))  # Still moderate because S1 peaks are close

  # For S2 trials: permuted should use other S2 template (nearby peak)
  s2_perm <- result$perm_bin_1[result$subject == "S2"]
  expect_true(all(s2_perm > 0.3))

  # Key test: S1 permutations should NOT use S2 templates
  # If they did, S1 fixations at (2,2)/(3,3) would get ~0 from S2's (7,7)/(8,8) peaks
  # The fact that perm values are moderate proves within-subject permutation
})

test_that("bin boundaries handle edge cases correctly", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5))
  )

  # Fixations at exact bin boundaries
  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(
      make_test_fixgroup(
        c(5, 5, 5, 5, 5),
        c(5, 5, 5, 5, 5),
        c(0, 100, 200, 300, 400)  # Boundaries: 0, 100, 200, 300, 400
      )
    )
  )

  # Bins: [0-200), [200-400), [400-500)
  # With times at boundaries, test which bin they fall into
  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0, 100, 200, 300, 400),
    time_bins = c(0, 200, 400, 500)
  )

  # Should have 3 bins
  expect_true("bin_1" %in% names(result))
  expect_true("bin_2" %in% names(result))
  expect_true("bin_3" %in% names(result))

  # All values should be valid (not NA) since all fixations are at (5,5) on peak
  expect_false(is.na(result$bin_1))
  expect_false(is.na(result$bin_2))
  expect_false(is.na(result$bin_3))
})

test_that("empty time bins return NA", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5))
  )

  # Fixations only in early time range
  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(
      make_test_fixgroup(
        c(5, 5),
        c(5, 5),
        c(0, 50)
      )
    )
  )

  # Sample times only in first 100ms, but bins extend to 300ms
  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = c(0, 50),  # Only early times
    time_bins = c(0, 100, 200, 300)  # Bins 2 and 3 have no samples
  )

  # Bin 1 should have values
  expect_false(is.na(result$bin_1))

  # Bins 2 and 3 should be NA (no time points in those ranges)
  expect_true(is.na(result$bin_2))
  expect_true(is.na(result$bin_3))
})

test_that("sampled time series has correct length", {
  template_tab <- tibble(
    trial_id = "A",
    density = list(make_test_density(5, 5))
  )

  source_tab <- tibble(
    trial_id = "A",
    fixgroup = list(make_test_fixgroup(c(5), c(5), c(0)))
  )

  times <- seq(0, 1000, by = 25)  # 41 time points

  result <- sample_density_time(
    template_tab = template_tab,
    source_tab = source_tab,
    match_on = "trial_id",
    times = times
  )

  # Sampled data frame should have same number of rows as time points
  expect_equal(nrow(result$sampled[[1]]), length(times))
  expect_equal(result$sampled[[1]]$time, times)
})
