# Regression tests for defects found while auditing eyesim against the eyes4s
# parity baseline. Each block names the audit item it covers.

audit_fg <- function() {
  fixation_group(x = c(10, 30, 80), y = c(10, 20, 40),
                 onset = c(0, 100, 200), duration = c(1, 2, 4))
}

audit_density <- function(fg, ...) {
  suppressMessages(eye_density(fg, sigma = 10, xbounds = c(0, 100), ybounds = c(0, 50),
                               outdim = c(5, 3), ...))
}

# Audit item 1 ---------------------------------------------------------------
test_that("eye_density honours explicit fixation weights", {
  fg <- audit_fg()
  d0 <- audit_density(fg)
  dw <- audit_density(fg, weights = c(4, 2, 1))
  expect_false(isTRUE(all.equal(d0$z, dw$z)))

  # Explicit weights equal to the durations reproduce duration weighting.
  expect_equal(audit_density(fg, weights = c(1, 2, 4))$z,
               audit_density(fg, duration_weighted = TRUE)$z)

  # Explicit weights take precedence over duration weighting.
  expect_equal(audit_density(fg, weights = c(4, 2, 1), duration_weighted = TRUE)$z, dw$z)

  # A zero weight removes a fixation: the map equals one built without it.
  expect_equal(audit_density(fg, weights = c(1, 1, 0))$z,
               audit_density(fg[1:2, ])$z)

  # Weights are aligned with the rows that survive the window filter.
  expect_equal(audit_density(fg, weights = c(4, 2, 1), window = c(50, 300))$z,
               audit_density(fg[2:3, ], weights = c(2, 1))$z)

  expect_error(audit_density(fg, weights = c(1, 2)), "one value per fixation")
  expect_error(audit_density(fg, weights = c(1, -1, 2)), "non-negative")
})
