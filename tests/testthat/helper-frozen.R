# Frozen GazeWeave reference fits were produced on aarch64-apple-darwin.
# Floating-point results reproduce to the last bit only on that platform
# (other CPUs and BLAS/libm builds differ in the final bits), so bit-identity
# is required there and a tight tolerance is used elsewhere.
frozen_platform <- function() {
  grepl("aarch64-apple-darwin", R.version$platform, fixed = TRUE)
}

expect_frozen <- function(object, expected, tolerance = 1e-10) {
  if (frozen_platform()) {
    testthat::expect_identical(object, expected)
  } else {
    testthat::expect_equal(object, expected, tolerance = tolerance)
  }
}
