# Evaluates `expr`, requiring at least one warning of `class` (or matching
# `regexp`) and muffling every such warning, however many there are. Small
# test fixtures trigger expected warnings (e.g. a calibration fallback), and
# leaving them unasserted would hide new ones in the test output.
expect_warnings <- function(expr, class = NULL, regexp = NULL) {
  seen <- 0L
  value <- withCallingHandlers(
    expr,
    warning = function(condition) {
      matches <- (is.null(class) || inherits(condition, class)) &&
        (is.null(regexp) || grepl(regexp, conditionMessage(condition)))
      if (matches) {
        seen <<- seen + 1L
        invokeRestart("muffleWarning")
      }
    }
  )
  testthat::expect_gt(seen, 0L, label = "number of expected warnings")
  value
}

expect_calibration_fallback <- function(expr) {
  expect_warnings(expr, class = "gaze_calibration_fallback")
}

expect_cca_group_fallback <- function(expr) {
  expect_warnings(expr, regexp = "no CCA model for group")
}
