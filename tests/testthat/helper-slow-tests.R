# Slow-test gate shared by every test file. EYESIM_SLOW_TESTS turns the slow
# checks on when it holds any non-empty value other than "false" or "0"
# (case-insensitive), so "true", "TRUE", "1" and "yes" all run them.
eyesim_slow_tests_enabled <- function(value = Sys.getenv("EYESIM_SLOW_TESTS")) {
  value <- trimws(value)
  nzchar(value) && !tolower(value) %in% c("false", "0")
}

skip_unless_slow_tests <- function(what = "slow check") {
  testthat::skip_if_not(
    eyesim_slow_tests_enabled(),
    paste0(what, " (set EYESIM_SLOW_TESTS=true)")
  )
}
