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

# CRAN tier. testthat and devtools::test() set NOT_CRAN=true, and so does the
# GitHub Actions check, so everything below runs in full there. On CRAN the
# heavy statistical and integration checks are skipped, and cross-platform
# parity checks run on a smaller sample so they still cover CRAN's platforms.
eyesim_on_cran <- function() {
  !identical(tolower(Sys.getenv("NOT_CRAN")), "true")
}

cran_sample_size <- function(full, cran) {
  if (eyesim_on_cran()) cran else full
}
