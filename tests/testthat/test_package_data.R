test_that("participant-level Wynn datasets are not distributed", {
  available <- utils::data(package = "eyesim")$results
  expect_false(any(c("wynn_study", "wynn_test") %in% available[, "Item"]))
})
