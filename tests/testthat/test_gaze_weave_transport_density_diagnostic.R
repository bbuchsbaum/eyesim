transport_density_source <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    gaze_weave_test_inst_path(
      "validation", "gaze-weave-transport-density-diagnostic.R"
    ),
    envir = environment
  )
  environment
}

test_that("graded response coding points toward oldness", {
  court <- transport_density_source()
  observed <- court$transport_density_oldness(c(0, 1, 2, 3, 4, NA))

  expect_s3_class(observed, "ordered")
  expect_identical(as.integer(observed), c(NA_integer_, 4:1, NA_integer_))
})

test_that("bare saliency fit recovers the 20-to-100 contrast", {
  court <- transport_density_source()
  fixture <- data.frame(
    saliency_z = rep(c(-2, -1, 0, 1, 2), 4),
    z_quality = rep(c(-1, 1), 10),
    z_duration = rep(c(-1, -1, 1, 1), 5)
  )
  fixture$gaze_info_bits <- 0.1 + 0.25 * fixture$saliency_z

  expect_equal(
    court$transport_density_fit_saliency(fixture), 1,
    tolerance = 1e-12
  )
  expect_equal(
    court$transport_density_fit_saliency(fixture, adjusted = TRUE), 1,
    tolerance = 1e-12
  )
})

test_that("ordinal loss selects the observed response probability", {
  court <- transport_density_source()
  probability <- matrix(
    c(0.1, 0.2, 0.3, 0.4, 0.4, 0.3, 0.2, 0.1),
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, as.character(1:4))
  )
  oldness <- ordered(c(4, 1), levels = 1:4)

  expect_equal(
    court$transport_density_ordinal_log_loss(probability, oldness),
    -log(c(0.4, 0.4))
  )
})

test_that("local Transport and density panels share frozen support", {
  result_dir <- gaze_weave_test_inst_path(
    "validation", "gaze-weave-transport-v3-retrieval-results"
  )
  comparator_dir <- gaze_weave_test_inst_path(
    "validation", "gaze-weave-recognition-full-cohort-results"
  )
  testthat::skip_if_not(
    file.exists(file.path(result_dir, "checkpoint-manifest.csv")) &&
      file.exists(file.path(comparator_dir, "manifest-md5.csv")),
    "Local-only full-cohort results are unavailable."
  )
  court <- transport_density_source()
  panel <- court$transport_density_read_panel(result_dir, comparator_dir)

  expect_equal(nrow(panel), 4L * 1295L)
  expect_identical(
    as.integer(table(panel$method)), rep(1295L, 4L)
  )
  expect_identical(
    sort(unique(panel$probe_type)), c("lure", "old")
  )
  expect_identical(
    as.integer(table(panel$probe_type[panel$method == "transport"])),
    c(655L, 640L)
  )
})
