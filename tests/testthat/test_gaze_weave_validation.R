validation_script <- system.file(
  "validation", "gaze-weave-comparison.R", package = "eyesim"
)
if (!nzchar(validation_script)) {
  validation_script <- testthat::test_path(
    "..", "..", "inst", "validation", "gaze-weave-comparison.R"
  )
}
source(validation_script, local = TRUE)

test_that("retrieval statistics give fractional credit under complete ties", {
  scores <- matrix(1, nrow = 4, ncol = 4)
  target <- seq_len(4)
  result <- validation_retrieval_statistics(scores, target)

  expect_equal(result[["pairwise_auc"]], 0.5)
  expect_equal(result[["top1_credit"]], 0.25)
  expect_equal(result[["reciprocal_rank"]], 0.4)
})

test_that("retrieval statistics treat roundoff-scale differences as ties", {
  scores <- matrix(1, nrow = 4, ncol = 4)
  scores[1, ] <- 1 + c(0, 1e-14, -1e-14, 2e-14)
  target <- seq_len(4)
  result <- validation_retrieval_statistics(scores, target)

  expect_equal(result[["pairwise_auc"]], 0.5)
  expect_equal(result[["top1_credit"]], 0.25)
  expect_equal(validation_row_statistics(scores, target), rep(0, 4))
})

test_that("duration-weighted density is invariant to fixation order", {
  coords <- rbind(c(20, 20), c(35, 70), c(75, 55), c(80, 15))
  forward <- validation_fixation_group(coords, duration = rep(1, 4))
  reverse <- validation_fixation_group(coords[4:1, ], duration = rep(1, 4))

  forward_signature <- validation_density_signature(forward)
  reverse_signature <- validation_density_signature(reverse)
  expect_equal(
    validation_density_similarity(forward_signature, reverse_signature),
    1,
    tolerance = 1e-12
  )
})

test_that("registered and raw baselines use the declared common warp", {
  reference <- validation_spatial_templates(4)
  scale <- 0.8
  translation <- c(5, -3)
  source <- lapply(reference, function(path) {
    coords <- sweep(cbind(path$x, path$y), 2, translation, FUN = "-") / scale
    validation_fixation_group(coords, path$duration)
  })
  spec <- gaze_transport_spec(
    spatial = gaze_gaussian_mixture(
      c(2.5, 5, 10), weights = c(0.45, 0.35, 0.2), unit = "px"
    ),
    warp = gaze_warp_contraction(
      center = c(50, 50), translation = TRUE
    ),
    screen = gaze_screen(100, 100, unit = "px"),
    backend = "reference",
    reliability = "none"
  )
  warp <- validation_fit_warp(reference, source, spec)
  registered <- lapply(source, validation_register_path, warp_model = warp)

  expect_equal(eyesim:::warp_parameters(warp)$scale, scale, tolerance = 1e-5)
  expect_equal(registered[[1]]$x, reference[[1]]$x, tolerance = 1e-5)
  expect_equal(registered[[1]]$y, reference[[1]]$y, tolerance = 1e-5)
})

test_that("participant statistic is zero when all candidate scores tie", {
  scores <- matrix(0.4, nrow = 5, ncol = 5)
  target <- seq_len(5)

  expect_equal(validation_row_statistics(scores, target), rep(0, 5))
  expect_equal(
    validation_row_statistics(scores, target, gaze = TRUE),
    rep(0, 5)
  )
})

test_that("conditional randomization is conservative for tied scores", {
  tied_case <- list(
    scores = list(
      gaze_weave = matrix(0, 4, 4),
      density_raw = matrix(1, 4, 4)
    ),
    target = seq_len(4)
  )
  operating <- validation_randomization_operating_characteristics(
    rep(list(tied_case), 4),
    sample_sizes = 4,
    repetitions = 20,
    randomization_draws = 39,
    signal_probabilities = 1,
    seed = 4
  )

  expect_equal(operating$estimated_power, c(0, 0))
  expect_equal(operating$estimated_type1, c(0, 0))
})
