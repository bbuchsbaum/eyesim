# Pre-CRAN stabilization of the GazeWeave interface.

test_that("fold assignment does not depend on RNGkind or collation", {
  tab <- data.frame(item = c("a_1", "B_2", "c-3", "D.4", "e5", "F6", "g_7", "H8"),
                    stratum = "s")
  reference <- make_gaze_weave_folds(tab, "item", "stratum", 2, seed = 3)$fold_id

  old_kind <- RNGkind()
  on.exit(suppressWarnings(RNGkind(old_kind[[1]], old_kind[[2]], old_kind[[3]])),
          add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(99)
  stream <- .Random.seed
  expect_identical(
    make_gaze_weave_folds(tab, "item", "stratum", 2, seed = 3)$fold_id, reference
  )
  # The caller's generator and stream are restored.
  expect_identical(RNGkind()[[1]], "L'Ecuyer-CMRG")
  expect_identical(.Random.seed, stream)

  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
  expect_identical(
    make_gaze_weave_folds(tab, "item", "stratum", 2, seed = 3)$fold_id, reference
  )
  expect_identical(RNGkind()[[3]], "Rounding")
  RNGkind("default", "default", "default")

  skip_if_not_installed("withr")
  expect_identical(
    withr::with_collate("C", make_gaze_weave_folds(tab, "item", "stratum", 2,
                                                   seed = 3)$fold_id),
    reference
  )
})

test_that("a missing seed is an error, not a random fold assignment", {
  tab <- data.frame(item = letters[1:6], stratum = "s")
  expect_error(make_gaze_weave_folds(tab, "item", "stratum", 2, seed = NULL),
               "seed must be a single number")
})

test_that("Replay warns when its calibration falls back to the prior", {
  tabs <- make_gaze_weave_cv_tables()
  expect_warning(
    fit <- suppressMessages(gaze_replay_cv(
      tabs$ref_tab, tabs$source_tab, match_on = c("participant", "image_id"),
      contrast_on = "participant", n_folds = 2
    )),
    class = "gaze_calibration_fallback"
  )
  expect_identical(fit$calibration_fallback_folds, 2L)
  expect_output(print(fit), "calibration fallback: 2 of 2 folds")
})

test_that("spatial and screen units must agree", {
  expect_error(
    gaze_transport_spec(spatial = gaze_gaussian_mixture(1, unit = "deg"),
                        screen = gaze_screen(800, 600)),
    "spatial model is in 'deg' but the screen is in 'px'"
  )
  expect_s3_class(
    gaze_transport_spec(spatial = gaze_gaussian_mixture(20, unit = "px"),
                        screen = gaze_screen(800, 600)),
    "gaze_transport_spec"
  )
})

test_that("pixel-scale paths with the degree-scale default kernel warn", {
  tabs <- make_gaze_weave_cv_tables()
  scale_paths <- function(paths) lapply(paths, function(path) {
    path$x <- path$x * 100
    path$y <- path$y * 100
    path
  })
  default_spatial <- gaze_transport_spec()$spatial
  expect_warning(
    check_gaze_spatial_scale(scale_paths(tabs$ref_tab$fixgroup),
                             scale_paths(tabs$source_tab$fixgroup),
                             default_spatial),
    class = "gaze_spatial_scale"
  )
  expect_silent(check_gaze_spatial_scale(tabs$ref_tab$fixgroup,
                                         tabs$source_tab$fixgroup,
                                         default_spatial))
})

test_that("a screen centre sets the screen limits", {
  centred <- gaze_screen(30, 20, unit = "deg", center = c(0, 0))
  expect_equal(centred$xlim, c(-15, 15))
  expect_equal(centred$ylim, c(-10, 10))
  corner <- gaze_screen(1024, 768)
  expect_equal(corner$xlim, c(0, 1024))
  expect_equal(corner$ylim, c(0, 768))
})

test_that("Replay warns when many fixations fall off the declared screen", {
  tabs <- make_gaze_weave_cv_tables()
  expect_warning(
    check_gaze_screen_clamping(tabs$ref_tab$fixgroup, tabs$source_tab$fixgroup,
                               gaze_screen(0.1, 0.1)),
    class = "gaze_screen_clamping"
  )
})

test_that("missing columns are named", {
  tabs <- make_gaze_weave_cv_tables()
  expect_error(
    gaze_transport_cv(tabs$ref_tab, tabs$source_tab,
                      match_on = c("participant", "image_idX")),
    "ref_tab \\(image_idX\\); source_tab \\(image_idX\\)"
  )
})

test_that("gaze_weave_cv routes template_on to Replay only", {
  tabs <- make_gaze_weave_cv_tables()
  expect_error(
    gaze_weave_cv(tabs$ref_tab, tabs$source_tab,
                  match_on = c("participant", "image_id"),
                  spec = gaze_transport_spec(), template_on = "image_id"),
    "template_on is available only for Replay"
  )
})

test_that("internal plumbing is not exported", {
  exports <- getNamespaceExports("eyesim")
  expect_false(any(c("gaze_local_order", "gaze_transport_prepare",
                     "gaze_transport_align_batch", "elastic_consensus_align")
                   %in% exports))
})
