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

  old_collate <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old_collate), add = TRUE)
  Sys.setlocale("LC_COLLATE", "C")
  expect_identical(
    make_gaze_weave_folds(tab, "item", "stratum", 2, seed = 3)$fold_id, reference
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

test_that("invalid fixation paths are reported by table, row and key", {
  tabs <- make_gaze_weave_cv_tables()
  bad <- tabs$source_tab
  bad$fixgroup[[3]]$x[[2]] <- NA
  bad$fixgroup[[5]]$onset[[2]] <- bad$fixgroup[[5]]$onset[[1]]
  match_on <- c("participant", "image_id")
  for (engine in list(
    function() gaze_transport_cv(tabs$ref_tab, bad, match_on = match_on),
    function() gaze_replay_cv(tabs$ref_tab, bad, match_on = match_on),
    function() gaze_baseline_cv(
      tabs$ref_tab, bad, match_on = match_on,
      spec = gaze_baseline_spec(screen = gaze_screen(20, 10, unit = "deg"),
                                density_sigmas = 1, methods = "density")
    )
  )) {
    error <- tryCatch(engine(), error = identity)
    expect_s3_class(error, "error")
    message <- conditionMessage(error)
    expect_match(message, "2 invalid fixation paths")
    expect_match(message, paste0(
      "source_tab row 3 \\(participant = ", bad$participant[[3]],
      ", image_id = ", bad$image_id[[3]], "\\), fixgroup: ",
      "GazeWeave inputs must contain only finite values"
    ))
    expect_match(message, "source_tab row 5 .*onsets must be unique")
  }

  bad_ref <- tabs$ref_tab
  bad_ref$fixgroup[[2]] <- bad_ref$fixgroup[[2]][0, ]
  expect_error(
    gaze_replay_cv(bad_ref, tabs$source_tab, match_on = match_on),
    "ref_tab row 2 .*at least one fixation"
  )
})

test_that("path validation caps the report and accepts valid paths", {
  tabs <- make_gaze_weave_cv_tables()
  bad <- tabs$source_tab
  for (row in 1:7) bad$fixgroup[[row]]$duration[] <- 0
  error <- tryCatch(
    stop_if_invalid_gaze_paths("Transport", bad, "fixgroup", "source_tab"),
    error = identity
  )
  expect_match(conditionMessage(error), "7 invalid fixation paths")
  expect_match(conditionMessage(error), "and 2 more")
  expect_null(stop_if_invalid_gaze_paths(
    "Transport", tabs$source_tab, "fixgroup", "source_tab"
  ))
  expect_identical(gaze_path_problem(list(x = 1)),
                   "GazeWeave inputs must be fixation_group objects or data frames.")
  expect_null(gaze_path_problem(tabs$source_tab$fixgroup[[1]]))
})

test_that("every engine returns the shared result columns", {
  tabs <- make_gaze_weave_cv_tables()
  match_on <- c("participant", "image_id")
  fits <- suppressWarnings(suppressMessages(list(
    transport = gaze_transport_cv(tabs$ref_tab, tabs$source_tab, match_on,
                                  contrast_on = "participant", n_folds = 2),
    replay = gaze_replay_cv(tabs$ref_tab, tabs$source_tab, match_on,
                            contrast_on = "participant", n_folds = 2),
    baseline = gaze_baseline_cv(
      tabs$ref_tab, tabs$source_tab, match_on, contrast_on = "participant",
      n_folds = 2,
      spec = gaze_baseline_spec(screen = gaze_screen(20, 10, unit = "deg"),
                                density_sigmas = 1, methods = "density")
    )
  )))
  for (engine in names(fits)) {
    fit <- fits[[engine]]
    results <- fit$results
    info <- paste("engine:", engine)
    expect_true(all(gaze_weave_core_columns %in% names(results)), info = info)
    expect_true(all(match_on %in% names(results)), info = info)
    for (column in setdiff(gaze_weave_core_columns, "candidates")) {
      expect_true(is.numeric(results[[column]]), info = paste(info, column))
    }
    expect_true(is.list(results$candidates), info = info)
    # Uncalibrated baseline methods have no calibrated probability.
    calibrated <- !is.na(results$gaze_info_bits)
    expect_true(any(calibrated) || engine == "baseline", info = info)
    expect_equal(results$gaze_info_bits[calibrated],
                 log2(results$posterior_true[calibrated] /
                        results$prior_true[calibrated]),
                 tolerance = 1e-10, info = info)
    expect_equal(results$prior_true, 1 / results$candidate_count,
                 tolerance = 1e-12, info = info)
    scored <- !is.na(results$template_rank)
    expect_true(any(scored), info = info)
    for (table in results$candidates[scored]) {
      expect_true(all(gaze_weave_core_candidate_columns %in% names(table)),
                  info = info)
    }
    expect_true(all(c("spec", "keys", "folds", "provenance") %in% names(fit)),
                info = info)
    expect_identical(fit$provenance$engine, c(
      transport = "transport", replay = "replay",
      baseline = "fair_baseline_court"
    )[[engine]], info = info)
  }
  expect_identical(fits$transport$provenance$revision, "2026.10")
  expect_identical(fits$replay$provenance$revision, "2026.10")
})

test_that("Transport-only arguments are refused for Replay", {
  tabs <- make_gaze_weave_cv_tables()
  expect_error(
    gaze_weave_cv(tabs$ref_tab, tabs$source_tab,
                  match_on = c("participant", "image_id"),
                  spec = gaze_replay_spec(), episode_on = "image_id"),
    "episode_on and priorvar are available only for Transport"
  )
})

test_that("input checks tolerate empty or tiny coordinate sets", {
  expect_identical(dim(gaze_path_coordinates(list())), c(0L, 2L))
  expect_identical(dim(gaze_path_coordinates(list(list(a = 1)))), c(0L, 2L))
  two <- list(data.frame(x = c(1, 2), y = c(1, 2)))
  expect_silent(check_gaze_spatial_scale(two, list(),
                                         gaze_transport_spec()$spatial))
  expect_silent(check_gaze_screen_clamping(list(), list(),
                                           gaze_screen(10, 10)))
})

test_that("print() on a Transport fit reports calibration fallback folds", {
  tabs <- make_gaze_weave_cv_tables()
  fit <- suppressWarnings(gaze_transport_cv(
    tabs$ref_tab, tabs$source_tab, match_on = c("participant", "image_id"),
    contrast_on = "participant", n_folds = 2,
    # Five inner folds need ten items per stratum; each has three.
    spec = gaze_transport_spec(calibration_folds = 5)
  ))
  expect_gt(fit$calibration_fallback_folds, 0L)
  expect_output(print(fit), "calibration fallback")
})

test_that("the optimized backend refuses a specification it cannot solve", {
  tabs <- make_gaze_weave_cv_tables()
  spec <- gaze_transport_spec(backend = "optimized", projection_method = "log")
  expect_error(
    gaze_transport_align(tabs$ref_tab$fixgroup[[1]],
                         tabs$source_tab$fixgroup[[1]], spec),
    "does not replace the log-domain projection"
  )
})

test_that("a prepared path must match the spec's whole chronology", {
  tabs <- make_gaze_weave_cv_tables()
  spec <- gaze_transport_spec(coverage_nodes = 2)
  coalescing <- gaze_transport_spec(
    coverage_nodes = 2,
    chronology = gaze_order_neighbours(2, coalesce_distance = 1)
  )
  source <- tabs$source_tab$fixgroup[[1]]
  references <- list(a = tabs$ref_tab$fixgroup[[1]],
                     b = tabs$ref_tab$fixgroup[[3]])
  prepared <- gaze_transport_prepare(source, spec)
  expect_error(gaze_transport_align_batch(references, prepared, coalescing),
               "different chronology")
  expect_error(gaze_transport_align_batch(
    lapply(references, gaze_transport_prepare, spec = spec), source,
    coalescing
  ), "different chronology")

  batch <- gaze_transport_align_batch(references, prepared, spec)
  expect_identical(names(batch), c("a", "b"))
  for (key in names(references)) {
    single <- gaze_transport_align(references[[key]], source, spec,
                                   candidate_key = key)
    expect_equal(batch[[key]]$log_score, single$log_score, tolerance = 1e-12)
  }
})

test_that("a screen-centred warp uses the screen centre", {
  tabs <- make_gaze_weave_cv_tables()
  screen <- gaze_screen(20, 10, unit = "deg")
  fit_with <- function(center) {
    fit <- suppressWarnings(gaze_transport_cv(
      tabs$ref_tab, tabs$source_tab, match_on = c("participant", "image_id"),
      contrast_on = "participant", n_folds = 2,
      spec = gaze_transport_spec(
        coverage_nodes = 2, screen = screen,
        warp = gaze_warp_contraction(center = center, translation = TRUE)
      )
    ))
    fit$results$gaze_info_bits
  }
  expect_identical(fit_with("screen"), fit_with(screen$center))
  expect_error(resolve_gaze_warp_center(gaze_warp_contraction(), NULL),
               "requires screen geometry")
})
