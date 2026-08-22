test_that("the common entry point exposes only Transport and Replay", {
  expect_named(
    eyesim:::gaze_engine_contract()$engines,
    c("transport", "replay")
  )
  expect_false(any(c("permute_on", "workers") %in%
                     names(formals(gaze_weave_cv))))

  expect_error(
    gaze_weave_cv(
      tibble::tibble(), tibble::tibble(), match_on = "item",
      spec = gaze_transport_spec()
    ),
    "matching, episode, prior, warp, and path columns"
  )
  expect_error(
    gaze_weave_cv(
      tibble::tibble(), tibble::tibble(), match_on = "item",
      engine = "transport", spec = make_small_gaze_replay_spec()
    ),
    "incompatible"
  )
  expect_error(
    gaze_weave_cv(
      tibble::tibble(), tibble::tibble(), match_on = "item",
      engine = "transport"
    ),
    "matching, episode, prior, warp, and path columns"
  )
})

test_that("the common entry point infers Replay from its specification", {
  tabs <- make_gaze_weave_cv_tables()
  replay <- gaze_weave_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = c("participant", "image_id"),
    contrast_on = "participant",
    n_folds = 2,
    spec = make_small_gaze_replay_spec(),
    seed = 29
  )

  expect_s3_class(replay, "gaze_replay_fit")
  expect_identical(replay$provenance$engine, "replay")
  expect_named(
    broom::tidy(replay),
    c("participant", "image_id", "gaze_info_bits")
  )
})

test_that("an engine specification is required when the engine is omitted", {
  expect_error(
    gaze_weave_cv(
      tibble::tibble(), tibble::tibble(), match_on = "item"
    ),
    "No GazeWeave engine is selected"
  )
})

test_that("Transport plots identify optimized correspondence", {
  reference <- make_gaze_fixations(rbind(c(0, 0), c(1, 1), c(2, 0)))
  source <- make_gaze_fixations(rbind(c(0.1, 0), c(1.1, 1), c(2.1, 0)))
  spec <- gaze_transport_spec(
    coverage_nodes = 2,
    entropy_schedule = 0.03,
    maxit = 40,
    tolerance = 1e-3,
    projection_maxit = 300,
    projection_tolerance = 1e-7,
    backend = "optimized",
    reliability = "none"
  )
  result <- gaze_transport_align(reference, source, spec)

  braid <- ggplot2::autoplot(result$alignment, type = "braid")
  diagnostics <- ggplot2::autoplot(result$alignment, type = "diagnostics")

  expect_s3_class(braid, "ggplot")
  expect_match(braid$labels$subtitle, "not a posterior")
  expect_match(diagnostics$labels$subtitle, "not a posterior")
})
