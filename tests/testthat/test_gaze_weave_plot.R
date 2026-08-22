test_that("Replay alignment plots identify posterior correspondence", {
  tabs <- make_gaze_weave_cv_tables()
  keep <- tabs$ref_tab$participant == "p1"
  training_ref <- tabs$ref_tab[keep, ]
  training_source <- tabs$source_tab[
    tabs$source_tab$participant == "p1", , drop = FALSE
  ]
  model <- fit_gaze_replay_model(
    training_ref,
    training_source,
    match_on = c("participant", "image_id"),
    spec = make_small_gaze_replay_spec()
  )
  result <- gaze_replay_align(
    training_ref$fixgroup[[1L]],
    training_source$fixgroup[[1L]],
    model
  )

  overlay <- ggplot2::autoplot(result$alignment, type = "overlay")
  braid <- ggplot2::autoplot(result$alignment, type = "braid")
  diagnostics <- ggplot2::autoplot(result$alignment, type = "diagnostics")

  expect_s3_class(overlay, "ggplot")
  expect_match(overlay$labels$subtitle, "posterior replay")
  expect_match(braid$labels$title, "Posterior")
  expect_match(diagnostics$labels$title, "posterior")
})

test_that("fitted-row selection requires one exact row", {
  fit <- list(results = data.frame(
    participant = c("p1", "p1", "p2"),
    item = c(1L, 2L, 1L)
  ))

  expect_identical(
    eyesim:::select_gaze_weave_result(
      fit, key = list(participant = "p1", item = 2L)
    ),
    2L
  )
  expect_error(
    eyesim:::select_gaze_weave_result(
      fit, key = list(participant = "p1")
    ),
    "exactly one"
  )
})
