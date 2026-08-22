test_that("the engine contract separates estimands and fixes one endpoint", {
  contract <- eyesim:::gaze_engine_contract()

  expect_equal(contract$version, 2L)
  expect_identical(contract$primary_endpoint, "gaze_info_bits")
  expect_identical(
    contract$estimands,
    c("coverage", "selection", "correspondence", "evidence")
  )
  expect_named(contract$engines, c("transport", "replay"))
  expect_identical(contract$engines$transport$directionality, "symmetric")
  expect_identical(contract$engines$replay$directionality, "encoding_to_recall")
})

test_that("candidate engine results require scientific provenance", {
  provenance <- list(
    engine_version = 2L,
    directionality = "symmetric",
    score_semantics = "negative_unregularized_profile_energy",
    duration_semantics = "unit_duration_mass",
    candidate_invariant = TRUE
  )
  result <- eyesim:::new_gaze_engine_result(
    engine = "transport",
    candidate_key = "item-1",
    log_score = -2.5,
    diagnostics = list(replay_coverage = 0.7),
    provenance = provenance
  )

  expect_s3_class(result, "gaze_engine_result")
  expect_identical(result$engine, "transport")
  expect_error(
    eyesim:::new_gaze_engine_result(
      "transport", "item-1", -2.5,
      diagnostics = list(replay_coverage = 0.7),
      provenance = within(provenance, candidate_invariant <- FALSE)
    ),
    "candidate-invariant"
  )
  expect_error(
    eyesim:::new_gaze_engine_result(
      "replay", "item-1", -2.5,
      diagnostics = list(replay_coverage = 0.7),
      provenance = provenance
    ),
    "directionality"
  )
})

test_that("candidate result sets reject mixed engines and degenerate scores", {
  make_result <- function(key, score, engine = "replay") {
    directionality <- if (engine == "replay") "encoding_to_recall" else "symmetric"
    eyesim:::new_gaze_engine_result(
      engine = engine,
      candidate_key = key,
      log_score = score,
      diagnostics = list(replay_coverage = 0.5),
      provenance = list(
        engine_version = 2L,
        directionality = directionality,
        score_semantics = "candidate_log_likelihood",
        duration_semantics = "normalized_gaze_time_grid",
        candidate_invariant = TRUE
      )
    )
  }

  expect_true(eyesim:::validate_gaze_engine_results(list(
    make_result("a", -1),
    make_result("b", -2)
  )))
  expect_error(
    eyesim:::validate_gaze_engine_results(list(
      make_result("a", -Inf),
      make_result("b", -Inf)
    )),
    "at least one"
  )
  expect_error(
    eyesim:::validate_gaze_engine_results(list(
      make_result("a", -1),
      make_result("b", -2, engine = "transport")
    )),
    "same engine"
  )
})
