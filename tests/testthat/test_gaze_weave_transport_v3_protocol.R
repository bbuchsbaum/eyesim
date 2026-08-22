transport_v3_validation_path <- function(...) {
  gaze_weave_test_inst_path("validation", ...)
}

test_that("the frozen development protocol preserves the blinded estimand", {
  protocol_path <- transport_v3_validation_path(
    "GAZEWEAVE-TRANSPORT-V3-PROTOCOL.md"
  )
  decision_path <- transport_v3_validation_path(
    "GAZEWEAVE-TRANSPORT-V3-DECISIONS.md"
  )
  manifest_path <- transport_v3_validation_path(
    "gaze-weave-transport-v3-manifest.json"
  )
  expect_true(file.exists(protocol_path))
  expect_true(file.exists(decision_path))
  expect_true(file.exists(manifest_path))

  protocol <- paste(readLines(protocol_path, warn = FALSE), collapse = "\n")
  decisions <- paste(readLines(decision_path, warn = FALSE), collapse = "\n")
  manifest <- paste(readLines(manifest_path, warn = FALSE), collapse = "\n")
  protocol_normalized <- gsub("[[:space:]]+", " ", protocol)

  required_protocol_terms <- c(
    "sole default inferential endpoint",
    "gaze_info_bits",
    "Pi = M Gamma",
    "selected reference and source mass",
    "conditional spatial residual",
    "conditional chronology",
    "correspondence smoothing",
    "warp penalty",
    "alignment ambiguity induced partly by correspondence regularization",
    "all retrieval subgroups remain sealed",
    "identical outer folds, candidate pools, priors, registration opportunities"
  )
  for (term in required_protocol_terms) {
    expect_match(protocol_normalized, term, fixed = TRUE)
  }
  expect_match(decisions, "reversal `>= 0.95`", fixed = TRUE)
  expect_match(decisions, "Median `<= 15` ms", fixed = TRUE)
  expect_match(manifest, '"version": "3.0.2"', fixed = TRUE)
  expect_match(manifest, '"primary_endpoint": "gaze_info_bits"', fixed = TRUE)
  expect_match(
    manifest,
    '"coupling_semantics": "optimized_correspondence_not_posterior"',
    fixed = TRUE
  )
  expect_match(manifest, '"reversal_length_absolute_band_max": 0.05', fixed = TRUE)
  expect_match(manifest, '"coverage_refinement_energy_max": 1e-3', fixed = TRUE)
  expect_match(manifest, '"native_reference_energy_max": 1e-8', fixed = TRUE)
  expect_match(manifest, '"pair_speedup_min": 4', fixed = TRUE)
})

test_that("estimator-lock files match their frozen hashes", {
  fixture_dir <- transport_v3_validation_path(
    "gaze-weave-transport-v3-golden"
  )
  manifest <- utils::read.csv(
    file.path(fixture_dir, "manifest-md5.csv"),
    stringsAsFactors = FALSE
  )
  paths <- file.path(fixture_dir, manifest$file)
  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), manifest$md5)

  summary <- utils::read.csv(
    file.path(fixture_dir, "expected-summary.csv"),
    stringsAsFactors = FALSE
  )
  expect_setequal(
    summary$fixture_id,
    c(
      "identical", "reversed", "local_swap", "offset",
      "partial_with_intrusions", "split_adjacent", "two_fixation_reversal"
    )
  )
  expect_true(all(summary$converged))
  expect_equal(
    summary[summary$fixture_id == "identical", "log_score"],
    summary[summary$fixture_id == "split_adjacent", "log_score"],
    tolerance = 1e-14
  )
})
