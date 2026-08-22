make_transport_v3_coverage_spec <- function(coverage_nodes = 2L,
                                            selection_weights = c(0.5, 0.5)) {
  gaze_transport_spec(
    coverage_nodes = coverage_nodes,
    selection_weights = selection_weights,
    entropy_schedule = 0.03,
    maxit = 50,
    tolerance = 1e-3,
    projection_maxit = 300,
    projection_tolerance = 1e-7,
    multistart = 1,
    backend = "reference"
  )
}

test_that("default coverage quadrature freezes the declared low-mass support", {
  spec <- gaze_transport_spec()
  expect_equal(spec$coverage$nodes, 12L)
  expect_equal(sum(spec$coverage$weight), 1, tolerance = 1e-14)
  expect_true(any(abs(spec$coverage$coverage - 0.206341023) < 1e-8))
  expect_true(any(spec$coverage$coverage < 0.25))
  expect_equal(spec$coverage_prior, c(2, 2))
  expect_error(
    gaze_transport_spec(selection_weights = c(0.2, 0.3)),
    "Symmetric"
  )
})

test_that("coverage integration is stable under quadrature refinement", {
  energy <- function(coverage, candidate) {
    if (candidate == 1L) {
      0.35 * coverage + 0.6 * (coverage - 0.32)^2
    } else {
      0.52 * coverage + 0.4 * (coverage - 0.68)^2
    }
  }
  coarse <- eyesim:::transport_v3_coverage_quadrature(12, c(2, 2))
  refined <- eyesim:::transport_v3_coverage_quadrature(24, c(2, 2))
  coarse_scores <- vapply(1:2, function(candidate) {
    eyesim:::transport_v3_integrate_coverage(
      energy(coarse$coverage, candidate), coarse
    )$log_score
  }, numeric(1))
  refined_scores <- vapply(1:2, function(candidate) {
    eyesim:::transport_v3_integrate_coverage(
      energy(refined$coverage, candidate), refined
    )$log_score
  }, numeric(1))
  coarse_evidence <- eyesim:::score_gaze_candidates(
    coarse_scores, true_index = 1, prior = c(0.5, 0.5)
  )
  refined_evidence <- eyesim:::score_gaze_candidates(
    refined_scores, true_index = 1, prior = c(0.5, 0.5)
  )

  expect_lte(max(abs(coarse_scores - refined_scores)), 1e-3)
  expect_lte(
    abs(coarse_evidence$gaze_info_bits - refined_evidence$gaze_info_bits),
    0.01
  )
})

test_that("actual pair and candidate evidence pass coverage refinement gates", {
  coords <- rbind(c(0, 0), c(1, 1.5), c(2.5, 0.5), c(3.5, 2))
  true_reference <- make_gaze_fixations(
    coords, duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  wrong_reference <- make_gaze_fixations(
    sweep(coords, 2, c(2.5, -1.5), FUN = "+"),
    duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  source <- make_gaze_fixations(
    sweep(coords, 2, c(0.1, -0.05), FUN = "+"),
    duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  make_spec <- function(nodes) {
    gaze_transport_spec(
      coverage_nodes = nodes,
      entropy_schedule = 0.03,
      maxit = 100,
      tolerance = 1e-4,
      projection_maxit = 300,
      projection_tolerance = 1e-7,
      backend = "reference"
    )
  }
  score_candidates <- function(spec) {
    c(
      gaze_transport_align(
        true_reference, source, spec, candidate_key = "true"
      )$log_score,
      gaze_transport_align(
        wrong_reference, source, spec, candidate_key = "wrong"
      )$log_score
    )
  }
  coarse <- score_candidates(make_spec(12))
  refined <- score_candidates(make_spec(24))
  coarse_evidence <- eyesim:::score_gaze_candidates(
    coarse, true_index = 1, prior = c(0.5, 0.5)
  )
  refined_evidence <- eyesim:::score_gaze_candidates(
    refined, true_index = 1, prior = c(0.5, 0.5)
  )

  expect_lte(max(abs(coarse - refined)), 1e-3)
  expect_lte(
    abs(coarse_evidence$gaze_info_bits - refined_evidence$gaze_info_bits),
    0.01
  )
})

test_that("Transport objective separates every scientific object", {
  set.seed(20260820)
  spec <- make_transport_v3_coverage_spec()
  reference_path <- make_gaze_fixations(
    rbind(c(0, 0), c(1, 2), c(3, 1), c(4, 3)),
    duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  source_path <- make_gaze_fixations(
    rbind(c(0.1, 0), c(1.2, 2), c(3.2, 1), c(4.1, 3)),
    duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  reference <- eyesim:::as_transport_v3_measure(
    reference_path, spec$chronology
  )
  source <- eyesim:::as_transport_v3_measure(source_path, spec$chronology)
  spatial_cost <- eyesim:::gaze_spatial_cost(
    reference$coords, source$coords, spec$spatial
  )
  correspondence <- matrix(stats::runif(16, 0.2, 1), 4, 4)
  correspondence <- correspondence / sum(correspondence)
  coverage <- 0.35
  coupling <- coverage * correspondence
  objective <- eyesim:::transport_v3_objective(
    coupling, reference, source, spatial_cost, spec, entropy = 0.03,
    warp_penalty = 0.07, gradient = TRUE
  )

  expect_equal(objective$coverage, coverage, tolerance = 1e-14)
  expect_equal(sum(objective$correspondence), 1, tolerance = 1e-14)
  expect_named(
    objective$components,
    c(
      "spatial", "chronology", "reference_selection", "source_selection",
      "warp_penalty", "correspondence_smoothing"
    )
  )
  expect_named(
    objective$conditional,
    c(
      "spatial", "chronology", "reference_selection", "source_selection",
      "correspondence_information"
    )
  )
  expect_equal(sum(objective$selected_reference_mass), coverage)
  expect_equal(sum(objective$selected_source_mass), coverage)
  expect_equal(
    objective$scientific,
    sum(objective$components[names(objective$components) !=
      "correspondence_smoothing"]),
    tolerance = 1e-14
  )
  expect_equal(
    objective$optimization,
    sum(objective$components),
    tolerance = 1e-14
  )
})

test_that("Transport fixed-mass objective has an analytic gradient", {
  set.seed(20260820)
  spec <- make_transport_v3_coverage_spec()
  path <- make_gaze_fixations(
    rbind(c(0, 0), c(1, 2), c(3, 1), c(4, 3)),
    duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  measure <- eyesim:::as_transport_v3_measure(path, spec$chronology)
  spatial_cost <- eyesim:::gaze_spatial_cost(
    measure$coords, measure$coords, spec$spatial
  )
  correspondence <- matrix(stats::runif(16, 0.5, 1), 4, 4)
  correspondence <- correspondence / sum(correspondence)
  coupling <- 0.3 * correspondence
  direction <- matrix(stats::rnorm(16), 4, 4)
  direction <- direction - mean(direction)
  direction <- direction / max(abs(direction))
  analytic <- eyesim:::transport_v3_objective(
    coupling, measure, measure, spatial_cost, spec, entropy = 0.03,
    gradient = TRUE
  )$gradient
  objective <- function(value) {
    eyesim:::transport_v3_objective(
      value, measure, measure, spatial_cost, spec, entropy = 0.03
    )$optimization
  }
  step <- 1e-7
  numerical <- (
    objective(coupling + step * direction) -
      objective(coupling - step * direction)
  ) / (2 * step)

  expect_equal(numerical, sum(analytic * direction), tolerance = 2e-7)
})

test_that("zero coverage is neutral and has no conditional singularity", {
  spec <- make_transport_v3_coverage_spec()
  path <- make_gaze_fixations(rbind(c(0, 0), c(1, 1)))
  measure <- eyesim:::as_transport_v3_measure(path, spec$chronology)
  zero <- matrix(0, 2, 2)
  result <- eyesim:::transport_v3_objective(
    zero,
    measure,
    measure,
    matrix(0, 2, 2),
    spec,
    entropy = 0.03,
    gradient = TRUE
  )

  expect_equal(result$coverage, 0)
  expect_equal(result$scientific, 0)
  expect_true(result$zero_coverage)
  expect_true(all(is.na(result$conditional)))
  expect_true(all(result$selected_reference_mass == 0))
  expect_true(all(is.finite(result$gradient)))
})

test_that("selecting one convenient pair pays both selection costs", {
  spec <- make_transport_v3_coverage_spec()
  path <- make_gaze_fixations(
    rbind(c(0, 0), c(2, 0), c(4, 0), c(6, 0))
  )
  measure <- eyesim:::as_transport_v3_measure(path, spec$chronology)
  coupling <- matrix(0, 4, 4)
  coupling[1, 1] <- 0.2
  objective <- eyesim:::transport_v3_objective(
    coupling,
    measure,
    measure,
    eyesim:::gaze_spatial_cost(
      measure$coords, measure$coords, spec$spatial
    ),
    spec,
    entropy = 0.03
  )

  expect_gt(objective$components[["reference_selection"]], 0)
  expect_gt(objective$components[["source_selection"]], 0)
  expect_gt(
    objective$scientific,
    objective$components[["spatial"]] +
      objective$components[["chronology"]]
  )
  expect_true(all(objective$selected_reference_mass <= measure$mass + 1e-14))
  expect_true(all(objective$selected_source_mass <= measure$mass + 1e-14))
})

test_that("reference alignment exports coverage, selection, and residuals", {
  spec <- make_transport_v3_coverage_spec()
  coords <- rbind(c(0, 0), c(1, 1.5), c(2.5, 0.5), c(3.5, 2))
  reference <- make_gaze_fixations(
    coords, duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  source <- make_gaze_fixations(
    sweep(coords, 2, c(0.1, -0.05), FUN = "+"),
    duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  split_source <- make_gaze_fixations(
    as.matrix(rbind(
      source[1, c("x", "y")],
      source[1, c("x", "y")],
      source[-1, c("x", "y")]
    )),
    duration = c(0.25, 0.75, 2, 1, 2),
    onset = c(0, 0.25, 1, 3, 4)
  )
  first <- gaze_transport_align(reference, source, spec, candidate_key = "a")
  split <- gaze_transport_align(
    reference, split_source, spec, candidate_key = "a"
  )

  expect_s3_class(first, "gaze_engine_result")
  expect_identical(first$engine, "transport")
  expect_equal(first$log_score, split$log_score, tolerance = 1e-8)
  expect_equal(
    first$diagnostics$matched_coverage,
    split$diagnostics$matched_coverage,
    tolerance = 1e-8
  )
  expect_equal(
    first$diagnostics$selected_reference_mass,
    split$diagnostics$selected_reference_mass,
    tolerance = 1e-8
  )
  expect_equal(
    first$diagnostics$selected_source_mass,
    split$diagnostics$selected_source_mass,
    tolerance = 1e-8
  )
  expect_named(
    first$diagnostics,
    c(
      "matched_coverage", "map_coverage", "selected_reference_mass",
      "selected_source_mass", "conditional_spatial_residual", "spatial_rmse",
      "conditional_chronology_residual", "local_order_preservation",
      "reference_selection_residual", "source_selection_residual",
      "correspondence_smoothing", "correspondence_information", "warp_penalty",
      "warp_scale", "warp_translation", "reference_dominance_error",
      "source_dominance_error", "start_spread", "start_scientific_spread",
      "coupling_change", "stationarity", "relative_objective_change"
    )
  )
  expect_true(first$convergence$converged)
})

test_that("one- and two-fixation pair alignments remain defined", {
  spec <- make_transport_v3_coverage_spec()
  one <- make_gaze_fixations(matrix(c(0, 0), 1, 2), duration = 1, onset = 0)
  two <- make_gaze_fixations(
    rbind(c(0, 0), c(1, 1)), duration = c(1, 1), onset = c(0, 1)
  )
  two_reversed <- make_gaze_fixations(
    rbind(c(1, 1), c(0, 0)), duration = c(1, 1), onset = c(0, 1)
  )
  one_fit <- gaze_transport_align(one, one, spec)
  two_fit <- gaze_transport_align(two, two_reversed, spec)

  expect_true(is.finite(one_fit$log_score))
  expect_true(is.finite(two_fit$log_score))
  expect_true(is.finite(one_fit$diagnostics$conditional_spatial_residual))
  expect_true(is.finite(two_fit$diagnostics$conditional_chronology_residual))
})
