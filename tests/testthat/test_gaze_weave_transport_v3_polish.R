make_transport_v3_polish_fixture <- function() {
  reference <- make_gaze_fixations(
    rbind(c(0, 0), c(1, 2), c(2.5, 0.8), c(4, 3)),
    duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  source <- make_gaze_fixations(
    rbind(c(0.1, 0), c(1.15, 2.1), c(2.4, 0.9), c(4.1, 2.9)),
    duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
  )
  list(reference = reference, source = source)
}

make_transport_v3_polish_spec <- function(backend = "optimized",
                                          polish = "audit") {
  gaze_transport_spec(
    coverage_nodes = 2,
    entropy_schedule = 0.03,
    maxit = 60,
    tolerance = 1e-3,
    projection_maxit = 300,
    projection_tolerance = 1e-7,
    backend = backend,
    polish = polish,
    polish_maxit = 30,
    polish_gap_tolerance = 1e-7,
    polish_relative_tolerance = 1e-10
  )
}

transport_v3_polish_barycentric <- function(coupling, reference_coords) {
  selected <- colSums(coupling)
  t(coupling) %*% reference_coords /
    pmax(selected, .Machine$double.xmin)
}

test_that("the exact linear oracle respects dominated marginals", {
  skip_if_not_installed("lpSolve")
  gradient <- matrix(c(1, 4, 3, 2), 2, 2)
  vertex <- eyesim:::transport_v3_linear_oracle(
    gradient,
    reference_mass = c(0.6, 0.4),
    source_mass = c(0.5, 0.5),
    coverage = 0.7
  )

  expect_equal(vertex, matrix(c(0.5, 0, 0, 0.2), 2, 2))
  expect_lte(max(eyesim:::transport_v3_feasibility(
    vertex, c(0.6, 0.4), c(0.5, 0.5), 0.7
  )), 1e-12)
  expect_equal(sum(gradient * vertex), 0.9, tolerance = 1e-12)
})

test_that("scientific polishing is feasible, monotone, and fully diagnosed", {
  skip_if_not_installed("lpSolve")
  fixture <- make_transport_v3_polish_fixture()
  entropic <- gaze_transport_align(
    fixture$reference,
    fixture$source,
    make_transport_v3_polish_spec(polish = "none")
  )
  polished <- gaze_transport_align(
    fixture$reference,
    fixture$source,
    make_transport_v3_polish_spec(polish = "audit")
  )
  diagnostic <- polished$diagnostics$polish

  expect_true(diagnostic$no_worse)
  expect_true(diagnostic$all_feasible)
  expect_lte(diagnostic$maximum_feasibility_error, 1e-8)
  expect_gte(polished$log_score, entropic$log_score - 1e-10)
  expect_named(
    diagnostic$trace,
    c(
      "iteration", "scientific", "gap", "step", "relative_improvement",
      "coverage_error", "reference_error", "source_error",
      "nonnegative_error"
    )
  )
  expect_true(all(diff(diagnostic$trace$scientific) <= 1e-10))
  expect_true(is.finite(diagnostic$map_gap))
  expect_match(polished$convergence$backend, "scientific_polish")
  expect_equal(
    polished$diagnostics$warp_scale,
    entropic$diagnostics$warp_scale,
    tolerance = 0
  )
  expect_equal(
    polished$diagnostics$warp_translation,
    entropic$diagnostics$warp_translation,
    tolerance = 0
  )
})

test_that("polishing is stable across independent structural starts", {
  skip_if_not_installed("lpSolve")
  fixture <- make_transport_v3_polish_fixture()
  spec <- make_transport_v3_polish_spec(backend = "reference")
  reference <- eyesim:::as_transport_v3_measure(
    fixture$reference, spec$chronology
  )
  source <- eyesim:::as_transport_v3_measure(fixture$source, spec$chronology)
  spatial_cost <- eyesim:::gaze_spatial_cost(
    reference$coords, source$coords, spec$spatial
  )
  coverage <- 0.5
  independent <- coverage * outer(reference$mass, source$mass)
  spatial_vertex <- eyesim:::transport_v3_linear_oracle(
    spatial_cost, reference$mass, source$mass, coverage
  )
  starts <- list(independent = independent, spatial = spatial_vertex)
  polished <- lapply(starts, function(start) {
    eyesim:::polish_transport_v3_mass(
      start, reference, source, spatial_cost, spec,
      maxit = 80, gap_tolerance = 1e-7, relative_tolerance = 1e-10
    )
  })
  energies <- vapply(polished, function(result) {
    result$objective$scientific
  }, numeric(1))
  barycentric <- lapply(polished, function(result) {
    transport_v3_polish_barycentric(result$coupling, reference$coords)
  })

  expect_lte(diff(range(energies)), 5e-4)
  expect_equal(
    sum(polished$independent$coupling),
    sum(polished$spatial$coupling),
    tolerance = 1e-12
  )
  expect_lte(max(polished$independent$feasibility), 1e-8)
  expect_lte(max(polished$spatial$feasibility), 1e-8)
  expect_lte(sqrt(mean((barycentric[[1]] - barycentric[[2]])^2)), 0.1)

  shifted_reference <- reference
  shifted_reference$coords <- sweep(
    shifted_reference$coords, 2, c(1.5, -0.8), FUN = "+"
  )
  candidate_energy <- function(candidate, start_name) {
    candidate_cost <- eyesim:::gaze_spatial_cost(
      candidate$coords, source$coords, spec$spatial
    )
    start <- if (identical(start_name, "independent")) {
      coverage * outer(candidate$mass, source$mass)
    } else {
      eyesim:::transport_v3_linear_oracle(
        candidate_cost, candidate$mass, source$mass, coverage
      )
    }
    eyesim:::polish_transport_v3_mass(
      start, candidate, source, candidate_cost, spec,
      maxit = 80, gap_tolerance = 1e-7, relative_tolerance = 1e-10
    )$objective$scientific
  }
  ranks <- lapply(names(starts), function(start_name) {
    sort(c(
      true = candidate_energy(reference, start_name),
      shifted = candidate_energy(shifted_reference, start_name)
    ), index.return = TRUE)$ix
  })
  expect_identical(ranks[[1]], ranks[[2]])
})

test_that("candidate order and backend do not change polished evidence", {
  skip_if_not_installed("lpSolve")
  skip_if_not(eyesim:::transport_v3_native_available())
  fixture <- make_transport_v3_polish_fixture()
  shifted <- make_gaze_fixations(
    cbind(fixture$reference$x + 1.5, fixture$reference$y - 0.8),
    duration = fixture$reference$duration,
    onset = fixture$reference$onset
  )
  references <- list(true = fixture$reference, shifted = shifted)
  optimized_spec <- make_transport_v3_polish_spec("optimized", "audit")
  reference_spec <- make_transport_v3_polish_spec("reference", "audit")
  forward <- gaze_transport_align_batch(
    references, fixture$source, optimized_spec, batch_size = 1
  )
  reverse <- gaze_transport_align_batch(
    rev(references), fixture$source, optimized_spec, batch_size = 2
  )
  oracle <- gaze_transport_align_batch(
    references, fixture$source, reference_spec, batch_size = 2
  )
  scores <- function(results) {
    stats::setNames(vapply(results, `[[`, numeric(1), "log_score"), names(results))
  }

  expect_equal(scores(forward), scores(reverse)[names(forward)], tolerance = 1e-10)
  expect_equal(scores(forward), scores(oracle), tolerance = 1e-7)
  expect_identical(
    names(sort(scores(forward), decreasing = TRUE)),
    names(sort(scores(oracle), decreasing = TRUE))
  )
  expect_true(all(vapply(forward, function(result) {
    result$diagnostics$polish$no_worse &&
      result$diagnostics$polish$all_feasible
  }, logical(1))))
})
