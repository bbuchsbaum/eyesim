make_transport_v3_backend_spec <- function(backend, spatial = NULL) {
  if (is.null(spatial)) {
    spatial <- gaze_gaussian_mixture(c(0.75, 1.5), weights = c(0.7, 0.3))
  }
  gaze_transport_spec(
    spatial = spatial,
    coverage_nodes = 2,
    entropy_schedule = 0.03,
    maxit = 50,
    tolerance = 1e-3,
    projection_maxit = 300,
    projection_tolerance = 1e-7,
    multistart = 1,
    backend = backend
  )
}

transport_v3_backend_barycentric <- function(result) {
  coupling <- result$alignment$coupling
  selected <- colSums(coupling)
  t(coupling) %*% result$alignment$reference$coords /
    pmax(selected, .Machine$double.xmin)
}

expect_transport_v3_backend_parity <- function(reference, source,
                                               spatial = NULL) {
  reference_fit <- gaze_transport_align(
    reference,
    source,
    make_transport_v3_backend_spec("reference", spatial),
    candidate_key = "fixture"
  )
  native_fit <- gaze_transport_align(
    reference,
    source,
    make_transport_v3_backend_spec("optimized", spatial),
    candidate_key = "fixture"
  )
  expect_equal(native_fit$log_score, reference_fit$log_score, tolerance = 1e-8)
  expect_equal(
    native_fit$alignment$profile$scientific_energy,
    reference_fit$alignment$profile$scientific_energy,
    tolerance = 1e-8
  )
  expect_equal(
    sum(native_fit$alignment$coupling),
    sum(reference_fit$alignment$coupling),
    tolerance = 1e-10
  )
  expect_equal(
    transport_v3_backend_barycentric(native_fit),
    transport_v3_backend_barycentric(reference_fit),
    tolerance = 1e-6
  )
  expect_equal(
    native_fit$diagnostics$warp_scale,
    reference_fit$diagnostics$warp_scale,
    tolerance = 1e-10
  )
  expect_equal(
    native_fit$diagnostics$warp_translation,
    reference_fit$diagnostics$warp_translation,
    tolerance = 1e-10
  )
  expect_true(native_fit$convergence$converged)
  expect_true(reference_fit$convergence$converged)
  invisible(list(reference = reference_fit, native = native_fit))
}

test_that("native backend matches the R oracle across fixture classes", {
  skip_if_not(eyesim:::transport_v3_native_available())
  set.seed(20260824)
  random_coords <- cbind(
    cumsum(stats::rnorm(6)),
    cumsum(stats::rnorm(6))
  )
  fixtures <- list(
    random = list(
      make_gaze_fixations(random_coords, duration = 1:6),
      make_gaze_fixations(random_coords + 0.08, duration = 6:1)
    ),
    sparse_duration = list(
      make_gaze_fixations(
        random_coords,
        duration = c(1e-4, 1, 1e-3, 4, 2e-4, 2)
      ),
      make_gaze_fixations(
        random_coords[c(1, 3, 2, 4, 6, 5), , drop = FALSE],
        duration = c(1e-4, 1e-3, 1, 4, 2, 2e-4)
      )
    ),
    one_fixation = list(
      make_gaze_fixations(matrix(c(0, 0), 1, 2)),
      make_gaze_fixations(matrix(c(0.1, -0.1), 1, 2))
    ),
    two_fixation = list(
      make_gaze_fixations(rbind(c(0, 0), c(1, 1))),
      make_gaze_fixations(rbind(c(1, 1), c(0, 0)))
    )
  )
  for (fixture in fixtures) {
    expect_transport_v3_backend_parity(fixture[[1]], fixture[[2]])
  }

  extreme_spatial <- gaze_gaussian_mixture(
    c(1e5, 2e5), weights = c(0.7, 0.3)
  )
  extreme <- random_coords * 1e6
  expect_transport_v3_backend_parity(
    make_gaze_fixations(extreme, duration = 1:6),
    make_gaze_fixations(extreme + 1e4, duration = 6:1),
    spatial = extreme_spatial
  )
})

test_that("native backend matches the public Transport-v2 path bank", {
  skip_if_not(eyesim:::transport_v3_native_available())
  input_path <- gaze_weave_test_inst_path(
    "validation", "gaze-weave-transport-v3-golden", "inputs.csv"
  )
  inputs <- utils::read.csv(input_path, stringsAsFactors = FALSE)
  fixture <- inputs[inputs$fixture_id == "partial_with_intrusions", ]
  make_role <- function(role) {
    rows <- fixture[fixture$role == role, ]
    fixation_group(
      x = rows$x,
      y = rows$y,
      onset = rows$onset,
      duration = rows$duration
    )
  }
  expect_transport_v3_backend_parity(
    make_role("reference"), make_role("source")
  )
})

test_that("batch size and candidate order cannot change native scores", {
  skip_if_not(eyesim:::transport_v3_native_available())
  coords <- rbind(c(0, 0), c(1, 2), c(3, 1), c(4, 3), c(6, 2))
  path <- function(offset) {
    make_gaze_fixations(
      sweep(coords, 2, offset, FUN = "+"),
      duration = c(1, 2, 1, 1, 2),
      onset = c(0, 1, 3, 4, 5)
    )
  }
  references <- list(
    item_a = path(c(0, 0)),
    item_b = path(c(1.2, -0.5)),
    item_c = path(c(-1.5, 0.8))
  )
  source <- path(c(0.1, -0.05))
  spec <- make_transport_v3_backend_spec("optimized")
  single <- gaze_transport_align_batch(
    references, source, spec, batch_size = 1
  )
  full <- gaze_transport_align_batch(
    references, source, spec, batch_size = 3
  )
  permuted <- gaze_transport_align_batch(
    references[c("item_c", "item_a", "item_b")],
    source,
    spec,
    batch_size = 2
  )
  scores <- function(results) {
    stats::setNames(
      vapply(results, `[[`, numeric(1), "log_score"),
      names(results)
    )
  }

  expect_equal(scores(single), scores(full), tolerance = 1e-12)
  expect_equal(
    scores(single), scores(permuted)[names(single)], tolerance = 1e-12
  )
  first_rank <- eyesim:::score_gaze_engine_results(
    single, true_key = "item_a"
  )$template_rank
  permuted_rank <- eyesim:::score_gaze_engine_results(
    permuted, true_key = "item_a"
  )$template_rank
  expect_equal(first_rank, permuted_rank, tolerance = 0)
  expect_true(all(vapply(single, function(result) {
    result$convergence$converged
  }, logical(1))))
})

test_that("auto backend falls back cleanly when native policy is unsupported", {
  path <- make_gaze_fixations(rbind(c(0, 0), c(1, 1), c(2, 0)))
  auto_spec <- gaze_transport_spec(
    coverage_nodes = 2,
    entropy_schedule = 0.03,
    maxit = 30,
    tolerance = 1e-3,
    projection_maxit = 300,
    projection_tolerance = 1e-7,
    multistart = 2,
    backend = "auto"
  )
  result <- gaze_transport_align(path, path, auto_spec)

  expect_identical(result$alignment$profile$backend, "reference_fallback")
  expect_true(result$alignment$profile$fallback)
  expect_match(
    result$alignment$profile$fallback_reason,
    "multistart",
    fixed = TRUE
  )
  expect_true(is.finite(result$log_score))
})
