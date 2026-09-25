# Transport solver revision 2026.10 ----------------------------------------

transport_revision_path <- function(x, y) {
  data.frame(
    x = x,
    y = y,
    onset = seq(0, by = 300, length.out = length(x)),
    duration = 250
  )
}

transport_revision_random_path <- function(n) {
  transport_revision_path(stats::runif(n, 1, 29), stats::runif(n, 1, 19))
}

transport_revision_noisy_path <- function(path, keep) {
  index <- sort(sample(nrow(path), keep))
  transport_revision_path(
    path$x[index] + stats::rnorm(length(index), 0, 1.5),
    path$y[index] + stats::rnorm(length(index), 0, 1.5)
  )
}

# A reproducible synthetic pair bank. Odd pairs are noisy partial replays of
# the reference; even pairs are independent (null) paths.
transport_revision_pairs <- function(n_pairs, sizes, seed) {
  set.seed(seed)
  lapply(seq_len(n_pairs), function(index) {
    reference <- transport_revision_random_path(sample(sizes, 1L))
    source <- if (index %% 2L == 1L) {
      transport_revision_noisy_path(
        reference, max(1L, nrow(reference) - 1L)
      )
    } else {
      transport_revision_random_path(sample(sizes, 1L))
    }
    list(reference = reference, source = source)
  })
}

test_that("specifications record the solver revision", {
  expect_identical(gaze_transport_spec()$revision, "2026.10")
  expect_identical(
    gaze_transport_spec(revision = "2026.08")$revision, "2026.08"
  )
  expect_error(gaze_transport_spec(revision = "2025.01"))
  legacy <- gaze_transport_spec()
  legacy$revision <- "edge_normalized_episode_transport"
  expect_identical(eyesim:::transport_v3_revision(legacy), "2026.08")
  legacy$revision <- NULL
  expect_identical(eyesim:::transport_v3_revision(legacy), "2026.08")
  unknown <- gaze_transport_spec()
  unknown$revision <- "2031.01"
  expect_error(
    eyesim:::transport_v3_revision(unknown),
    "Unknown Transport solver revision"
  )
  path <- transport_revision_path(c(2, 5, 9), c(3, 7, 4))
  expect_error(
    gaze_transport_align(path, path, unknown),
    "Unknown Transport solver revision"
  )
})

test_that("a specification saved at the base commit scores exactly as legacy", {
  # fixtures/transport-spec-04a5b25.rds is gaze_transport_spec() serialised by
  # commit 04a5b25, before solver revisions existed; its revision slot holds
  # the estimand name. The expected scores were computed by that commit.
  frozen <- readRDS(test_path("fixtures", "transport-spec-04a5b25.rds"))
  expect_identical(frozen$revision, "edge_normalized_episode_transport")
  pinned <- gaze_transport_spec(revision = "2026.08")
  simple <- list(
    transport_revision_path(c(2, 5, 9), c(3, 7, 4)),
    transport_revision_path(c(2.4, 5.3, 8.6), c(3.2, 6.5, 4.4))
  )
  stationary <- list(
    transport_revision_path(
      c(23.7, 7.7, 22.43, 6.4), c(11.57, 5.02, 15.13, 13.19)
    ),
    transport_revision_path(c(21.2, 7.15), c(15.37, 14))
  )
  base_scores <- list(
    simple = -0.048261237015642555,
    stationary = -0.58519445915782109
  )
  pairs <- list(simple = simple, stationary = stationary)
  for (name in names(pairs)) {
    pair <- pairs[[name]]
    from_frozen <- gaze_transport_align(pair[[1]], pair[[2]], frozen)
    from_pinned <- gaze_transport_align(pair[[1]], pair[[2]], pinned)
    expect_identical(from_frozen$log_score, base_scores[[name]])
    expect_identical(from_pinned$log_score, base_scores[[name]])
    expect_identical(from_frozen$convergence$revision, "2026.08")
  }
})

test_that("a stationary stage converges in the optimized backend", {
  skip_if_not(eyesim:::transport_v3_native_available())
  reference <- transport_revision_path(
    c(23.7, 7.7, 22.43, 6.4), c(11.57, 5.02, 15.13, 13.19)
  )
  source <- transport_revision_path(c(21.2, 7.15), c(15.37, 14))

  optimized <- gaze_transport_align(
    reference, source, gaze_transport_spec(backend = "optimized")
  )
  oracle <- gaze_transport_align(
    reference, source, gaze_transport_spec(backend = "reference")
  )

  expect_true(optimized$convergence$converged)
  expect_true(oracle$convergence$converged)
  expect_identical(optimized$convergence$backend, "native_rcpparmadillo")
  expect_false(optimized$convergence$fallback)
  expect_lt(abs(optimized$log_score - oracle$log_score), 1e-6)
})

test_that("the stop rule does not stop on one heavily backtracked step", {
  skip_if_not(eyesim:::transport_v3_native_available())
  # Pair 2 of the stop-rule probe (set.seed(21)): two independent paths.
  reference <- transport_revision_path(
    c(24.474, 25.334, 6.407, 7.057, 19.212, 10.385, 15.214, 19.28, 28.036),
    c(10.264, 2.11, 3.718, 12.44, 2.853, 14.908, 8.384, 16.664, 14.897)
  )
  source <- transport_revision_path(
    c(18.65, 24.246, 4.654, 25.692, 7.721, 19.149, 7.537),
    c(2.249, 1.599, 17.605, 7.152, 17.89, 10.193, 3.376)
  )

  # Oracle: reference backend, multistart = 2, tolerance = 1e-8,
  # maxit = 3000 (revision 2026.10). It takes several seconds, so its score
  # is frozen here; inst/validation is not involved. Backend agreement alone
  # is not sufficient, because both backends could share a premature stop.
  # Revision 2026.08 scores this pair -1.3634 (a 0.14-nat stop-rule error).
  oracle_score <- -1.2199479011

  default <- gaze_transport_align(
    reference, source, gaze_transport_spec(backend = "optimized")
  )
  tight <- gaze_transport_align(
    reference, source,
    gaze_transport_spec(
      backend = "optimized", tolerance = 1e-8, maxit = 3000L
    )
  )

  expect_true(default$convergence$converged)
  expect_lt(abs(default$log_score - oracle_score), 0.01)
  expect_lt(abs(default$log_score - tight$log_score), 0.01)
})

test_that("an uncertifiable residual is recorded as a projection-limited stall", {
  skip_if_not(eyesim:::transport_v3_native_available())
  reference <- transport_revision_path(c(2, 5, 9), c(3, 7, 4))
  source <- transport_revision_path(c(2.4, 5.3, 8.6), c(3.2, 6.5, 4.4))
  # A residual tolerance far below what 1e-8 projections can resolve.
  spec <- gaze_transport_spec(
    backend = "optimized", tolerance = 1e-12, coverage_nodes = 4L
  )
  expect_no_warning(result <- gaze_transport_align(reference, source, spec))
  expect_identical(result$convergence$status, "stalled_projection_limited")
  expect_false(result$convergence$converged)
  expect_identical(result$convergence$backend, "native_rcpparmadillo")
  expect_true(is.finite(result$log_score))
  terminations <- unlist(lapply(result$convergence$masses, function(stages) {
    vapply(stages, `[[`, character(1), "termination")
  }))
  expect_true("stalled_projection_limited" %in% terminations)
})

test_that("auto records the backend and never falls back silently", {
  skip_if_not(eyesim:::transport_v3_native_available())
  reference <- transport_revision_path(c(2, 5, 9), c(3, 7, 4))
  source <- transport_revision_path(c(2.4, 5.3, 8.6), c(3.2, 6.5, 4.4))
  spec <- gaze_transport_spec(backend = "auto")

  native <- gaze_transport_align(reference, source, spec)
  expect_identical(native$convergence$backend, "native_rcpparmadillo")
  expect_false(native$convergence$fallback)
  expect_true(is.na(native$convergence$fallback_reason))

  local_mocked_bindings(
    solve_transport_v3_profile_native = function(...) {
      stop("synthetic native failure")
    }
  )
  expect_warning(
    fallback <- gaze_transport_align(reference, source, spec),
    class = "gaze_transport_backend_fallback"
  )
  expect_identical(fallback$convergence$backend, "reference_fallback")
  expect_true(fallback$convergence$fallback)
  expect_match(fallback$convergence$fallback_reason, "synthetic native failure")

  legacy <- gaze_transport_spec(backend = "auto", revision = "2026.08")
  expect_no_warning(
    silent <- gaze_transport_align(reference, source, legacy)
  )
  expect_identical(silent$convergence$backend, "reference_fallback")
})

test_that("unsupported native policies route to the reference backend", {
  path <- transport_revision_path(c(0, 1, 2), c(0, 1, 0))
  spec <- gaze_transport_spec(
    coverage_nodes = 2, entropy_schedule = 0.03, maxit = 30,
    tolerance = 1e-3, projection_maxit = 300, projection_tolerance = 1e-7,
    multistart = 2, backend = "auto"
  )
  expect_no_warning(result <- gaze_transport_align(path, path, spec))
  expect_identical(result$convergence$backend, "reference")
  expect_false(result$convergence$fallback)
  expect_match(result$convergence$backend_route, "multistart", fixed = TRUE)
})

test_that("cross-validation reports backend fallback counts", {
  skip_if_not(eyesim:::transport_v3_native_available())
  make_path <- function(anchor, source = FALSE) {
    offset <- if (source) c(0.1, 0.05) else c(0, 0)
    coords <- rbind(
      anchor + offset, anchor + offset + c(1, 0.3),
      anchor + offset + c(2, -0.2)
    )
    make_gaze_fixations(coords, duration = c(1, 2, 1), onset = c(0, 1, 3))
  }
  references <- expand.grid(
    participant = c("p1", "p2"), item = 1:3, stringsAsFactors = FALSE
  )
  references$fixgroup <- lapply(seq_len(nrow(references)), function(i) {
    make_path(c(references$item[[i]] * 3, 0))
  })
  sources <- references[c("participant", "item")]
  sources$fixgroup <- lapply(seq_len(nrow(sources)), function(i) {
    make_path(c(sources$item[[i]] * 3, 0), source = TRUE)
  })
  spec <- gaze_transport_spec(
    coverage_nodes = 2, entropy_schedule = 0.03, maxit = 40,
    tolerance = 1e-3, projection_maxit = 300, projection_tolerance = 1e-7,
    backend = "auto", reliability = "none"
  )
  run_cv <- function() {
    gaze_transport_cv(
      references, sources, match_on = c("participant", "item"),
      contrast_on = "participant", spec = spec, n_folds = 2
    )
  }

  clean <- run_cv()
  expect_true(all(clean$results$backend_fallbacks == 0L))
  expect_identical(clean$solver$backend_fallback_count, 0L)
  expect_identical(clean$solver$heldout_backends, "native_rcpparmadillo")
  expect_identical(clean$solver$revision, "2026.10")
  expect_true(is.integer(clean$results$solver_stalled))
  expect_identical(
    sum(clean$results$solver_stalled), clean$solver$heldout_stalled_count
  )

  local_mocked_bindings(
    solve_transport_v3_profile_native = function(...) {
      stop("synthetic native failure")
    }
  )
  expect_warning(
    fallen <- run_cv(),
    class = "gaze_transport_backend_fallback"
  )
  expect_true(all(fallen$results$backend_fallbacks > 0L))
  expect_gt(fallen$solver$backend_fallback_count, 0L)
  expect_identical(
    sum(fallen$results$backend_fallbacks),
    fallen$solver$heldout_backend_fallback_count
  )
})

test_that("native and reference backends agree on random synthetic pairs", {
  skip_if_not(eyesim:::transport_v3_native_available())
  pairs <- transport_revision_pairs(10L, 2:5, seed = 20260924L)
  optimized <- gaze_transport_spec(backend = "optimized")
  reference <- gaze_transport_spec(backend = "reference")
  for (pair in pairs) {
    native <- tryCatch(
      gaze_transport_align(pair$reference, pair$source, optimized),
      error = function(condition) condition
    )
    expect_false(inherits(native, "error"))
    if (inherits(native, "error")) next
    oracle <- gaze_transport_align(pair$reference, pair$source, reference)
    expect_true(native$convergence$converged)
    expect_lt(abs(native$log_score - oracle$log_score), 1e-4)
  }
})

test_that("native backend never fails on 150 random synthetic pairs", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.getenv("EYESIM_SLOW_TESTS")))
  skip_if_not(eyesim:::transport_v3_native_available())
  pairs <- transport_revision_pairs(150L, 2:12, seed = 11L)
  optimized <- gaze_transport_spec(backend = "optimized")
  reference <- gaze_transport_spec(backend = "reference")
  failures <- 0L
  differences <- numeric(0)
  for (pair in pairs) {
    native <- tryCatch(
      gaze_transport_align(pair$reference, pair$source, optimized),
      error = function(condition) condition
    )
    if (inherits(native, "error")) {
      failures <- failures + 1L
      next
    }
    oracle <- gaze_transport_align(pair$reference, pair$source, reference)
    differences <- c(differences, abs(native$log_score - oracle$log_score))
  }
  expect_identical(failures, 0L)
  expect_lte(max(differences), 1e-4)
})

test_that("the polish gap is finite and small at a converged solution", {
  skip_if_not_installed("lpSolve")
  reference <- transport_revision_path(c(2, 5, 9, 12), c(3, 7, 4, 8))
  source <- transport_revision_path(c(2.4, 5.3, 8.6), c(3.2, 6.5, 4.4))
  result <- gaze_transport_align(
    reference, source,
    gaze_transport_spec(backend = "reference", polish = "audit")
  )
  polish <- result$diagnostics$polish
  expect_true(polish$all_converged)
  gaps <- vapply(result$alignment$profile$fits, function(fit) {
    fit$polish$gap
  }, numeric(1))
  reasons <- vapply(result$alignment$profile$fits, function(fit) {
    fit$polish$gap_reason
  }, character(1))
  # Every gap is either finite or NA with a stated reason; never a number
  # produced by the floored boundary gradient.
  expect_true(all(is.finite(gaps) | (is.na(gaps) & nzchar(reasons))))
  expect_true(all(is.na(reasons[is.finite(gaps)])))
  # The objective is non-convex, so the gap is a stationarity measure only.
  # At this converged, interior solution it is small on the nats scale.
  expect_true(is.finite(polish$map_gap))
  expect_lte(polish$maximum_gap, 0.01)
})

test_that("the polish gap is NA with a reason on a selection boundary", {
  skip_if_not_installed("lpSolve")
  coupling <- matrix(c(0.3, 0, 0, 0.2), 2, 2)
  expect_true(eyesim:::transport_v3_selection_boundary(
    matrix(c(0.5, 0, 0, 0), 2, 2), c(0.5, 0.5), c(0.5, 0.5)
  ))
  expect_false(eyesim:::transport_v3_selection_boundary(
    coupling, c(0.5, 0.5), c(0.5, 0.5)
  ))
  # A pair whose low-coverage nodes select a single fixation pair.
  reference <- transport_revision_path(
    c(24.474, 25.334, 6.407, 7.057, 19.212, 10.385, 15.214, 19.28, 28.036),
    c(10.264, 2.11, 3.718, 12.44, 2.853, 14.908, 8.384, 16.664, 14.897)
  )
  source <- transport_revision_path(
    c(18.65, 24.246, 4.654, 25.692, 7.721, 19.149, 7.537),
    c(2.249, 1.599, 17.605, 7.152, 17.89, 10.193, 3.376)
  )
  result <- gaze_transport_align(
    reference, source,
    gaze_transport_spec(backend = "optimized", polish = "audit")
  )
  fits <- result$alignment$profile$fits
  gaps <- vapply(fits, function(fit) fit$polish$gap, numeric(1))
  reasons <- vapply(fits, function(fit) fit$polish$gap_reason, character(1))
  expect_true(any(is.na(gaps)))
  expect_match(reasons[is.na(gaps)], "^selection_boundary")
  expect_identical(
    result$diagnostics$polish$undefined_gap_nodes, sum(is.na(gaps))
  )
})
