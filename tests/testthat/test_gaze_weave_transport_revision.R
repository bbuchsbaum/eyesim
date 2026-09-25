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
  oracle_score <- -1.2068974004

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

# Largest objective decrease that a small projected mirror step achieves from
# each converged coverage node of an alignment (tightly projected, final
# entropy). Zero means no small step descends.
transport_revision_small_step_descent <- function(result, reference, source,
                                                  spec) {
  reference <- eyesim:::as_transport_v3_measure(reference, spec$chronology)
  source <- eyesim:::as_transport_v3_measure(source, spec$chronology)
  cost <- eyesim:::gaze_spatial_cost(
    reference$coords, source$coords, spec$spatial
  )
  entropy <- utils::tail(spec$entropy_schedule, 1)
  rows <- seq_along(reference$mass)
  columns <- seq_along(source$mass)
  tight <- spec$control
  tight$projection_tolerance <- 1e-11
  vapply(result$alignment$profile$fits, function(fit) {
    if (!identical(fit$status, "converged")) return(0)
    augmented <- fit$augmented
    current <- eyesim:::transport_v3_objective(
      augmented[rows, columns, drop = FALSE], reference, source, cost, spec,
      entropy, gradient = TRUE
    )
    gradient <- current$gradient - stats::median(current$gradient)
    best <- 0
    for (step in c(0.1, 0.03, 0.01, 0.003)) {
      kernel <- augmented
      kernel[rows, columns] <- pmax(augmented[rows, columns], 1e-300) *
        exp(pmax(pmin(-step * gradient, 50), -50))
      projection <- eyesim:::project_partial_coupling_revised(
        kernel, reference$mass, source$mass, fit$coverage, tight
      )
      if (!isTRUE(projection$converged)) next
      value <- eyesim:::transport_v3_objective(
        projection$plan[rows, columns, drop = FALSE], reference, source,
        cost, spec, entropy
      )$optimization
      best <- max(best, current$optimization - value)
    }
    best
  }, numeric(1))
}

transport_revision_pair134 <- function() {
  list(
    reference = transport_revision_path(
      c(1.59957729466259, 28.085163069889, 2.61879675090313,
        3.98616542108357, 11.4314303006977),
      c(5.67472292575985, 9.4715854995884, 4.62438578438014,
        8.86520146718249, 5.57847462547943)
    ),
    source = transport_revision_path(
      c(7.00039075873792, 10.6210063286126, 18.2622979432344,
        23.798659957014, 2.2444723052904),
      c(13.8152038375847, 18.9951667808928, 4.35126664955169,
        17.9943125946447, 14.0034753344953)
    )
  )
}

transport_revision_pair118 <- function() {
  list(
    reference = transport_revision_path(
      c(19.4531975416467, 7.80971233732998, 17.2545288726687,
        25.749476599507, 15.7219495503232),
      c(4.1737774903886, 7.78590210527182, 8.19896909128875,
        15.0142093077302, 1.23240248532966)
    ),
    source = transport_revision_path(
      c(19.0306649431586, 5.78759774472564, 4.19821098353714,
        8.76898297201842, 11.7428692597896, 1.33927661459893,
        26.1695548668504, 16.166871888563, 23.1607949361205,
        2.313935123384, 17.7435952061787),
      c(15.0350786959752, 9.20990427210927, 9.78326073940843,
        18.8875433998182, 3.77925990754738, 5.44136211648583,
        2.79010437708348, 1.51434923103079, 9.80730763170868,
        13.8354574074037, 16.7554191015661)
    )
  )
}

test_that("a saturated mirror step never certifies a non-stationary node", {
  skip_if_not(eyesim:::transport_v3_native_available())
  pair <- transport_revision_pair134()
  for (backend in c("optimized", "reference")) {
    spec <- gaze_transport_spec(backend = backend)
    result <- gaze_transport_align(pair$reference, pair$source, spec)
    descent <- transport_revision_small_step_descent(
      result, pair$reference, pair$source, spec
    )
    expect_lte(max(descent), 1e-6)
  }
})

test_that("review pair 118 scores agree across mirror step sizes", {
  # A regression check on one pair only; in general scores still depend on
  # step_size through the local optimum reached (see NEWS).
  skip_if_not(eyesim:::transport_v3_native_available())
  pair <- transport_revision_pair118()
  scores <- vapply(c(2, 0.25, 0.05), function(step) {
    gaze_transport_align(
      pair$reference, pair$source,
      gaze_transport_spec(
        backend = "optimized", step_size = step, maxit = 20000L
      )
    )$log_score
  }, numeric(1))
  expect_lte(diff(range(scores)), 1e-4)
})

test_that("the chronology term is continuous as edge mass vanishes", {
  reference_relation <- matrix(0, 3, 3)
  reference_relation[1, 2] <- reference_relation[2, 3] <- 1
  source_relation <- reference_relation
  plan <- function(epsilon) {
    # All mass on the edge-free diagonal cell (1, 3), plus leakage onto the
    # edge-carrying cells.
    correspondence <- matrix(epsilon, 3, 3)
    correspondence[1, 3] <- 1
    correspondence / sum(correspondence)
  }
  exact <- eyesim:::transport_v3_edge_terms(
    plan(0), reference_relation, source_relation,
    pseudo_count = eyesim:::transport_v3_chronology_pseudo_count
  )$residual
  leaky <- eyesim:::transport_v3_edge_terms(
    plan(1e-9), reference_relation, source_relation,
    pseudo_count = eyesim:::transport_v3_chronology_pseudo_count
  )$residual
  expect_equal(exact, 1)
  expect_lt(abs(leaky - exact), 1e-6)

  # The pseudo-count gradient matches finite differences.
  set.seed(4)
  correspondence <- matrix(stats::runif(9), 3, 3)
  correspondence <- correspondence / sum(correspondence)
  terms <- eyesim:::transport_v3_edge_terms(
    correspondence, reference_relation, source_relation, gradient = TRUE,
    pseudo_count = eyesim:::transport_v3_chronology_pseudo_count
  )
  # Mass-preserving directional derivatives: move h from cell 1 to each cell.
  for (cell in 2:9) {
    shifted <- correspondence
    shifted[cell] <- shifted[cell] + 1e-7
    shifted[1] <- shifted[1] - 1e-7
    numeric <- (eyesim:::transport_v3_edge_terms(
      shifted, reference_relation, source_relation,
      pseudo_count = eyesim:::transport_v3_chronology_pseudo_count
    )$residual - terms$residual) / 1e-7
    expect_equal(
      numeric, terms$gradient[cell] - terms$gradient[1], tolerance = 1e-4
    )
  }
})

test_that("a one-fixation source gets equal chronology from reordered candidates", {
  source <- transport_revision_path(10.2, 6.1)
  coords <- list(x = c(4, 10, 16, 22), y = c(5, 6, 8, 5))
  forward <- transport_revision_path(coords$x, coords$y)
  shuffled <- transport_revision_path(coords$x[c(3, 1, 4, 2)],
                                      coords$y[c(3, 1, 4, 2)])
  for (backend in c("optimized", "reference")) {
    spec <- gaze_transport_spec(backend = backend)
    a <- gaze_transport_align(forward, source, spec)
    b <- gaze_transport_align(shuffled, source, spec)
    expect_equal(
      a$alignment$profile$conditional_chronology,
      b$alignment$profile$conditional_chronology,
      tolerance = 1e-10
    )
    expect_equal(a$log_score, b$log_score, tolerance = 1e-8)
  }
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
    projection_method = "log", backend = "auto"
  )
  expect_no_warning(result <- gaze_transport_align(path, path, spec))
  expect_identical(result$convergence$backend, "reference")
  expect_false(result$convergence$fallback)
  expect_match(result$convergence$backend_route, "log-domain", fixed = TRUE)
})

test_that("native multistart matches the reference oracle", {
  skip_if_not(eyesim:::transport_v3_native_available())
  reference <- transport_revision_path(c(2, 5, 9, 12), c(3, 7, 4, 8))
  source <- transport_revision_path(c(2.4, 5.3, 8.6), c(3.2, 6.5, 4.4))
  native <- gaze_transport_align(
    reference, source,
    gaze_transport_spec(backend = "optimized", multistart = 2L)
  )
  oracle <- gaze_transport_align(
    reference, source,
    gaze_transport_spec(backend = "reference", multistart = 2L)
  )
  expect_identical(native$convergence$backend, "native_rcpparmadillo")
  expect_true(all(vapply(native$alignment$profile$fits, function(fit) {
    "spatial_native" %in% names(fit$start_objectives)
  }, logical(1))))
  expect_lt(abs(native$log_score - oracle$log_score), 1e-6)
})

test_that("both backends read the chronology pseudo-count from one constant", {
  skip_if_not(eyesim:::transport_v3_native_available())
  reference <- transport_revision_path(c(2, 5, 9, 12), c(3, 7, 4, 8))
  source <- transport_revision_path(c(2.4, 5.3, 8.6), c(3.2, 6.5, 4.4))
  score <- function(backend) {
    gaze_transport_align(
      reference, source, gaze_transport_spec(backend = backend)
    )$log_score
  }
  default <- score("optimized")
  local_mocked_bindings(transport_v3_chronology_pseudo_count = 0.2)
  native <- score("optimized")
  oracle <- score("reference")
  expect_gt(abs(native - default), 1e-4)
  expect_lt(abs(native - oracle), 1e-6)
})

test_that("reaching maxit is recorded and scored, not an error", {
  skip_if_not(eyesim:::transport_v3_native_available())
  reference <- transport_revision_path(c(2, 5, 9, 12), c(3, 7, 4, 8))
  source <- transport_revision_path(c(2.4, 5.3, 8.6), c(3.2, 6.5, 4.4))
  spec <- gaze_transport_spec(backend = "optimized", step_size = 0.05, maxit = 5L)
  expect_no_error(result <- gaze_transport_align(reference, source, spec))
  expect_identical(result$convergence$status, "not_converged")
  expect_identical(result$convergence$backend, "native_rcpparmadillo")
  expect_false(result$convergence$fallback)
  expect_true(is.finite(result$log_score))
  expect_no_warning(
    auto <- gaze_transport_align(
      reference, source,
      gaze_transport_spec(backend = "auto", step_size = 0.05, maxit = 5L)
    )
  )
  expect_false(auto$convergence$fallback)
})

# Pair 94 of the review probe set (set.seed(11) generator).
transport_revision_pair94 <- function() {
  list(
    reference = transport_revision_path(
      c(20.1050154371187, 2.82494349312037, 7.02438600268215),
      c(8.91253510117531, 15.691644763574, 4.01146831410006)
    ),
    source = transport_revision_path(
      c(22.3944645132869, 8.49569070152938, 26.5495868250728),
      c(7.96401504566893, 1.03310775477439, 11.2457623830996)
    )
  )
}

test_that("a plan with a zeroed support cell is not certified in place", {
  pairs <- list(
    transport_revision_pair94(),
    list(
      reference = transport_revision_path(
        c(24.116132248193, 15.0186096066609, 2.27322669886053,
          21.8344945404679),
        c(7.28982275770977, 18.4285079068504, 14.6237265886739,
          4.59919406240806)
      ),
      source = transport_revision_path(
        c(23.551057927151, 14.9630431319283, 24.4763657430114),
        c(6.9930906277275, 16.9068108187189, 6.88592172687604)
      )
    )
  )
  spec <- gaze_transport_spec(backend = "reference")
  entropy <- utils::tail(spec$entropy_schedule, 1)
  for (pair in pairs) {
    reference <- eyesim:::as_transport_v3_measure(pair$reference, spec$chronology)
    source <- eyesim:::as_transport_v3_measure(pair$source, spec$chronology)
    swapped <- eyesim:::gaze_measure_order_key(reference) <
      eyesim:::gaze_measure_order_key(source)
    if (swapped) {
      held <- reference
      reference <- source
      source <- held
    }
    cost <- eyesim:::gaze_spatial_cost(
      reference$coords, source$coords, spec$spatial
    )
    rows <- seq_along(reference$mass)
    columns <- seq_along(source$mass)
    objective <- function(plan) {
      eyesim:::transport_v3_objective(
        plan[rows, columns, drop = FALSE], reference, source, cost, spec,
        entropy
      )$optimization
    }
    result <- gaze_transport_align(pair$reference, pair$source, spec)
    fits <- result$alignment$profile$fits
    converged <- which(vapply(fits, function(fit) {
      identical(fit$status, "converged")
    }, logical(1)))
    fit <- fits[[converged[[ceiling(length(converged) / 2)]]]]
    plan <- fit$augmented
    if (swapped) plan <- t(plan)
    for (what in c("real", "slack")) {
      kernel <- plan
      if (what == "real") {
        cell <- which(plan[rows, columns] == max(plan[rows, columns]),
                      arr.ind = TRUE)[1, ]
        kernel[cell[[1]], cell[[2]]] <- 1e-250
      } else {
        kernel[which.max(plan[rows, ncol(plan)]), ncol(plan)] <- 1e-250
      }
      control <- spec$control
      control$projection_tolerance <- 1e-12
      trapped <- eyesim:::project_partial_coupling_revised(
        kernel, reference$mass, source$mass, fit$coverage, control
      )$plan
      stage <- eyesim:::transport_v3_reference_stage_revised(
        trapped, reference, source, cost, fit$coverage, spec, entropy
      )
      if (isTRUE(stage$converged)) {
        expect_lte(
          objective(stage$augmented), objective(plan) + 1e-4,
          label = paste("trapped", what, "cell certified")
        )
      } else {
        succeed()
      }
    }
  }
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
  expect_identical(
    sum(clean$results$solver_not_converged),
    clean$solver$heldout_not_converged_count
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
    # A projection-limited stall or maxit is a recorded outcome, not a
    # failure. Near the projection-noise floor the two backends can end the
    # same node with different such statuses, so only failures are excluded.
    expect_false(identical(native$convergence$status, "numerical_failure"))
    expect_false(identical(oracle$convergence$status, "numerical_failure"))
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
