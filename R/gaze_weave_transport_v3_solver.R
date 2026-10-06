# Pure-R Transport oracle -----------------------------------------------

as_transport_v3_measure <- function(x, chronology) {
  measure <- if (inherits(x, "gaze_measure")) x else as_gaze_measure(x, chronology)
  if (!inherits(measure, "gaze_measure")) {
    stop("x must produce a gaze_measure.")
  }
  if (identical(chronology$type, "order_neighbours")) {
    measure$relation <- 1 * (measure$relation > 0)
  }
  measure
}

transport_v3_jensen_shannon <- function(observed, target) {
  if (!is.numeric(observed) || !is.numeric(target) ||
      length(observed) != length(target) || any(!is.finite(observed)) ||
      any(!is.finite(target)) || any(observed < 0) || any(target < 0) ||
      abs(sum(observed) - 1) > 1e-8 || abs(sum(target) - 1) > 1e-8) {
    stop("observed and target must be finite non-negative unit masses.")
  }
  midpoint <- (observed + target) / 2
  positive_observed <- observed > 0
  positive_target <- target > 0
  0.5 * (
    sum(observed[positive_observed] * log(
      observed[positive_observed] / midpoint[positive_observed]
    )) +
      sum(target[positive_target] * log(
        target[positive_target] / midpoint[positive_target]
      ))
  )
}

transport_v3_jensen_shannon_gradient <- function(observed, target) {
  midpoint <- (observed + target) / 2
  0.5 * log(
    pmax(observed, .Machine$double.xmin) /
      pmax(midpoint, .Machine$double.xmin)
  )
}

# Revision 2026.10 mutual information. Identical to gaze_mutual_information()
# except that a positive cell whose independent reference mass underflows to
# zero is evaluated in the log domain instead of contributing Inf.
transport_v3_mutual_information <- function(correspondence) {
  if (!is.matrix(correspondence) || any(!is.finite(correspondence)) ||
      any(correspondence < 0) || abs(sum(correspondence) - 1) > 1e-8) {
    stop("correspondence must be a finite non-negative unit-mass matrix.")
  }
  alpha <- rowSums(correspondence)
  beta <- colSums(correspondence)
  reference <- outer(alpha, beta)
  positive <- correspondence > 0
  value <- correspondence[positive]
  ratio <- value / reference[positive]
  unstable <- !(reference[positive] > 0) | !is.finite(ratio)
  terms <- value * log(ratio)
  if (any(unstable)) {
    cells <- which(positive, arr.ind = TRUE)[unstable, , drop = FALSE]
    terms[unstable] <- value[unstable] * (
      log(value[unstable]) - log(alpha[cells[, 1L]]) - log(beta[cells[, 2L]])
    )
  }
  sum(terms)
}

transport_v3_objective <- function(coupling, reference, source, spatial_cost,
                                   spec, entropy, warp_penalty = 0,
                                   gradient = FALSE) {
  if (!is.matrix(coupling) || any(dim(coupling) != dim(spatial_cost)) ||
      any(!is.finite(coupling)) || any(coupling < 0)) {
    stop("coupling must be a compatible finite non-negative matrix.")
  }
  if (!is.numeric(warp_penalty) || length(warp_penalty) != 1L ||
      !is.finite(warp_penalty) || warp_penalty < 0) {
    stop("warp_penalty must be one finite non-negative value.")
  }
  coverage <- sum(coupling)
  if (coverage < 0 || coverage > 1 + 1e-8) {
    stop("coupling coverage must lie in [0, 1].")
  }
  if (coverage <= .Machine$double.eps) {
    zero_components <- c(
      spatial = 0,
      chronology = 0,
      reference_selection = 0,
      source_selection = 0,
      warp_penalty = warp_penalty,
      correspondence_smoothing = 0
    )
    result <- list(
      scientific = warp_penalty,
      optimization = warp_penalty,
      coverage = 0,
      correspondence = NULL,
      selected_reference_mass = rep(0, length(reference$mass)),
      selected_source_mass = rep(0, length(source$mass)),
      components = zero_components,
      conditional = c(
        spatial = NA_real_,
        chronology = NA_real_,
        reference_selection = NA_real_,
        source_selection = NA_real_,
        correspondence_information = NA_real_
      ),
      zero_coverage = TRUE
    )
    if (gradient) {
      result$gradient <- matrix(0, nrow(coupling), ncol(coupling))
    }
    return(result)
  }

  correspondence <- coupling / coverage
  reference_selected_unit <- rowSums(correspondence)
  source_selected_unit <- colSums(correspondence)
  spatial <- sum(correspondence * spatial_cost)
  edge <- transport_v3_edge_terms(
    correspondence,
    reference$relation,
    source$relation,
    gradient = gradient,
    pseudo_count = if (transport_v3_revised(spec)) {
      transport_v3_chronology_pseudo_count
    } else {
      0
    }
  )
  reference_selection <- transport_v3_jensen_shannon(
    reference_selected_unit, reference$mass
  )
  source_selection <- transport_v3_jensen_shannon(
    source_selected_unit, source$mass
  )
  correspondence_information <- if (transport_v3_revised(spec)) {
    transport_v3_mutual_information(correspondence)
  } else {
    gaze_mutual_information(correspondence)
  }
  scientific_components <- c(
    spatial = coverage * spatial,
    chronology = coverage * spec$temporal_weight * edge$residual,
    reference_selection = coverage * spec$selection_weights[[1]] *
      reference_selection,
    source_selection = coverage * spec$selection_weights[[2]] * source_selection,
    warp_penalty = warp_penalty
  )
  smoothing <- entropy * coverage * correspondence_information
  components <- c(
    scientific_components,
    correspondence_smoothing = smoothing
  )
  result <- list(
    scientific = sum(scientific_components),
    optimization = sum(scientific_components) + smoothing,
    coverage = coverage,
    correspondence = correspondence,
    selected_reference_mass = rowSums(coupling),
    selected_source_mass = colSums(coupling),
    components = components,
    conditional = c(
      spatial = spatial,
      chronology = edge$residual,
      reference_selection = reference_selection,
      source_selection = source_selection,
      correspondence_information = correspondence_information
    ),
    edge = edge,
    zero_coverage = FALSE
  )
  if (!gradient) return(result)

  reference_selection_gradient <- transport_v3_jensen_shannon_gradient(
    reference_selected_unit, reference$mass
  )
  source_selection_gradient <- transport_v3_jensen_shannon_gradient(
    source_selected_unit, source$mass
  )
  result$gradient <- spatial_cost +
    spec$temporal_weight * edge$gradient +
    outer(
      spec$selection_weights[[1]] * reference_selection_gradient,
      rep(1, ncol(coupling))
    ) +
    outer(
      rep(1, nrow(coupling)),
      spec$selection_weights[[2]] * source_selection_gradient
    ) +
    entropy * gaze_mutual_information_gradient(correspondence)
  result
}

# Revision 2026.10 first-order stopping rule. It is a heuristic model of the
# decrease still available at a node, not a proof of optimality (the
# objective is non-convex). With dual prices a, b fitted by P-weighted least
# squares to T (T = centred gradient G on real cells, 0 on slack cells),
# delta_ij = a_i + b_j - T_ij is the first-order gain per unit of mass moved
# into cell ij. The predicted remaining decrease is
#   real cells:  sum entropy * P * (exp(x) - 1 - x),  x = delta / entropy,
#                the exact gain of relaxing one cell against its entropic
#                barrier (equal to P delta^2 / (2 entropy) for small x);
#   slack cells: sum P delta^2 / (2 entropy) (no barrier; quadratic model).
# A slack cell with negligible mass (< 1e-12) and delta > 1e-3 is a trap:
# mirror descent grows it only multiplicatively, so it cannot re-enter the
# support in useful time. The same holds for a real cell whose relaxed mass
# P exp(delta / entropy) (uncapped, with P = 0 counted as the smallest
# normal double) exceeds 1e-9 while P < 1e-12. Traps block certification and
# trigger a reseed. Mirrors stationarity_check() in src/transport_v3.cpp.
transport_v3_stationarity <- function(augmented, centered_gradient, entropy) {
  n_rows <- nrow(augmented)
  n_columns <- ncol(augmented)
  weight <- augmented
  weight[n_rows, n_columns] <- 0
  target <- matrix(0, n_rows, n_columns)
  target[-n_rows, -n_columns] <- centered_gradient
  row_weight <- rowSums(weight)
  column_weight <- colSums(weight)
  normal <- rbind(
    cbind(diag(row_weight, n_rows), weight),
    cbind(t(weight), diag(column_weight, n_columns))
  )
  right <- c(rowSums(weight * target), colSums(weight * target))
  free <- seq_len(n_rows + n_columns - 1L)
  reduced <- normal[free, free, drop = FALSE]
  diag(reduced) <- diag(reduced) + 1e-14 * max(diag(reduced))
  solution <- tryCatch(
    solve(reduced, right[free]),
    error = function(condition) NULL
  )
  if (is.null(solution) || any(!is.finite(solution))) {
    return(list(predicted = Inf, trap = FALSE, s0 = Inf))
  }
  solution <- c(solution, 0)
  delta <- outer(
    solution[seq_len(n_rows)], solution[n_rows + seq_len(n_columns)],
    FUN = "+"
  ) - target
  real <- matrix(FALSE, n_rows, n_columns)
  real[-n_rows, -n_columns] <- TRUE
  slack <- !real
  slack[n_rows, n_columns] <- FALSE
  x <- pmin(delta[real] / entropy, 700)
  mass <- weight[real]
  gain_real <- sum(entropy * mass * (exp(x) - 1 - x))
  gain_slack <- sum(weight[slack] * delta[slack]^2) / (2 * entropy)
  # Uncapped relaxed log-mass; an exactly zero cell is floored at the
  # smallest normal double so that it is not invisible to the test.
  relaxed_log_mass <- log(pmax(mass, .Machine$double.xmin)) +
    delta[real] / entropy
  real_trap <- any(mass < 1e-12 & relaxed_log_mass > log(1e-9))
  slack_trap <- any(weight[slack] < 1e-12 & delta[slack] > 1e-3)
  list(
    predicted = gain_real + gain_slack,
    trap = real_trap || slack_trap,
    s0 = sum(weight * delta^2)
  )
}

# Revision 2026.10 reseed: mix the plan with the independent start and
# re-project, lifting every cell to a non-negligible mass so that a trapped
# cell can re-enter the support.
transport_v3_reseed_fraction <- 1e-3

# One entropy stage of the revision 2026.10 mirror-descent solver. The native
# backend (src/transport_v3.cpp, revision = 1) implements the same rules:
#
# * Objective: the chronology residual is 1 - 2A / (R + S + kappa), with
#   kappa = transport_v3_chronology_pseudo_count passed from R.
# * Stopping rule (first-order, heuristic): a stage stops as "stationary"
#   when the plan is feasible, no cell is trapped, and the first-order model
#   of the decrease still available at the plan (transport_v3_stationarity:
#   entropic single-cell relaxation gains on real cells plus a quadratic
#   model on slack cells, from P-weighted least-squares dual prices) is at
#   most `tolerance` nats. This bounds neither the distance to a local
#   optimum along slow, low-curvature directions nor escape from saddles; on
#   the review set (60 pairs, 720 nodes, default tolerance 1e-6) the q95
#   decrease that a tight continuation still found at certified nodes was
#   9.7e-7, but 9 of 692 nodes still decreased by 1e-5 to 1.4e-3.
# * Traps: a negligible-mass cell whose relaxed mass exceeds 1e-9 (real) or
#   whose dual gain exceeds 1e-3 (slack) blocks certification and triggers a
#   reseed (mix 1e-3 of the independent start, re-project), at most three
#   times per stage.
# * Step limit: every trial step is capped at
#   step_limit = min(step_size, 50 / max|G|), so the +/-50 exponent clamp
#   never binds. After an accepted trial the next start step doubles.
# * Projection: standard Sinkhorn, then a damped log-domain dual Newton
#   finisher, then log-domain Sinkhorn.
# * Noise-limited iterations: both revisions accept a trial that rises by at
#   most 1e-12 relative, so a line search "fails" only when every trial rose
#   by more than that, which at 1e-8 projection accuracy is noise-dominated;
#   a failure establishes nothing. When no trial is accepted, or the accepted
#   decrease is at most 10 * projection_tolerance * max(1, |f|) at a step
#   backtracked below step_limit / 8, the stage ends
#   "stalled_projection_limited": not converged, still scored. A stage that
#   reaches maxit ends "maxit" and its node is "not_converged": recorded and
#   scored, not an error.
# * Backtracking continues past 21 trials while the mirror exponent is large.
# * Starts: every coverage node is solved from the independent start, the
#   spatial start when multistart = 2, and the adjacent-coverage
#   continuation. Fits whose stages all ended converged or stalled are
#   preferred over fits that hit maxit; within that tier the lowest
#   regularized objective is kept (ties keep the earlier start). The result
#   is still a local optimum of a non-convex objective.
transport_v3_reference_stage_revised <- function(augmented, reference, source,
                                                 spatial_cost, coverage, spec,
                                                 entropy) {
  real_rows <- seq_len(length(reference$mass))
  real_columns <- seq_len(length(source$mass))
  row_target <- c(reference$mass, 1 - coverage)
  column_target <- c(source$mass, 1 - coverage)
  control <- spec$control
  step <- control$step_size
  stage_converged <- FALSE
  stalled <- FALSE
  termination <- "maxit"
  final_change <- Inf
  final_objective_change <- Inf
  final_step <- NA_real_
  reduced_gradient <- Inf
  reseeds <- 0L
  projection <- list(error = Inf, converged = FALSE, method = NA_character_)
  iteration <- 0L
  for (iteration in seq_len(control$maxit)) {
    coupling <- augmented[real_rows, real_columns, drop = FALSE]
    current <- transport_v3_objective(
      coupling, reference, source, spatial_cost, spec, entropy,
      gradient = TRUE
    )
    centered_gradient <- current$gradient - stats::median(current$gradient)
    feasibility <- max(
      abs(rowSums(augmented) - row_target),
      abs(colSums(augmented) - column_target)
    )
    stationarity <- transport_v3_stationarity(
      augmented, centered_gradient, entropy
    )
    reduced_gradient <- stationarity$predicted
    if (stationarity$trap && reseeds < 3L) {
      # A cell outside the support should re-enter it: mix in the
      # independent start and re-project, then keep iterating.
      reseeds <- reseeds + 1L
      mixed <- (1 - transport_v3_reseed_fraction) * augmented +
        transport_v3_reseed_fraction *
          transport_augmented_start(reference$mass, source$mass, coverage)
      reseeded <- project_partial_coupling_revised(
        mixed, reference$mass, source$mass, coverage, control
      )
      if (isTRUE(reseeded$converged)) {
        augmented <- reseeded$plan
        projection <- reseeded
        step <- control$step_size
        next
      }
    }
    if (feasibility <= control$projection_tolerance && !stationarity$trap &&
        stationarity$predicted <= control$tolerance) {
      stage_converged <- TRUE
      termination <- "stationary"
      # The certificate concerns the current plan: report its own marginal
      # error, and the method of the projection that produced it.
      projection <- list(
        error = feasibility,
        converged = TRUE,
        method = if (is.na(projection$method)) "input" else projection$method,
        fallback_from_standard = projection$fallback_from_standard
      )
      break
    }
    gradient_scale <- max(abs(centered_gradient))
    step_limit <- if (gradient_scale > 0) {
      min(control$step_size, 50 / gradient_scale)
    } else {
      control$step_size
    }
    scale <- max(1, abs(current$optimization))
    accepted <- FALSE
    trial_step <- min(step, step_limit)
    proposal <- augmented
    proposal_objective <- current
    trial_projection <- projection
    backtrack <- 0L
    repeat {
      if (backtrack > 20L &&
          (trial_step * gradient_scale <= 1e-4 || backtrack > 80L)) {
        break
      }
      backtrack <- backtrack + 1L
      kernel <- augmented
      update <- exp(pmax(pmin(-trial_step * centered_gradient, 50), -50))
      kernel[real_rows, real_columns] <- pmax(coupling, 1e-300) * update
      trial_projection <- project_partial_coupling_revised(
        kernel, reference$mass, source$mass, coverage, control
      )
      if (!trial_projection$converged) {
        trial_step <- trial_step / 2
        next
      }
      proposal <- trial_projection$plan
      proposal_objective <- transport_v3_objective(
        proposal[real_rows, real_columns, drop = FALSE],
        reference, source, spatial_cost, spec, entropy
      )
      if (is.finite(proposal_objective$optimization) &&
          proposal_objective$optimization <= current$optimization +
            1e-12 * max(1, abs(current$optimization))) {
        accepted <- TRUE
        break
      }
      trial_step <- trial_step / 2
    }
    decrease <- current$optimization - proposal_objective$optimization
    noise_limited <- !accepted || (
      trial_step < step_limit / 8 &&
        decrease <= 10 * control$projection_tolerance * scale
    )
    if (noise_limited) {
      stalled <- TRUE
      termination <- "stalled_projection_limited"
      final_objective_change <- if (accepted) abs(decrease) / scale else Inf
      break
    }
    final_change <- max(abs(proposal - augmented))
    final_objective_change <- abs(decrease) / scale
    final_step <- trial_step
    projection <- trial_projection
    augmented <- proposal
    # Doubling (2026.08 used 1.1) lets a backtracked step recover quickly.
    step <- min(trial_step * 2, control$step_size)
  }
  list(
    augmented = augmented,
    converged = stage_converged,
    stalled = stalled,
    final_change = final_change,
    final_objective_change = final_objective_change,
    final_step = final_step,
    history = list(
      entropy = entropy,
      iterations = iteration,
      converged = stage_converged,
      coupling_change = final_change,
      relative_objective_change = final_objective_change,
      projected_update_residual = final_change /
        max(final_step, .Machine$double.eps),
      projection_error = projection$error,
      projection_converged = isTRUE(projection$converged),
      projection_method = projection$method,
      projection_fallback = isTRUE(projection$fallback_from_standard),
      termination = termination,
      step = final_step,
      predicted_decrease = reduced_gradient,
      reseeds = reseeds
    )
  )
}

solve_transport_v3_mass_reference <- function(reference, source, spatial_cost,
                                              coverage, spec,
                                              initial_augmented = NULL) {
  revised <- transport_v3_revised(spec)
  n_reference <- length(reference$mass)
  n_source <- length(source$mass)
  real_rows <- seq_len(n_reference)
  real_columns <- seq_len(n_source)
  starts <- transport_initial_plans(
    reference, source, spatial_cost, coverage, spec
  )
  if (!is.null(initial_augmented)) {
    project <- if (revised) {
      project_partial_coupling_revised
    } else {
      project_partial_coupling
    }
    projected <- project(
      initial_augmented,
      reference$mass,
      source$mass,
      coverage,
      spec$control
    )
    if (projected$converged) {
      # Coverage continuation is within one candidate and follows the same
      # ascending-node policy for every candidate. It never carries state
      # across candidates. Revision 2026.10 solves from the continuation in
      # addition to the structural starts, so the node's solution does not
      # depend on how far the previous node was optimized.
      if (revised) {
        starts$adjacent_coverage <- projected$plan
      } else {
        starts <- list(adjacent_coverage = projected$plan)
      }
    }
  }

  fits <- lapply(names(starts), function(start_name) {
    augmented <- starts[[start_name]]
    history <- list()
    all_stages_converged <- TRUE
    projection_converged <- TRUE
    final_change <- Inf
    final_objective_change <- Inf
    final_step <- NA_real_
    terminations <- character(0)
    for (entropy in spec$entropy_schedule) {
      if (revised) {
        stage <- transport_v3_reference_stage_revised(
          augmented, reference, source, spatial_cost, coverage, spec, entropy
        )
        augmented <- stage$augmented
        # Keep the last accepted step across stages, as the native backend
        # does: a stage that takes no step (stationary at its input plan)
        # leaves final_step NA, which would make the residual NA.
        if (is.finite(stage$final_step)) {
          final_change <- stage$final_change
          final_objective_change <- stage$final_objective_change
          final_step <- stage$final_step
        }
        history[[length(history) + 1L]] <- c(
          stage$history, list(start = start_name)
        )
        projection_converged <- projection_converged &&
          stage$history$projection_converged
        all_stages_converged <- all_stages_converged && stage$converged
        terminations <- c(terminations, stage$history$termination)
        next
      }
      step <- spec$control$step_size
      stage_converged <- FALSE
      projection <- list(error = Inf, converged = FALSE, method = NA_character_)
      for (iteration in seq_len(spec$control$maxit)) {
        coupling <- augmented[real_rows, real_columns, drop = FALSE]
        current <- transport_v3_objective(
          coupling, reference, source, spatial_cost, spec, entropy,
          gradient = TRUE
        )
        centered_gradient <- current$gradient - stats::median(current$gradient)
        accepted <- FALSE
        trial_step <- step
        proposal <- augmented
        proposal_objective <- current
        for (backtrack in 0:20) {
          kernel <- augmented
          update <- exp(pmax(pmin(-trial_step * centered_gradient, 50), -50))
          kernel[real_rows, real_columns] <-
            pmax(coupling, 1e-300) * update
          projection <- project_partial_coupling(
            kernel,
            reference$mass,
            source$mass,
            coverage,
            spec$control
          )
          if (!projection$converged) {
            trial_step <- trial_step / 2
            next
          }
          proposal <- projection$plan
          proposal_coupling <- proposal[real_rows, real_columns, drop = FALSE]
          proposal_objective <- transport_v3_objective(
            proposal_coupling,
            reference,
            source,
            spatial_cost,
            spec,
            entropy
          )
          if (is.finite(proposal_objective$optimization) &&
              proposal_objective$optimization <= current$optimization +
                1e-12 * max(1, abs(current$optimization))) {
            accepted <- TRUE
            break
          }
          trial_step <- trial_step / 2
        }
        if (!accepted) break
        final_change <- max(abs(proposal - augmented))
        final_objective_change <- abs(
          proposal_objective$optimization - current$optimization
        ) / max(1, abs(current$optimization))
        final_step <- trial_step
        augmented <- proposal
        step <- min(trial_step * 1.1, spec$control$step_size)
        if (isTRUE(
          final_objective_change <= spec$control$tolerance &&
            projection$error <= spec$control$projection_tolerance
        )) {
          stage_converged <- TRUE
          break
        }
      }
      history[[length(history) + 1L]] <- list(
        entropy = entropy,
        iterations = iteration,
        converged = stage_converged,
        coupling_change = final_change,
        relative_objective_change = final_objective_change,
        projected_update_residual = final_change /
          max(final_step, .Machine$double.eps),
        projection_error = projection$error,
        projection_converged = projection$converged,
        projection_method = projection$method,
        projection_fallback = isTRUE(projection$fallback_from_standard)
      )
      projection_converged <- projection_converged && projection$converged
      all_stages_converged <- all_stages_converged && stage_converged
    }
    coupling <- augmented[real_rows, real_columns, drop = FALSE]
    final <- transport_v3_objective(
      coupling,
      reference,
      source,
      spatial_cost,
      spec,
      utils::tail(spec$entropy_schedule, 1)
    )
    fit <- list(
      start = start_name,
      coupling = coupling,
      augmented = augmented,
      objective = final,
      history = history,
      converged = projection_converged && all_stages_converged,
      projection_converged = projection_converged,
      stages_converged = all_stages_converged,
      final_change = final_change,
      final_objective_change = final_objective_change,
      projected_update_residual = final_change /
        max(final_step, .Machine$double.eps)
    )
    if (revised) {
      fit$status <- if (isTRUE(fit$converged)) {
        "converged"
      } else if (any(terminations == "maxit")) {
        "not_converged"
      } else {
        "stalled_projection_limited"
      }
    }
    fit
  })
  optimization <- vapply(fits, function(fit) {
    fit$objective$optimization
  }, numeric(1))
  scientific <- vapply(fits, function(fit) {
    fit$objective$scientific
  }, numeric(1))
  converged <- vapply(fits, function(fit) isTRUE(fit$converged), logical(1))
  # Revision 2026.08 silently restarted a non-converged warm-started node from
  # the cold start (the native backend never did). Revision 2026.10 removes
  # that restart: every node is solved from both the independent and the
  # continuation start in both backends, and the node status is reported.
  if (!revised && !any(converged) && !is.null(initial_augmented)) {
    fallback <- solve_transport_v3_mass_reference(
      reference, source, spatial_cost, coverage, spec,
      initial_augmented = NULL
    )
    fallback$continuation_fallback <- TRUE
    fallback$continuation_start_objectives <- stats::setNames(
      optimization, names(starts)
    )
    return(fallback)
  }
  best <- if (revised) {
    # Prefer fits whose stages all ended certified or projection-limited
    # over fits that hit maxit; within the preferred tier keep the lowest
    # finite regularized objective (ties keep the first start).
    preferred <- vapply(fits, function(fit) {
      fit$status %in% c("converged", "stalled_projection_limited")
    }, logical(1))
    tier <- if (any(preferred)) which(preferred) else seq_along(fits)
    finite <- tier[is.finite(optimization[tier])]
    if (length(finite) == 0L) {
      fits[[tier[[1L]]]]
    } else {
      fits[[finite[[which.min(optimization[finite])]]]]
    }
  } else {
    fits[[which.min(optimization)]]
  }
  best$coverage <- coverage
  best$start_objectives <- stats::setNames(optimization, names(starts))
  best$start_scientific_objectives <- stats::setNames(scientific, names(starts))
  best$start_spread <- diff(range(optimization))
  best$start_scientific_spread <- diff(range(scientific))
  best
}

solve_transport_v3_profile_reference <- function(reference, source, spec,
                                                 quadrature = spec$coverage) {
  spatial_cost <- gaze_spatial_cost(
    reference$coords, source$coords, spec$spatial
  )
  fits <- vector("list", length(quadrature$coverage))
  previous_augmented <- NULL
  for (index in seq_along(quadrature$coverage)) {
    fits[[index]] <- solve_transport_v3_mass_reference(
      reference,
      source,
      spatial_cost,
      quadrature$coverage[[index]],
      spec,
      initial_augmented = previous_augmented
    )
    previous_augmented <- fits[[index]]$augmented
  }
  structure(
    list(
      coverage = quadrature$coverage,
      quadrature = quadrature,
      scientific_energy = vapply(fits, function(fit) {
        fit$objective$scientific
      }, numeric(1)),
      regularized_energy = vapply(fits, function(fit) {
        fit$objective$optimization
      }, numeric(1)),
      conditional_spatial = vapply(fits, function(fit) {
        fit$objective$conditional[["spatial"]]
      }, numeric(1)),
      conditional_chronology = vapply(fits, function(fit) {
        fit$objective$conditional[["chronology"]]
      }, numeric(1)),
      reference_selection = vapply(fits, function(fit) {
        fit$objective$conditional[["reference_selection"]]
      }, numeric(1)),
      source_selection = vapply(fits, function(fit) {
        fit$objective$conditional[["source_selection"]]
      }, numeric(1)),
      correspondence_information = vapply(fits, function(fit) {
        fit$objective$conditional[["correspondence_information"]]
      }, numeric(1)),
      fits = fits,
      spatial_cost = spatial_cost,
      backend = "reference"
    ),
    class = c("gaze_transport_profile", "list")
  )
}

transport_v3_integrate_coverage <- function(scientific_energy, quadrature) {
  if (!is.numeric(scientific_energy) ||
      length(scientific_energy) != length(quadrature$coverage) ||
      any(!is.finite(scientific_energy))) {
    stop("scientific_energy must match the finite quadrature nodes.")
  }
  log_component <- quadrature$log_weight - scientific_energy
  log_score <- gaze_log_sum_exp(log_component)
  posterior <- exp(log_component - log_score)
  list(
    log_score = log_score,
    log_component = log_component,
    posterior = posterior,
    expected_coverage = sum(quadrature$coverage * posterior),
    map_index = which.max(posterior),
    map_coverage = quadrature$coverage[[which.max(posterior)]]
  )
}

transport_v3_profile_evidence <- function(profile) {
  transport_v3_integrate_coverage(
    profile$scientific_energy, profile$quadrature
  )
}

transpose_transport_v3_profile <- function(profile) {
  profile$spatial_cost <- t(profile$spatial_cost)
  profile$fits <- lapply(profile$fits, function(fit) {
    fit$coupling <- t(fit$coupling)
    fit$augmented <- t(fit$augmented)
    selected_reference <- fit$objective$selected_reference_mass
    fit$objective$selected_reference_mass <-
      fit$objective$selected_source_mass
    fit$objective$selected_source_mass <- selected_reference
    reference_selection <- fit$objective$components[["reference_selection"]]
    fit$objective$components[["reference_selection"]] <-
      fit$objective$components[["source_selection"]]
    fit$objective$components[["source_selection"]] <- reference_selection
    fit
  })
  selection <- profile$reference_selection
  profile$reference_selection <- profile$source_selection
  profile$source_selection <- selection
  profile
}

gaze_measure_order_key <- function(measure) {
  values <- c(
    nrow(measure$coords),
    as.numeric(measure$coords),
    as.numeric(measure$mass),
    as.numeric(measure$relation)
  )
  paste(
    formatC(values, digits = 17L, format = "fg", flag = "#"),
    collapse = "|"
  )
}

solve_transport_v3_profile_symmetric <- function(reference, source, spec,
                                                 quadrature = spec$coverage) {
  reference_key <- gaze_measure_order_key(reference)
  source_key <- gaze_measure_order_key(source)
  if (reference_key >= source_key) {
    return(solve_transport_v3_profile_ordered(
      reference, source, spec, quadrature
    ))
  }
  transpose_transport_v3_profile(
    solve_transport_v3_profile_ordered(source, reference, spec, quadrature)
  )
}

# Per-coverage-node solver status. Revision 2026.10 fits carry an explicit
# status; 2026.08 fits are either converged or not.
transport_v3_node_status <- function(fits) {
  vapply(fits, function(fit) {
    if (!is.null(fit$status)) return(fit$status)
    if (isTRUE(fit$converged)) "converged" else "not_converged"
  }, character(1))
}

# Worst node status of one alignment.
transport_v3_alignment_status <- function(node_status) {
  severity <- c(
    "converged", "stalled_projection_limited", "not_converged",
    "numerical_failure"
  )
  severity[[max(match(node_status, severity), 1L, na.rm = TRUE)]]
}

transport_v3_alignment_from_profile <- function(profile, reference, source,
                                                registered_source, warp, spec) {
  evidence <- transport_v3_profile_evidence(profile)
  map_fit <- profile$fits[[evidence$map_index]]
  coupling <- map_fit$coupling
  correspondence <- coupling / sum(coupling)
  dx <- outer(reference$coords[, 1], registered_source$coords[, 1], FUN = "-")
  dy <- outer(reference$coords[, 2], registered_source$coords[, 2], FUN = "-")
  spatial_rmse <- sqrt(sum(correspondence * (dx^2 + dy^2)))
  warp_summary <- warp_parameters(warp)
  diagnostics <- list(
    matched_coverage = evidence$expected_coverage,
    map_coverage = evidence$map_coverage,
    selected_reference_mass = map_fit$objective$selected_reference_mass,
    selected_source_mass = map_fit$objective$selected_source_mass,
    conditional_spatial_residual = map_fit$objective$conditional[["spatial"]],
    spatial_rmse = spatial_rmse,
    conditional_chronology_residual =
      map_fit$objective$conditional[["chronology"]],
    local_order_preservation =
      1 - map_fit$objective$conditional[["chronology"]],
    reference_selection_residual =
      map_fit$objective$conditional[["reference_selection"]],
    source_selection_residual =
      map_fit$objective$conditional[["source_selection"]],
    correspondence_smoothing =
      map_fit$objective$components[["correspondence_smoothing"]],
    correspondence_information =
      map_fit$objective$conditional[["correspondence_information"]],
    warp_penalty = map_fit$objective$components[["warp_penalty"]],
    warp_scale = warp_summary$scale,
    warp_translation = warp_summary$translation,
    reference_dominance_error = max(
      map_fit$objective$selected_reference_mass - reference$mass
    ),
    source_dominance_error = max(
      map_fit$objective$selected_source_mass - registered_source$mass
    ),
    start_spread = map_fit$start_spread,
    start_scientific_spread = map_fit$start_scientific_spread,
    coupling_change = map_fit$final_change,
    stationarity = map_fit$projected_update_residual,
    relative_objective_change = map_fit$final_objective_change
  )
  if (!is.null(profile$polish)) {
    polish_gaps <- vapply(profile$fits, function(fit) {
      fit$polish$gap
    }, numeric(1))
    polish_improvements <- vapply(profile$fits, function(fit) {
      fit$polish$improvement
    }, numeric(1))
    polish_feasibility <- vapply(profile$fits, function(fit) {
      max(fit$polish$feasibility)
    }, numeric(1))
    diagnostics$polish <- list(
      enabled = TRUE,
      oracle = profile$polish$oracle,
      map_gap = map_fit$polish$gap,
      maximum_gap = if (!transport_v3_revised(spec)) {
        max(polish_gaps)
      } else if (any(is.finite(polish_gaps))) {
        max(polish_gaps[is.finite(polish_gaps)])
      } else {
        NA_real_
      },
      map_improvement = map_fit$polish$improvement,
      total_improvement = sum(polish_improvements),
      maximum_feasibility_error = max(polish_feasibility),
      no_worse = profile$polish$no_worse,
      all_feasible = profile$polish$all_feasible,
      all_converged = profile$polish$all_converged,
      stopping_reason = map_fit$polish$stopping_reason,
      trace = map_fit$polish$trace
    )
    if (transport_v3_revised(spec)) {
      diagnostics$polish$gap_reason <- map_fit$polish$gap_reason
      diagnostics$polish$undefined_gap_nodes <- sum(is.na(polish_gaps))
    }
  }
  entropic_converged <- all(vapply(
    profile$fits, `[[`, logical(1), "converged"
  ))
  node_status <- transport_v3_node_status(profile$fits)
  structure(
    list(
      log_score = evidence$log_score,
      profile = profile,
      coverage_posterior = evidence$posterior,
      expected_coverage = evidence$expected_coverage,
      map_coverage = evidence$map_coverage,
      coupling = coupling,
      correspondence = correspondence,
      reference = reference,
      source = source,
      registered_source = registered_source,
      warp = warp,
      spec = spec,
      diagnostics = diagnostics,
      convergence = list(
        converged = entropic_converged &&
          (is.null(profile$polish) || profile$polish$all_feasible),
        method = "fixed_mass_mirror_descent_reference",
        backend = profile$backend,
        fallback = isTRUE(profile$fallback),
        fallback_reason = if (is.null(profile$fallback_reason)) {
          NA_character_
        } else {
          profile$fallback_reason
        },
        backend_route = if (is.null(profile$backend_route)) {
          NA_character_
        } else {
          profile$backend_route
        },
        revision = transport_v3_revision(spec),
        status = transport_v3_alignment_status(node_status),
        node_status = node_status,
        masses = lapply(profile$fits, `[[`, "history"),
        polish = profile$polish
      )
    ),
    class = c("gaze_transport_alignment", "list")
  )
}

#' Align gaze paths with edge-normalized Transport
#'
#' @param reference,source Fixation paths or prepared gaze measures.
#' @param spec A [gaze_transport_spec()].
#' @param warp_model Candidate-invariant fitted warp. A non-identity warp must
#'   be learned outside the scored pair.
#' @param candidate_key Candidate identifier.
#'
#' @return A symmetric `gaze_engine_result` with a coverage-integrated
#'   Transport alignment.
#' @export
gaze_transport_align <- function(reference, source, spec,
                                    warp_model = NULL,
                                    candidate_key = "candidate") {
  if (!inherits(spec, "gaze_transport_spec")) {
    stop("spec must be created by gaze_transport_spec().")
  }
  if (is.null(warp_model)) {
    if (!identical(spec$warp$type, "none")) {
      stop("A non-identity warp must be learned from independent trials.")
    }
    warp_model <- identity_gaze_warp_model()
  }
  reference_measure <- as_transport_v3_measure(reference, spec$chronology)
  source_measure <- as_transport_v3_measure(source, spec$chronology)
  registered_source <- apply_gaze_warp_model(source_measure, warp_model)
  profile <- solve_transport_v3_profile_symmetric(
    reference_measure, registered_source, spec
  )
  if (identical(spec$control$polish, "audit")) {
    profile <- polish_transport_v3_profile(
      profile, reference_measure, registered_source, spec
    )
  }
  alignment <- transport_v3_alignment_from_profile(
    profile,
    reference_measure,
    source_measure,
    registered_source,
    warp_model,
    spec
  )
  new_gaze_engine_result(
    engine = "transport",
    candidate_key = candidate_key,
    log_score = alignment$log_score,
    diagnostics = alignment$diagnostics,
    alignment = alignment,
    convergence = alignment$convergence,
    provenance = list(
      engine_version = 3L,
      solver_revision = transport_v3_revision(spec),
      directionality = "symmetric",
      score_semantics = "coverage_integrated_negative_scientific_energy",
      duration_semantics = "unit_duration_mass_with_fixed_ordinal_neighbours",
      candidate_invariant = TRUE,
      coverage_prior = spec$coverage_prior,
      coverage_nodes = spec$coverage$nodes,
      correspondence_semantics = "optimized_alignment_not_posterior"
    )
  )
}
