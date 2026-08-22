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
    gradient = gradient
  )
  reference_selection <- transport_v3_jensen_shannon(
    reference_selected_unit, reference$mass
  )
  source_selection <- transport_v3_jensen_shannon(
    source_selected_unit, source$mass
  )
  correspondence_information <- gaze_mutual_information(correspondence)
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

solve_transport_v3_mass_reference <- function(reference, source, spatial_cost,
                                              coverage, spec,
                                              initial_augmented = NULL) {
  n_reference <- length(reference$mass)
  n_source <- length(source$mass)
  real_rows <- seq_len(n_reference)
  real_columns <- seq_len(n_source)
  starts <- transport_initial_plans(
    reference, source, spatial_cost, coverage, spec
  )
  if (!is.null(initial_augmented)) {
    projected <- project_partial_coupling(
      initial_augmented,
      reference$mass,
      source$mass,
      coverage,
      spec$control
    )
    if (projected$converged) {
      # Coverage continuation is within one candidate and follows the same
      # ascending-node policy for every candidate. It never carries state
      # across candidates.
      starts <- list(adjacent_coverage = projected$plan)
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
    for (entropy in spec$entropy_schedule) {
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
    list(
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
  })
  optimization <- vapply(fits, function(fit) {
    fit$objective$optimization
  }, numeric(1))
  scientific <- vapply(fits, function(fit) {
    fit$objective$scientific
  }, numeric(1))
  converged <- vapply(fits, function(fit) isTRUE(fit$converged), logical(1))
  if (!any(converged) && !is.null(initial_augmented)) {
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
  best <- fits[[which.min(optimization)]]
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
      maximum_gap = max(polish_gaps),
      map_improvement = map_fit$polish$improvement,
      total_improvement = sum(polish_improvements),
      maximum_feasibility_error = max(polish_feasibility),
      no_worse = profile$polish$no_worse,
      all_feasible = profile$polish$all_feasible,
      all_converged = profile$polish$all_converged,
      stopping_reason = map_fit$polish$stopping_reason,
      trace = map_fit$polish$trace
    )
  }
  entropic_converged <- all(vapply(
    profile$fits, `[[`, logical(1), "converged"
  ))
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
