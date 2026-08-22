# Optimized and batched Transport backend --------------------------------

transport_v3_native_available <- function() {
  is.loaded("_eyesim_transport_v3_profile_native_cpp", PACKAGE = "eyesim")
}

solve_transport_v3_profile_native <- function(reference, source, spec,
                                              quadrature = spec$coverage) {
  if (!transport_v3_native_available()) {
    stop("Transport native backend is unavailable on this platform.")
  }
  if (spec$control$multistart != 1L) {
    stop("The native backend currently requires multistart = 1.")
  }
  if (identical(spec$control$projection_method, "log")) {
    stop("The native backend does not replace the log-domain projection oracle.")
  }
  spatial_cost <- gaze_spatial_cost(
    reference$coords, source$coords, spec$spatial
  )
  native <- transport_v3_profile_native_cpp(
    reference_mass = reference$mass,
    source_mass = source$mass,
    reference_relation = reference$relation,
    source_relation = source$relation,
    spatial_cost = spatial_cost,
    coverage_values = quadrature$coverage,
    entropy_schedule = spec$entropy_schedule,
    temporal_weight = spec$temporal_weight,
    selection_weight = spec$selection_weights[[1]],
    maxit = spec$control$maxit,
    step_size = spec$control$step_size,
    tolerance = spec$control$tolerance,
    projection_maxit = spec$control$projection_maxit,
    projection_tolerance = spec$control$projection_tolerance
  )
  fits <- native$fits
  if (length(fits) != length(quadrature$coverage) ||
      any(!vapply(fits, function(fit) isTRUE(fit$native_ok), logical(1)))) {
    stop("Transport native solver reported a numerical failure.")
  }
  final_entropy <- utils::tail(spec$entropy_schedule, 1)
  fits <- lapply(fits, function(fit) {
    # Re-evaluate every native plan with the readable R oracle. Candidate
    # evidence therefore never depends on a duplicated native formula.
    fit$objective <- transport_v3_objective(
      fit$coupling,
      reference,
      source,
      spatial_cost,
      spec,
      final_entropy
    )
    fit$start_objectives[[1]] <- fit$objective$optimization
    fit$start_scientific_objectives[[1]] <- fit$objective$scientific
    fit
  })
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
      backend = "native_rcpparmadillo",
      fallback = FALSE,
      native_provenance = list(
        kernel = "batched_fixed_mass_mirror_descent_and_standard_sinkhorn",
        scientific_recheck = "pure_R_oracle",
        candidate_warm_start = FALSE,
        coverage_warm_start = TRUE
      )
    ),
    class = c("gaze_transport_profile", "list")
  )
}

solve_transport_v3_profile_ordered <- function(reference, source, spec,
                                               quadrature = spec$coverage) {
  backend <- spec$control$backend
  if (identical(backend, "reference")) {
    return(solve_transport_v3_profile_reference(
      reference, source, spec, quadrature
    ))
  }
  native <- tryCatch(
    solve_transport_v3_profile_native(reference, source, spec, quadrature),
    error = function(condition) condition
  )
  native_success <- !inherits(native, "error") &&
    all(vapply(native$fits, function(fit) isTRUE(fit$converged), logical(1)))
  if (native_success) return(native)
  if (identical(backend, "optimized")) {
    if (inherits(native, "error")) stop(conditionMessage(native))
    stop("Transport optimized backend did not converge.")
  }
  reference <- solve_transport_v3_profile_reference(
    reference, source, spec, quadrature
  )
  reference$fallback <- TRUE
  reference$fallback_reason <- if (inherits(native, "error")) {
    conditionMessage(native)
  } else {
    "native_nonconvergence"
  }
  reference$backend <- "reference_fallback"
  reference
}

#' Prepare one immutable gaze path for repeated Transport scoring
#'
#' @param path A fixation path or prepared gaze measure.
#' @param spec A [gaze_transport_spec()].
#'
#' @return A prepared path with duration mass and local-edge terms.
#' @export
gaze_transport_prepare <- function(path, spec) {
  if (!inherits(spec, "gaze_transport_spec")) {
    stop("spec must be created by gaze_transport_spec().")
  }
  measure <- as_transport_v3_measure(path, spec$chronology)
  structure(
    list(
      measure = measure,
      fixation_count = nrow(measure$coords),
      effective_fixation_count = 1 / sum(measure$mass^2),
      duration_concentration = sum(measure$mass^2),
      immutable = TRUE,
      chronology = list(
        clock = spec$chronology$clock,
        neighbours = spec$chronology$neighbours,
        edge_count = sum(measure$relation > 0)
      )
    ),
    class = c("gaze_transport_prepared", "list")
  )
}

#' Batch candidate alignments with Transport
#'
#' @param references Named list of candidate fixation paths or prepared paths.
#' @param source One source fixation path or prepared path.
#' @param spec A [gaze_transport_spec()].
#' @param batch_size Positive number of candidates scheduled per batch.
#'
#' @return Named candidate `gaze_engine_result` objects in input order.
#' @export
gaze_transport_align_batch <- function(references, source, spec,
                                          batch_size = length(references)) {
  if (!is.list(references) || length(references) == 0L) {
    stop("references must contain at least one candidate path.")
  }
  candidate_keys <- names(references)
  if (is.null(candidate_keys) || any(!nzchar(candidate_keys)) ||
      anyDuplicated(candidate_keys)) {
    stop("references must have complete unique candidate names.")
  }
  batch_size <- as.integer(batch_size)
  if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L) {
    stop("batch_size must be a positive integer.")
  }
  prepared_source <- if (inherits(source, "gaze_transport_prepared")) {
    source
  } else {
    gaze_transport_prepare(source, spec)
  }
  prepared_references <- lapply(references, function(reference) {
    if (inherits(reference, "gaze_transport_prepared")) {
      reference
    } else {
      gaze_transport_prepare(reference, spec)
    }
  })
  batches <- split(
    seq_along(prepared_references),
    ceiling(seq_along(prepared_references) / batch_size)
  )
  results <- vector("list", length(prepared_references))
  for (batch in batches) {
    for (index in batch) {
      results[[index]] <- gaze_transport_align(
        prepared_references[[index]]$measure,
        prepared_source$measure,
        spec,
        candidate_key = candidate_keys[[index]]
      )
    }
  }
  stats::setNames(results, candidate_keys)
}
