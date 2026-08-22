# Conditional-gradient Transport audit ----------------------------------

transport_v3_feasibility <- function(coupling, reference_mass, source_mass,
                                     coverage) {
  c(
    coverage = abs(sum(coupling) - coverage),
    reference = max(c(0, rowSums(coupling) - reference_mass)),
    source = max(c(0, colSums(coupling) - source_mass)),
    nonnegative = max(c(0, -coupling))
  )
}

transport_v3_linear_oracle <- function(gradient, reference_mass, source_mass,
                                       coverage) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Transport polish requires the optional lpSolve package.")
  }
  n_reference <- length(reference_mass)
  n_source <- length(source_mass)
  n_variables <- n_reference * n_source
  constraints <- matrix(0, n_reference + n_source + 1L, n_variables)
  for (reference_index in seq_len(n_reference)) {
    indices <- reference_index + (seq_len(n_source) - 1L) * n_reference
    constraints[reference_index, indices] <- 1
  }
  for (source_index in seq_len(n_source)) {
    indices <- seq_len(n_reference) + (source_index - 1L) * n_reference
    constraints[n_reference + source_index, indices] <- 1
  }
  constraints[n_reference + n_source + 1L, ] <- 1
  solution <- lpSolve::lp(
    direction = "min",
    objective.in = as.vector(gradient),
    const.mat = constraints,
    const.dir = c(
      rep("<=", n_reference + n_source), "="
    ),
    const.rhs = c(reference_mass, source_mass, coverage),
    transpose.constraints = TRUE
  )
  if (solution$status != 0L) {
    stop("Transport linear minimization oracle failed with status ",
         solution$status, ".")
  }
  matrix(solution$solution, n_reference, n_source)
}

transport_v3_line_search <- function(coupling, direction, reference, source,
                                     spatial_cost, spec) {
  objective <- function(step) {
    transport_v3_objective(
      coupling + step * direction,
      reference,
      source,
      spatial_cost,
      spec,
      entropy = 0
    )$scientific
  }
  grid <- seq(0, 1, length.out = 21L)
  values <- vapply(grid, objective, numeric(1))
  best <- which.min(values)
  candidates <- data.frame(step = grid, value = values)
  lower <- grid[[max(1L, best - 1L)]]
  upper <- grid[[min(length(grid), best + 1L)]]
  refined <- stats::optimize(
    objective,
    interval = c(lower, upper),
    tol = 1e-10
  )
  candidates <- rbind(
    candidates,
    data.frame(step = refined$minimum, value = refined$objective)
  )
  candidates[which.min(candidates$value), , drop = FALSE]
}

polish_transport_v3_mass <- function(coupling, reference, source,
                                     spatial_cost, spec,
                                     maxit = spec$control$polish_maxit,
                                     gap_tolerance =
                                       spec$control$polish_gap_tolerance,
                                     relative_tolerance =
                                       spec$control$polish_relative_tolerance) {
  coverage <- sum(coupling)
  current <- transport_v3_objective(
    coupling,
    reference,
    source,
    spatial_cost,
    spec,
    entropy = 0,
    gradient = TRUE
  )
  initial_scientific <- current$scientific
  trace <- vector("list", maxit)
  converged <- FALSE
  stopping_reason <- "maximum_iterations"
  final_iteration <- 0L

  for (iteration in seq_len(maxit)) {
    vertex <- transport_v3_linear_oracle(
      current$gradient, reference$mass, source$mass, coverage
    )
    direction <- vertex - coupling
    gap <- sum(current$gradient * (coupling - vertex))
    feasibility <- transport_v3_feasibility(
      coupling, reference$mass, source$mass, coverage
    )
    trace[[iteration]] <- data.frame(
      iteration = iteration - 1L,
      scientific = current$scientific,
      gap = gap,
      step = 0,
      relative_improvement = 0,
      coverage_error = feasibility[["coverage"]],
      reference_error = feasibility[["reference"]],
      source_error = feasibility[["source"]],
      nonnegative_error = feasibility[["nonnegative"]],
      stringsAsFactors = FALSE
    )
    final_iteration <- iteration - 1L
    if (gap <= gap_tolerance) {
      converged <- TRUE
      stopping_reason <- "conditional_gradient_gap"
      break
    }
    line <- transport_v3_line_search(
      coupling,
      direction,
      reference,
      source,
      spatial_cost,
      spec
    )
    proposal <- coupling + line$step[[1]] * direction
    proposal_objective <- transport_v3_objective(
      proposal,
      reference,
      source,
      spatial_cost,
      spec,
      entropy = 0,
      gradient = TRUE
    )
    relative_improvement <- (
      current$scientific - proposal_objective$scientific
    ) / max(1, abs(current$scientific))
    if (proposal_objective$scientific > current$scientific + 1e-10) {
      stop("Transport polish line search worsened the scientific objective.")
    }
    trace[[iteration]]$step <- line$step[[1]]
    trace[[iteration]]$relative_improvement <- relative_improvement
    coupling <- proposal
    current <- proposal_objective
    final_iteration <- iteration
    if (relative_improvement <= relative_tolerance) {
      converged <- TRUE
      stopping_reason <- "relative_scientific_energy"
      break
    }
  }
  trace <- do.call(rbind, Filter(Negate(is.null), trace))
  feasibility <- transport_v3_feasibility(
    coupling, reference$mass, source$mass, coverage
  )
  final_vertex <- transport_v3_linear_oracle(
    current$gradient, reference$mass, source$mass, coverage
  )
  final_gap <- sum(current$gradient * (coupling - final_vertex))
  list(
    coupling = coupling,
    objective = transport_v3_objective(
      coupling, reference, source, spatial_cost, spec, entropy = 0
    ),
    initial_scientific = initial_scientific,
    improvement = initial_scientific - current$scientific,
    gap = final_gap,
    feasibility = feasibility,
    trace = trace,
    converged = converged,
    stopping_reason = stopping_reason,
    iterations = final_iteration,
    oracle = "lpSolve_exact_transport_linear_program"
  )
}

polish_transport_v3_profile <- function(profile, reference, source, spec) {
  polished <- lapply(profile$fits, function(fit) {
    result <- polish_transport_v3_mass(
      fit$coupling,
      reference,
      source,
      profile$spatial_cost,
      spec
    )
    fit$entropic_objective <- fit$objective
    fit$entropic_coupling <- fit$coupling
    fit$coupling <- result$coupling
    fit$objective <- result$objective
    fit$polish <- result
    fit
  })
  profile$entropic_scientific_energy <- profile$scientific_energy
  profile$scientific_energy <- vapply(polished, function(fit) {
    fit$objective$scientific
  }, numeric(1))
  profile$fits <- polished
  profile$polish <- list(
    enabled = TRUE,
    oracle = "lpSolve_exact_transport_linear_program",
    no_worse = all(profile$scientific_energy <=
      profile$entropic_scientific_energy + 1e-10),
    all_feasible = all(vapply(polished, function(fit) {
      max(fit$polish$feasibility) <= 1e-8
    }, logical(1))),
    all_converged = all(vapply(polished, function(fit) {
      fit$polish$converged
    }, logical(1)))
  )
  profile$backend <- paste0(profile$backend, "+scientific_polish")
  profile
}
