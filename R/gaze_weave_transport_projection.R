# Fixed-mass transport projection primitives -------------------------------

gaze_mutual_information <- function(correspondence) {
  if (!is.matrix(correspondence) || any(!is.finite(correspondence)) ||
      any(correspondence < 0) || abs(sum(correspondence) - 1) > 1e-8) {
    stop("correspondence must be a finite non-negative unit-mass matrix.")
  }
  alpha <- rowSums(correspondence)
  beta <- colSums(correspondence)
  reference <- outer(alpha, beta)
  positive <- correspondence > 0
  sum(correspondence[positive] *
        log(correspondence[positive] / reference[positive]))
}

gaze_mutual_information_gradient <- function(correspondence) {
  alpha <- rowSums(correspondence)
  beta <- colSums(correspondence)
  safe <- pmax(correspondence, .Machine$double.xmin)
  safe_alpha <- pmax(alpha, .Machine$double.xmin)
  safe_beta <- pmax(beta, .Machine$double.xmin)
  log(safe) - outer(log(safe_alpha), log(safe_beta), FUN = "+")
}


transport_augmented_start <- function(reference_mass, source_mass, coverage) {
  n_reference <- length(reference_mass)
  n_source <- length(source_mass)
  augmented <- matrix(0, n_reference + 1L, n_source + 1L)
  augmented[seq_len(n_reference), seq_len(n_source)] <-
    coverage * outer(reference_mass, source_mass)
  augmented[seq_len(n_reference), n_source + 1L] <-
    (1 - coverage) * reference_mass
  augmented[n_reference + 1L, seq_len(n_source)] <-
    (1 - coverage) * source_mass
  augmented
}

masked_sinkhorn_projection_log <- function(kernel, row_target, column_target,
                                           mask, maxit, tolerance) {
  active_rows <- row_target > 0
  active_columns <- column_target > 0
  active_mask <- mask[active_rows, active_columns, drop = FALSE]
  if (any(rowSums(active_mask) == 0L) || any(colSums(active_mask) == 0L)) {
    stop("Masked projection has no support for a positive marginal target.")
  }
  active_kernel <- kernel[active_rows, active_columns, drop = FALSE]
  log_kernel <- matrix(-Inf, nrow(active_kernel), ncol(active_kernel))
  log_kernel[active_mask] <- log(pmax(
    active_kernel[active_mask], .Machine$double.xmin
  ))
  active_row_target <- row_target[active_rows]
  active_column_target <- column_target[active_columns]
  log_u <- rep(0, length(active_row_target))
  log_v <- rep(0, length(active_column_target))
  converged <- FALSE
  error <- Inf
  iteration <- 0L
  for (iteration in seq_len(maxit)) {
    log_u <- log(active_row_target) - row_log_sum_exp(
      sweep(log_kernel, 2L, log_v, FUN = "+")
    )
    log_v <- log(active_column_target) - row_log_sum_exp(
      sweep(t(log_kernel), 2L, log_u, FUN = "+")
    )
    if (iteration == 1L || iteration %% 10L == 0L || iteration == maxit) {
      active_plan <- exp(
        sweep(sweep(log_kernel, 1L, log_u, FUN = "+"), 2L, log_v, FUN = "+")
      )
      error <- max(
        abs(rowSums(active_plan) - active_row_target),
        abs(colSums(active_plan) - active_column_target)
      )
      if (is.finite(error) && error <= tolerance) {
        converged <- TRUE
        break
      }
    }
  }
  active_plan <- exp(
    sweep(sweep(log_kernel, 1L, log_u, FUN = "+"), 2L, log_v, FUN = "+")
  )
  plan <- matrix(0, nrow(kernel), ncol(kernel))
  plan[active_rows, active_columns] <- active_plan
  plan[!mask] <- 0
  list(
    plan = plan,
    converged = converged,
    error = error,
    iterations = iteration,
    method = "log"
  )
}

masked_sinkhorn_projection_standard <- function(
    kernel, row_target, column_target, mask, maxit, tolerance) {
  active_rows <- row_target > 0
  active_columns <- column_target > 0
  active_mask <- mask[active_rows, active_columns, drop = FALSE]
  if (any(rowSums(active_mask) == 0L) || any(colSums(active_mask) == 0L)) {
    stop("Masked projection has no support for a positive marginal target.")
  }
  active_kernel <- kernel[active_rows, active_columns, drop = FALSE]
  positive_kernel <- active_kernel[active_mask]
  kernel_scale <- max(positive_kernel)
  scaled_kernel <- active_kernel / kernel_scale
  scaled_kernel[!active_mask] <- 0
  scaled_positive <- scaled_kernel[active_mask]
  numerical_failure <-
    !is.finite(kernel_scale) || kernel_scale <= 0 ||
    any(!is.finite(scaled_positive)) || any(scaled_positive <= 0)
  if (numerical_failure) {
    return(list(
      plan = NULL, converged = FALSE, error = Inf, iterations = 0L,
      method = "standard", numerical_failure = TRUE
    ))
  }

  active_row_target <- row_target[active_rows]
  active_column_target <- column_target[active_columns]
  u <- rep(1, length(active_row_target))
  v <- rep(1, length(active_column_target))
  converged <- FALSE
  error <- Inf
  iteration <- 0L
  for (iteration in seq_len(maxit)) {
    kernel_v <- as.vector(scaled_kernel %*% v)
    if (any(!is.finite(kernel_v)) || any(kernel_v <= 0)) {
      numerical_failure <- TRUE
      break
    }
    u <- active_row_target / kernel_v
    kernel_u <- as.vector(crossprod(scaled_kernel, u))
    if (any(!is.finite(kernel_u)) || any(kernel_u <= 0)) {
      numerical_failure <- TRUE
      break
    }
    v <- active_column_target / kernel_u
    if (any(!is.finite(u)) || any(!is.finite(v))) {
      numerical_failure <- TRUE
      break
    }
    if (iteration == 1L || iteration %% 10L == 0L || iteration == maxit) {
      row_current <- u * as.vector(scaled_kernel %*% v)
      column_current <- v * as.vector(crossprod(scaled_kernel, u))
      error <- max(
        abs(row_current - active_row_target),
        abs(column_current - active_column_target)
      )
      if (is.finite(error) && error <= tolerance) {
        converged <- TRUE
        break
      }
    }
  }
  if (numerical_failure) {
    return(list(
      plan = NULL, converged = FALSE, error = Inf, iterations = iteration,
      method = "standard", numerical_failure = TRUE
    ))
  }
  active_plan <- scaled_kernel * tcrossprod(u, v)
  plan <- matrix(0, nrow(kernel), ncol(kernel))
  plan[active_rows, active_columns] <- active_plan
  plan[!mask] <- 0
  list(
    plan = plan,
    converged = converged,
    error = error,
    iterations = iteration,
    method = "standard",
    numerical_failure = FALSE
  )
}

masked_sinkhorn_projection <- function(kernel, row_target, column_target, mask,
                                       maxit, tolerance,
                                       method = c("auto", "standard", "log")) {
  method <- match.arg(method)
  if (identical(method, "log")) {
    return(masked_sinkhorn_projection_log(
      kernel, row_target, column_target, mask, maxit, tolerance
    ))
  }
  standard <- masked_sinkhorn_projection_standard(
    kernel, row_target, column_target, mask, maxit, tolerance
  )
  if (identical(method, "standard") || standard$converged) {
    return(standard)
  }
  fallback <- masked_sinkhorn_projection_log(
    kernel, row_target, column_target, mask, maxit, tolerance
  )
  fallback$fallback_from_standard <- TRUE
  fallback
}

# Damped Newton ascent on the dual of the masked KL projection. Used by the
# Transport revision 2026.10 solver as a finisher when standard Sinkhorn
# exhausts its iterations: when the plan's support splits into blocks joined
# only by small entries, Sinkhorn's slowest mode (the relative scaling of the
# blocks) contracts arbitrarily slowly, whereas Newton resolves it directly.
# The native backend implements the identical iteration.
masked_newton_projection <- function(kernel, row_target, column_target, mask,
                                     tolerance, maxit = 100L) {
  failure <- list(
    plan = NULL, converged = FALSE, error = Inf, iterations = 0L,
    method = "newton", numerical_failure = TRUE
  )
  if (any(is.na(kernel)) || any(kernel[mask] < 0) ||
      any(!is.finite(kernel[mask]))) {
    return(failure)
  }
  # Work with the log kernel so that entries spanning hundreds of orders of
  # magnitude neither overflow nor underflow.
  log_kernel <- suppressWarnings(log(kernel))
  log_kernel[!mask] <- -Inf
  finite <- is.finite(log_kernel)
  if (!any(finite)) return(failure)
  log_kernel <- log_kernel - max(log_kernel[finite])
  n_rows <- nrow(log_kernel)
  n_columns <- ncol(log_kernel)
  free <- seq_len(n_rows + n_columns - 1L)
  log_sum_exp <- function(values) {
    top <- max(values)
    if (!is.finite(top)) return(top)
    top + log(sum(exp(values - top)))
  }
  row_potential <- rep(0, n_rows)
  column_potential <- rep(0, n_columns)
  # Three log-domain Sinkhorn sweeps place the potentials near the solution.
  for (pass in seq_len(3L)) {
    row_potential <- log(row_target) - apply(
      sweep(log_kernel, 2L, column_potential, FUN = "+"), 1L, log_sum_exp
    )
    column_potential <- log(column_target) - apply(
      sweep(log_kernel, 1L, row_potential, FUN = "+"), 2L, log_sum_exp
    )
  }
  if (any(!is.finite(row_potential)) || any(!is.finite(column_potential))) {
    return(failure)
  }
  plan_at <- function(row_potential, column_potential) {
    exp(log_kernel + outer(row_potential, column_potential, FUN = "+"))
  }
  dual_at <- function(row_potential, column_potential, plan) {
    sum(row_target * row_potential) + sum(column_target * column_potential) -
      sum(plan)
  }
  plan <- plan_at(row_potential, column_potential)
  error <- Inf
  iteration <- 0L
  for (iteration in seq_len(maxit)) {
    row_sums <- rowSums(plan)
    column_sums <- colSums(plan)
    gradient <- c(row_target - row_sums, column_target - column_sums)
    error <- max(abs(gradient))
    if (!is.finite(error)) return(failure)
    if (error <= tolerance) break
    hessian <- rbind(
      cbind(diag(row_sums, n_rows), plan),
      cbind(t(plan), diag(column_sums, n_columns))
    )[free, free, drop = FALSE]
    diag(hessian) <- diag(hessian) + 1e-14 * max(diag(hessian))
    direction <- tryCatch(
      solve(hessian, gradient[free]),
      error = function(condition) NULL
    )
    if (is.null(direction) || any(!is.finite(direction))) return(failure)
    direction <- c(direction, 0)
    row_direction <- direction[seq_len(n_rows)]
    column_direction <- direction[n_rows + seq_len(n_columns)]
    current_dual <- dual_at(row_potential, column_potential, plan)
    slope <- sum(gradient * direction)
    step <- 1
    accepted <- FALSE
    while (step >= 1e-12) {
      trial_rows <- row_potential + step * row_direction
      trial_columns <- column_potential + step * column_direction
      trial_plan <- plan_at(trial_rows, trial_columns)
      trial_dual <- dual_at(trial_rows, trial_columns, trial_plan)
      # Armijo ascent on the dual, or (once the dual increase falls below
      # rounding) a decrease of the marginal error.
      trial_error <- max(abs(c(
        row_target - rowSums(trial_plan),
        column_target - colSums(trial_plan)
      )))
      if (is.finite(trial_dual) &&
          (trial_dual >= current_dual + 1e-4 * step * slope ||
           trial_error < error)) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) break
    row_potential <- trial_rows
    column_potential <- trial_columns
    plan <- trial_plan
  }
  error <- max(abs(c(row_target - rowSums(plan), column_target - colSums(plan))))
  if (!is.finite(error)) return(failure)
  plan[!mask] <- 0
  list(
    plan = plan,
    converged = is.finite(error) && error <= tolerance,
    error = error,
    iterations = iteration,
    method = "newton",
    numerical_failure = !is.finite(error)
  )
}

# Revision 2026.10 fixed-mass projection: standard Sinkhorn, then the Newton
# finisher when Sinkhorn exhausts its iterations, then (for the "auto"
# method only) the log-domain oracle.
project_partial_coupling_revised <- function(augmented_kernel, reference_mass,
                                             source_mass, coverage, control) {
  method <- if (is.null(control$projection_method)) {
    "auto"
  } else {
    control$projection_method
  }
  if (identical(method, "log")) {
    return(project_partial_coupling(
      augmented_kernel, reference_mass, source_mass, coverage, control
    ))
  }
  n_reference <- length(reference_mass)
  n_source <- length(source_mass)
  row_target <- c(reference_mass, 1 - coverage)
  column_target <- c(source_mass, 1 - coverage)
  mask <- matrix(TRUE, n_reference + 1L, n_source + 1L)
  mask[n_reference + 1L, n_source + 1L] <- FALSE
  standard <- masked_sinkhorn_projection_standard(
    augmented_kernel, row_target, column_target, mask,
    control$projection_maxit, control$projection_tolerance
  )
  if (standard$converged) return(standard)
  newton <- masked_newton_projection(
    augmented_kernel, row_target, column_target, mask,
    control$projection_tolerance
  )
  newton$fallback_from_standard <- TRUE
  if (newton$converged || identical(method, "standard")) return(newton)
  fallback <- masked_sinkhorn_projection_log(
    augmented_kernel, row_target, column_target, mask,
    control$projection_maxit, control$projection_tolerance
  )
  fallback$fallback_from_standard <- TRUE
  fallback
}

project_partial_coupling <- function(augmented_kernel, reference_mass,
                                     source_mass, coverage, control) {
  n_reference <- length(reference_mass)
  n_source <- length(source_mass)
  row_target <- c(reference_mass, 1 - coverage)
  column_target <- c(source_mass, 1 - coverage)
  mask <- matrix(TRUE, n_reference + 1L, n_source + 1L)
  mask[n_reference + 1L, n_source + 1L] <- FALSE
  masked_sinkhorn_projection(
    augmented_kernel,
    row_target,
    column_target,
    mask,
    maxit = control$projection_maxit,
    tolerance = control$projection_tolerance,
    method = if (is.null(control$projection_method)) {
      "auto"
    } else {
      control$projection_method
    }
  )
}

transport_initial_plans <- function(reference, source, spatial_cost,
                                       coverage, spec) {
  base <- transport_augmented_start(reference$mass, source$mass, coverage)
  starts <- list(independent = base)
  if (spec$control$multistart >= 2L) {
    n_reference <- length(reference$mass)
    n_source <- length(source$mass)
    spatial <- base
    scale <- max(stats::median(spatial_cost), 0.05)
    spatial[seq_len(n_reference), seq_len(n_source)] <-
      pmax(base[seq_len(n_reference), seq_len(n_source)], 1e-300) *
      exp(-pmin(spatial_cost / scale, 700))
    project <- if (transport_v3_revised(spec)) {
      project_partial_coupling_revised
    } else {
      project_partial_coupling
    }
    projected <- project(
      spatial, reference$mass, source$mass, coverage, spec$control
    )
    if (projected$converged) {
      starts$spatial <- projected$plan
    }
  }
  starts
}
