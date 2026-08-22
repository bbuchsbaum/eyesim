# Transport specification ------------------------------------------------

gauss_legendre_rule <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("n must be an integer of at least two.")
  }
  index <- seq_len(n - 1L)
  off_diagonal <- index / sqrt(4 * index^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(index, index + 1L)] <- off_diagonal
  jacobi[cbind(index + 1L, index)] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  order_index <- order(decomposition$values)
  list(
    nodes = decomposition$values[order_index],
    weights = 2 * decomposition$vectors[1L, order_index]^2
  )
}

transport_v3_coverage_quadrature <- function(n = 8L,
                                             prior_shape = c(2, 2)) {
  if (!is.numeric(prior_shape) || length(prior_shape) != 2L ||
      any(!is.finite(prior_shape)) || any(prior_shape <= 0)) {
    stop("prior_shape must contain two finite positive values.")
  }
  rule <- gauss_legendre_rule(n)
  coverage <- (rule$nodes + 1) / 2
  integration_weight <- rule$weights / 2
  log_prior_density <- stats::dbeta(
    coverage, prior_shape[[1]], prior_shape[[2]], log = TRUE
  )
  log_weight <- log(integration_weight) + log_prior_density
  log_weight <- log_weight - gaze_log_sum_exp(log_weight)
  structure(
    list(
      coverage = coverage,
      weight = exp(log_weight),
      log_weight = log_weight,
      integration_weight = integration_weight,
      prior_shape = as.numeric(prior_shape),
      nodes = as.integer(n),
      method = "gauss_legendre_beta_prior"
    ),
    class = c("gaze_coverage_quadrature", "list")
  )
}

#' Specify edge-normalized GazeWeave Transport
#'
#' Transport separates matched coverage, normalized correspondence,
#' selected source and target mass, spatial fidelity, directed local order,
#' correspondence smoothing, and candidate-invariant registration. Coverage
#' is integrated by Gauss-Legendre quadrature under a fixed beta prior.
#'
#' @param spatial A [gaze_gaussian_mixture()] specification.
#' @param chronology Ordinal chronology from [gaze_order_neighbours()].
#' @param coverage_prior Two positive beta-prior shape parameters.
#' @param coverage_nodes Number of default coverage quadrature nodes.
#' @param temporal_weight Non-negative directed local-order weight.
#' @param selection_weights Non-negative reference and source selection weights.
#' @param entropy_schedule Positive correspondence-smoothing continuation values.
#' @param temperature_bounds Calibration temperature bounds.
#' @param warp Candidate-invariant cross-fitted warp specification.
#' @param screen Optional screen geometry.
#' @param maxit,step_size,tolerance Reference solver controls.
#' @param projection_maxit,projection_tolerance,projection_method Fixed-mass
#'   projection controls.
#' @param multistart One or two common structural starts.
#' @param backend Pair solver backend. `"reference"` is the readable R oracle;
#'   `"optimized"` uses the estimator-preserving batched implementation;
#'   `"auto"` uses the optimized backend with a clean reference fallback.
#' @param polish Optional scientific-objective Frank-Wolfe audit. The default
#'   keeps the entropic solution; `"audit"` polishes every coverage node.
#' @param polish_maxit,polish_gap_tolerance,polish_relative_tolerance
#'   Conditional-gradient stopping controls.
#' @param reliability Response-blind calibration policy. The default learns
#'   effective-fixation shrinkage on inner out-of-fold predictions and includes
#'   the temperature-only solution as an exact boundary.
#' @param reliability_kappa_bounds Non-negative bounds for the shrinkage scale.
#' @param calibration_folds,calibration_seed Inner calibration-fold policy.
#' @param log_temperature_prior_sd,log1p_kappa_prior_sd Calibration penalties.
#'
#' @return A frozen `gaze_transport_spec`.
#' @export
gaze_transport_spec <- function(
    spatial = gaze_gaussian_mixture(c(0.75, 1.5), weights = c(0.7, 0.3)),
    chronology = gaze_order_neighbours(neighbours = 2),
    coverage_prior = c(2, 2), coverage_nodes = 12L,
    temporal_weight = 2, selection_weights = c(0.5, 0.5),
    entropy_schedule = c(0.05, 0.015),
    temperature_bounds = c(0.05, 20),
    warp = gaze_warp_none(), screen = NULL,
    maxit = 1000L, step_size = 2, tolerance = 5e-5,
    projection_maxit = 1000L, projection_tolerance = 1e-8,
    projection_method = c("auto", "standard", "log"),
    multistart = 1L, backend = c("auto", "optimized", "reference"),
    polish = c("none", "audit"), polish_maxit = 200L,
    polish_gap_tolerance = 1e-7, polish_relative_tolerance = 1e-8,
    reliability = c("effective_fixations", "none"),
    reliability_kappa_bounds = c(0, 100),
    calibration_folds = 2L, calibration_seed = 20260822L,
    log_temperature_prior_sd = 1, log1p_kappa_prior_sd = 1) {
  if (!inherits(spatial, "gaze_spatial_spec")) {
    stop("spatial must be created by gaze_gaussian_mixture().")
  }
  if (!inherits(chronology, "gaze_chronology_spec") ||
      !identical(chronology$type, "order_neighbours")) {
    stop("Transport chronology must use gaze_order_neighbours().")
  }
  if (!is.numeric(coverage_prior) || length(coverage_prior) != 2L ||
      any(!is.finite(coverage_prior)) || any(coverage_prior <= 0)) {
    stop("coverage_prior must contain two finite positive beta shapes.")
  }
  coverage_nodes <- as.integer(coverage_nodes)
  if (length(coverage_nodes) != 1L || is.na(coverage_nodes) ||
      coverage_nodes < 2L) {
    stop("coverage_nodes must be an integer of at least two.")
  }
  if (!is.numeric(temporal_weight) || length(temporal_weight) != 1L ||
      !is.finite(temporal_weight) || temporal_weight < 0) {
    stop("temporal_weight must be one finite non-negative value.")
  }
  if (!is.numeric(selection_weights) || length(selection_weights) != 2L ||
      any(!is.finite(selection_weights)) || any(selection_weights < 0)) {
    stop("selection_weights must contain two finite non-negative values.")
  }
  if (abs(selection_weights[[1]] - selection_weights[[2]]) > 1e-14) {
    stop("Symmetric Transport requires equal reference and source selection weights.")
  }
  if (!is.numeric(entropy_schedule) || length(entropy_schedule) == 0L ||
      any(!is.finite(entropy_schedule)) || any(entropy_schedule <= 0)) {
    stop("entropy_schedule must contain finite positive values.")
  }
  if (!is.numeric(temperature_bounds) || length(temperature_bounds) != 2L ||
      any(!is.finite(temperature_bounds)) || any(temperature_bounds <= 0) ||
      temperature_bounds[[1]] >= temperature_bounds[[2]]) {
    stop("temperature_bounds must be increasing finite positive values.")
  }
  if (!inherits(warp, "gaze_warp_spec")) {
    stop("warp must be a GazeWeave warp specification.")
  }
  if (!is.null(screen) && !inherits(screen, "gaze_screen")) {
    stop("screen must be NULL or created by gaze_screen().")
  }
  if (identical(warp$type, "contraction") &&
      is.character(warp$center) && is.null(screen)) {
    stop("screen-centred contraction requires a screen specification.")
  }
  integer_controls <- list(
    maxit = as.integer(maxit),
    projection_maxit = as.integer(projection_maxit),
    multistart = as.integer(multistart)
  )
  if (is.na(integer_controls$maxit) || integer_controls$maxit < 1L ||
      is.na(integer_controls$projection_maxit) ||
      integer_controls$projection_maxit < 10L ||
      is.na(integer_controls$multistart) ||
      !integer_controls$multistart %in% 1:2) {
    stop("Invalid Transport integer solver controls.")
  }
  numeric_controls <- c(
    step_size = step_size,
    tolerance = tolerance,
    projection_tolerance = projection_tolerance
  )
  if (any(!is.finite(numeric_controls)) || any(numeric_controls <= 0)) {
    stop("Transport numeric solver controls must be finite and positive.")
  }

  projection_method <- match.arg(projection_method)
  backend <- match.arg(backend)
  polish <- match.arg(polish)
  reliability <- match.arg(reliability)
  polish_maxit <- as.integer(polish_maxit)
  if (length(polish_maxit) != 1L || is.na(polish_maxit) || polish_maxit < 1L) {
    stop("polish_maxit must be a positive integer.")
  }
  polish_tolerances <- c(
    gap = polish_gap_tolerance,
    relative = polish_relative_tolerance
  )
  if (any(!is.finite(polish_tolerances)) || any(polish_tolerances <= 0)) {
    stop("polish tolerances must be finite and positive.")
  }
  if (!is.numeric(reliability_kappa_bounds) ||
      length(reliability_kappa_bounds) != 2L ||
      any(!is.finite(reliability_kappa_bounds)) ||
      reliability_kappa_bounds[[1L]] != 0 ||
      reliability_kappa_bounds[[2L]] <= 0) {
    stop("reliability_kappa_bounds must run from zero to a positive value.")
  }
  calibration_folds <- as.integer(calibration_folds)
  calibration_seed <- as.integer(calibration_seed)
  if (length(calibration_folds) != 1L || is.na(calibration_folds) ||
      calibration_folds < 2L || length(calibration_seed) != 1L ||
      is.na(calibration_seed)) {
    stop("calibration_folds and calibration_seed must be valid scalar integers.")
  }
  calibration_prior_sd <- c(
    temperature = log_temperature_prior_sd,
    reliability = log1p_kappa_prior_sd
  )
  if (any(is.na(calibration_prior_sd)) ||
      any(calibration_prior_sd <= 0)) {
    stop("calibration prior standard deviations must be positive.")
  }
  quadrature <- transport_v3_coverage_quadrature(
    coverage_nodes, coverage_prior
  )
  structure(
    list(
      spatial = spatial,
      chronology = chronology,
      coverage_prior = as.numeric(coverage_prior),
      coverage = quadrature,
      temporal_weight = as.numeric(temporal_weight),
      selection_weights = as.numeric(selection_weights),
      entropy_schedule = as.numeric(entropy_schedule),
      temperature_bounds = as.numeric(temperature_bounds),
      reliability = reliability,
      reliability_kappa_bounds = as.numeric(reliability_kappa_bounds),
      calibration = list(
        folds = calibration_folds,
        seed = calibration_seed,
        log_temperature_prior_sd = as.numeric(log_temperature_prior_sd),
        log1p_kappa_prior_sd = as.numeric(log1p_kappa_prior_sd)
      ),
      warp = warp,
      screen = screen,
      control = c(
        integer_controls,
        as.list(numeric_controls),
        list(
          projection_method = projection_method,
          backend = backend,
          polish = polish,
          polish_maxit = polish_maxit,
          polish_gap_tolerance = as.numeric(polish_gap_tolerance),
          polish_relative_tolerance = as.numeric(polish_relative_tolerance)
        )
      ),
      revision = "edge_normalized_episode_transport",
      estimand = "edge_normalized_episode_transport"
    ),
    class = c("gaze_transport_spec", "list")
  )
}

#' @export
print.gaze_transport_spec <- function(x, ...) {
  cat("GazeWeave Transport specification\n")
  cat("  coverage: ", x$coverage$nodes, " Gauss-Legendre nodes; Beta(",
      paste(x$coverage_prior, collapse = ", "), ") prior\n", sep = "")
  cat("  chronology: next", x$chronology$neighbours, "ordinal neighbours\n")
  cat("  backend:", x$control$backend, "\n")
  invisible(x)
}
