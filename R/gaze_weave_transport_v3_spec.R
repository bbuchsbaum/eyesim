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

# Solver revision of a Transport specification. Specifications serialised
# before revisions existed carry no revision or, from 2026-08 onwards, the
# estimand name in the `revision` slot; both are solved with the frozen
# 2026.08 rules. Any other value is rejected.
transport_v3_revision <- function(spec) {
  revision <- spec$revision
  if (is.null(revision) ||
      identical(revision, "edge_normalized_episode_transport")) {
    return("2026.08")
  }
  if (identical(revision, "2026.08") || identical(revision, "2026.10")) {
    return(revision)
  }
  stop(
    "Unknown Transport solver revision: ",
    paste(format(revision), collapse = ", "),
    ". Supported revisions are \"2026.10\" and \"2026.08\"."
  )
}

# Default revision 2026.10 stopping tolerance (predicted remaining decrease,
# nats). On the 60-pair review set, a tight continuation of certified nodes
# found a further decrease of at most 1e-6 at the 95th percentile, but more
# than 1e-5 on 9 of 692 nodes (up to 1.4e-3, slow directions and saddles
# that a first-order rule cannot detect). Tightening to 1e-7..1e-9 did not
# remove those nodes and multiplied projection-limited stalls.
transport_v3_default_tolerance <- 1e-6

transport_v3_revised <- function(spec) {
  identical(transport_v3_revision(spec), "2026.10")
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
#' @param entropy_schedule Positive correspondence-smoothing continuation
#'   values; the last value defines the optimized objective. `NULL` (the
#'   default) uses `c(0.15, 0.05, 0.015)` under revision `"2026.10"` and
#'   `c(0.05, 0.015)` under `"2026.08"`. The extra, smoother first stage
#'   reduces, but does not remove, the dependence of the local optimum
#'   reached on `step_size`.
#' @param temperature_bounds Calibration temperature bounds.
#' @param warp Candidate-invariant cross-fitted warp specification.
#' @param screen Optional screen geometry.
#' @param maxit,step_size Mirror-descent iteration limit and largest step.
#' @param tolerance Stopping tolerance. Under revision `"2026.10"` its unit is
#'   the regularized objective (nats): a stage stops when a first-order model
#'   predicts at most `tolerance` of remaining decrease at the current plan
#'   (see `revision`). `NULL` (the default) uses `1e-6` under `"2026.10"` and
#'   the relative objective change `5e-5` under `"2026.08"`.
#' @param projection_maxit,projection_tolerance,projection_method Fixed-mass
#'   projection controls.
#' @param multistart One or two common structural starts. Under revision
#'   `"2026.10"` both backends always add the adjacent-coverage continuation,
#'   and `2` adds a spatial start (natively supported). On the review set,
#'   `2` lowered the 95th-percentile gap to a multistart oracle only from
#'   0.0056 to 0.0053 nats at 1.48 times the runtime, and on a second review
#'   set it never moved a score by more than 2e-5, so it is opt-in.
#' @param backend Pair solver backend. `"reference"` is the readable R oracle;
#'   `"optimized"` uses the estimator-preserving batched implementation;
#'   `"auto"` uses the optimized backend with a reference fallback. Under
#'   revision `"2026.10"`, a specification the native backend cannot solve
#'   (`projection_method = "log"`) is routed to the reference backend for
#'   every pair, a node that reaches `maxit` is recorded as
#'   `"not_converged"` and scored rather than raising an error, and a
#'   per-pair numerical fallback raises a warning of class
#'   `gaze_transport_backend_fallback`. The backend used is recorded in every
#'   alignment's `convergence` element.
#' @param polish Optional scientific-objective Frank-Wolfe audit. The default
#'   keeps the entropic solution; `"audit"` polishes every coverage node.
#' @param polish_maxit,polish_gap_tolerance,polish_relative_tolerance
#'   Conditional-gradient stopping controls.
#' @param revision Solver revision. `"2026.10"` (the default) is the
#'   corrected solver and estimand. The chronology residual is
#'   `1 - 2A / (R + S + 1e-3)`, continuous as the selected edge mass
#'   vanishes. The stopping rule is first-order: a stage stops when the plan
#'   is feasible, no cell outside the support is trapped, and a first-order
#'   model predicts at most `tolerance` of remaining decrease. It does not
#'   certify a local optimum: slow directions and saddles can remain (on the
#'   review set about 1% of stopped nodes could still be lowered by 1e-5 to
#'   1.4e-3). Trapped cells are reseeded. Mirror steps are capped so the
#'   exponent never saturates. A line search limited by projection noise
#'   records the stage as `"stalled_projection_limited"` and `maxit` as
#'   `"not_converged"`; both are scored. Standard Sinkhorn is finished by a
#'   dual Newton projection and then log-domain Sinkhorn. Every coverage node
#'   is solved from the independent and the continuation start, replacing
#'   the reference backend's silent cold restart. The result is a local
#'   optimum of a non-convex objective. Scores fall short of a multistart
#'   oracle (one-sided, mostly in null and sparse-source pairs): by up to
#'   0.12 nats with 95th percentile 0.0056 on one 60-pair review set, and
#'   up to 0.117 with 95th percentile 0.023 (null pairs 0.033, matched 0.006)
#'   on a second 88-pair, half-null set. `backend = "auto"` warns about and records
#'   every numerical fallback. The polish Frank-Wolfe gap is `NA` with a
#'   reason on a selection boundary. `"2026.08"` reproduces the frozen August
#'   2026 solver exactly; the frozen validation courts pin it.
#'   Specifications saved before revisions existed are solved as
#'   `"2026.08"`.
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
    entropy_schedule = NULL,
    temperature_bounds = c(0.05, 20),
    warp = gaze_warp_none(), screen = NULL,
    maxit = 1000L, step_size = 2, tolerance = NULL,
    projection_maxit = 1000L, projection_tolerance = 1e-8,
    projection_method = c("auto", "standard", "log"),
    multistart = 1L, backend = c("auto", "optimized", "reference"),
    polish = c("none", "audit"), polish_maxit = 200L,
    polish_gap_tolerance = 1e-7, polish_relative_tolerance = 1e-8,
    reliability = c("effective_fixations", "none"),
    reliability_kappa_bounds = c(0, 100),
    calibration_folds = 2L, calibration_seed = 20260822L,
    log_temperature_prior_sd = 1, log1p_kappa_prior_sd = 1,
    revision = c("2026.10", "2026.08")) {
  revision <- match.arg(revision)
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
  if (is.null(entropy_schedule)) {
    entropy_schedule <- if (identical(revision, "2026.10")) {
      c(0.15, 0.05, 0.015)
    } else {
      c(0.05, 0.015)
    }
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
  if (is.null(tolerance)) {
    tolerance <- if (identical(revision, "2026.10")) {
      transport_v3_default_tolerance
    } else {
      5e-5
    }
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
      revision = revision,
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
  cat("  solver revision:", transport_v3_revision(x), "\n")
  invisible(x)
}
