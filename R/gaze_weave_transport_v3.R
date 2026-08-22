# Edge-normalized GazeWeave Transport ------------------------------------

validate_transport_v3_edge_inputs <- function(correspondence,
                                               reference_relation,
                                               source_relation) {
  if (!is.matrix(correspondence) || !is.numeric(correspondence) ||
      any(!is.finite(correspondence)) || any(correspondence < 0)) {
    stop("correspondence must be a finite non-negative numeric matrix.")
  }
  if (!is.matrix(reference_relation) || !is.numeric(reference_relation) ||
      nrow(reference_relation) != ncol(reference_relation) ||
      nrow(reference_relation) != nrow(correspondence) ||
      any(!is.finite(reference_relation)) || any(reference_relation < 0)) {
    stop("reference_relation must be a compatible finite non-negative square matrix.")
  }
  if (!is.matrix(source_relation) || !is.numeric(source_relation) ||
      nrow(source_relation) != ncol(source_relation) ||
      nrow(source_relation) != ncol(correspondence) ||
      any(!is.finite(source_relation)) || any(source_relation < 0)) {
    stop("source_relation must be a compatible finite non-negative square matrix.")
  }
  mass <- sum(correspondence)
  if (!is.finite(mass) || mass <= 0 || abs(mass - 1) > 1e-8) {
    stop("correspondence must have unit mass.")
  }
  invisible(TRUE)
}

transport_v3_edge_terms <- function(correspondence, reference_relation,
                                    source_relation, gradient = FALSE,
                                    edge_tolerance = 1e-15) {
  validate_transport_v3_edge_inputs(
    correspondence, reference_relation, source_relation
  )
  if (length(edge_tolerance) != 1L || !is.finite(edge_tolerance) ||
      edge_tolerance < 0) {
    stop("edge_tolerance must be one finite non-negative value.")
  }

  reference_selected <- rowSums(correspondence)
  source_selected <- colSums(correspondence)
  reference_edge_mass <- as.numeric(crossprod(
    reference_selected,
    reference_relation %*% reference_selected
  ))
  source_edge_mass <- as.numeric(crossprod(
    source_selected,
    source_relation %*% source_selected
  ))
  forward_agreement <- reference_relation %*% correspondence %*%
    t(source_relation)
  agreement <- sum(correspondence * forward_agreement)
  denominator <- reference_edge_mass + source_edge_mass

  if (denominator <= edge_tolerance) {
    residual <- 0
    derivative <- if (gradient) {
      matrix(0, nrow(correspondence), ncol(correspondence))
    } else {
      NULL
    }
  } else {
    raw_residual <- 1 - 2 * agreement / denominator
    numerical_slack <- 64 * .Machine$double.eps
    if (raw_residual < -numerical_slack ||
        raw_residual > 1 + numerical_slack) {
      stop("edge-conditioned Dice residual fell outside [0, 1].")
    }
    residual <- min(1, max(0, raw_residual))
    derivative <- NULL
    if (gradient) {
      agreement_gradient <- forward_agreement +
        t(reference_relation) %*% correspondence %*% source_relation
      reference_gradient <- as.vector(
        (reference_relation + t(reference_relation)) %*%
          reference_selected
      )
      source_gradient <- as.vector(
        (source_relation + t(source_relation)) %*% source_selected
      )
      denominator_gradient <- outer(reference_gradient, rep(1, ncol(correspondence))) +
        outer(rep(1, nrow(correspondence)), source_gradient)
      derivative <- -2 * (
        agreement_gradient * denominator - agreement * denominator_gradient
      ) / denominator^2
    }
  }

  list(
    residual = residual,
    agreement = agreement,
    reference_edge_mass = reference_edge_mass,
    source_edge_mass = source_edge_mass,
    denominator = denominator,
    reference_selected = reference_selected,
    source_selected = source_selected,
    gradient = derivative,
    semantics = "directed_edge_conditioned_dice"
  )
}

transport_v3_edge_residual <- function(correspondence, reference_relation,
                                       source_relation,
                                       edge_tolerance = 1e-15) {
  transport_v3_edge_terms(
    correspondence,
    reference_relation,
    source_relation,
    gradient = FALSE,
    edge_tolerance = edge_tolerance
  )$residual
}

transport_v3_edge_gradient <- function(correspondence, reference_relation,
                                       source_relation,
                                       edge_tolerance = 1e-15) {
  transport_v3_edge_terms(
    correspondence,
    reference_relation,
    source_relation,
    gradient = TRUE,
    edge_tolerance = edge_tolerance
  )$gradient
}
