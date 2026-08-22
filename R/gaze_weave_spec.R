#' Declare screen geometry for GazeWeave
#'
#' @param width,height Positive screen dimensions.
#' @param unit Coordinate unit. One of `"px"`, `"deg"`, or `"normalized"`.
#' @param center Optional two-element screen centre. Defaults to half the screen
#'   dimensions.
#'
#' @return A `gaze_screen` specification.
#' @export
gaze_screen <- function(width, height, unit = c("px", "deg", "normalized"),
                        center = NULL) {
  unit <- match.arg(unit)
  if (length(width) != 1L || !is.finite(width) || width <= 0 ||
      length(height) != 1L || !is.finite(height) || height <= 0) {
    stop("width and height must be finite positive scalars.")
  }
  if (is.null(center)) {
    center <- c(width / 2, height / 2)
  }
  if (!is.numeric(center) || length(center) != 2L || any(!is.finite(center))) {
    stop("center must be a finite numeric vector of length two.")
  }

  structure(
    list(
      width = as.numeric(width),
      height = as.numeric(height),
      unit = unit,
      center = as.numeric(center),
      xlim = c(0, as.numeric(width)),
      ylim = c(0, as.numeric(height))
    ),
    class = c("gaze_screen", "list")
  )
}

#' Declare a multiscale Gaussian spatial model for GazeWeave
#'
#' The resulting spatial cost is the negative log of a normalized mixture of
#' Gaussian overlap kernels. Mixture weights are normalized to sum to one.
#'
#' @param sigmas Positive Gaussian bandwidths.
#' @param weights Optional non-negative mixture weights.
#' @param unit Coordinate unit for `sigmas`.
#'
#' @return A `gaze_spatial_spec` object.
#' @export
gaze_gaussian_mixture <- function(sigmas, weights = NULL,
                                  unit = c("native", "px", "deg", "normalized")) {
  unit <- match.arg(unit)
  if (!is.numeric(sigmas) || length(sigmas) == 0L ||
      any(!is.finite(sigmas)) || any(sigmas <= 0)) {
    stop("sigmas must contain finite positive values.")
  }
  if (is.null(weights)) {
    weights <- rep(1 / length(sigmas), length(sigmas))
  }
  if (!is.numeric(weights) || length(weights) != length(sigmas) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("weights must be finite, non-negative, and match sigmas.")
  }
  weights <- weights / sum(weights)

  structure(
    list(sigmas = as.numeric(sigmas), weights = as.numeric(weights), unit = unit),
    class = c("gaze_spatial_spec", "list")
  )
}

#' Declare the local chronology model for GazeWeave
#'
#' @param horizon Positive temporal horizon on normalized trial time.
#'
#' @return A `gaze_chronology_spec` object.
#' @export
gaze_local_order <- function(horizon = 0.2) {
  if (length(horizon) != 1L || !is.finite(horizon) || horizon <= 0) {
    stop("horizon must be a finite positive scalar.")
  }
  structure(
    list(type = "local_order", horizon = as.numeric(horizon)),
    class = c("gaze_chronology_spec", "list")
  )
}

#' Declare fixed-neighbour ordinal chronology for Transport
#'
#' @param neighbours Number of forward fixation neighbours represented in the
#'   directed relation graph.
#' @param coalesce_distance Maximum distance for merging immediately adjacent
#'   fixation atoms. The default zero merges only exactly identical locations.
#' @param unit Coordinate unit for `coalesce_distance`.
#'
#' @return A `gaze_chronology_spec` object.
#' @export
gaze_order_neighbours <- function(neighbours = 2L, coalesce_distance = 0,
                                  unit = c("native", "px", "deg", "normalized")) {
  neighbours <- as.integer(neighbours)
  unit <- match.arg(unit)
  if (length(neighbours) != 1L || is.na(neighbours) || neighbours < 1L) {
    stop("neighbours must be a positive integer.")
  }
  if (!is.numeric(coalesce_distance) || length(coalesce_distance) != 1L ||
      !is.finite(coalesce_distance) || coalesce_distance < 0) {
    stop("coalesce_distance must be a finite non-negative scalar.")
  }
  structure(
    list(
      type = "order_neighbours",
      clock = "order",
      neighbours = neighbours,
      coalesce_distance = as.numeric(coalesce_distance),
      unit = unit
    ),
    class = c("gaze_chronology_spec", "list")
  )
}


#' Declare no geometric registration for GazeWeave
#'
#' @return A `gaze_warp_spec` object.
#' @export
gaze_warp_none <- function() {
  structure(list(type = "none", fit_by = NULL), class = c("gaze_warp_spec", "list"))
}

#' Declare cross-fitted contraction registration for GazeWeave
#'
#' @param center Either `"screen"` or a finite two-element numeric vector.
#' @param translation Whether to estimate translation in addition to isotropic
#'   scale.
#' @param fit_by Optional columns defining separate warp strata, such as
#'   `"participant"`.
#' @param shrink Positive stabilization term for moment ratios.
#'
#' @return A `gaze_warp_spec` object.
#' @export
gaze_warp_contraction <- function(center = "screen", translation = TRUE,
                                  fit_by = NULL, shrink = 1e-6) {
  valid_center <- (is.character(center) && length(center) == 1L &&
                     identical(center, "screen")) ||
    (is.numeric(center) && length(center) == 2L && all(is.finite(center)))
  if (!valid_center) {
    stop("center must be 'screen' or a finite numeric vector of length two.")
  }
  if (!is.logical(translation) || length(translation) != 1L || is.na(translation)) {
    stop("translation must be TRUE or FALSE.")
  }
  if (!is.null(fit_by) && (!is.character(fit_by) || length(fit_by) == 0L)) {
    stop("fit_by must be NULL or a non-empty character vector.")
  }
  if (length(shrink) != 1L || !is.finite(shrink) || shrink <= 0) {
    stop("shrink must be a finite positive scalar.")
  }

  structure(
    list(
      type = "contraction",
      center = center,
      translation = translation,
      fit_by = fit_by,
      shrink = as.numeric(shrink)
    ),
    class = c("gaze_warp_spec", "list")
  )
}
