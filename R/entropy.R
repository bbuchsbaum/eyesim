#' Entropy of fixation patterns
#'
#' These methods quantify how concentrated or diffuse a fixation pattern is
#' using Shannon entropy. For density maps, entropy is computed from the
#' normalized density surface. For fixation groups, entropy can be computed from
#' either a derived density map or a discrete occupancy grid.
#'
#' @param x The input object.
#' @param normalize Logical; if `TRUE` (default), divide entropy by the maximum
#'   possible entropy for the number of valid bins so results lie in `[0, 1]`.
#' @param base Logarithm base used in the Shannon entropy. Defaults to
#'   `exp(1)`; use `2` for bits.
#' @param method For `fixation_group` objects, one of `"density"` (default) or
#'   `"grid"`.
#' @param sigma Optional bandwidth for density-based entropy on fixation groups.
#'   If `NULL`, `suggest_sigma()` is used.
#' @param xbounds,ybounds Optional display bounds for fixation groups. If not
#'   supplied, the observed fixation ranges are used with a small padding.
#' @param outdim Grid dimensions for density-based entropy from fixation groups.
#' @param grid Grid dimensions for occupancy-grid entropy from fixation groups.
#' @param duration_weighted Logical; if `TRUE`, duration-weighted KDE is used
#'   for density-based entropy.
#' @param aggregate For multiscale density objects, one of `"mean"` (default)
#'   or `"none"`.
#' @param ... Additional arguments passed to `eye_density()` when
#'   `method = "density"`.
#'
#' @return A numeric entropy value. Multiscale density objects return either the
#'   mean entropy across scales or a named numeric vector, depending on
#'   `aggregate`.
#' @export
fixation_entropy.default <- function(x, ...) {
  stop("No fixation_entropy() method for objects of class ", paste(class(x), collapse = "/"), ".")
}

#' @rdname fixation_entropy.default
#' @export
fixation_entropy.eye_density <- function(x, normalize = TRUE, base = exp(1), ...) {
  entropy_from_mass(x$z, normalize = normalize, base = base)
}

#' @rdname fixation_entropy.default
#' @export
fixation_entropy.density <- function(x, normalize = TRUE, base = exp(1), ...) {
  entropy_from_mass(x$z, normalize = normalize, base = base)
}

#' @rdname fixation_entropy.default
#' @export
fixation_entropy.eye_density_multiscale <- function(x, normalize = TRUE, base = exp(1),
                                                    aggregate = c("mean", "none"), ...) {
  aggregate <- match.arg(aggregate)
  ent <- vapply(x, fixation_entropy.eye_density, numeric(1), normalize = normalize, base = base, ...)
  sigmas <- attr(x, "sigmas_vector")
  if (is.null(sigmas)) {
    sigmas <- vapply(x, function(el) el$sigma %||% NA_real_, numeric(1))
  }
  names(ent) <- paste0("sigma_", sigmas)

  if (aggregate == "mean") {
    return(mean(ent, na.rm = TRUE))
  }
  ent
}

#' @rdname fixation_entropy.default
#' @export
fixation_entropy.fixation_group <- function(x, normalize = TRUE, base = exp(1),
                                            method = c("density", "grid"),
                                            sigma = NULL,
                                            xbounds = NULL, ybounds = NULL,
                                            outdim = c(50, 50),
                                            grid = c(10, 10),
                                            duration_weighted = FALSE,
                                            ...) {
  method <- match.arg(method)

  if (nrow(x) == 0L) {
    return(NA_real_)
  }

  if (method == "density") {
    bounds <- resolve_fixation_entropy_bounds(x, xbounds = xbounds, ybounds = ybounds)
    if (is.null(sigma)) {
      sigma <- suggest_sigma(x, xbounds = bounds$xbounds, ybounds = bounds$ybounds)
    }
    dens <- eye_density(
      x,
      sigma = sigma,
      xbounds = bounds$xbounds,
      ybounds = bounds$ybounds,
      outdim = outdim,
      duration_weighted = duration_weighted,
      ...
    )
    if (is.null(dens)) {
      return(NA_real_)
    }
    return(fixation_entropy(dens, normalize = normalize, base = base))
  }

  counts <- grid_fixation_counts(x, grid = grid, xbounds = xbounds, ybounds = ybounds)
  entropy_from_mass(counts, normalize = normalize, base = base)
}

resolve_fixation_entropy_bounds <- function(x, xbounds = NULL, ybounds = NULL) {
  pad_range <- function(v) {
    rng <- range(v, na.rm = TRUE)
    span <- diff(rng)
    if (!is.finite(span) || span <= .Machine$double.eps) {
      span <- 1
    }
    rng + c(-0.05, 0.05) * span
  }

  list(
    xbounds = if (is.null(xbounds)) pad_range(x$x) else xbounds,
    ybounds = if (is.null(ybounds)) pad_range(x$y) else ybounds
  )
}

grid_fixation_counts <- function(x, grid = c(10, 10), xbounds = NULL, ybounds = NULL) {
  if (length(grid) != 2L || any(grid < 1L)) {
    stop("grid must be a length-2 vector of positive integers.")
  }

  bounds <- resolve_fixation_entropy_bounds(x, xbounds = xbounds, ybounds = ybounds)
  gx <- seq(bounds$xbounds[1], bounds$xbounds[2], length.out = grid[1] + 1L)
  gy <- seq(bounds$ybounds[1], bounds$ybounds[2], length.out = grid[2] + 1L)

  ix <- findInterval(x$x, gx, rightmost.closed = TRUE, all.inside = TRUE)
  iy <- findInterval(x$y, gy, rightmost.closed = TRUE, all.inside = TRUE)

  counts <- matrix(0, nrow = grid[1], ncol = grid[2])
  for (k in seq_along(ix)) {
    counts[ix[[k]], iy[[k]]] <- counts[ix[[k]], iy[[k]]] + 1
  }
  counts
}

entropy_from_mass <- function(mass, normalize = TRUE, base = exp(1)) {
  vals <- as.numeric(mass)
  vals <- vals[is.finite(vals)]

  if (length(vals) == 0L) {
    return(NA_real_)
  }

  total <- sum(vals)
  if (!is.finite(total) || total <= .Machine$double.eps) {
    return(NA_real_)
  }

  probs <- vals / total
  probs <- probs[probs > 0]
  if (length(probs) == 0L) {
    return(NA_real_)
  }

  entropy <- -sum(probs * log(probs, base = base))
  if (!normalize) {
    return(entropy)
  }

  support_n <- length(vals)
  if (support_n <= 1L) {
    return(0)
  }

  entropy / log(support_n, base = base)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
