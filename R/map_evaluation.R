#' Sample a smooth fixation density map with a set of discrete fixations.
#'
#' This function samples a smooth fixation density map represented by the object \code{x} with a set of discrete fixations provided in \code{fix}.
#'
#' @param times A vector of numeric values representing the time points at which the density map should be sampled (default is NULL).
#' 
#' @details The function first checks if the \code{times} parameter is NULL. If so, it directly samples the density map using the coordinates of the fixations in the \code{fix} argument. If the \code{times} parameter is provided, the function first calls the \code{sample_fixations} function to generate a new fixation sequence with the specified time points, and then samples the density map using the coordinates of the new fixation sequence. The result is a data frame containing the sampled density values and the corresponding time points.
#'
#' Each point is looked up at the nearest lattice point, separately in x and y;
#' there is no interpolation. The nearest index is found with \code{round()}, so a
#' coordinate exactly midway between two grid points is resolved half to even on
#' the index scale (on the grid 0, 25, 50, 75, 100, the value 12.5 maps to 25 and
#' 37.5 also maps to 25). Coordinates outside the lattice are clamped to the
#' nearest edge point, so an off-grid fixation takes the edge value rather than
#' \code{NA}.
#'
#' @param normalize A character string specifying how to normalize the density
#'   map before sampling. One of:
#'   \describe{
#'     \item{\code{"none"}}{No normalization (default). Returns raw density values.}
#'     \item{\code{"max"}}{Divide by the maximum density value, yielding values between zero and one.}
#'     \item{\code{"sum"}}{Divide by the sum of all density values, yielding a
#'       probability distribution.}
#'     \item{\code{"zscore"}}{Z-score the density map (subtract mean, divide by SD).
#'       Useful for comparing across maps with different scales.}
#'   }
#' @param transform Map value transformation, applied after normalization.
#'   `log` computes `log(z + floor)`; `none` leaves values unchanged.
#' @param floor Explicit positive additive offset required for `log`. For maps
#'   normalized to unit sum, `1 / (10 * length(x$z))` is a possible choice.
#'   This is an offset, not a hard lower bound.
#' @rdname sample_density
#' @importFrom stats approx
#' @md
#' @export
sample_density.density <- function(x, fix, times = NULL, normalize = c("none", "max", "sum", "zscore"), ...,
                                   transform = c("none", "log"), floor = NULL) {
  normalize <- match.arg(normalize)
  transform <- match.arg(transform)
  zmat <- prepare_map_values(x$z, normalize, transform, floor)

  if (is.null(times)) {
    idx <- density_grid_index(fix, x$x, x$y)
    data.frame(z = zmat[idx], time = fix$onset)
  } else {
    fg <- sample_fixations(fix, times)
    idx <- density_grid_index(fg, x$x, x$y)
    data.frame(z = zmat[idx], time = times)
  }
}

# Normalize a density matrix as described in sample_density.density().
normalize_density_values <- function(zmat, normalize) {
  if (normalize == "max") {
    mx <- max(zmat, na.rm = TRUE)
    if (mx > 0) zmat <- zmat / mx
  } else if (normalize == "sum") {
    s <- sum(zmat, na.rm = TRUE)
    if (s > 0) zmat <- zmat / s
  } else if (normalize == "zscore") {
    mu <- mean(zmat, na.rm = TRUE)
    sd_val <- stats::sd(as.vector(zmat), na.rm = TRUE)
    if (sd_val > 0) zmat <- (zmat - mu) / sd_val
  }
  zmat
}

# Two-column (ix, iy) matrix index of the grid cells nearest to the x/y
# coordinates of `fix`, clamped to the grid.
density_grid_index <- function(fix, grid_x, grid_y) {
  nearest_index <- function(coord, grid) {
    ind <- round(approx(grid, seq_along(grid), coord, rule = 2)$y)
    ind[ind < 1L] <- 1L
    ind[ind > length(grid)] <- length(grid)
    ind
  }

  cds <- cbind(fix$x, fix$y)
  ix <- nearest_index(cds[, 1], grid_x)
  iy <- nearest_index(cds[, 2], grid_y)
  cbind(ix, iy)
}

# Pre-normalized maps for the templates in `dens_list` (restricted to `which`)
# that sample_density() would send to sample_density.density() and that lie on
# the lattice of the first such template. Returns list(x, y, z): the shared
# grid vectors and a list with the normalized matrix for each eligible
# template, NULL elsewhere. x and y are NULL when no template is eligible.
density_lookup_maps <- function(dens_list, normalize, which = seq_along(dens_list), enabled = TRUE,
                                transform = "none", floor = NULL) {
  ret <- list(x = NULL, y = NULL, z = vector("list", length(dens_list)))
  if (!enabled) {
    return(ret)
  }

  uses_density_method <- function(obj) {
    for (cl in class(obj)) {
      if (!is.null(utils::getS3method("sample_density", cl, optional = TRUE))) {
        return(identical(cl, "density"))
      }
    }
    FALSE
  }

  for (j in unique(which)) {
    obj <- dens_list[[j]]
    if (!is.list(obj) || !uses_density_method(obj) ||
        !is.numeric(obj$x) || !is.numeric(obj$y) || !is.matrix(obj$z) ||
        !identical(dim(obj$z), c(length(obj$x), length(obj$y)))) {
      next
    }
    if (is.null(ret$x)) {
      ret$x <- obj$x
      ret$y <- obj$y
    } else if (!identical(obj$x, ret$x) || !identical(obj$y, ret$y)) {
      next
    }
    # Maps whose normalization warns or fails keep the per-call path, which
    # reports the condition as before.
    ret$z[j] <- list(tryCatch(prepare_map_values(obj$z, normalize, transform, floor),
                              warning = function(w) NULL,
                              error = function(e) NULL))
  }
  ret
}

# Shared value preparation. Map z-scoring is distinct from fitting score scaling.
prepare_map_values <- function(z, normalize = "none", transform = "none", floor = NULL) {
  transform <- match.arg(transform, c("none", "log"))
  if (transform == "log") {
    check_log_floor(floor)
    if (normalize == "zscore") stop("Log transformation cannot follow zscore normalization.")
    if (any(z < 0, na.rm = TRUE)) stop("Log transformation requires non-negative map values.")
  }
  z <- normalize_density_values(z, normalize)
  if (transform == "log") log(z + floor) else z
}

check_log_floor <- function(floor) {
  if (!is.numeric(floor) || length(floor) != 1L || !is.finite(floor) || floor <= 0) {
    stop("floor must be an explicit positive finite number for log transformation.")
  }
}

validate_density_grid <- function(map) {
  valid_axis <- function(x) is.numeric(x) && length(x) >= 2L &&
    all(is.finite(x)) && all(diff(x) > 0)
  if (!is.list(map) || !valid_axis(map$x) || !valid_axis(map$y) ||
      !is.matrix(map$z) || !is.numeric(map$z) ||
      !identical(dim(map$z), c(length(map$x), length(map$y)))) {
    stop("Maps require increasing finite x/y grids and a matching numeric z matrix.")
  }
  invisible(map)
}

# Legacy extension methods need not accept newly added arguments. Only forward
# transformation arguments when requested, retaining the original default call.
sample_map_dispatch <- function(map, fix, times, normalize, transform, floor) {
  args <- list(x = map, fix = fix, times = times, normalize = normalize)
  if (transform != "none" || !is.null(floor)) {
    args$transform <- transform
    args$floor <- floor
  }
  do.call(sample_density, args)
}

#' @rdname sample_density
#' @param format For a named list of maps, return long rows or one column per map.
#' @md
#' @export
sample_density.list <- function(x, fix, times = NULL, ..., format = c("long", "wide")) {
  format <- match.arg(format)
  if (!length(x) || is.null(names(x)) || anyNA(names(x)) || any(names(x) == "") ||
      anyDuplicated(names(x))) stop("Supply a nonempty, uniquely named list of maps.")
  if (any(names(x) %in% c("time", "fixation"))) stop("Map names cannot be time or fixation.")
  values <- lapply(x, sample_density, fix = fix, times = times, ...)
  n <- nrow(values[[1]])
  if (!all(vapply(values, function(v) nrow(v) == n &&
                  identical(v$time, values[[1]]$time), logical(1)))) {
    stop("Sampled maps must return aligned rows and times.")
  }
  if (format == "wide") {
    out <- data.frame(fixation = seq_len(n), time = values[[1]]$time)
    for (nm in names(values)) out[[nm]] <- values[[nm]]$z
    return(out)
  }
  dplyr::bind_rows(lapply(seq_along(x), function(j) {
    data.frame(fixation = seq_len(n), reference = rep(names(x)[j], n),
               z = values[[j]]$z, time = values[[j]]$time)
  }))
}

#' @rdname sample_density
#' @param map_var,name_var For a map table, columns containing maps and unique names.
#' @md
#' @export
sample_density.data.frame <- function(x, fix, times = NULL, ...,
                                      map_var = "density", name_var = "reference") {
  if (!all(c(map_var, name_var) %in% names(x))) stop("Missing map or name column.")
  sample_density(stats::setNames(x[[map_var]], x[[name_var]]), fix, times, ...)
}


# Two density maps are comparable cell by cell only when they share a lattice.
density_lattices_equal <- function(x, y) {
  gx <- as.numeric(x$x)
  gy <- as.numeric(x$y)
  hx <- as.numeric(y$x)
  hy <- as.numeric(y$y)
  length(gx) == length(hx) && length(gy) == length(hy) &&
    isTRUE(all.equal(gx, hx)) && isTRUE(all.equal(gy, hy)) &&
    identical(dim(x$z), dim(y$z))
}

check_same_density_lattice <- function(x, y) {
  if (!density_lattices_equal(x, y)) {
    # Classed so callers that turn per-pair errors into NA can let it through.
    stop(structure(
      class = c("eyesim_lattice_mismatch", "error", "condition"),
      list(
        message = paste0(
          "similarity(): density maps are on different lattices. Both maps must ",
          "have the same x and y grid coordinates (same bounds and outdim)."
        ),
        call = NULL
      )
    ))
  }
  invisible(TRUE)
}

