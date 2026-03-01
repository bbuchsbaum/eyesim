#' Calculate Eye Density by Groups
#'
#' This function calculates the eye density for fixations, grouped by specified variables.
#'
#' @param x A data frame containing fixations and additional grouping variables.
#' @param groups A character vector specifying the grouping variables to use.
#' @param sigma A numeric value or a numeric vector specifying the bandwidth(s) for the kernel density estimation. If a vector is provided, multiscale densities are computed. Default is 50.
#' @param xbounds A numeric vector of length 2 specifying the x-axis bounds for the density calculation (default is c(0, 1000)).
#' @param ybounds A numeric vector of length 2 specifying the y-axis bounds for the density calculation (default is c(0, 1000)).
#' @param outdim A numeric vector of length 2 specifying the dimensions of the output density matrix (default is c(100, 100)).
#' @param duration_weighted A logical value indicating whether the density should be weighted by fixation duration (default is TRUE).
#' @param window A numeric vector of length 2 specifying the time window for selecting fixations (default is NULL).
#' @param min_fixations Minimum number of fixations required for computing a density map.
#'   Rows with fewer fixations after optional filtering will receive `NULL` in the
#'   result column. Default is 2.
#' @param keep_vars A character vector specifying additional variables to keep in the output (default is NULL).
#' @param fixvar A character string specifying the name of the column containing fixation groups (default is "fixgroup").
#' @param result_name A character string specifying the name for the density result variable (default is "density").
#' @param ... Additional arguments passed to the `eye_density.fixation_group` function.
#'
#' @return A data frame containing the original grouping variables, the fixations, and the density result. If `sigma` is a vector, the column specified by `result_name` will contain `eye_density_multiscale` objects.
#'
#' @examples
#' # Create a data frame with fixations and a grouping variable
#' fixations <- data.frame(
#'   subject = rep(c("A", "B"), each = 25),
#'   x = runif(50, 0, 1000),
#'   y = runif(50, 0, 1000),
#'   duration = runif(50, 1, 5),
#'   onset = seq(1, 50)
#' )
#' eyetab <- eye_table("x", "y", "duration", "onset", groupvar=c("subject"), data=fixations)
#'
#' # Calculate eye density by subject
#' result <- density_by(eyetab, groups = "subject")
#'
#' @export
#' @import rlang
#' @importFrom dplyr group_by group_split bind_rows
#' @importFrom tibble as_tibble
density_by <- function(x, groups, sigma=50, xbounds=c(0, 1000), ybounds=c(0, 1000), outdim=c(100,100),
                       duration_weighted=TRUE, window=NULL, min_fixations=2,
                       keep_vars=NULL, fixvar="fixgroup", result_name="density", ...) {

  ## Window filtering that leaves fewer than min_fixations fixations is handled

  ## by eye_density.fixation_group(), which returns NULL. NULL results are
  ## removed below with a warning.

  rname <- rlang::sym(result_name)
  vars <- c(groups, keep_vars)

  if (!missing(groups) && !is.null(groups) ) {
    ret_list <- x %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(groups))) %>%
      dplyr::group_split() %>%
      lapply(function(df) {
        g <- do.call(rbind, df[[fixvar]])
        d <- eye_density(g, sigma,
                         xbounds = xbounds, ybounds = ybounds, outdim = outdim,
                         duration_weighted = duration_weighted, window = window,
                         min_fixations = min_fixations,
                         origin = attr(x, "origin"), ...)
        tibble::as_tibble(df[1, vars]) %>%
          dplyr::mutate(!!fixvar := list(g), !!rname := list(d))
      })
    ret <- dplyr::bind_rows(ret_list)
  } else {
    fx <- do.call(rbind, x[[fixvar]])
    d <- eye_density(fx, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim,
                     duration_weighted=duration_weighted, window=window,
                     min_fixations=min_fixations,
                     origin=attr(x, "origin"), ...)
    ret <- tibble(!!fixvar := list(fx), !!rname := list(d))

  }

  # Remove rows where density computation failed (NULL results)
  if (any(vapply(ret[[result_name]], is.null, logical(1)))) {
    warning("Removing rows with NULL density results in density_by().")
    ret <- ret %>%
      dplyr::ungroup() %>%
      dplyr::filter(!vapply(.data[[result_name]], is.null, logical(1)))
  }

  ret

}


#' Suggest Kernel Bandwidth for Density Estimation
#'
#' Estimates a reasonable kernel bandwidth (sigma) for fixation density maps
#' using a 2D variant of Silverman's rule of thumb, adapted for eye-tracking
#' data. The estimate accounts for spatial spread of fixations and sample size.
#'
#' @param x A \code{fixation_group} object, or a numeric vector of x-coordinates.
#' @param y A numeric vector of y-coordinates (only used when \code{x} is not a
#'   \code{fixation_group}).
#' @param xbounds Optional numeric vector of length 2 for the x-axis display
#'   bounds. Used to scale the estimate relative to display size.
#' @param ybounds Optional numeric vector of length 2 for the y-axis display
#'   bounds. Used to scale the estimate relative to display size.
#'
#' @return A single numeric value representing the suggested sigma
#'   (kernel standard deviation in coordinate units).
#'
#' @details
#' The function uses a 2D Silverman rule: \eqn{\sigma = n^{-1/6} \cdot
#' \sqrt{(IQR_x^2 + IQR_y^2)/2} / 1.349}. When display bounds are provided,
#' the result is clamped to between 1\% and 15\% of the mean display dimension
#' to avoid extreme values.
#'
#' @examples
#' fg <- fixation_group(x = runif(30, 0, 1280), y = runif(30, 0, 1024),
#'                      onset = cumsum(rep(200, 30)), duration = rep(200, 30))
#' suggest_sigma(fg, xbounds = c(0, 1280), ybounds = c(0, 1024))
#'
#' @export
#' @importFrom stats IQR
suggest_sigma <- function(x, y = NULL, xbounds = NULL, ybounds = NULL) {
  if (inherits(x, "fixation_group")) {
    yvals <- x$y
    xvals <- x$x
  } else {
    if (is.null(y)) stop("'y' must be provided when 'x' is not a fixation_group.")
    xvals <- x
    yvals <- y
  }

  n <- length(xvals)
  if (n < 2) return(NA_real_)

  iqr_x <- IQR(xvals)
  iqr_y <- IQR(yvals)

  # 2D Silverman rule: use geometric mean of IQRs, scaled by n^{-1/6}
  spread <- sqrt((iqr_x^2 + iqr_y^2) / 2) / 1.349
  sigma <- spread * n^(-1/6)

  # Clamp to reasonable fraction of display if bounds provided
  if (!is.null(xbounds) && !is.null(ybounds)) {
    display_scale <- mean(c(diff(range(xbounds)), diff(range(ybounds))))
    sigma <- max(sigma, display_scale * 0.01)  # at least 1% of display
    sigma <- min(sigma, display_scale * 0.15)  # at most 15% of display
  }

  sigma
}
