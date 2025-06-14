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

  ## TODO what happens if window produces fixations < 0?

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
    #browser()
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


#' @keywords internal
#' @noRd
rank_trans <- scales::trans_new(name="rank",
                                transform=function(x) { rank(x) },
                                inverse=function(x) (length(x)+1) - rank(x))

#' @keywords internal
#' @noRd
cuberoot_trans <- scales::trans_new(name="rank",
                                    transform=function(x) { x^(1/3) },
                                    inverse=function(x) x^3)
