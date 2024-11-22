
#' Calculate Eye Density by Groups
#'
#' This function calculates the eye density for fixations, grouped by specified variables.
#'
#' @param x A data frame containing fixations and additional grouping variables.
#' @param groups A character vector specifying the grouping variables to use.
#' @param sigma A numeric value specifying the bandwidth for the kernel density estimation (default is 50).
#' @param xbounds A numeric vector of length 2 specifying the x-axis bounds for the density calculation (default is c(0, 1000)).
#' @param ybounds A numeric vector of length 2 specifying the y-axis bounds for the density calculation (default is c(0, 1000)).
#' @param outdim A numeric vector of length 2 specifying the dimensions of the output density matrix (default is c(100, 100)).
#' @param duration_weighted A logical value indicating whether the density should be weighted by fixation duration (default is TRUE).
#' @param window A numeric vector of length 2 specifying the time window for selecting fixations (default is NULL).
#' @param keep_vars A character vector specifying additional variables to keep in the output (default is NULL).
#' @param fixvar A character string specifying the name of the fixation variable (default is "fixgroup").
#' @param result_name A character string specifying the name for the density result variable (default is "density").
#' @param ... Additional arguments passed to the `eye_density` function.
#'
#' @return A data frame containing the original grouping variables, the fixations, and the density result.
#'
#' @examples
#' # Create a data frame with fixations and a grouping variable
#' fixations <- data.frame(
#'   subject = rep(c("A", "B"), each = 25),
#'   x = runif(50, 0, 1000),
#'   y = runif(50, 0, 1000),
#'   duration = rep(1, 50),
#'   onset = seq(1, 50)
#' )
#'
#' # Calculate eye density by subject
#' result <- density_by(fixations, groups = "subject")
#'
#' @export
#' @family density_by
#' @import rlang
#' @importFrom dplyr group_by do rowwise
#' @importFrom tibble as_tibble
#' Calculate Eye Density by Groups
#'
#' This function calculates the eye density for fixations, grouped by specified variables.
#'
#' @param x A data frame containing fixations and additional grouping variables.
#' @param groups A character vector specifying the grouping variables to use.
#' @param sigma A numeric value specifying the bandwidth for the kernel density estimation (default is 50).
#' @param xbounds A numeric vector of length 2 specifying the x-axis bounds for the density calculation (default is c(0, 1000)).
#' @param ybounds A numeric vector of length 2 specifying the y-axis bounds for the density calculation (default is c(0, 1000)).
#' @param outdim A numeric vector of length 2 specifying the dimensions of the output density matrix (default is c(100, 100)).
#' @param duration_weighted A logical value indicating whether the density should be weighted by fixation duration (default is TRUE).
#' @param window A numeric vector of length 2 specifying the time window for selecting fixations (default is NULL).
#' @param keep_vars A character vector specifying additional variables to keep in the output (default is NULL).
#' @param result_name A character string specifying the name for the density result variable (default is "density").
#' @param ... Additional arguments passed to the `eye_density.fixation_group` function.
#'
#' @return A data frame containing the original grouping variables, the fixations, and the density result.
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
#'
#' # Calculate eye density by subject
#' result <- density_by(fixations, groups = "subject")
#'
#' @export
#' @import rlang
#' @importFrom dplyr group_by do rowwise
#' @importFrom tibble as_tibble
density_by <- function(x, groups, sigma = 50, xbounds = c(0, 1000), ybounds = c(0, 1000), outdim = c(100, 100),
                       duration_weighted = TRUE, window = NULL, keep_vars = NULL, fixvar = "fixgroup", result_name = "density", ...) {
  ## TODO what happens if window produces fixations < 0?

  rname <- rlang::sym(result_name)
  vars <- c(groups, keep_vars)

  if (!missing(groups) && !is.null(groups)) {
    ret <- x %>%
      group_by(across(all_of(groups))) %>%
      do({
        g <- do.call(rbind, .[[fixvar]])
        cbind(.[1, vars], tibble(!!fixvar := list(g)))
      }) %>%
      rowwise() %>%
      do({
        d <- eye_density(.[[fixvar]], sigma, xbounds = xbounds, ybounds = ybounds, outdim = outdim,
                         duration_weighted = duration_weighted, window = window, origin = attr(x, "origin"), ...)
        cbind(as_tibble(.[vars]), tibble(!!fixvar := list(.[[fixvar]]), !!rname := list(d)))
      })
  } else {
    fx <- do.call(rbind, x[[fixvar]])
    d <- eye_density(fx, sigma, xbounds = xbounds, ybounds = ybounds, outdim = outdim,
                     duration_weighted = duration_weighted, window = window, origin = attr(x, "origin"), ...)
    ret <- tibble(!!fixvar := list(fx), !!rname := list(d))
  }

  ret
}


# density_by <- function(x, groups, sigma=50, xbounds=c(0, 1000), ybounds=c(0, 1000), outdim=c(100,100),
#                        duration_weighted=TRUE, window=NULL, keep_vars=NULL, fixvar="fixgroup", result_name="density", ...) {
#
#   ## TODO what happens if window produces fixations < 0?
#
#   rname <- rlang::sym(result_name)
#   vars <- c(groups, keep_vars)
#
#   if (!missing(groups) && !is.null(groups) ) {
#     ret <- x %>% group_by(.dots=groups) %>% do( {
#       g <- do.call(rbind, .[[fixvar]])
#       cbind(.[1,vars],tibble(!!fixvar := list(g)))
#     }) %>% rowwise() %>% do( {
#       d <- eye_density(.[[fixvar]], sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim,
#                        duration_weighted=duration_weighted, window=window, origin=attr(x, "origin"), ...)
#       cbind(as_tibble(.[vars]), tibble(!!fixvar :=list(.[[fixvar]]), !!rname := list(d)))
#     })
#   } else {
#     #browser()
#     fx <- do.call(rbind, x[[fixvar]])
#     d <- eye_density(fx, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim,
#                      duration_weighted=duration_weighted, window=window,origin=attr(x, "origin"), ...)
#     ret <- tibble(!!fixvar := list(fx), !!rname := list(d))
#
#   }
#
#   ret
#
# }
#
# density_by <- function(x, groups, sigma = 50,
#                        xbounds = c(0, 1000), ybounds = c(0, 1000),
#                        outdim = c(100, 100),
#                        duration_weighted = TRUE,
#                        window = NULL,
#                        keep_vars = NULL,
#                        fixvar = "fixgroup",
#                        result_name = "density", ...) {
#
#   require(dplyr)
#   require(rlang)
#   require(tibble)
#
#   rname <- sym(result_name)
#   vars <- keep_vars  # Exclude groups from vars to avoid conflicts
#
#   # Group and compute density
#   if (!missing(groups) && !is.null(groups)) {
#     ret <- x %>%
#       group_by(across(all_of(groups))) %>%
#       group_modify(~ {
#         .x <- .x  # The data for the current group
#
#         # Extract fixations from the list-column if fixvar is specified
#         if (!is.null(fixvar) && fixvar %in% names(.x)) {
#           fixations_list <- .x[[fixvar]]
#           # Assuming fixations_list is a list of fixation data frames
#           fixations <- do.call(rbind, fixations_list)
#         } else {
#           fixations <- .x
#         }
#
#         density_map <- eye_density.fixation_group(
#           x = fixations,
#           sigma = sigma,
#           xbounds = xbounds,
#           ybounds = ybounds,
#           outdim = outdim,
#           duration_weighted = duration_weighted,
#           window = window,
#           origin = attr(x, "origin"),
#           ...
#         )
#
#         # Prepare the return data frame without grouping variables
#         ret_df <- tibble(!!rname := list(density_map))
#
#         # Include additional variables if specified
#         if (!is.null(vars)) {
#           additional_vars <- .x[1, vars, drop = FALSE]
#           ret_df <- bind_cols(additional_vars, ret_df)
#         }
#
#         ret_df
#       }) %>%
#       ungroup()
#   } else {
#     # No grouping, compute density for all fixations
#     if (!is.null(fixvar) && fixvar %in% names(x)) {
#       fixations_list <- x[[fixvar]]
#       fixations <- do.call(rbind, fixations_list)
#     } else {
#       fixations <- x
#     }
#
#     density_map <- eye_density.fixation_group(
#       x = fixations,
#       sigma = sigma,
#       xbounds = xbounds,
#       ybounds = ybounds,
#       outdim = outdim,
#       duration_weighted = duration_weighted,
#       window = window,
#       origin = attr(x, "origin"),
#       ...
#     )
#
#     ret <- tibble(!!rname := list(density_map))
#   }
#
#   return(ret)
# }


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
