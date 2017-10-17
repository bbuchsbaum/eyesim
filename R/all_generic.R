
#' eye_density
#'
#' compute a density map for a set of eye fixations
#'
#' @param x
#' @param sigma
#' @param xbounds
#' @param ybounds
#' @param outdim
#' @param weights
#' @param normalize
#' @export
eye_density <- function(x, sigma, xbounds, ybounds, outdim, weights, normalize,...) {
  UseMethod("eye_density", x)
}


#' coords
#'
#' @param x
#' @export
coords <- function(x) {
  UseMethod("coords", x)
}


#' rep_fixations
#'
#' @param x
#' @param resolution
#' @export
rep_fixations <- function(x, resolution) {
  UseMethod("rep_fixations", x)
}

#' similarity
#'
#' @param x
#' @param y
#' @param method
#' @export
similarity <- function(x, y, method) {
  UseMethod("similarity", x)
}

