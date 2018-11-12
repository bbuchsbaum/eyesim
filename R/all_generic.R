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


#' @param x the object to sample
#' @param fix the fixations
#' @export
sample_density <- function(x,fix,...) {
  UseMethod("sample_density", x)
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

#' sample_fixations
#'
#' @param x the fixation group
#' @param time the continuous time points to sample
#' @export
sample_fixations <- function(x, time) {
  UseMethod("sample_fixations", x)
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

density_matrix <- function(x, groups,...) {
  UseMethod("density_matrix", x)
}
