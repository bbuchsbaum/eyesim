#' eye_density
#'
#' compute a density map for a set of eye fixations
#'
#' @param x the fixations
#' @param sigma the standard deviation of the kernel
#' @param xbounds the x bounds
#' @param ybounds the y bounds
#' @param outdim the output dimensionality of the density map
#' @param weights fixation weights
#' @param normalize whether to normalize the output map
#' @export
eye_density <- function(x, sigma, xbounds, ybounds, outdim, weights, normalize,...) {
  UseMethod("eye_density", x)
}

#' @export
get_density <- function(x,...) {
  UseMethod("get_density", x)
}

#' sample a smooth fixation density map with a set of discrete fixations
#'
#' @param x the object to sample
#' @param fix the fixations
#' @export
sample_density <- function(x,fix,...) {
  UseMethod("sample_density", x)
}



#' extract coordinates
#'
#' @param x the object
#' @export
coords <- function(x) {
  UseMethod("coords", x)
}

#' rescale spatial coordinats
#'
#' @param x the object to rescale
#' @param sx the x scale factor
#' @param sy the y scale factor
rescale <- function(x, sx, sy) {
  UseMethod("rescale", x)
}


#' replicate a fixation sequence
#'
#' @param x the object
#' @param resolution the temporal resolution of the replicated fixations
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
#' @param x the first object to compare
#' @param y the second object to compare
#' @param method the similarity metric
#' @param ... extra args
#' @export
similarity <- function(x, y, method, ...) {
  UseMethod("similarity", x)
}

#' density matrix
#'
#' @param x the object
#' @param groups grouping variable
#' @param ... extra args
#' @export
density_matrix <- function(x, groups,...) {
  UseMethod("density_matrix", x)
}

#' center eye-movements
#'
#' @param x the object
#' @param origin the origin of the new coordinate system
#' @param ... extra args
#' @export
center <- function(x, origin, ...) {
  UseMethod("center", x)
}


#' normalize eye-movements to unit range
#'
#' @param x the object
#' @param xbounds the x bounds
#' @param ybounds the y bounds
#' @param ... extra args
#' @export
normalize <- function(x, xbounds, ybounds, ...) {
  UseMethod("normalize", x)
}

#' construct a scanpath of a fixation_group of related object
#'
#' a scanpath contains polor coordinates (rho, theta) along with absolute x and y spatial coordinates.
#'
#'
#'
#' @param x the fixations
#' @param ... extra args
#' @export
#' @examples
#'
#' fg <- fixation_group(x=c(.1,.5,1), y=c(1,.5,1), onset=1:3, duration=rep(1,3))
#' sp <- scanpath(fg)
#'
scanpath <- function(x, ...) {
  UseMethod("scanpath", x)
}

#' add scanpath to dataset
#'
#' @param x the input dataset
#' @param ... extra args
#' @export
add_scanpath <- function(x, ...) {
  UseMethod("add_scanpath", x)
}
