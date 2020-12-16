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

#' @export
get_density <- function(x,...) {
  UseMethod("get_density", x)
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

#' rescale spatial coordinats
#'
#' @param x the object to rescale
#' @param sx the x scale factor
#' @param sy the y scale factor
rescale <- function(x, sx, sy) {
  UseMethod("rescale", x)
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
#' @param x the first object to compare
#' @param y the second object to compare
#' @param method the similarity metric
#' @export
similarity <- function(x, y, method, ...) {
  UseMethod("similarity", x)
}


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
#' @param xbounds
#' @param ybounds
#' @param ... extra args
#' @export
normalize <- function(x, xbounds, ybounds, ...) {
  UseMethod("normalize", x)
}

#' compute scanpath of a fixation_group of related object
#'
#' @param x the fixations
#' @param ... extra args
#' @export
scanpath <- function(x, ...) {
  UseMethod("scanpath", x)
}

#' add scanpath to dataset
#'
#' @param x the inpuy dataset
#' @param ... extra args
#' @export
add_scanpath <- function(x, ...) {
  UseMethod("add_scanpath", x)
}
