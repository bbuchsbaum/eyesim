

eye_density <- function(x, sigma, xbounds, ybounds, outdim, weights, normalize) {
  UseMethod("eye_density", x)
}


coords <- function(x) {
  UseMethod("coords", x)
}


rep_fixations <- function(x, resolution) {
  UseMethod("rep_fixations", x)
}

similarity <- function(x, y, method) {
  UseMethod("similarity", x)
}

