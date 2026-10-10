#' eye_density
#'
#'
#' Compute a density map for a set of eye fixations.
#'
#' This function uses kernel density estimation to compute a density map of eye fixations.
#' It takes various parameters to control the computation and output of the density map.
#'
#' @param x A data frame containing the eye fixations.
#' @param sigma Numeric, the standard deviation of the kernel used for density estimation.
#' @param xbounds Numeric vector of length 2, defining the minimum and maximum x-axis bounds of the density map.
#' @param ybounds Numeric vector of length 2, defining the minimum and maximum y-axis bounds of the density map.
#' @param outdim Numeric vector of length 2, specifying the number of rows and columns in the output density map.
#' @param weights Optional, a numeric vector of fixation weights to be used in the density estimation.
#' @param normalize Logical, whether to normalize the output density map such that its values sum to 1.
#' @param ... Additional arguments passed on to the method.
#'
#' @return An object of class "eye_density", "density", and "list" containing the computed density map and other relevant information.
#' @seealso \code{\link{eye_density.fixation_group}} for an example of a specific method implementation for this generic.
#' @export
#' @family eye_density
eye_density <- function(x, sigma, xbounds, ybounds, outdim, weights, normalize,...) {
  UseMethod("eye_density", x)
}


#' get_density
#'
#' Get the density matrix from an object.
#'
#' @export
#' @family get_density
#' @param x the object
#' @param ... additional arguments
#' @return The density values as a numeric matrix (for an \code{eye_density},
#'   its \code{z} component, with rows indexing \code{x} and columns \code{y}).
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' ed <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300),
#'                   outdim = c(40, 30))
#' dim(get_density(ed))
get_density <- function(x,...) {
  UseMethod("get_density", x)
}

#' sample_density
#'
#' Sample a smooth fixation density map with a set of discrete fixations.
#'
#' This function samples a given smooth fixation density map with a set of discrete fixations
#' to estimate the density at the locations of those fixations.
#'
#' @param x An object representing the smooth fixation density map to be sampled.
#' @param fix A data frame containing the discrete fixations used for sampling the density map.
#' @param ... Additional arguments passed on to the method, including
#'   \code{normalize} (one of \code{"none"}, \code{"max"}, \code{"sum"}, or
#'   \code{"zscore"}) to control density map normalization before sampling.
#'
#' @return A data frame with columns 'z' (density estimates at fixation locations) and 'time' (onset time of fixations).
#' @seealso \code{\link{sample_density.density}} for an example of a specific method implementation for this generic.
#' @export
sample_density <- function(x,fix,...) {
  UseMethod("sample_density", x)
}



#' extract coordinates
#'
#' @param x the object
#' @return A two-column numeric matrix with columns \code{x} and \code{y},
#'   one row per fixation.
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' coords(fg)
#' @export
coords <- function(x) {
  UseMethod("coords", x)
}

#' rescale
#'
#' Rescale spatial coordinates.
#'
#' This function rescales the spatial coordinates of an object by the given scale factors.
#'
#' @param x An object containing spatial coordinates to be rescaled.
#' @param sx A numeric value representing the x-axis scale factor.
#' @param sy A numeric value representing the y-axis scale factor.
#'
#' @family rescale
#' @return An object with rescaled spatial coordinates.
#' @seealso \code{\link{rescale.fixation_group}} for an example of a specific method implementation for this generic.
#' @export
rescale <- function(x, sx, sy) {
  UseMethod("rescale", x)
}


#' rep_fixations
#'
#' Replicate a fixation sequence.
#'
#' This function replicates a fixation sequence with a specified temporal resolution.
#' It can be useful when working with fixation data that needs to be resampled or when
#' creating fixation sequences with consistent temporal spacing.
#'
#' For fixation groups, each fixation is repeated
#' \code{floor(duration * resolution)} times, and at least once. The floor
#' tolerates floating-point error, so a duration of 0.29 at resolution 100 gives
#' 29 copies rather than 28.
#'
#' @param x An object representing a fixation sequence.
#' @param resolution A numeric value representing the temporal resolution of the replicated fixations.
#'
#' @return An object containing the replicated fixation sequence with the specified temporal resolution.
#' @seealso \code{\link{rep_fixations.fixation_group}} for an example of a specific method implementation for this generic.
#' @export
rep_fixations <- function(x, resolution) {
  UseMethod("rep_fixations", x)
}

#' sample_fixations
#'
#' Sample fixation coordinates at arbitrary time points.
#'
#' Each time point takes the coordinates of the most recent fixation whose onset
#' is at or before it. Fixation durations are not used: gaze during a saccade is
#' attributed to the preceding fixation, and the last fixation is held for all
#' later time points. Time points before the first onset give \code{NA}.
#'
#' @param x the fixation group
#' @param time the continuous time points to sample
#' @param ... Additional arguments passed to methods.
#' @return A data frame of class \code{sampled_fixation_group} with one row per
#'   time point and columns \code{x}, \code{y} (\code{NA} before the first
#'   onset), \code{onset} (the time point) and \code{duration} (1).
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' sample_fixations(fg, time = c(-50, 0, 250, 500))
#' @export
sample_fixations <- function(x, time, ...) {
  UseMethod("sample_fixations", x)
}

#' Compute Similarity Between Two Objects
#'
#' This function computes the similarity between two objects using a specified similarity metric.
#'
#' @param x The first object to compare.
#' @param y The second object to compare.
#' @param method A character string specifying the similarity metric to be used.
#' @param ... Additional arguments passed to the similarity computation method.
#'
#' @details Density maps are compared cell by cell, so two maps must share a
#' lattice: the same \code{x} and \code{y} grid coordinates, as produced by the
#' same bounds and \code{outdim}. Comparing maps on different lattices is an
#' error. A plain numeric \code{y} must have one value per grid cell of \code{x}.
#'
#' @return A numeric value representing the similarity between the two input objects.
#'
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' fg2 <- fixation_group(x = c(110, 190, 320), y = c(90, 160, 210),
#'                       onset = c(0, 200, 400), duration = c(200, 200, 200))
#' d1 <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300),
#'                   outdim = c(40, 30))
#' d2 <- eye_density(fg2, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300),
#'                   outdim = c(40, 30))
#' similarity(d1, d2, method = "pearson")
#' similarity(d1, d2, method = "cosine")
#'
#' @family similarity
#' @export
similarity <- function(x, y, method, ...) {
  UseMethod("similarity", x)
}

#' Compute Fixation Entropy
#'
#' This generic computes entropy summaries for fixation patterns represented as
#' fixation groups or density maps.
#'
#' @param x The input object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric entropy value, or a named numeric vector for multiscale
#'   objects when requested.
#' @export
fixation_entropy <- function(x, ...) {
  UseMethod("fixation_entropy", x)
}

#' Center Eye-Movements in a New Coordinate System
#'
#' This function centers the eye-movements in a new coordinate system with the specified origin.
#'
#' @param x The input object containing the eye-movements to be centered.
#' @param origin A vector containing the x and y coordinates of the new coordinate system's origin.
#' @param ... Additional arguments to be passed to the centering method.
#'
#' @return An object containing the centered eye-movements in the new coordinate system.
#'
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' # Centre on the mean fixation position (the default) ...
#' center(fg)
#' # ... or on a chosen origin, such as the centre of an 800 x 600 screen.
#' center(fg, origin = c(400, 300))
#'
#' @export
center <- function(x, origin, ...) {
  UseMethod("center", x)
}


#' Normalize Eye-Movements to Unit Range
#'
#' This function normalizes the eye-movements to a unit range based on the specified x and y bounds.
#'
#' @param x The input object containing the eye-movements to be normalized.
#' @param xbounds A vector containing the minimum and maximum x bounds for normalization.
#' @param ybounds A vector containing the minimum and maximum y bounds for normalization.
#' @param ... Additional arguments to be passed to the normalization method.
#'
#' @return An object containing the normalized eye-movements in the unit range.
#'
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' normalize(fg, xbounds = c(0, 400), ybounds = c(0, 300))
#'
#' @export
normalize <- function(x, xbounds, ybounds, ...) {
  UseMethod("normalize", x)
}

#' Construct a Scanpath of a Fixation Group of Related Objects
#'
#' This function creates a scanpath containing polar coordinates (rho, theta) along with
#' absolute x and y spatial coordinates for a given fixation group.
#'
#' @param x The fixations.
#' @param ... Extra arguments.
#'
#' @return A scanpath object.
#' @export
#' @examples
#' # Create a fixation group
#' fg <- fixation_group(x=c(.1,.5,1), y=c(1,.5,1), onset=1:3, duration=rep(1,3))
#' # Create a scanpath for the fixation group
#' sp <- scanpath(fg)
scanpath <- function(x, ...) {
  UseMethod("scanpath", x)
}

#' Add Scanpath to Dataset
#'
#' This function adds a scanpath to the input dataset.
#'
#' @param x The input dataset to which the scanpath will be added.
#' @param ... Additional arguments to be passed to the method that adds the scanpath.
#'
#' @return A dataset with the added scanpath.
#'
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' tab <- tibble::tibble(trial = 1, fixgroup = list(fg))
#' sp <- add_scanpath(tab)$scanpath[[1]]
#' # Saccade amplitude (rho) and direction (theta) from each fixation to the next
#' as.data.frame(sp)[, c("x", "y", "rho", "theta")]
#'
#' @export
add_scanpath <- function(x, ...) {
  UseMethod("add_scanpath", x)
}
