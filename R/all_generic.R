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

#' @export
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
#' @param ... Additional arguments passed on to the method.
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
#' @seealso \code{\link{rescale.fixation_group}} and \code{\link{rescale.eye_density}} for examples of specific method implementations for this generic.
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
#' @param x the fixation group
#' @param time the continuous time points to sample
#' @export
sample_fixations <- function(x, time) {
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
#' @return A numeric value representing the similarity between the two input objects.
#'
#' @examples
#' # Example usage of the similarity function
#' object1 <- # first object data
#' object2 <- # second object data
#' similarity_value <- similarity(object1, object2, method = "some_method")
#'
#' @family similarity
#' @export
similarity <- function(x, y, method, ...) {
  UseMethod("similarity", x)
}

#' Compute Density Matrix for a Given Object
#'
#' This function computes the density matrix for a given object, optionally taking into account a grouping variable.
#'
#' @param x The input object for which the density matrix should be computed.
#' @param groups An optional grouping variable to consider when computing the density matrix.
#' @param ... Additional arguments to be passed to the density matrix computation method.
#'
#' @return A density matrix representing the input object, taking into account the specified grouping variable if provided.
#'
#' @examples
#' # Example usage of the density_matrix function
#' input_object <- # input object data
#' grouping_variable <- # optional grouping variable
#' result_density_matrix <- density_matrix(input_object, groups = grouping_variable)
#'
#' @export
density_matrix <- function(x, groups,...) {
  UseMethod("density_matrix", x)
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
#' # Example usage of the center function
#' input_object <- # input object data containing eye-movements
#' new_origin <- c(500, 500) # New coordinate system origin
#' centered_object <- center(input_object, origin = new_origin)
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
#' # Example usage of the normalize function
#' input_object <- # input object data containing eye-movements
#' x_bounds <- c(0, 1000) # X bounds for normalization
#' y_bounds <- c(0, 1000) # Y bounds for normalization
#' normalized_object <- normalize(input_object, xbounds = x_bounds, ybounds = y_bounds)
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
#' # Example usage of the add_scanpath function
#' input_dataset <- # input dataset data
#' # Additional arguments required for the specific method that adds the scanpath
#' updated_dataset <- add_scanpath(input_dataset, ...)
#'
#' @export
add_scanpath <- function(x, ...) {
  UseMethod("add_scanpath", x)
}
