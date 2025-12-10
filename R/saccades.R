
#' Convert Cartesian Coordinates to Polar Coordinates
#'
#' This function converts Cartesian coordinates (x, y) to polar coordinates (rho, theta).
#'
#' @param x A numeric vector representing the x-coordinates.
#' @param y A numeric vector representing the y-coordinates.
#'
#' @return A matrix with two columns, where the first column is rho (the radial coordinate)
#'   and the second column is theta (the angular coordinate).
#' @examples
#' cart2pol(c(1, 2), c(2, 2))
cart2pol <- function(x, y) {
  rho <- sqrt(x ^ 2 + y ^ 2)
  theta = atan2(y, x)
  cbind(rho, theta)
}

rad2deg <- function(rad) {(rad * 180) / (pi)}



#' Calculate the Angle Between Two Vectors
#'
#' This function calculates the angle (in degrees) between two vectors, x1 and x2.
#'
#' @param x1 A numeric vector.
#' @param x2 A numeric vector.
#'
#' @return A numeric value representing the angle between the two vectors in degrees.
#' @export
#' @examples
#' calcangle(c(1, 2), c(2, 2))
calcangle <- function(x1, x2) {
  rad2deg(acos(sum(x1 * x2) / (norm(x1, "2") * norm(x2, "2"))))
}


#' Add Scanpath to a Data Frame
#'
#' This function adds a scanpath to a data frame.
#'
#' @param x A data frame.
#' @param outvar The output variable name for the scanpath. Defaults to "scanpath".
#' @param fixvar The fixation group variable name. Defaults to "fixgroup".
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with the added scanpath.
#' @export
#' @examples
#' # Create a data frame with a fixation group
#' df <- data.frame(x = 1:5, y = 6:10, fixgroup = rep(1, 5))
#' # Add a scanpath to the data frame
#' df <- add_scanpath.data.frame(df)
add_scanpath.data.frame <- function(x, outvar="scanpath", fixvar="fixgroup", ...) {
  x %>% mutate(!!outvar := list(scanpath(.data[[fixvar]][[1]])))
}

#' Add Scanpath to an Eye Table
#'
#' This function adds a scanpath to an eye table.
#'
#' @param x An eye table object.
#' @param outvar The output variable name for the scanpath. Defaults to "scanpath".
#' @param fixvar The fixation group variable name. Defaults to "fixgroup".
#' @param ... Additional arguments (currently unused).
#'
#' @return An eye table object with the added scanpath.
#' @export
#' @examples
#' # Create an eye table with a fixation group
#' df <- data.frame(x = 1:5, y = 6:10, fixgroup = rep(1, 5))
#' eye_table_df <- as_eye_table(df)
#' # Add a scanpath to the eye table
#' eye_table_df <- add_scanpath.eye_table(eye_table_df)
add_scanpath.eye_table <- function(x, outvar="scanpath", fixvar="fixgroup", ...) {
  x %>% mutate(!!outvar := list(scanpath(.data[[fixvar]][[1]])))
}



#' Create a Scanpath for a Fixation Group
#'
#' This function creates a scanpath for a fixation group.
#'
#' @param x A fixation group object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A scanpath object.
#' @export
#' @examples
#' # Create a fixation group
#' fixgroup <- data.frame(x = 1:5, y = 6:10)
#' # Create a scanpath for the fixation group
#' scanpath_obj <- scanpath.fixation_group(fixgroup)
scanpath.fixation_group <- function(x,...) {

  lenx <- diff(x$x)
  leny <- diff(x$y)

  polar <- cart2pol(lenx, leny)
  x <- x %>%
    mutate(lenx=c(lenx,0), leny=c(leny,0), rho=c(polar[,1],0), theta=c(polar[,2], 0))
  class(x) <- c("scanpath", class(x))
  x
}



