#' Resolve transformation
#'
#' Returns the appropriate transformation object for a given name.
#'
#' @param transform Character string specifying the transformation. Supported
#'   values are "identity", "sqroot", "curoot", and "rank".
#' @return A transformation object or the string "identity" when no matching
#'   transformation is found.
#' @keywords internal
#' @noRd
resolve_transform <- function(transform) {
  if (transform == "identity") {
    "identity"
  } else if (transform == "sqroot") {
    squareroot_trans
  } else if (transform == "curoot") {
    cuberoot_trans
  } else if (transform == "rank") {
    rank_trans
  } else {
    "identity"
  }
}
