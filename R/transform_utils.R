#' Internal transformation utilities for visualization
#'
#' These scales transformations are used throughout the package when
#' plotting eye density data. They are not exported.
#'
#' @keywords internal
#' @noRd
rank_trans <- scales::trans_new(
  name = "rank",
  transform = function(x) rank(x),
  inverse = function(x) (length(x) + 1) - rank(x)
)

#' @keywords internal
#' @noRd
cuberoot_trans <- scales::trans_new(
  name = "curoot",
  transform = function(x) x^(1/3),
  inverse = function(x) x^3
)

#' @keywords internal
#' @noRd
squareroot_trans <- scales::trans_new(
  name = "sqroot",
  transform = function(x) x^(1/2),
  inverse = function(x) x^2
)

