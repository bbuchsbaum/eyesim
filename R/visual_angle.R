#' Convert between visual angle and screen pixels
#'
#' Converts a visual angle to the on-screen extent in pixels, and back, for a
#' display viewed at a given distance. The extent is centred on the line of
#' sight: an angle \eqn{\theta} spans \eqn{2 d \tan(\theta / 2)} on the screen,
#' where \eqn{d} is the viewing distance. Use it to set distance thresholds such
#' as the \code{radius} of \code{\link{fixation_revisits}} in degrees.
#'
#' @param deg Visual angle in degrees (vectorised).
#' @param px Extent in pixels (vectorised).
#' @param distance Eye-to-screen distance.
#' @param screen_size Physical width (or height) of the screen, in the same units
#'   as \code{distance}.
#' @param screen_res Screen resolution in pixels along the same dimension as
#'   \code{screen_size}.
#'
#' @return \code{visual_angle_to_px()} returns pixels; \code{px_to_visual_angle()}
#'   returns degrees.
#'
#' @examples
#' # 2 degrees on a 52 cm wide, 1920 px display viewed from 60 cm
#' visual_angle_to_px(2, distance = 60, screen_size = 52, screen_res = 1920)
#' px_to_visual_angle(77, distance = 60, screen_size = 52, screen_res = 1920)
#' @export
visual_angle_to_px <- function(deg, distance, screen_size, screen_res) {
  check_display_geometry(distance, screen_size, screen_res)
  if (any(!is.finite(deg) | deg < 0 | deg >= 180)) stop("deg must be finite and in [0, 180).")
  2 * distance * tan(deg * pi / 360) * screen_res / screen_size
}

#' @rdname visual_angle_to_px
#' @export
px_to_visual_angle <- function(px, distance, screen_size, screen_res) {
  check_display_geometry(distance, screen_size, screen_res)
  if (any(!is.finite(px) | px < 0)) stop("px must be finite and non-negative.")
  2 * atan(px * screen_size / screen_res / (2 * distance)) * 180 / pi
}

check_display_geometry <- function(distance, screen_size, screen_res) {
  ok <- function(v) length(v) == 1L && is.numeric(v) && is.finite(v) && v > 0
  if (!ok(distance) || !ok(screen_size) || !ok(screen_res)) {
    stop("distance, screen_size and screen_res must each be a single positive number.")
  }
}
