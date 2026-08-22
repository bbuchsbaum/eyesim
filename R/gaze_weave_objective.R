# GazeWeave objective --------------------------------------------------------

row_log_sum_exp <- function(x) {
  row_max <- apply(x, 1L, max)
  row_max + log(rowSums(exp(sweep(x, 1L, row_max, FUN = "-"))))
}

gaze_spatial_cost <- function(reference_coords, source_coords, spatial) {
  if (!is.matrix(reference_coords) || ncol(reference_coords) != 2L ||
      !is.matrix(source_coords) || ncol(source_coords) != 2L) {
    stop("Gaze coordinates must be two-column matrices.")
  }
  dx <- outer(reference_coords[, 1], source_coords[, 1], FUN = "-")
  dy <- outer(reference_coords[, 2], source_coords[, 2], FUN = "-")

  log_terms <- matrix(
    unlist(lapply(seq_along(spatial$sigmas), function(i) {
      # Form dimensionless displacements before squaring. Besides avoiding
      # unnecessary overflow, this makes equivalent pixel/degree problems
      # follow the same numerical path when coordinates and bandwidths share
      # a conversion factor.
      scaled_distance_sq <-
        (dx / spatial$sigmas[[i]])^2 + (dy / spatial$sigmas[[i]])^2
      log(spatial$weights[[i]]) - as.vector(scaled_distance_sq) / 4
    }), use.names = FALSE),
    nrow = length(dx),
    ncol = length(spatial$sigmas)
  )
  cost <- -row_log_sum_exp(log_terms)
  # Remove sub-ULP conversion noise before the non-convex continuation path.
  # Fourteen significant digits retain far more precision than the declared
  # solver tolerances while making equivalent coordinate units reproducible.
  matrix(signif(cost, 14L),
         nrow = nrow(reference_coords), ncol = nrow(source_coords))
}
