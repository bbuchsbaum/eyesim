#' @noRd
emd_position_similarity <- function(fg1, fg2, screensize) {
  # Extract x and y coordinates
  points1 <- as.matrix(fg1[, c("x", "y")])
  points2 <- as.matrix(fg2[, c("x", "y")])

  # Calculate the Earth Mover's Distance
  emd_value <- emdw(points1, fg1$duration, points2, fg2$duration)

  # Normalize EMD to convert it into a similarity metric
  max_dist <- sqrt(screensize[1]^2 + screensize[2]^2)
  normalized_emd <- emd_value / max_dist

  # Convert to similarity
  similarity <- 1 - normalized_emd
  return(similarity)
}


#' @noRd
duration_diff_1d <- function(dur1, dur2, cds) {
  dur1 <- dur1[cds[,1]]
  dur2 <- dur2[cds[,2]]

  adiff <- abs(dur1-dur2)
  denom <- pmax(dur1, dur2)
  adiff/denom
}


#' @noRd
angle_diff_1d <- function(theta1, theta2, cds) {
  theta1 <- theta1[cds[,1]]
  theta2 <- theta2[cds[,2]]

  theta1 <- ifelse(theta1 < 0, pi + (pi + theta1), theta1)
  theta2 <- ifelse(theta2 < 0, pi + (pi + theta2), theta2)

  adiff <- abs(theta1-theta2)
  adiff <- ifelse(adiff > pi, 2 * pi - adiff, adiff)
  adiff
}

#' @noRd
vector_diff_1d <- function(x,y, v1, metric=c("l1", "l2"), cds) {
  metric <- match.arg(metric)
  x1 = x[[v1]]
  x2 = y[[v1]]

  #cds <- arrayInd(as.integer(path), c(nrow(x),nrow(y)))
  #cds <- cbind(cds[,2], cds[,1])

  if (metric == "l2") {
    abs(x1[cds[,1]] - x2[cds[,2]])
  } else {
    x1[cds[,1]] - x2[cds[,2]]
  }

}

#' @noRd
vector_diff_2d <- function(x,y, v1,v2, cds) {
  x1 = x[[v1]]
  x2 = y[[v1]]
  y1 = x[[v2]]
  y2 = y[[v2]]

  #cds <- arrayInd(as.integer(path), c(nrow(x),nrow(y)))
  #cds <- cbind(cds[,2], cds[,1])

  sqrt((x1[cds[,1]] - x2[cds[,2]])^2 + (y1[cds[,1]] - y2[cds[,2]])^2)

}


#' @noRd
create_graph <- function(x, y) {
  M <- proxy::dist(cbind(x$lenx, x$leny), cbind(y$lenx, y$leny))
  M_assignment <- matrix(seq(nrow(M) * ncol(M)), nrow(M), ncol(M), byrow = TRUE)

  nr <- nrow(M)
  nc <- ncol(M)

  edges_right <- if (nc > 1) {
    start <- M_assignment[, -nc, drop = FALSE]
    cbind(as.vector(start), as.vector(start + 1), as.vector(M[, -1, drop = FALSE]))
  } else {
    NULL
  }

  edges_down <- if (nr > 1) {
    start <- M_assignment[-nr, , drop = FALSE]
    cbind(as.vector(start), as.vector(start + nc), as.vector(M[-1, , drop = FALSE]))
  } else {
    NULL
  }

  edges_diag <- if (nr > 1 && nc > 1) {
    start <- M_assignment[-nr, -nc, drop = FALSE]
    cbind(as.vector(start), as.vector(start + nc + 1), as.vector(M[-1, -1, drop = FALSE]))
  } else {
    NULL
  }

  last_node <- M_assignment[nr, nc]
  edges_last <- matrix(c(last_node, last_node, 0), nrow = 1)

  out <- rbind(edges_right, edges_down, edges_diag, edges_last)

  g <- igraph::graph_from_data_frame(out[, 1:2])
  igraph::E(g)$weight <- out[, 3]
  spath <- igraph::shortest_paths(g, 1, max(igraph::V(g)), weights = igraph::E(g)$weight, output = "both", predecessors = TRUE)
  list(g = g, vpath = spath$vpath[[1]], epath = spath$epath[[1]], M = M, M_assignment = M_assignment, pred = spath$predecessors)
}



#' @export
install_multimatch <- function() {
  reticulate::py_install("multimatch_gaze", pip=TRUE)
}


#' Compute MultiMatch Metrics for Scanpath Similarity
#'
#' This function computes multiple similarity metrics between two scanpaths, including vector, direction, length, position, duration, and EMD-based position similarity.
#'
#' @param x A data frame representing the first scanpath. Must contain at least three columns: \code{x}, \code{y}, and \code{onset}, and at least three rows.
#' @param y A data frame representing the second scanpath. Must contain at least three columns: \code{x}, \code{y}, and \code{onset}, and at least three rows.
#' @param screensize A numeric vector of length 2 indicating the width and height of the screen in pixels.
#'
#' @details
#' The function computes six different similarity metrics between the scanpaths \code{x} and \code{y}:
#' \itemize{
#'   \item \code{mm_vector}: Similarity based on the 2D vectors between fixations.
#'   \item \code{mm_direction}: Similarity based on the direction (angle) of saccades between fixations.
#'   \item \code{mm_length}: Similarity based on the length of saccades between fixations.
#'   \item \code{mm_position}: Similarity based on the spatial position of fixations.
#'   \item \code{mm_duration}: Similarity based on the duration of fixations.
#'   \item \code{mm_position_emd}: Order-insensitive similarity based on the Earth Mover's Distance (EMD) between the spatial positions of fixations.
#' }
#'
#' The function ensures that both scanpaths have strictly increasing onset times and contain at least three fixations. It also normalizes the similarity scores to lie between 0 and 1, with higher values indicating greater similarity.
#'
#' @return A named numeric vector with the following elements:
#' \describe{
#'   \item{mm_vector}{Similarity based on the 2D vectors between fixations.}
#'   \item{mm_direction}{Similarity based on the direction (angle) of saccades between fixations.}
#'   \item{mm_length}{Similarity based on the length of saccades between fixations.}
#'   \item{mm_position}{Similarity based on the spatial position of fixations.}
#'   \item{mm_duration}{Similarity based on the duration of fixations.}
#'   \item{mm_position_emd}{Order-insensitive similarity based on the Earth Mover's Distance (EMD) between the spatial positions of fixations.}
#' }
#'
#' @references
#' Dewhurst, R., Nyström, M., Jarodzka, H., Foulsham, T., Johansson, R., & Holmqvist, K. (2012).
#' It depends on how you look at it: Scanpath comparison in multiple dimensions with MultiMatch, a vector-based approach.
#' Behavior research methods, 44, 1079-1100.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' scanpath1 <- data.frame(x = runif(10, 0, 500), y = runif(10, 0, 500), onset = cumsum(runif(10, 1, 5)))
#' scanpath2 <- data.frame(x = runif(10, 0, 500), y = runif(10, 0, 500), onset = cumsum(runif(10, 1, 5)))
#' screensize <- c(500, 500)
#' similarity_scores <- multi_match(scanpath1, scanpath2, screensize)
#' print(similarity_scores)
#' }
#'
#' @importFrom dplyr arrange
#' @importFrom stats median
#' @export
multi_match <- function(x,y, screensize) {

  if (any(diff(x$onset) <= 0)) {
    stop("multi_match: x `onset` vector must be strictly increasing")
  }

  if (any(diff(y$onset) <= 0)) {
    stop("multi_match: y `onset` vector must be strictly increasing")
  }

  if (nrow(x) < 3 || nrow(y) < 3) {
    warning("multi_match requires 3 or more coordinates in each scanpath, returning NAs")
    return(c(mm_vector=NA, mm_direction=NA,
             mm_length=NA, mm_position=NA,
             mm_duration=NA))

  }

  sacx <- x[1:(nrow(x)-1),]
  sacy <- y[1:(nrow(y)-1),]

  gout <- create_graph(sacx,sacy)
  p <- as.integer(gout$vpath)
  rnum <- ceiling(p / ncol(gout$M))
  cnum <-  p %% ncol(gout$M)
  cnum[cnum==0] <- ncol(gout$M)

  cds <- cbind(rnum, cnum)

  vector_d <- vector_diff_2d(sacx, sacy, "lenx", "leny", cds)
  vector_sim <-  1- median(vector_d) / (2 * sqrt(screensize[1]^2 + (screensize[2]^2)))

  direction_d <- angle_diff_1d(sacx$theta, sacy$theta, cds)
  direction_sim <- 1 - median(direction_d) / pi

  duration_d <- duration_diff_1d(sacx$duration, sacy$duration, cds)
  duration_sim <- 1 - median(duration_d)

  length_d <- abs(vector_diff_1d(sacx, sacy, "rho", metric="l1", cds))
  length_sim <- 1 - (median(length_d)) / (sqrt(screensize[1]^2 + (screensize[2]^2)))

  position_d <- vector_diff_2d(sacx, sacy, "x", "y", cds)
  position_sim <- 1 - (median(position_d)) / (sqrt(screensize[1]^2 + (screensize[2]^2)))

  emd_position_sim <- emd_position_similarity(sacx, sacy, screensize)

  c(mm_vector=vector_sim, mm_direction=direction_sim,
    mm_length=length_sim, mm_position=position_sim,
    mm_duration=duration_sim,
    mm_position_emd=emd_position_sim)
}


#' @noRd
py_multi_match <- function(fg1, fg2,
                        screensize,
                        grouping=FALSE,
                        tdir=25,
                        tdur=.05,
                        tamp=100) {

  if (!requireNamespace("reticulate")) {
    stop("multi_match requires access to python library `multimatch_gaze` via `reticulate`, please install")
  }

  if (!exists("mmgaze")) {
    mmgaze <<- try(reticulate::import("multimatch_gaze"))
    if (inherits(mmgaze, "try-error")) {
      stop("cannot load python module `multimatch_gaze`")
    }
  }

  fg1 <- fg1 %>% arrange(onset)
  fg2 <- fg2 %>% arrange(onset)

  fix1 <- fg1[, c("x", "y", "duration")]
  fix2 <- fg2[, c("x", "y", "duration")]

  colnames(fix1) <- c("start_x", "start_y", "duration")
  colnames(fix2) <- c("start_x", "start_y", "duration")

  ret <- mmgaze$docomparison(fix1, fix2, as.integer(screensize), grouping=grouping, TDir=tdir,TDur=tdur,TAmp=tamp)
  names(ret) <- c("mm_vector", "mm_direction", "mm_length", "mm_position", "mm_duration")
  ret
}

#' @noRd
#' Compute weighted Earth Mover's Distance (Wasserstein-1) between two 2-D point clouds.
#'
#' This helper tries to use the T4transport package (preferred) and falls back to the
#' transport package if available. Points are supplied as two-column matrices with
#' corresponding non-negative weights that need not sum to one (they will be
#' normalised internally).
#'
#' @param x Matrix of coordinates (n \times 2).
#' @param wx Numeric vector of weights for `x` (length n).
#' @param y Matrix of coordinates (m \times 2).
#' @param wy Numeric vector of weights for `y` (length m).
#' @return A single numeric value – the Earth Mover's Distance.
#' @keywords internal
emdw <- function(x, wx, y, wy, lambda = 0.01) {
  # Prefer emdist::emdw if available (package already in Imports)
  if (requireNamespace("emdist", quietly = TRUE)) {
    return(emdist::emdw(x, wx, y, wy))
  }

  # Fall back to T4transport if available
  if (requireNamespace("T4transport", quietly = TRUE)) {
    dmat <- proxy::dist(x, y)
    wx <- wx / sum(wx)
    wy <- wy / sum(wy)
    return(T4transport::sinkhornD(dmat, wx = wx, wy = wy, lambda = lambda)$distance)
  }

  # Finally, attempt with transport package
  if (requireNamespace("transport", quietly = TRUE)) {
    df1 <- data.frame(x = x[, 1], y = x[, 2], mass = wx / sum(wx))
    df2 <- data.frame(x = y[, 1], y = y[, 2], mass = wy / sum(wy))
    res <- transport::transport(df1, df2, p = 1)
    return(sum(res$dist * res$mass))
  }

  stop("Could not compute EMD: please install the 'emdist', 'T4transport', or 'transport' package.")
}

