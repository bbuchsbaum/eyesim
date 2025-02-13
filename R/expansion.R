
#' Estimate Scaling Parameters for Fixation Data
#'
#' This function estimates the scaling parameters for two sets of fixation data using the Hausdorff distance as an optimization objective.
#'
#' @param x A data frame or matrix containing the x coordinates of the first set of fixations.
#' @param y A data frame or matrix containing the y coordinates of the second set of fixations.
#' @param lower A numeric vector of length 2 specifying the lower bounds for the scaling parameters (default is c(0.1, 0.1)).
#' @param upper A numeric vector of length 2 specifying the upper bounds for the scaling parameters (default is c(10, 10)).
#' @param window A numeric vector of length 2 specifying the time window to restrict the fixations in `y` (default is NULL, which considers all fixations).
#'
#' @return A list containing the estimated scaling parameters.
#'
#' @keywords internal
#' @importFrom pracma hausdorff_dist
estimate_scale <- function(x, y, lower=c(.1,.1), upper=c(10,10), window) {

  if (!is.null(window)) {
    y <- subset(x, onset >= window[1] & onset < window[2])
  }

  if (nrow(x) == 0 || nrow(y) == 0) {
    return(list(par=c(1,1)))
  }
  #if (!is.null(window_y)) {
  #  x <- subset(y, onset >= window_y[1] & onset < window_y[2])
  #}

  cx <- as.matrix(x[,1:2])
  cy <- as.matrix(y[,1:2])
  par <- c(1,1)
  f <- function(p) {
    #browser()
    newy <- cy %*% diag(c(p[1],p[2]))
    pracma::hausdorff_dist(newy,cx)
  }

  ret <- optim(par, f, lower=lower, upper=upper, method="L-BFGS")
}


#' Match Scaling Parameters for Fixation Data
#'
#' This function matches the scaling parameters of fixation data between a reference table and a source table based on a common matching variable.
#'
#' @param ref_tab A data frame containing the reference fixation data.
#' @param source_tab A data frame containing the source fixation data.
#' @param match_on A string specifying the variable name in both `ref_tab` and `source_tab` to match on.
#' @param refvar A string specifying the variable name in `ref_tab` containing the reference fixations (default is "fixgroup").
#' @param sourcevar A string specifying the variable name in `source_tab` containing the source fixations (default is "fixgroup").
#' @param window A numeric vector of length 2 specifying the time window to restrict the fixations in the source fixation data (default is NULL, which considers all fixations).
#' @param ... Additional arguments passed to the `estimate_scale` function.
#'
#' @return A data frame containing the original source fixation data with additional columns for the matched scaling parameters.
#'
#' @examples
#' # Example usage of the match_scale function
#' ref_tab <- # reference fixation data
#' source_tab <- # source fixation data
#' matched_data <- match_scale(ref_tab, source_tab, match_on = "subject_id")
#'
#' @importFrom dplyr ungroup mutate filter bind_rows
#' @importFrom purrr pmap
match_scale <- function(ref_tab, source_tab, match_on,
                        refvar="fixgroup",sourcevar="fixgroup",
                        window,...) {
  if (!is.null(window) ) {
    assertthat::assert_that(window[2] > window[1])
  }

  matchind <- match(source_tab[[match_on]], ref_tab[[match_on]])

  source_tab <- source_tab %>% ungroup() %>% mutate(matchind=matchind)
  if (any(is.na(matchind))) {
    warning("did not find matching template map for all source maps. Removing non-matching elements.")
    source_tab <- source_tab %>% filter(!is.na(matchind))
    matchind <- matchind[!is.na(matchind)]
  }

  ret <- source_tab %>% purrr::pmap(function(...) {
    . <- list(...)
    d1 <- ref_tab[[refvar]][[.$matchind]]
    d2 <- .[[sourcevar]]
    res=estimate_scale(d1, d2, window=window)
    tibble(scale_x=res$par[1], scale_y=res$par[2])
  }) %>% bind_rows()

  source_tab %>% mutate(scale_x=ret$scale_x, scale_y=ret$scale_y)
}


#' @export
#' @family rescale
rescale.fixation_group <- function(x, sx,sy) {
  x %>% mutate(x=x*sx, y=y*sy)
}




