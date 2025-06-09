library(microbenchmark)
library(eyesim)

# Use the wynn_study dataset bundled with the package
# This dataset contains fixation groups for many subjects and images
# and is large enough to observe performance differences

data(wynn_study)

# Preserve the previous implementation for comparison

#' Original density_by using do()/rowwise()
density_by_old <- function(x, groups, sigma=50, xbounds=c(0,1000), ybounds=c(0,1000),
                           outdim=c(100,100), duration_weighted=TRUE, window=NULL,
                           min_fixations=2, keep_vars=NULL, fixvar="fixgroup",
                           result_name="density", ...) {
  rname <- rlang::sym(result_name)
  vars <- c(groups, keep_vars)
  if (!missing(groups) && !is.null(groups)) {
    ret <- x %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(groups))) %>%
      do({
        g <- do.call(rbind, .[[fixvar]])
        cbind(.[1, vars], tibble::tibble(!!fixvar := list(g)))
      }) %>%
      dplyr::rowwise() %>%
      do({
        d <- eye_density(.[[fixvar]], sigma,
                         xbounds = xbounds, ybounds = ybounds, outdim = outdim,
                         duration_weighted = duration_weighted, window = window,
                         min_fixations = min_fixations,
                         origin = attr(x, "origin"), ...)
        cbind(tibble::as_tibble(.[vars]),
              tibble::tibble(!!fixvar := list(.[[fixvar]]), !!rname := list(d)))
      })
  } else {
    fx <- do.call(rbind, x[[fixvar]])
    d <- eye_density(fx, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim,
                     duration_weighted=duration_weighted, window=window,
                     min_fixations=min_fixations,
                     origin=attr(x, "origin"), ...)
    ret <- tibble::tibble(!!fixvar := list(fx), !!rname := list(d))
  }
  if (any(vapply(ret[[result_name]], is.null, logical(1)))) {
    ret <- ret %>%
      dplyr::ungroup() %>%
      dplyr::filter(!vapply(.data[[result_name]], is.null, logical(1)))
  }
  ret
}

# Benchmark old vs new implementation
bench <- microbenchmark(
  old = density_by_old(wynn_study, groups = c("ImageNumber", "Subject")),
  new = density_by(wynn_study, groups = c("ImageNumber", "Subject")),
  times = 3
)
print(bench)
