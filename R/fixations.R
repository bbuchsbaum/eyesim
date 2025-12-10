
#' @export
summary.fixation_group <- function(object) {
  cen_x <- mean(object$x)
  cen_y <- mean(object$y)
  nfix <- nrow(object)
  sd_x <- sd(object$x)
  sd_y <- sd(object$y)
  tibble(cen_x=cen_x, cen_y=cen_y, sd_x=sd_x, sd_y=sd_y, nfix=nfix)
}

#' Create a Fixation Group Object
#'
#' This function creates a fixation group object containing fixation data
#' with x/y coordinates, duration, onset times, and an optional group index.
#'
#' @param x A numeric vector of x-coordinates for each fixation.
#' @param y A numeric vector of y-coordinates for each fixation.
#' @param duration A numeric vector of fixation durations. If missing, computed from onset differences.
#' @param onset A numeric vector of fixation onset times.
#' @param group An optional group identifier (default is 0).
#'
#' @return A tibble of class "fixation_group" with columns: index, x, y, duration, onset, group_index.
#'
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#'
#' @importFrom tibble tibble
#' @export
fixation_group <- function(x, y, duration, onset, group=0) {
  assert_that(length(x) == length(y))
  assert_that(length(x) == length(onset))

  if (missing(duration)) {
    duration <- c(diff(onset), 0)
  }

  assert_that(length(x) == length(duration))

  ret <- tibble(index=1:length(x),
                x=x,y=y, duration=duration,
                onset=onset, group_index=group)
  class(ret) <- c("fixation_group", class(ret))
  ret
}


#' @rdname rep_fixations
#' @export
rep_fixations.fixation_group <- function(x, resolution=100) {
  nreps <- as.integer(x$duration/ (1/resolution))
  nreps[nreps < 1] <- 1
  x <- x[rep(1:nrow(x), nreps),]
  x
}

#' @rdname sample_fixations
#' @param fast Logical. If TRUE (default), uses faster approximation method.
#' @export
sample_fixations.fixation_group <- function(x, time, fast=TRUE, ...) {


  ret <- if (fast) {
    x1 <- approx(x$onset, x$x, xout=time, method="constant", f=0)
    y1 <- approx(x$onset, x$y, xout=time, method="constant", f=0)
    data.frame(x=x1$y, y=y1$y, onset=time, duration=rep(1,length(time)))

  } else {
    purrr::map(time, function(t) {
      if (t < x$onset[1]) {
        c(x=NA,y=NA, onset=t, duration=NA)
      } else {
        delta <- t - x$onset
        valid <- which(delta >= 0)
        len <- length(valid)
        if (len == 0) {
          c(x=NA,y=NA, onset=t, duration=NA)
        } else {
          c(x=x$x[len], y=x$y[len], onset=t, duration=0)
        }
      }
    }) %>% map_dfr(bind_rows)
  }

  class(ret) <- c("sampled_fixation_group", "fixation_group", class(ret))
  ret
}

#' @export
coords.fixation_group <- function(x) {
  res <- cbind(x$x, x$y)
  colnames(res) <- c("x", "y")
  res
}

#' @rdname center
#' @export
center.fixation_group <- function(x, origin=NULL, ...) {
  if (is.null(origin)) {
    origin <- c(mean(x$x), mean(x$y))
  }
  out <- x %>% mutate(x = x - origin[1], y= y - origin[2])
  out
}

#' @rdname normalize
#' @export
normalize.fixation_group <- function(x, xbounds, ybounds, ...) {
  x %>% mutate(x=(x - xbounds[1])/xbounds[2], y=(y-ybounds[1])/ybounds[2])
}
