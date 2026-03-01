
#' @importFrom stats sd
#' @export
summary.fixation_group <- function(object, ...) {
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
#' @importFrom stats approx
#' @importFrom purrr map_dfr
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
  x %>% mutate(x=(x - xbounds[1])/(xbounds[2] - xbounds[1]),
               y=(y - ybounds[1])/(ybounds[2] - ybounds[1]))
}


#' @export
print.fixation_group <- function(x, ...) {
  cat("Fixation group:", nrow(x), "fixations\n")
  cat("  x range: [", round(min(x$x), 1), ", ", round(max(x$x), 1), "]\n", sep = "")
  cat("  y range: [", round(min(x$y), 1), ", ", round(max(x$y), 1), "]\n", sep = "")
  cat("  duration range: [", round(min(x$duration), 1), ", ", round(max(x$duration), 1), "]\n", sep = "")
  cat("  onset range: [", round(min(x$onset), 1), ", ", round(max(x$onset), 1), "]\n", sep = "")
  invisible(x)
}


#' Concatenate Fixation Groups
#'
#' Combines multiple fixation_group objects into one by row-binding. Onset
#' times are shifted so that each subsequent group continues from where the
#' previous one ended.
#'
#' @param ... \code{fixation_group} objects to concatenate.
#' @param recursive Ignored (present for S3 method compatibility).
#'
#' @return A single \code{fixation_group} object.
#' @export
#'
#' @examples
#' fg1 <- fixation_group(x = c(1, 2), y = c(3, 4),
#'                        onset = c(0, 100), duration = c(100, 100))
#' fg2 <- fixation_group(x = c(5, 6), y = c(7, 8),
#'                        onset = c(0, 100), duration = c(100, 100))
#' combined <- c(fg1, fg2)
c.fixation_group <- function(..., recursive = FALSE) {
  args <- list(...)
  # Filter NULLs
  args <- args[!vapply(args, is.null, logical(1))]
  if (length(args) == 0) return(NULL)

  # Validate all are fixation_group
  is_fg <- vapply(args, inherits, logical(1), what = "fixation_group")
  if (!all(is_fg)) {
    stop("All arguments to c.fixation_group() must be fixation_group objects.")
  }

  if (length(args) == 1) return(args[[1]])

  # Shift onsets sequentially
  result <- args[[1]]
  for (i in 2:length(args)) {
    fg <- args[[i]]
    offset <- max(result$onset) + max(result$duration[nrow(result)])
    fg$onset <- fg$onset + offset
    result <- rbind(result, fg)
  }

  # Reindex
  result$index <- seq_len(nrow(result))
  class(result) <- c("fixation_group", setdiff(class(result), "fixation_group"))
  result
}
