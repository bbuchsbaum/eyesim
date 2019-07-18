

#' @importFrom dplyr do group_by select filter
#' @import magrittr
#' @importFrom assertthat assert_that
#' @export
eye_table <- function(x, y, duration, onset, groupvar, vars=NULL, data, clip_bounds=c(0,1280, 0,1280), relative_coords=TRUE) {

  data <- data %>% rename_(x=x, y=y, duration=duration, onset=onset)

  data <- if (is.null(vars)) {
    data %>% select_("x","y","duration", "onset", .dots=c(groupvar))
  } else {
    data %>% select_("x","y","duration", "onset", .dots=c(vars, groupvar))
  }


  xdir <- sign(clip_bounds[2] - clip_bounds[1])
  ydir <- sign(clip_bounds[4] - clip_bounds[3])
  xr <- sort(c(clip_bounds[1], clip_bounds[2]))
  yr <- sort(c(clip_bounds[3], clip_bounds[4]))

  data <- data %>% filter(x >= xr[1] & x <= xr[2]
                        & y >= yr[1] & y <= yr[2])


  if (relative_coords) {
    data <- data %>% mutate(x = (x - clip_bounds[1]) * xdir, y = (y - clip_bounds[3]) * ydir)
  }



  res <- data %>%
    group_by_(.dots=groupvar) %>%
    do({
        cbind(.[1,], tibble(fixgroup=list(fixation_group(.[["x"]], .[["y"]], .[["duration"]], .[["onset"]]))))
    }) %>% select_(.dots=c("fixgroup", vars))

  class(res) <- c("eye_table", class(res))

  if (relative_coords) {
    xr1 <- (xr[1] - clip_bounds[1]) * xdir
    xr2 <- (xr[2] - clip_bounds[1]) * xdir
    yr1 <- (yr[1] - clip_bounds[3]) * ydir
    yr2 <- (yr[2] - clip_bounds[3]) * ydir
    attr(res, "origin") <- c((xr1+xr2)/2, (yr1+yr2)/2)
  } else {
    attr(res, "origin") <- c((xr[1]+xr[2])/2, (yr[1]+yr[2])/2)
  }
  res

}








#' @export
coords.fixation_group <- function(x) {
  res <- cbind(x$x, x$y)
  colnames(res) <- c("x", "y")
  res
}


#' @importFrom tibble tibble
#' @export
fixation_group <- function(x, y, duration, onset) {
  ret <- tibble(x=x,y=y, duration=duration, onset=onset)
  class(ret) <- c("fixation_group", class(ret))
  ret
}


#' @export
rep_fixations.fixation_group <- function(x, resolution=100) {
  nreps <- as.integer(x$duration/resolution)
  nreps[nreps < 1] <- 1
  x <- x[rep(1:nrow(x), nreps),]
  x
}

#' @export
sample_fixations.fixation_group <- function(x, time) {

  ret <- purrr::map(time, function(t) {
    if (t < x$onset[1]) {
      c(x=NA,y=NA, onset=t, duration=NA)
    } else {
      delta <- t- x$onset
      valid <- which(delta >= 0)
      len <- length(valid)
      if (len == 0) {
        c(x=NA,y=NA, onset=t, duration=NA, delta=NA)
      } else {
        c(x=x$x[len], y=x$y[len], onset=t, duration=t-x$onset[len])
      }
    }
  }) %>% map_dfr(bind_rows)

  class(ret) <- c("sampled_fixation_group", "fixation_group", class(ret))
  ret
}


