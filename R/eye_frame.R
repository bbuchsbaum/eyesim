
#' @importFrom dplyr do group_by select filter
#' @import magrittr
#' @importFrom assertthat assert_that
#' @param x
#' @param y
#' @param duration
#' @param onset
#' @param groupvar
#' @param vars
#' @param data
#' @param clip_bounds
#' @param relaive_coords
#' @export
eye_table <- function(x, y, duration, onset, groupvar, vars=NULL, data,
                      clip_bounds=c(0,1280, 0,1280), relative_coords=TRUE) {

  assertthat::assert_that(inherits(data, "data.frame"))

  colmapping <- c("x","y","duration", "onset")
  names(colmapping) <- c(x,y,duration,onset)


  data <- data %>% rename_with(.cols=all_of(c(x,y,duration,onset)),
                               .fn = function(x){colmapping[x]}) %>% as_tibble()

  data <- if (is.null(vars)) {
    #data %>% select_at("x","y","duration", "onset", .dots=c(groupvar))
    data %>% select_at(c("x","y","duration", "onset", groupvar))
  } else {
    data %>% select_at(c("x","y","duration", "onset", vars, groupvar))
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
    group_by_at(groupvar) %>%
    do({
        cbind(.[1,], tibble(fixgroup=list(fixation_group(.[["x"]], .[["y"]], .[["duration"]], .[["onset"]]))))
    }) %>% select_at(c("fixgroup", vars))

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







