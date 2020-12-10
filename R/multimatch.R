
#' @export
install_multimatch <- function() {
  reticulate::py_install("multimatch_gaze", pip=TRUE)
}

# vector_diff <- function(x,y)  {
#   lenx1 <- x$lenx
#   lenx2 <- y$lenx
#
#   leny1 <- x$leny
#   leny2 <- y$leny
#
#   M <- do.call(rbind, lapply(1:length(lenx), function(i) {
#     x_diff <- abs(rep(lenx1[i], length(lenx2)) - lenx2)
#     y_diff <- abs(rep(leny1[i], length(leny2)) - leny2)
#     sqrt(x_diff^2 + y_diff^2)
#   }) )
#
#
#
# }


#' multi_match
#'
#' @param fg1 the first \code{fixation_group}
#' @param fg2 the second \code{fixation_group}
#' @param screensize a two element vector indicating screen size
#' @param tdir
multi_match <- function(fg1, fg2,
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

