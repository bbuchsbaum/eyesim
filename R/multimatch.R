

mmgaze <- NULL


install_multimatch <- function() {
  reticulate::py_install("multimatch_gaze", pip=TRUE)
}


multi_match <- function(fg1, fg2,
                        screensize,
                        grouping=FALSE,
                        tdir=25,
                        tdur=.05,
                        tamp=100) {

  if (!requireNamespace("reticulate")) {
    stop("multi_match requires access to python library `multimatch_gaze` via `reticulate`, please install")
  }

  if (is.null(mmgaze)) {
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

  ret <- mmgaze$docomparison(fix1, fix2, as.integer(screensize), grouping, tdir,tdur,tamp)
  names(ret) <- c("mm_vector", "mm_direction", "mm_length", "mm_position", "mm_duration")
  ret
}

