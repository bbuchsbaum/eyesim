
#' @noRd
#' @keywords internal
crqa <- function(fg1, fg2, radius=60, delay=1, embed=1, rescale=0, metric=c("euclidean", "manhattan")) {
  if (!requireNamespace("crqa", quietly = TRUE)) {
    stop("Package 'crqa' is required for this function. Install it with install.packages('crqa').")
  }
  metric <- match.arg(metric)
  nr1 <- nrow(fg1)
  nr2 <- nrow(fg2)
  nr <- min(nr1,nr2)
  if (nr == 0L) {
    stop("crqa requires non-empty fixation groups.")
  }

  ts1 <- as.matrix(fg1[seq_len(nr), c("x", "y")])
  ts2 <- as.matrix(fg2[seq_len(nr), c("x", "y")])
  ret <- crqa::crqa(ts1, ts2, method="mdcrqa", radius=radius, delay=delay,
                    embed=embed, rescale=rescale, metric=metric)
  ret
}
