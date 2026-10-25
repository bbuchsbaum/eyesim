
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

  n_emb <- nr - (embed - 1) * delay
  if (n_emb < 1L) {
    stop("crqa: embed/delay leave no embedded points for ", nr, " fixations.")
  }
  # crqa's own embedding drops a dimension when only one embedded point is
  # left (it subsets without drop = FALSE), producing a spurious 2x2 plot.
  # Embed that single point here and hand crqa an already-embedded series.
  if (embed > 1 && n_emb == 1L) {
    lagged <- function(ts) do.call(cbind, lapply(seq_len(embed) - 1L, function(j)
      ts[1L + j * delay, , drop = FALSE]))
    ts1 <- lagged(ts1)
    ts2 <- lagged(ts2)
    embed <- 1
  }
  ret <- crqa::crqa(ts1, ts2, method="mdcrqa", radius=radius, delay=delay,
                    embed=embed, rescale=rescale, metric=metric)
  ret
}
