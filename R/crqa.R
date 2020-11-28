
#' @import crqa
crqa <- function(fg1, fg2, radius=60, delay=1, embed=1, rescale=0, metric=c("euclidean", "manhattan")) {
  nr1 <- nrow(fg1)
  nr2 <- nrow(fg2)
  nr <- min(nr1,nr2)

  ts1 <- as.matrix(fg1[1:nr, 1:2])
  ts2 <- as.matrix(fg1[1:nr, 1:2])
  ret <- crqa(ts1,ts2, method="mdcrqa", radius=radius)
  ret
}
