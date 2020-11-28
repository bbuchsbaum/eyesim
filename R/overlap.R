

fixation_overlap <- function(x, y, dthresh=60, time_samples=seq(0,max(x$onset), by=20), dist_method=c("euclidean", "manhattan")) {
  method <- match.arg(dist_method)

  fx1 <- sample_fixations(x, time_samples)
  fx2 <- sample_fixations(y, time_samples)

  d <- proxy::dist(fx1[,1:2], fx2[,1:2], pairwise=TRUE, method=method)
  overlap <- sum(d[!is.na(d)] < dthresh)
  perc <- overlap/length(d)
}
