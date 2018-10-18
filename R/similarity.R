#' template_similarity
#'
#' @param ref_tab
#' @param source_tab
#' @param match_on
#' @param method
#' @export
template_similarity <- function(ref_tab, source_tab, match_on,
                                method=c("spearman", "pearson", "cosine", "l1", "jaccard", "dcov"),
                                permutations=10) {


  method <- match.arg(method)
  message("template_similarity: similarity metric is ", method)



  matchind <- match(source_tab[[match_on]], ref_tab[[match_on]])

  source_tab <- source_tab %>% ungroup() %>% mutate(matchind=matchind)

  if (any(is.na(matchind))) {
    warning("did not find matching template map for all source maps. Removing non-matching elements.")
    source_tab <- source_tab %>% filter(!is.na(matchind))
    matchind <- matchind[!is.na(matchind)]
  }

  ret <- source_tab %>% rowwise() %>% do( {

    d1 <- ref_tab$density[[.$matchind]]
    d2 <- .$density

    sim <- similarity(d1,d2, method=method)

    if (permutations > 0) {
      mind <- sample(matchind, permutations)
      psim <- mean(sapply(mind, function(i) {
        similarity(ref_tab$density[[i]], d2, method=method)
      }))

      data.frame(eye_sim=sim, perm_sim=psim, eye_sim_diff=sim-psim)
    } else {
      data.frame(eye_sim=sim)
    }
  })

  if (permutations > 0) {
    source_tab %>% mutate(eye_sim=ret$eye_sim, perm_sim=ret$perm_sim, eye_sim_diff=ret$eye_sim_diff)
  } else {
    source_tab %>% mutate(eye_sim=ret$eye_sim)
  }
}

#' @export
#' @importFrom MASS kde2d
eye_density.fixation_group <- function(x, sigma=50, xbounds=c(min(x$x), max(x$x)), ybounds=c(min(x$y), max(x$y)),
                                       outdim=c(xbounds[2] - xbounds[1], ybounds[2] - ybounds[1]),
                                       normalize=TRUE, duration_weighted=FALSE, angular=FALSE, angle_bins=25, origin=c(0,0)) {
  wts <- if (duration_weighted) {
    x$duration
  } else {
    rep(1, nrow(x))
  }

  print(angular)
  out <- if (angular) {
    xrep <- rep_fixations(x, 50)
    theta <- to_angle(xrep$x - origin[1], xrep$y - origin[2])
    group <- cut(theta, seq(-pi/2, pi/2, length.out=angle_bins))
    ang <- as.vector(table(group))
    list(z=ang)
  } else if (duration_weighted) {
    xrep <- rep_fixations(x, 50)
    kde2d(xrep$x, xrep$y, n=outdim, h=sigma, lims=c(xbounds, ybounds))
  } else {
    kde2d(x$x, x$y, n=outdim, h=sigma, lims=c(xbounds, ybounds))
  }

  if (normalize) {
    out$z <- out$z/sum(out$z)
  }

  out$fixgroup <- x

  class(out) <- c("eye_density", "list")
  out

  #ret <- blurpoints(cbind(x$x, x$y), sigma=sigma, xbounds=xbounds, ybounds=ybounds, resolution=resolution, weights=wts, normalize=normalize)
}


to_angle <- function(x, y) {
  r <- sqrt(x^2 + y^2)
  asin(x/r)
}


#' @importFrom proxy simil
#' @export
similarity.eye_density <- function(x, y, method=c("pearson", "spearman", "cosine", "l1", "jaccard", "dcov")) {
  method=match.arg(method)

  if (inherits(y, "eye_density")) {
    y <- y$z
  }

  if (method=="pearson" || method == "spearman") {
    cor(as.vector(x$z), as.vector(y), method=method)
  } else if (method == "cosine") {
    proxy::simil(as.vector(x$z), as.vector(y), method="cosine", by_rows=FALSE)[,]
  } else if (method == "l1") {
    z <- as.vector(x$z)
    y <- as.vector(y)
    x1 <- z/sum(z)
    x2 <- y/sum(y)
    1-(proxy::dist(x1, x2, method="Manhattan", by_rows=FALSE)[,])
  } else if (method == "jaccard") {
    proxy::simil(as.vector(x$z), as.vector(y), method="eJaccard", by_rows=FALSE)[,]
  } else if (method=="dcov") {
    z <- as.vector(x$z)
    y <- as.vector(y)
    x1 <- z/sum(z)
    x2 <- y/sum(y)
    energy::dcor(x1,x2)
  }

}
