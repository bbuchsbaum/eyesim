


fixation_similarity <- function(ref_tab, source_tab, match_on, permutations=0, method="sinkhorn", window=NULL, ...) {

  matchind <- match(source_tab[[match_on]], ref_tab[[match_on]])

  source_tab <- source_tab %>% ungroup() %>% mutate(matchind=matchind)

  if (any(is.na(matchind))) {
    warning("did not find matching template map for all source maps. Removing non-matching elements.")
    source_tab <- source_tab %>% filter(!is.na(matchind))
    matchind <- matchind[!is.na(matchind)]
  }

  if ( !is.null(window) ) {
    assertthat::assert_that(window[2] > window[1])
  }

  ret <- source_tab %>% rowwise() %>% do( {

    d1 <- ref_tab[["fixgroup"]][[.$matchind]]
    d2 <- .[["fixgroup"]]


    sim <- similarity(d1,d2, method=method, window=window)

    if (permutations > 0) {
      #browser()
      mind <- if (permutations < length(mind)) {
        mind <- sample(matchind, permutations)
      } else {
        matchind
      }

      mind <- mind[!mind %in% .$matchind]

      psim <- mean(sapply(mind, function(i) {
        similarity(ref_tab[[refvar]][[i]], d2, method=method)
      }))

      data.frame(eye_sim=sim, perm_sim=psim, eye_sim_diff=sim-psim)
    } else {
      sim <- similarity(d1,d2, method=method, window=window, ...)
      data.frame(eye_sim=sim)
    }
  })

  out <- source_tab %>% mutate(eye_sim=ret$eye_sim)


}


#' template_similarity
#'
#' compute similarity betwen each density map in a \code{source_tab} with a matching density map in \code{ref_tab}
#'
#' @param ref_tab the table of reference density maps
#' @param source_tab the table of source density maps
#' @param match_on the variable to match on
#' @param permute_on the variable used to stratify permutations
#' @param refvar the name of the variable containing density maps in reference table
#' @param sourcvar the name of the variable containing density maps in source table
#' @param method the similarity method
#' @param permutations the number of permutations for the baseline map
#' @export
template_similarity <- function(ref_tab, source_tab, match_on, permute_on = NULL, refvar="density", sourcevar="density",
                                method=c("spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov"),
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

  if (!is.null(permute_on)) {
    assertthat::assert_that(permute_on %in% names(source_tab) && assertthat::assert_that(permute_on %in% names(ref_tab)))
    match_split <- split(matchind, source_tab[[permute_on]])
  }

  ret <- source_tab %>% rowwise() %>% do( {

    d1 <- ref_tab[[refvar]][[.$matchind]]
    d2 <- .[[sourcevar]]

    sim <- similarity(d1,d2, method=method)

    if (permutations > 0) {

      mind <- if (!is.null(permute_on)) {
        match_split[[as.character(.[[permute_on]])]]
      } else {
        matchind
      }

      mind <- if (permutations < length(mind)) {
        mind <- sample(matchind, permutations)
      } else {
        matchind
      }

      mind <- mind[!mind %in% .$matchind]

      psim <- mean(sapply(mind, function(i) {
        similarity(ref_tab[[refvar]][[i]], d2, method=method)
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
sample_density.density <- function(x, fix, times=NULL) {
  if (is.null(times)) {
    cds <- as.matrix(cbind(fix$x, fix$y))
    ix <- sapply(cds[,1], function(p) which.min(abs(x$x - p)))
    iy <- sapply(cds[,2], function(p) which.min(abs(x$y - p)))
    data.frame(z=x$z[cbind(ix,iy)], time=fix$onset)
  } else {
    fg <- sample_fixations(fix, times)
    cds <- as.matrix(cbind(fg$x, fg$y))
    res <- lapply(1:nrow(cds), function(i) {
      i1 <- which.min(abs(x$x - cds[i,1]))
      i2 <- which.min(abs(x$y - cds[i,2]))
      if (length(i1) > 0) {
        x$z[i1,i2]
      } else {
        NA
      }
    })

    data.frame(z=unlist(res), time=times)
  }
}


#' @export
gen_density <- function(x,y,z) {
  if (!all(dim(x) == c(length(x), length(y)))) {
    stop("length of x and y must equal nrow(z) and ncol(z)")
  }

  out <- list(
    x=x,
    y=y,
    z=z)

  class(out) <- c("density", "list")
  out
}

#' @export
get_density.eye_density <- function(x, ...) {
  x$z
}

#' @export
#' @importFrom purrr cross_df
#' @importFrom dplyr mutate
as.data.frame.eye_density <- function(x, ...) {
  z <- x$z
  kde_df <- x %>%
    .[c("x", "y")] %>%
    purrr::cross_df() %>%
    dplyr::mutate(z = as.vector(z))

  kde_df
}

#' @export
print.eye_density <- function(x) {
  cat("fixation density map", "\n")
  cat("xlim: ", range(x$x), "\n")
  cat("ylim: ", range(x$y), "\n")
  cat("z range: ", range(x$z), "\n")
}



#' @export
#' @importFrom MASS kde2d
eye_density.fixation_group <- function(x, sigma=50, xbounds=c(min(x$x), max(x$x)), ybounds=c(min(x$y), max(x$y)),
                                       outdim=c(xbounds[2] - xbounds[1], ybounds[2] - ybounds[1]),
                                       normalize=TRUE, duration_weighted=FALSE, window=NULL, angular=FALSE, angle_bins=25, origin=c(0,0)) {

  if (!is.null(window)) {
    assertthat::assert_that(window[2] > window[1])
    assert_that("onset" %in% colnames(x))
    x <- filter(x, onset >= window[1] & onset < window[2])
    assertthat::assert_that(nrow(x) > 0)
  }


  wts <- if (duration_weighted) {
    x$duration
  } else {
    rep(1, nrow(x))
  }



  out <- if (angular) {
    xrep <- rep_fixations(x, 50)
    theta <- to_angle(xrep$x - origin[1], xrep$y - origin[2])
    group <- cut(theta, seq(-pi/2, pi/2, length.out=angle_bins))
    ang <- as.vector(table(group))
    list(z=ang)
  } else if (duration_weighted || !is.null(window)) {
    xrep <- rep_fixations(x, 50)
    xrep <- x
    kde2d(xrep$x, xrep$y, h=sigma, n=outdim, lims=c(xbounds, ybounds))
  } else {
    kde2d(x$x, x$y, n=outdim, h=sigma, lims=c(xbounds, ybounds))
  }

  if (normalize) {
    out$z <- out$z/sum(out$z)
  }

  out$z <- zapsmall(out$z)

  out$fixgroup <- x

  class(out) <- c("eye_density", "density", "list")
  out

  #ret <- blurpoints(cbind(x$x, x$y), sigma=sigma, xbounds=xbounds, ybounds=ybounds, resolution=resolution, weights=wts, normalize=normalize)
}

Ops.eye_density <- function(e1,e2) {
    op = .Generic[[1]]
    switch(op,
           `-` = {
             delta <- e1$z - e2$z
             #delta <- (delta - min(delta))/(max(delta) - min(delta))
             structure(list(x=e1$x, y=e2$y, z=delta, fixgorup=rbind(e1$fixgroup, e2$fixgroup)),
                       class=c("eye_density_delta", "eye_density", "density", "list"))

           },
           `+` = {
             add = (e1$z + e2$z)/2
             structure(list(x=e1$x, y=e2$y, z=add, fixgorup=rbind(e1$fixgroup, e2$fixgroup)),
                       class=c("eye_density_add", "eye_density", "density", "list"))

           },
           `/` = {
             div = log(e1$z/e2$z)
             structure(list(x=e1$x, y=e2$y, z=div, fixgorup=rbind(e1$fixgroup, e2$fixgroup)),
                       class=c("eye_density_div", "eye_density", "density", "list"))

           },
           stop("undefined operation")
    )
}




to_angle <- function(x, y) {
  r <- sqrt(x^2 + y^2)
  asin(x/r)
}


sigmoid <- function (x, a = 1, b = 0)  {
  if (length(x) == 0)
    return(c())
  stopifnot(is.numeric(x), is.numeric(a), is.numeric(b))
  a <- a[1]
  b <- b[1]
  1/(1 + exp(-a * (x - b)))
}

#' @export
similarity.fixation_group <- function(x, y, method="sinkhorn",window=NULL,
                                xdenom=1000, ydenom=1000, tdenom=3000,
                                tweight=.8,  lambda=200) {

  if (!is.null(window)) {
      x <- filter(x, onset >= window[1] & onset < window[2])
      y <- filter(y, onset >= window[1] & onset < window[2])
  }

  if (nrow(x) == 0) {
    warning("no observations in 'x'")
    return(NA)
  }

  if (nrow(y) == 0) {
    warning("no observations in 'y'")
    return(NA)
  }

  xy1 <- cbind(x$x/xdenom, x$y/ydenom)
  xy2 <- cbind(y$x/xdenom, y$y/ydenom)

  #spd <- proxy::dist(xy1,xy2)
  #td <- proxy::dist(x$onset/tdenom, y$onset/tdenom)
  #d <- spd + tweight*td

  xyt1 <- cbind(x$x/xdenom, x$y/ydenom, x$onset/tdenom * tweight)
  xyt2 <- cbind(y$x/xdenom, y$y/ydenom, y$onset/tdenom * tweight)
  d <- proxy::dist(xyt1, xyt2)

  #stw1 <- sigmoid(x$onset, a=a, b=b)
  #stw2 <- sigmoid(y$onset, a=a, b=b)

  d0 <- T4transport::sinkhornD(d,wx=x$duration, wy=y$duration, lambda=lambda)$distance
  -log(d0)
}




#' @importFrom proxy simil
#' @export
similarity.density <- function(x, y, method=c("pearson", "spearman", "fisherz", "cosine", "l1", "jaccard", "dcov")) {
  method=match.arg(method)

  if (inherits(y, "density")) {
    y <- y$z
  }

  compute_similarity(x$z, as.vector(y), method)

}

compute_similarity <- function(x,y, method=c("pearson", "spearman", "fisherz", "cosine", "l1", "jaccard", "dcov")) {
  method=match.arg(method)
  if (method=="pearson" || method == "spearman") {
    cor(as.vector(x), as.vector(y), method=method)
  } else if (method == "fisherz") {
    atanh(cor(as.vector(x), as.vector(y)))
  } else if (method == "cosine") {
    proxy::simil(as.vector(x), as.vector(y), method="cosine", by_rows=FALSE)[,]
  } else if (method == "l1") {
    z <- as.vector(x)
    y <- as.vector(y)
    x1 <- z/sum(z)
    x2 <- y/sum(y)
    1-(proxy::dist(x1, x2, method="Manhattan", by_rows=FALSE)[,])
  } else if (method == "jaccard") {
    proxy::simil(as.vector(x), as.vector(y), method="eJaccard", by_rows=FALSE)[,]
  } else if (method=="dcov") {
    z <- as.vector(x)
    y <- as.vector(y)
    x1 <- z/sum(z)
    x2 <- y/sum(y)
    energy::dcor(x1,x2)
  }

}



kde2d_weighted <- function (x, y, h, n = 25, lims = c(range(x), range(y)), w)
{
  nx <- length(x)
  if (length(y) != nx)
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  #browser()
  n <- rep(n, length.out = 2L)
  gx <- seq(lims[1], lims[2], length = n[1])
  gy <- seq(lims[3], lims[4], length = n[2])
  h <- if (missing(h))
    c(bandwidth.nrd(x), bandwidth.nrd(y))
  else rep(h, length.out = 2L)
  if (any(h <= 0))
    stop("bandwidths must be strictly positive")
  if (missing(w))
    w <- numeric(nx) + 1
  h <- h/4
  ax <- outer(gx, x, "-")/h[1]
  ay <- outer(gy, y, "-")/h[2]
  z <- (matrix(rep(w, n), nrow = n, ncol = nx, byrow = TRUE) *
          matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n,
                                                 nx))/(sum(w) * h[1] * h[2])
  return(list(x = gx, y = gy, z = z))
}
