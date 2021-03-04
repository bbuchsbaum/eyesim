

run_similarity_analysis <- function(ref_tab, source_tab, match_on, permutations, permute_on=NULL, method,
                                    refvar, sourcevar, window=NULL, ...) {
  args <- list(...)

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

  do_sim <- function(d1,d2,method, window=NULL) {
    if (!is.null(window)) {
      p <- purrr::partial(similarity, d1, d2, method=method, window=window)
      sim <- unlist(do.call(p, args))
    } else {
      p <- purrr::partial(similarity, d1, d2, method=method)
      sim <- unlist(do.call(p, args))
    }
  }

  ##browser()
  ret <- source_tab %>% furrr::future_pmap(function(...) {
    . <- list(...)
    d1 <- ref_tab[[refvar]][[.$matchind]]
    d2 <- .[[sourcevar]]

    #p <- purrr::partial(similarity, d1, d2, method=method)

    sim <- do_sim(d1,d2,method, window)

    if (permutations > 0) {

      mind <- if (!is.null(permute_on)) {
        ## limit matching indices to permute variable
        match_split[[as.character(.[[permute_on]])]]
      } else {
        matchind
      }

      if (permutations < length(mind)) {
        mind <- sample(mind, permutations)
      }

      elnum <- match(.$matchind, mind)

      if (!is.na(elnum)) {
        mind <- mind[-elnum]
      }

      psim <- do.call(rbind, lapply(mind, function(i) {
        d1p <- ref_tab[[refvar]][[i]]
        do_sim(d1p, d2, method, window)

      }))

      if (ncol(psim) > 1) {
        psim <- colMeans(psim)
        cnames <- c(names(sim), paste0("perm_", names(sim)), paste0(names(sim), "_diff"))
        c(unlist(sim), psim, unlist(sim)-psim) %>% set_names(cnames) %>% bind_rows()
      } else {
        tibble(eye_sim=sim, perm_sim=mean(psim), eye_sim_diff=sim - mean(psim))
      }
    } else {
        if (length(sim) == 1) {
          tibble(eye_sim=sim)
        } else {
          sim %>% bind_rows()
        }
    }

  }) %>% bind_rows()

  source_tab %>% bind_cols(ret)


}



#' fixation_similarity
#'
#' compute similarity between each fixation group in a \code{source_tab} with a matching fixation group in \code{ref_tab}
#'
#' @inheritsParam template_similarity
#' @param window
#' @export
fixation_similarity <- function(ref_tab, source_tab, match_on, permutations=0, permute_on=NULL,
                                method=c("sinkhorn", "overlap"),
                                refvar="fixgroup", sourcevar="fixgroup", window=NULL, ...) {
  if (!is.null(window) ) {
    assertthat::assert_that(window[2] > window[1])
  }
  message("fixation_similarity: similarity metric is ", method)

  method <- match.arg(method)
  run_similarity_analysis(ref_tab,source_tab, match_on, permutations, permute_on, method, refvar, sourcevar, window, ...)

}

#' @inheritParams template_similarity
#' @export
scanpath_similarity <- function(ref_tab, source_tab, match_on, permutations=0, permute_on=NULL,
                                method=c("multimatch"),
                                refvar="scanpath", sourcevar="scanpath", window=NULL, ...) {

  if (!is.null(window) ) {
    assertthat::assert_that(window[2] > window[1])
  }

  message("scan_similarity: similarity metric is ", method)
  run_similarity_analysis(ref_tab,source_tab, match_on, permutations, permute_on,
                          method, refvar, sourcevar, window, ...)

}



#' template_similarity
#'
#' compute similarity between each density map in a \code{source_tab} with a matching density map in \code{ref_tab}
#'
#' @param ref_tab the table of reference density maps
#' @param source_tab the table of source density maps
#' @param match_on the variable to match on
#' @param permute_on the variable used to stratify permutations
#' @param refvar the name of the variable containing density maps in reference table
#' @param sourcevar the name of the variable containing density maps in source table
#' @param method the similarity method
#' @param permutations the number of permutations for the baseline map
#' @param ... extra args to pass to `similarity` function
#' @export
template_similarity <- function(ref_tab, source_tab, match_on, permute_on = NULL, refvar="density", sourcevar="density",
                                method=c("spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov"),
                                permutations=10, ...) {


  method <- match.arg(method)
  message("template_similarity: similarity metric is ", method)
  run_similarity_analysis(ref_tab,source_tab, match_on, permutations, permute_on, method, refvar, sourcevar,...)
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
  if (!all(dim(z) == c(length(x), length(y)))) {
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
                                       normalize=TRUE, duration_weighted=FALSE, window=NULL,  origin=c(0,0)) {

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

  out <- if (duration_weighted || !is.null(window)) {
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
similarity.scanpath <- function(x, y, method=c("multimatch"),
                                      window=NULL,
                                      screensize=NULL,...) {

  if (!inherits(y, "scanpath")) {
    stop("`y` must be of type `scanpath`")
  }

  if (!is.null(window)) {
    #print(paste("window", window))
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

  if (is.null(screensize)) {
    stop("method `multi_match` requires a `screensize` argument (e.g. c(1000,1000)")
  }

  multi_match(x,y,screensize)

}


#' @export
similarity.fixation_group <- function(x, y, method=c("sinkhorn", "overlap"),
                                      window=NULL,
                                xdenom=1000, ydenom=1000, tdenom=3000,
                                tweight=.8,  lambda=.1, dthresh=40,
                                time_samples=NULL, screensize=NULL,...) {
  method <- match.arg(method)

  if (!inherits(y, "fixation_group")) {
    stop("`y` must be of type `fixation_group`")
  }

  if (!is.null(window)) {
      #print(paste("window", window))
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

  if (method == "sinkhorn") {
    xy1 <- cbind(x$x/xdenom, x$y/ydenom)
    xy2 <- cbind(y$x/xdenom, y$y/ydenom)

    xyt1 <- cbind(x$x/xdenom, x$y/ydenom, x$onset/tdenom * tweight)
    xyt2 <- cbind(y$x/xdenom, y$y/ydenom, y$onset/tdenom * tweight)

    d <- proxy::dist(xyt1, xyt2)

    #stw1 <- sigmoid(x$onset, a=a, b=b)
    #stw2 <- sigmoid(y$onset, a=a, b=b)
    xdur <- x$duration/sum(x$duration)
    ydur <- y$duration/sum(y$duration)

    d0 <- T4transport::sinkhornD(d,wx=xdur, wy=ydur, lambda=lambda)$distance
    1/(1+d0)
  } else if (method == "overlap") {
    if (is.null(time_samples)) {
      stop("method `overlap` requires a vector of `time_samples`")
    }
    fixation_overlap(x, y, dthresh=dthresh, time_samples=time_samples)
  }

}




#' @importFrom proxy simil
#' @export
similarity.density <- function(x, y, method=c("pearson", "spearman", "fisherz", "cosine",
                                              "l1", "jaccard", "dcov")) {
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
