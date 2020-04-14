#' template_similarity
#'
#' @param ref_tab
#' @param source_tab
#' @param match_on
#' @param method
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
  kde_df <- x %>%
    .[c("x", "y")] %>%
    purrr::cross_df() %>%
    dplyr::mutate(density = as.vector(x$z))
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
                                       normalize=TRUE, duration_weighted=FALSE, angular=FALSE, angle_bins=25, origin=c(0,0)) {
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
  } else if (duration_weighted) {
    xrep <- rep_fixations(x, 50)
    kde2d(xrep$x, xrep$y, n=outdim, h=sigma, lims=c(xbounds, ybounds))
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


to_angle <- function(x, y) {
  r <- sqrt(x^2 + y^2)
  asin(x/r)
}




#' @importFrom proxy simil
#' @export
similarity.density <- function(x, y, method=c("pearson", "spearman", "fisherz", "cosine", "l1", "jaccard", "dcov")) {
  method=match.arg(method)

  if (inherits(y, "density")) {
    y <- y$z
  }

  compute_similarity(x$z, as.vector(y))

}

compute_similarity <- function(x,y, method=c("pearson", "spearman", "fisherz", "cosine", "l1", "jaccard", "dcov")) {
  method=match.arg(method)
  if (method=="pearson" || method == "spearman") {
    cor(as.vector(x), as.vector(y), method=method)
  } else if (method == "fisherz") {
    atanh(cor(as.vector(x), as.vector(y), method=method))
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
