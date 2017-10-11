
#' @export
#' @importFrom dplyr do group_by select filter
#' @importFrom assertthat assert_that
create_eye_table <- function(x, y, duration, onset, groupvar, vars=NULL, data) {

  data <- if (is.null(vars)) {
    data %>% select_(.dots=c(x,y,duration, onset, groupvar))
  } else {

    data %>% select_(.dots=c(x,y,duration, onset, vars, groupvar))
  }


  res <- data %>%
    group_by_(.dots=groupvar) %>%
    do(
      {

        cbind(.[1,], tibble(fixgroup=list(fixation_group(.[[x]], .[[y]], .[[duration]], .[[onset]]))))
      }) %>% select_(.dots=c("fixgroup", vars))

  class(res) <- c("eye_table", class(res))
  res

}



#' @importFrom ggplot2 ggplot
plot.fixation_group <- function(x, type=c("contour", "density", "density_alpha"), bandwidth=100,
                                xlim=range(x$x),
                                ylim=range(x$y), size_points=TRUE) {
  type <- match.arg(type)

  if (size_points) {
    ps <- (x$duration - min(x$duration)) / (max(x$duration) - min(x$duration))
    x$psize <- ps*2 + 1
  }

  p <- ggplot(data=x, aes(x=x, y=y)) +
    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2])

  if (size_points) {
    p <- p + geom_point(aes(size=psize))
  } else {
    p <- p + geom_point()
  }

  p <- if (type== "contour") {
    p + stat_density_2d(aes(colour = ..level..), h=bandwidth)
  } else if (type == "density") {
    p + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE, h=bandwidth)
  } else {
    p + stat_density2d(aes(alpha = ..density..), geom = "raster", contour = FALSE, h=bandwidth)
  }

  p
}

#' @export
coords.fixation_group <- function(x) {
  res <- cbind(x$x, x$y)
  colnames(res) <- c("x", "y")
  res
}


#' @export
density_by <- function(x, groups, sigma=50, xbounds=c(0, 1000), ybounds=c(0, 1000), outdim=c(100,100), ...) {
  ret <- x %>% group_by_(.dots=groups) %>% do( {
    g <- do.call(rbind, .$fixgroup)
    cbind(.[1,groups],tibble(fixgroup=list(g)))
  }) %>% rowwise() %>% do( {
    d <- eye_density(.$fixgroup, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim, ...)
    cbind(as_tibble(.[groups]), tibble( fixgroup=list(.$fixgroup), density=list(d)))
  })

  ret

}

#' @param ref_tab
#' @param source_tab
#' @param match_on
#' @param method
#' @export
target_similarity <- function(ref_tab, source_tab, match_on, method="pearson") {
  matchind <- match(delay_dens[[match_on]], study_dens[[match_on]])
  source_tab <- source_tab %>% ungroup() %>% mutate(matchind=matchind)
  ret <- source_tab %>% rowwise() %>% do( {
    d1 <- ref_tab$density[[.$matchind]]
    d2 <- .$density
    sim <- similarity(d1,d2, method=method)
    data.frame(eye_sim=sim)
  })

  source_tab %>% mutate(eye_sim=ret$eye_sim)

}


eye_density.fixation_group <- function(x, sigma=50, xbounds=c(min(x$x), max(x$x)), ybounds=c(min(x$y), max(x$y)),
                                       outdim=c(xbounds[2] - xbounds[1], ybounds[2] - ybounds[1]),
                                       duration_weighted=FALSE, normalize=TRUE) {
  wts <- if (duration_weighted) {
    x$duration
  } else {
    rep(1, nrow(x))
  }

  out <- if (duration_weighted) {
    xrep <- rep_fixations(x, 50)
    kde2d(xrep$x, xrep$y, n=outdim, h=sigma, lims=c(xbounds, ybounds))
  } else {
    kde2d(x$x, x$y, n=outdim, h=sigma, lims=c(xbounds, ybounds))
  }

  if (normalize) {
    out$z <- out$z/sum(out$z)
  }

  class(out) <- c("eye_density", "list")
  out

  #ret <- blurpoints(cbind(x$x, x$y), sigma=sigma, xbounds=xbounds, ybounds=ybounds, resolution=resolution, weights=wts, normalize=normalize)
}


cosine_sim <- function(x,y){
  dot.prod <- x%*%y
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

#' @importFrom proxy
similarity.eye_density <- function(x, y, method=c("pearson", "spearman", "cosine")) {
  method=match.arg(method)
  if (method=="pearson" || method == "spearman") {
    cor(as.vector(x$z), as.vector(y$z), method=method)
  } else if (method == "cosine") {
    cosine_sim(as.vector(x$z), as.vector(y$z))
  } else {
    stop()
  }

}


#' @importFrom tibble tibble
fixation_group <- function(x, y, duration, onset) {
  ret <- tibble(x=x,y=y, duration=duration, onset=onset)
  class(ret) <- c("fixation_group", class(ret))
  ret
}


#' @export
rep_fixations.fixation_group <- function(x, resolution=100) {
  nreps <- as.integer(x$duration/resolution)
  nreps[nreps < 1] <- 1
  x <- x[rep(1:nrow(x), nreps),]
  x
}


