


#' @importFrom ggplot2 ggplot
plot.fixation_group <- function(x, type=c("contour", "density", "density_alpha"), bandwidth=100, xlim=range(x$x), ylim=range(x$y), size_points=TRUE) {
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

density_by <- function(x, groups, sigma=50, xbounds=c(min(x$x), max(x$x)), ybounds=c(min(x$y), max(x$y)), outdim=c(100,100), ...) {
  ret <- x %>% group_by_(.dots=groups) %>% do( {
    g <- do.call(rbind, .$fixgroup)
    cbind(.[1,groups],tibble(fixgroup=list(g)))
  }) %>% rowwise() %>% do( {

    d <- eye_density(.$fixgroup, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim, ...)
    cbind(as_tibble(.[groups]), tibble( fixgroup=list(.$fixgroup), density=list(d)))
  })

  ret

  #cbind(., tibble(density=list(eye_density(.$fixgroup, sigma, ...))))

}



eye_density.fixation_group <- function(x, sigma=50, xbounds=c(min(x$x), max(x$x)), ybounds=c(min(x$y), max(x$y)),
                                       outdim=c(xbounds[2] - xbounds[1], ybounds[2] - ybounds[1]), duration_weighted=FALSE, normalize=TRUE) {
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

rep_fixations.fixation_group <- function(x, resolution=100) {
  nreps <- as.integer(x$duration/resolution)
  nreps[nreps < 1] <- 1
  x <- x[rep(1:nrow(x), nreps),]
  x
}


#' @export
#' @importFrom dplyr do group_by select filter
create_eye_table <- function(x="CURRENT_FIX_X", y="CURRENT_FIX_Y", duration="CURRENT_FIX_DURATION", onset="CURRENT_FIX_START", groupvar, vars, data) {

  if (missing(vars)) {
    other <- names(data)
    vars <- other[!(other %in% c(x,y,duration, onset, groupvar))]
  }

  data %>% select_(c(x,y,duration, onset, vars, groupvar))

  res <- data %>%
    group_by_(.dots=groupvar) %>%
    do(cbind(.[1,], tibble(fixgroup=list(fixation_group(.[[x]], .[[y]], .[[duration]], .[[onset]]))))) %>%
    select_(.dots=c("fixgroup", vars))

  class(res) <- c("eye_table", class(res))
  res

}
