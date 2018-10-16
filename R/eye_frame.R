
#' @export
#' @importFrom dplyr do group_by select filter
#' @import magrittr
#' @importFrom assertthat assert_that
eye_table <- function(x, y, duration, onset, groupvar, vars=NULL, data, clip_bounds=c(0,1280, 0,1280), relative_coords=TRUE) {

  data <- data %>% rename_(x=x, y=y, duration=duration, onset=onset)

  data <- if (is.null(vars)) {
    data %>% select_("x","y","duration", "onset", .dots=c(groupvar))
  } else {
    data %>% select_("x","y","duration", "onset", .dots=c(vars, groupvar))
  }


  xdir <- sign(clip_bounds[2] - clip_bounds[1])
  ydir <- sign(clip_bounds[4] - clip_bounds[3])
  xr <- sort(c(clip_bounds[1], clip_bounds[2]))
  yr <- sort(c(clip_bounds[3], clip_bounds[4]))

  data <- data %>% filter(x >= xr[1] & x <= xr[2]
                        & y >= yr[1] & y <= yr[2])


  if (relative_coords) {
    data <- data %>% mutate(x = (x - clip_bounds[1]) * xdir, y = (y - clip_bounds[3]) * ydir)
  }



  res <- data %>%
    group_by_(.dots=groupvar) %>%
    do({
        cbind(.[1,], tibble(fixgroup=list(fixation_group(.[["x"]], .[["y"]], .[["duration"]], .[["onset"]]))))
    }) %>% select_(.dots=c("fixgroup", vars))

  class(res) <- c("eye_table", class(res))
  res

}


rank_trans <- scales::trans_new(name="rank",
                                transform=function(x) { browser(); rank(x) },
                                inverse=function(x) (length(x)+1) - rank(x))

cuberoot_trans <- scales::trans_new(name="rank",
                                transform=function(x) { x^(1/3) },
                                inverse=function(x) x^3)

#' @import ggplot2
#' @importFrom ggplot2 ggplot aes annotation_raster geom_point
#' @importFrom imager load.image
#' @importFrom RColorBrewer brewer.pal
#' @export
plot.fixation_group <- function(x, type=c("contour", "density", "raster"), bandwidth=100,
                                xlim=range(x$x),
                                ylim=range(x$y),
                                size_points=TRUE,
                                show_points=TRUE,
                                bins=max(as.integer(length(x$x)/10),4),
                                bg_image=NULL,
                                colours=rev(RColorBrewer.brewer.pal(n=10, "Spectral")),
                                alpha_range=c(.5,1)) {
  type <- match.arg(type)

  if (size_points) {
    ps <- (x$duration - min(x$duration)) / (max(x$duration) - min(x$duration))
    x$psize <- ps*2 + 1
  }

  p <- ggplot(data=x, aes(x=x, y=y)) +
    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2])

  if (!is.null(bg_image)) {
    im <- imager::load.image(bg_image)
    p <- p + annotation_raster(as.raster(im),
                               xmin=xlim[1],
                               xmax=xlim[2],
                               ymin=ylim[1],
                               ymax=ylim[2])
  }

  if (show_points) {
    if (size_points) {
      p <- p + geom_point(aes(size=psize))
    } else {
      p <- p + geom_point()
    }
  }


  p <- if (type== "contour") {
    p + stat_density_2d(aes(colour = ..level..), h=bandwidth) +
      theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) + guides(size = "none")
  } else if (type == "density") {
    p + stat_density2d(aes(fill = ..level.., alpha=..level..), geom = "polygon", bins=bins, h=bandwidth)   +
      scale_fill_gradientn(colours=rev(brewer.pal(n=10, "Spectral")), guide=FALSE, trans=cuberoot_trans) +
      scale_alpha_continuous(range=alpha_range, trans=cuberoot_trans, guide=FALSE) +
      theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) + guides(size = "none")
  } else if (type == "raster") {
      p + stat_density_2d(aes(fill = ..density.., alpha=..density..), geom="raster", bins=bins, h=bandwidth, contour = FALSE, interpolate=TRUE) +
      scale_fill_gradientn(colours=rev(brewer.pal(n=10, "Spectral")), trans=cuberoot_trans, guide = FALSE)+
      scale_alpha_continuous(range=alpha_range, guide = FALSE, trans=cuberoot_trans) +
        theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank()) + guides(size = "none")
  } else {
    stop(paste("unrecognized 'type' ", type))

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
density_by <- function(x, groups, sigma=50, xbounds=c(0, 1000), ybounds=c(0, 1000), outdim=c(100,100),
                       duration_weighted=TRUE, ...) {
  ret <- x %>% group_by_(.dots=groups) %>% do( {
    g <- do.call(rbind, .$fixgroup)
    cbind(.[1,groups],tibble(fixgroup=list(g)))
  }) %>% rowwise() %>% do( {
    d <- eye_density(.$fixgroup, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim,
                     duration_weighted=duration_weighted,...)
    cbind(as_tibble(.[groups]), tibble( fixgroup=list(.$fixgroup), density=list(d)))
  })

  ret

}

#' @export
#' @importFrom quantreg rq
template_regression <- function(ref_tab, source_tab, match_on,
                                baseline_tab, baseline_key,
                                method=c("lm", "rlm", "rank")) {
  method <- match.arg(method)
  matchind <- match(source_tab[[match_on]], ref_tab[[match_on]])
  source_tab <- source_tab %>% ungroup() %>% mutate(matchind=matchind)

  if (any(is.na(matchind))) {
    warning("did not find matching template map for all source maps. Removing non-matching elements.")
    source_tab <- source_tab %>% filter(!is.na(matchind))
    matchind <- matchind[!is.na(matchind)]
  }


  ret <- source_tab %>% rowwise() %>% do( {
    id <- which(study_dens_subj_avg[[baseline_key]] == .[[baseline_key]][1])
    bdens <- study_dens_subj_avg$density[[id]]
    d1 <- ref_tab$density[[.$matchind]]
    d2 <- .$density

    df1 <- data.frame(y=as.vector(d2$z), baseline=as.vector(bdens$z), x2=as.vector(d1$z))

    est <- if (method == "lm") {
      res <- lm(y ~ baseline + x2, data=df1)
      coef(res)[2:3]
    } else if (method == "rlm") {
      res <- rlm(y ~ baseline + x2, data=df1, maxit=100)
      coef(res)[2:3]
    } else if (method == "rank") {
      #browser()
      #res <- rfit(y ~ baseline + x2, data=df1)
      #coef(res)[2:3]
      res <- pcor(df1, method="spearman")
      res$estimate[2:3,1]
    } else {
      stop()
    }

    data.frame(b0=est[1], b1=est[2])
  })


  source_tab %>% mutate(beta_baseline=ret$b0, beta_source=ret$b1)

}

#' template_similarity
#'
#' @param ref_tab
#' @param source_tab
#' @param match_on
#' @param method
#' @export
template_similarity <- function(ref_tab, source_tab, match_on,
                                method=c("spearman", "pearson", "cosine", "l1", "jaccard"),
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
                                       normalize=TRUE, duration_weighted=FALSE) {
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

  out$fixgroup <- x

  class(out) <- c("eye_density", "list")
  out

  #ret <- blurpoints(cbind(x$x, x$y), sigma=sigma, xbounds=xbounds, ybounds=ybounds, resolution=resolution, weights=wts, normalize=normalize)
}


# cosine_sim <- function(x,y){
#   dot.prod <- x%*%y
#   norm.x <- norm(x,type="2")
#   norm.y <- norm(y,type="2")
#   theta <- acos(dot.prod / (norm.x * norm.y))
#   as.numeric(theta)
# }


#' @importFrom proxy simil
#' @export
similarity.eye_density <- function(x, y, method=c("pearson", "spearman", "cosine", "l1", "jaccard")) {
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
  }

}


#' @importFrom tibble tibble
#' @export
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


