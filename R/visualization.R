

#' @import colorplane
#' @export
plot.eye_density <- function(x,show_points=TRUE) {
  stop("not implemented")
  xlim <- range(x$x)
  ylim <- range(x$y)

  dfx <- data.frame(x=xlim, y=ylim)

  clrs <- map_colors(IntensityColorPlane(as.numeric(x$z)))
  z <- matrix(clrs@clrs, nrow(x$z), ncol(x$z))

  p <- ggplot(data=dfx, aes(x=x, y=y))
  p <- p + annotation_raster(as.raster(t(z)),
                               xmin=xlim[1],
                               xmax=xlim[2],
                               ymin=ylim[1],
                               ymax=ylim[2])

}



#' @import ggplot2
#' @importFrom ggplot2 ggplot aes annotation_raster geom_point
#' @importFrom imager load.image
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
#' @examples
#'
#' fg <- fixation_group(x=runif(50, 0, 100), y=runif(50, 0, 100), duration=rep(1,50), onset=seq(1,50))
plot.fixation_group <- function(x, type=c("contour", "density", "raster"), bandwidth=100,
                                xlim=range(x$x),
                                ylim=range(x$y),
                                size_points=TRUE,
                                show_points=TRUE,
                                bins=max(as.integer(length(x$x)/10),4),
                                bg_image=NULL,
                                colours=rev(RColorBrewer.brewer.pal(n=10, "Spectral")),
                                alpha_range=c(.5,1),
                                transform=c("identity", "sqroot", "curoot", "rank")) {
  type <- match.arg(type)
  transform <- match.arg(transform)

  if (size_points) {
    ps <- (x$duration - min(x$duration)) / (max(x$duration) - min(x$duration))
    x$psize <- ps*2 + 1
  }

  p <- ggplot(data=x, aes(x=x, y=y)) +
    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))

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

  trans <- if (transform == "identity") {
    "identity"
  } else if (transform == "sqroot") {
    squareroot_trans
  } else if (transform == "curoot") {
    cuberoot_trans
  } else if (transform == "rank") {
    rank_trans
  } else {
    "identity"
  }


  p <- if (type== "contour") {
    #dens <- as.data.frame.eye_density(eye_density(x, sigma=bandwidth))
    #p + geom_contour_filled(aes(x, y, z = density), dens, alpha=alpha) +
    #  guides(size = "none") +
    #  theme_void()

    p + stat_density_2d(aes(fill=stat(level)), h=bandwidth, geom="polygon", alpha=alpha) +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
      guides(size = "none")
  } else if (type == "density") {
    p + stat_density2d(aes(fill = ..level.., alpha=..level..), geom = "polygon", bins=bins, h=bandwidth)   +
      scale_fill_gradientn(colours=rev(brewer.pal(n=10, "Spectral")), guide=FALSE, trans=cuberoot_trans) +
      scale_alpha_continuous(range=alpha_range, trans=trans, guide=FALSE) +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) + guides(size = "none") +
      scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  } else if (type == "raster") {
    p + stat_density_2d(aes(fill = ..density.., alpha=..density..), geom="raster", bins=bins,
                        h=bandwidth, contour = FALSE, interpolate=TRUE) +
      scale_fill_gradientn(colours=rev(brewer.pal(n=10, "Spectral")), trans=cuberoot_trans, guide = FALSE) +
      scale_alpha_continuous(range=alpha_range, guide = FALSE, trans=trans) +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) + guides(size = "none") +
      scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  } else {
    stop(paste("unrecognized 'type' ", type))

  }

  p
}


rank_trans <- scales::trans_new(name="rank",
                                transform=function(x) { rank(x) },
                                inverse=function(x) (length(x)+1) - rank(x))

cuberoot_trans <- scales::trans_new(name="curoot",
                                    transform=function(x) { x^(1/3) },
                                    inverse=function(x) x^3)

squareroot_trans <- scales::trans_new(name="sqroot",
                                    transform=function(x) { x^(1/2) },
                                    inverse=function(x) x^2)


