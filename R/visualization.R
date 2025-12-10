#' Animate a Fixation Scanpath with gganimate
#'
#' This function creates an animated visualization of a fixation scanpath using gganimate.
#'
#' @param x A `fixation_group` object.
#' @param bg_image An optional image file name to use as the background.
#' @param xlim The range in x coordinates (default: range of x values in the fixation group).
#' @param ylim The range in y coordinates (default: range of y values in the fixation group).
#' @param alpha The opacity of each dot (default: 1).
#' @param anim_over Animate over index (ordered) or onset (real time) (default: c("index", "onset")).
#' @param type Display as points or a raster (default: c("points", "raster")).
#' @param time_bin The size of the time bins (default: 1).
#'
#' @return A gganimate object representing the animated scanpath.
#' @import gganimate
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous geom_point scale_colour_gradientn theme_void labs theme scale_alpha_continuous
#' @importFrom ggplot2 annotation_raster stat_density_2d scale_fill_gradientn guides
#' @examples
#' # Create a fixation group
#' fg <- fixation_group(x=c(.1,.5,1), y=c(1,.5,1), onset=1:3, duration=rep(1,3))
#' # Animate the scanpath for the fixation group
#' anim_sp <- anim_scanpath(fg)
#' @export
#' @family visualization
anim_scanpath <- function(x, bg_image=NULL, xlim=range(x$x),
                          ylim=range(x$y), alpha=1,
                          anim_over=c("index", "onset"),
                          type=c("points", "raster"),
                          time_bin=1) {

  anim_over <- match.arg(anim_over)
  type <- match.arg(type)

  if (time_bin > 1) {
    x <- x %>% mutate(time_bin = round(onset/time_bin))
    anim_over = "time_bin"
  }

  p <- ggplot(data=x, aes(x=x, y=y)) +
    scale_x_continuous(expand=expansion(mult = c(.1, .1)), limits=c(xlim[1], xlim[2])) +
    scale_y_continuous(expand=expansion(mult = c(.1, .1)), limits=c(ylim[1], ylim[2]))

  if (!is.null(bg_image)) {
    im <- imager::load.image(bg_image)
    p <- p + annotation_raster(as.raster(im),
                               xmin=xlim[1],
                               xmax=xlim[2],
                               ymin=ylim[1],
                               ymax=ylim[2])
  }




  p <- if (type == "points") {
    p + geom_point(aes(colour=onset, size=15), alpha=alpha, show.legend=FALSE) +
      scale_colour_gradientn(colours=rev(brewer.pal(n=10, "Spectral"))) +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
      labs(title = 'Onset: {frame_time}') +
     gganimate::transition_time(.data[[anim_over]])
  } else {
    p + stat_density_2d(aes(fill = ..density.., alpha=..density..), geom="raster", bins=20,
                        h=100, contour = FALSE, interpolate=TRUE) +
      scale_fill_gradientn(colours=rev(brewer.pal(n=10, "Spectral")), guide = "none") +
      scale_alpha_continuous(range=c(.5,1), guide = "none") +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) + guides(size = "none") +
      labs(title = 'Onset: {frame_time}') +
      gganimate::transition_time(.data[[anim_over]])

  }


  p


}

#' Plot Eye Density
#'
#' This function creates a plot of the eye density using ggplot2.
#'
#' @param x An "eye_density" object.
#' @param alpha The opacity level for the raster layer (default: 0.8).
#' @param bg_image An optional image file name to use as the background.
#' @param transform The transformation to apply to the density values (default: c("identity", "sqroot", "curoot", "rank")).
#' @param ... Additional args
#' @return A ggplot object representing the eye density plot.
#' @import colorplane
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradientn theme annotation_raster
#' @importFrom ggplot2 element_blank
#' @family visualization
#' @examples
#' # Assume `ed` is an "eye_density" object
#' # Plot the eye density
#' plot_eye_density <- plot.eye_density(ed)
#' @export
plot.eye_density <- function(x, alpha=.8, bg_image=NULL,
                             transform=c("identity", "sqroot", "curoot", "rank"), ...) {
  transform <- match.arg(transform)
  trans <- resolve_transform(transform)

  xlim <- range(x$x)
  ylim <- range(x$y)

  dfx <- expand.grid(x=x$x, y=x$y)
  dfx$z <- as.vector(x$z)

  p <- ggplot(data=dfx, aes(x=x, y=y, fill=z))

  if (!is.null(bg_image)) {
    im <- imager::load.image(bg_image)
    p <- p + annotation_raster(as.raster(im),
                               xmin=xlim[1],
                               xmax=xlim[2],
                               ymin=ylim[1],
                               ymax=ylim[2])
  }

  p + geom_raster(alpha=alpha) +
    scale_fill_gradientn(colours=rev(brewer.pal(n=10, "Spectral")), trans=trans, guide = "none") +
    theme(axis.text        = element_blank(),
                          axis.ticks       = element_blank(),
                          axis.title       = element_blank(),
                          panel.background = element_blank())



}


# fix_to_graph <- function(fixgroup) {
#   nodes <- data.frame(name=paste0("N", 1:nrow(fixgroup)), onset=fixgroup$onset,
#                       duration=fixgroup$duration)
#   edges <- data.frame(from=1:(nrow(fixgroup)-1),
#                       to=2:nrow(fixgroup))
#
#   tbl_graph(nodes = nodes, edges = edges)
#
# }

#' Plot a fixation_group object
#'
#' This function creates a plot of the fixation_group object using ggplot2.
#' Different plot types and options can be specified to customize the output.
#'
#' @param x A fixation_group object.
#' @param type The type of plot to display (default: c("points", "contour", "filled_contour", "density", "raster")).
#' @param bandwidth The bandwidth for the kernel density estimator (default: 60).
#' @param xlim The x-axis limits (default: range of x values in the fixation_group object).
#' @param ylim The y-axis limits (default: range of y values in the fixation_group object).
#' @param size_points Whether to size points according to fixation duration (default: TRUE).
#' @param show_points Whether to show the fixations as points (default: TRUE).
#' @param show_path Whether to show the fixation path (default: TRUE).
#' @param bins Number of bins for density calculations (default: max(as.integer(length(x$x)/10), 4)).
#' @param bg_image An optional background image file name.
#' @param colours Color palette to use for the plot (default: rev(RColorBrewer.brewer.pal(n=10, "Spectral"))).
#' @param alpha_range A vector specifying the minimum and maximum alpha values for density plots (default: c(0.5, 1)).
#' @param alpha The opacity level for the points (default: 0.8).
#' @param window A vector specifying the time window for selecting fixations (default: NULL).
#' @param transform The transformation to apply to the density values (default: c("identity", "sqroot", "curoot", "rank")).
#' @param ... Additional arguments (currently unused).
#' @return A ggplot object representing the fixation group plot.
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes annotation_raster geom_point
#' @importFrom imager load.image
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr filter
#' @family visualization
#' @examples
#' # Create a fixation_group object
#' fg <- fixation_group(x=runif(50, 0, 100), y=runif(50, 0, 100), duration=rep(1,50), onset=seq(1,50))
#' # Plot the fixation group
#' plot_fixation_group(fg)
#' @export
plot.fixation_group <- function(x, type=c("points", "contour", "filled_contour", "density", "raster"),
                                bandwidth=60,
                                xlim=range(x$x),
                                ylim=range(x$y),
                                size_points=TRUE,
                                show_points=TRUE,
                                show_path=TRUE,
                                bins=max(as.integer(length(x$x)/10),4),
                                bg_image=NULL,
                                colours=rev(RColorBrewer.brewer.pal(n=10, "Spectral")),
                                alpha_range=c(.5,1),
                                alpha=.8,
                                window=NULL,
                                transform=c("identity", "sqroot", "curoot", "rank"), ...) {
  type <- match.arg(type)
  transform <- match.arg(transform)

  if (!is.null(window)) {
    assertthat::assert_that(length(window)==2)
    assertthat::assert_that(window[2] > window[1])
    x <- filter(x, onset >= window[1] & onset < window[2])
  }

  if (size_points) {
    d_range <- max(x$duration) - min(x$duration)
    if (d_range == 0) {
      x$psize <- 1
    } else {
      ps <- (x$duration - min(x$duration)) / d_range
      x$psize <- ps * 2 + 1
    }
  }

  p <- ggplot(data=x, aes(x=x, y=y)) +
    #xlim(xlim[1], xlim[2]) +
    #ylim(ylim[1], ylim[2]) +
    scale_x_continuous(expand=expansion(mult = c(.1, .1)), limits=c(xlim[1], xlim[2])) +
    scale_y_continuous(expand=expansion(mult = c(.1, .1)), limits=c(ylim[1], ylim[2]))

  if (!is.null(bg_image)) {
    im <- imager::load.image(bg_image)
    p <- p + annotation_raster(as.raster(im),
                               xmin=xlim[1],
                               xmax=xlim[2],
                               ymin=ylim[1],
                               ymax=ylim[2])
  }



  trans <- resolve_transform(transform)


  p <- if (type== "contour") {
    #dens <- as.data.frame.eye_density(eye_density(x, sigma=bandwidth))
    #p + geom_contour_filled(aes(x, y, z = density), dens, alpha=alpha) +
    #  guides(size = "none") +
    #  theme_void()

    p + stat_density_2d(aes(colour=..level..), h=bandwidth) +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
      guides(size = "none")
  } else if (type == "filled_contour") {
    ##dens <- as.data.frame.eye_density(eye_density(x, sigma=bandwidth))
    p + geom_density_2d_filled(alpha=alpha, h=bandwidth) +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
      #scale_alpha_continuous(range=alpha_range, trans=trans, guide=FALSE) +
      theme(legend.position = "none") +
      guides(size = "none")

  } else if (type == "density") {
    p + stat_density2d(aes(fill = ..level.., alpha=..level..), geom = "polygon", bins=bins, h=bandwidth)   +
      scale_fill_gradientn(colours=rev(brewer.pal(n=10, "Spectral")), guide="none", trans=cuberoot_trans) +
      scale_alpha_continuous(range=alpha_range, trans=trans, guide="none") +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
      guides(size = "none")
      #+
      #scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  } else if (type == "raster") {
    p + stat_density_2d(aes(fill = ..density.., alpha=..density..), geom="raster", bins=bins,
                        h=bandwidth, contour = FALSE, interpolate=TRUE) +
      scale_fill_gradientn(colours=rev(brewer.pal(n=10, "Spectral")), trans=cuberoot_trans, guide = "none") +
      scale_alpha_continuous(range=alpha_range, guide = "none", trans=trans) +
      theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
      guides(size = "none") #+
      #scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))
  } else if (type == "points") {
    if (show_path) {
      p <- p +  geom_path(aes(x,y, colour=onset), show.legend=FALSE)
    }

    p <- p + theme_void() + theme(panel.grid = element_blank(), panel.border = element_blank()) +
      guides(size = "none")

  } else {
    stop(paste("unrecognized 'type' ", type))

  }

  if (show_points) {
    if (size_points) {
      p <- p + geom_point(aes(size=psize, colour=onset), alpha=alpha, show.legend=FALSE) +
        scale_colour_gradient(low = "yellow", high = "red", na.value = NA)
      if (nrow(x) < 50) {
        p <- p + geom_text(aes(x,y, label=index))
      }
    } else {
      p <- p + geom_point(aes(colour=onset), alpha=alpha, show.legend=FALSE) +
        scale_colour_gradient(low = "yellow", high = "red", na.value = NA)
      if (nrow(x) < 50) {
        p <- p + geom_text(aes(x,y, label=index))
      }
    }
  }

  p
}

#' @noRd
rank_trans <- scales::trans_new(name="rank",
                                transform=function(x) { rank(x) },
                                inverse=function(x) (length(x)+1) - rank(x))

#' @noRd
cuberoot_trans <- scales::trans_new(name="curoot",
                                    transform=function(x) { x^(1/3) },
                                    inverse=function(x) x^3)

#' @noRd
squareroot_trans <- scales::trans_new(name="sqroot",
                                    transform=function(x) { x^(1/2) },
                                    inverse=function(x) x^2)


