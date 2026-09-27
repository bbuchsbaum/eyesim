# Background stimulus layer shared by the stimulus-space plots. The image is
# drawn with its stored colours (no contrast stretching) and stretched to
# the xlim/ylim extent, top row at ylim[2].
eyesim_bg_layer <- function(bg_image, xlim, ylim) {
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop("Package 'imager' is required for bg_image support. Install it with install.packages('imager').")
  }
  im <- if (is.character(bg_image)) imager::load.image(bg_image) else bg_image
  ras <- if (inherits(im, "cimg")) {
    grDevices::as.raster(im, rescale = max(im, na.rm = TRUE) > 1)
  } else {
    grDevices::as.raster(im)
  }
  annotation_raster(ras, xmin = xlim[1], xmax = xlim[2],
                    ymin = ylim[1], ymax = ylim[2])
}

# Fixed-aspect stimulus coordinates. Limits are applied by the coordinate
# system, so no data are dropped and densities are not truncated.
eyesim_spatial_coord <- function(xlim, ylim, bg_image = NULL) {
  coord_fixed(ratio = 1, xlim = xlim, ylim = ylim,
              expand = is.null(bg_image), clip = "on")
}


# Peak of the kernel density on the grid ggplot2's stat_density_2d uses.
kde_peak <- function(x, y, h, xlim, ylim, n = 100L) {
  xr <- range(c(x, xlim))
  yr <- range(c(y, ylim))
  max(MASS::kde2d(x, y, h = h, n = n, lims = c(xr, yr))$z)
}

# Thin white contour lines that keep density edges visible over a stimulus.
density_halo <- function(h, breaks = NULL, bins = NULL) {
  args <- list(h = h, colour = "white", linewidth = 0.3, alpha = 0.8)
  if (!is.null(breaks)) args$breaks <- breaks else args$bins <- bins
  do.call(ggplot2::geom_density_2d, args)
}

#' Animate a Fixation Scanpath with gganimate
#'
#' Creates an animated scanpath: fixations appear one at a time (by order or
#' by onset), earlier fixations stay visible as faded marks, and colour
#' encodes onset on the shared eyesim time scale.
#'
#' @param x A `fixation_group` object.
#' @param bg_image An optional image file name (or `cimg`) to use as the background.
#' @param xlim The range in x coordinates (default: range of x values in the fixation group).
#' @param ylim The range in y coordinates (default: range of y values in the fixation group).
#' @param alpha The opacity of each dot (default: 1).
#' @param anim_over Animate over index (ordered) or onset (real time) (default: c("index", "onset")).
#' @param type Display as points or a raster (default: c("points", "raster")).
#' @param time_bin The size of the time bins (default: 1).
#'
#' @return A gganimate object representing the animated scanpath.
#' @importFrom ggplot2 ggplot aes geom_point geom_path labs theme annotation_raster
#' @importFrom ggplot2 stat_density_2d coord_fixed
#' @importFrom grDevices as.raster
#' @examples
#' # Create a fixation group
#' fg <- fixation_group(x=c(.1,.5,1), y=c(1,.5,1), onset=1:3, duration=rep(1,3))
#' # Animate the scanpath for the fixation group
#' if (requireNamespace("gganimate", quietly = TRUE)) {
#'   anim_sp <- anim_scanpath(fg)
#' }
#' @export
#' @family visualization
anim_scanpath <- function(x, bg_image=NULL, xlim=range(x$x),
                          ylim=range(x$y), alpha=1,
                          anim_over=c("index", "onset"),
                          type=c("points", "raster"),
                          time_bin=1) {

  if (!requireNamespace("gganimate", quietly = TRUE)) {
    stop("Package 'gganimate' is required for anim_scanpath(). Install it with install.packages('gganimate').")
  }

  anim_over <- match.arg(anim_over)
  type <- match.arg(type)

  if (time_bin > 1) {
    x <- x %>% mutate(time_bin = round(onset/time_bin))
    anim_over = "time_bin"
  }
  frame_title <- switch(anim_over,
    index = "Fixation {round(frame_time)}",
    onset = "Onset {round(frame_time)}",
    time_bin = "Time bin {round(frame_time)}"
  )

  p <- ggplot(data=x, aes(x=x, y=y))
  if (!is.null(bg_image)) {
    p <- p + eyesim_bg_layer(bg_image, xlim, ylim)
  }

  p <- if (type == "points") {
    p + geom_point(aes(fill = onset), shape = 21, size = 4.5, stroke = 0.5,
                   colour = "white", alpha = alpha, show.legend = FALSE) +
      scale_colour_eyesim_time(aesthetics = "fill") +
      gganimate::shadow_mark(alpha = 0.35, size = 2.2)
  } else {
    p + stat_density_2d(aes(fill = after_stat(density), alpha = after_stat(density)),
                        geom = "raster", h = 100, contour = FALSE, interpolate = TRUE) +
      scale_fill_eyesim_density(transparent = FALSE, guide = "none") +
      scale_alpha_continuous(range = c(0, 0.9), guide = "none")
  }

  p + eyesim_spatial_coord(xlim, ylim, bg_image) +
    theme_eyesim_spatial() +
    labs(title = frame_title) +
    gganimate::transition_time(.data[[anim_over]])
}

# Legend wording for an eye_density map: a normalized map holds the
# proportion of fixation density in each grid cell.
eye_density_quantity <- function(x) {
  total <- sum(x$z, na.rm = TRUE)
  if (isTRUE(abs(total - 1) < 1e-6)) "Fixation\nproportion\nper cell" else "Fixation\ndensity"
}

#' Plot Eye Density
#'
#' Draws a fixation density map on the shared eyesim density scale. Zero
#' density is transparent, so a background stimulus shows through, and the
#' map keeps the aspect ratio of its coordinate grid. Over a `bg_image`,
#' thin white contour lines at five evenly spaced levels (from zero to the
#' upper limit or the map's maximum) keep density edges visible.
#'
#' @param x An "eye_density" object.
#' @param alpha Maximum opacity of the density layer (default: 0.8); lower
#'   densities are progressively more transparent.
#' @param bg_image An optional image file name (or `cimg`) to use as the background.
#' @param transform The transformation to apply to the density values (default: c("identity", "sqroot", "curoot", "rank")).
#' @param colours Optional vector of colours for the density ramp
#'   (default: the eyesim density ramp).
#' @param legend Whether to show the density colour bar (default: `TRUE`).
#' @param limits Optional density limits for the colour scale (default: from
#'   zero to the map's maximum). Give several maps the same `limits` to
#'   compare them on one scale.
#' @param ... Additional args
#' @return A ggplot object representing the eye density plot.
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradientn theme annotation_raster
#' @importFrom ggplot2 element_blank
#' @importFrom grDevices as.raster
#' @family visualization
#' @examples
#' # Create a fixation group and compute eye density
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' ed <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300))
#' # Plot the eye density
#' plot(ed)
#' @export
plot.eye_density <- function(x, alpha=.8, bg_image=NULL,
                             transform=c("identity", "sqroot", "curoot", "rank"),
                             colours = NULL, legend = TRUE, limits = NULL, ...) {
  transform <- match.arg(transform)

  xlim <- range(x$x)
  ylim <- range(x$y)

  dfx <- expand.grid(x=x$x, y=x$y)
  dfx$z <- as.vector(x$z)

  p <- ggplot(data=dfx, aes(x=x, y=y, fill=z))

  if (!is.null(bg_image)) {
    p <- p + eyesim_bg_layer(bg_image, xlim, ylim)
  }

  # The frame is the grid extent (the stimulus extent), so the image fills
  # it; the outer half of the edge raster cells is clipped.
  p <- p + geom_raster() +
    density_scale(colours, transform, alpha_range = c(0, alpha), limits = limits,
                  title = density_legend_title(transform, eye_density_quantity(x))) +
    coord_fixed(ratio = 1, expand = FALSE, xlim = xlim, ylim = ylim) +
    theme_eyesim_spatial()
  if (!is.null(bg_image)) {
    # Thin white contours keep density edges visible over the stimulus.
    top <- if (!is.null(limits) && is.finite(limits[[2]])) limits[[2]] else max(dfx$z, na.rm = TRUE)
    p <- p + geom_contour(data = dfx, mapping = aes(x = x, y = y, z = z),
                          inherit.aes = FALSE, breaks = seq(0, top, length.out = 6L)[-1],
                          colour = "white", linewidth = 0.3, alpha = 0.8)
  }
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  p
}


#' Plot a fixation_group object
#'
#' Plots a fixation group in stimulus coordinates with the shared eyesim
#' style. `type = "points"` draws the scanpath: fixations coloured by onset
#' and sized by duration, joined in order. The density types draw a kernel
#' density of the fixations on the eyesim density scale, with the fixations
#' overlaid in grey. The plot keeps a 1:1 aspect ratio; `xlim`/`ylim` set the
#' visible extent (and the extent of `bg_image`) without dropping data.
#' A ring marks the first fixation.
#'
#' With fewer than 50 fixations, fixation numbers are placed when the plot
#' is drawn, using its actual size: next to their own fixation, or further
#' out with a leader line, never covering another number or marker. Numbers
#' that cannot be placed are omitted, and a note in the panel says how many.
#'
#' `type = "density"` and `type = "filled_contour"` draw non-overlapping
#' density bands between evenly spaced edges from zero, each coloured at its
#' midpoint; the band containing zero is transparent, and the stepped colour
#' bar shows the same edges. Over a `bg_image`, thin white contour lines
#' separate the density from the stimulus colours.
#'
#' @param x A fixation_group object.
#' @param type The type of plot to display (default: c("points", "contour", "filled_contour", "density", "raster")).
#' @param bandwidth The bandwidth for the kernel density estimator (default: 60).
#' @param xlim The x-axis limits (default: range of x values in the fixation_group object).
#' @param ylim The y-axis limits (default: range of y values in the fixation_group object).
#' @param size_points Whether to size points according to fixation duration;
#'   point area is proportional to duration (default: TRUE).
#' @param show_points Whether to show the fixations as points (default: TRUE).
#' @param show_path Whether to show the fixation path (default: TRUE).
#' @param bins Number of density bands for `type = "density"`
#'   (default: max(as.integer(length(x$x)/10), 4)); `type = "filled_contour"`
#'   uses 10 bands.
#' @param bg_image An optional background image file name (or `cimg`).
#' @param colours Optional colour ramp. For `type = "points"` it replaces the
#'   onset (time) ramp; for the density types it replaces the density ramp.
#'   Default `NULL` uses the eyesim ramps.
#' @param alpha_range Opacity at zero and at maximum density for the density
#'   layers (default: c(0, 0.9), so zero density is transparent).
#' @param alpha The opacity level for the points (default: 0.8).
#' @param window A vector specifying the time window for selecting fixations (default: NULL).
#' @param transform The transformation applied to the density colour scale
#'   (default: c("identity", "sqroot", "curoot", "rank")).
#' @param legend Whether to show legends: onset (points) or density (density
#'   types), plus duration when `size_points = TRUE` (default: `TRUE`).
#' @param limits Optional density limits for the colour scale of the density
#'   types (default: from zero to the maximum). Use the same `limits` across
#'   plots to compare them on one scale; for the banded types they also fix
#'   the band edges. Densities above the upper limit take the top colour.
#' @param aspect `"equal"` (default) keeps one data unit the same length on
#'   both axes, so saccade angles and cluster shapes are true; wide scanpaths
#'   then leave blank space. `"free"` fills the panel and distorts the
#'   geometry.
#' @param ... Additional arguments (currently unused).
#' @return A ggplot object representing the fixation group plot.
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes annotation_raster geom_point
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices as.raster
#' @importFrom dplyr filter
#' @family visualization
#' @examples
#' # Create a fixation_group object
#' fg <- fixation_group(x=runif(50, 0, 100), y=runif(50, 0, 100), duration=rep(1,50), onset=seq(1,50))
#' # Plot the fixation group using the S3 method
#' plot(fg)
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
                                colours=NULL,
                                alpha_range=c(0, .9),
                                alpha=.8,
                                window=NULL,
                                transform=c("identity", "sqroot", "curoot", "rank"),
                                legend=TRUE, limits=NULL,
                                aspect=c("equal", "free"), ...) {
  type <- match.arg(type)
  transform <- match.arg(transform)
  aspect <- match.arg(aspect)
  density_type <- type != "points"
  col <- eyesim_colours()

  if (!is.null(window)) {
    assertthat::assert_that(length(window)==2)
    assertthat::assert_that(window[2] > window[1])
    x <- filter(x, onset >= window[1] & onset < window[2])
  }

  # Without a stimulus, pad default limits so density contours can close.
  if (density_type && is.null(bg_image)) {
    if (missing(xlim)) xlim <- xlim + c(-1, 1) * bandwidth
    if (missing(ylim)) ylim <- ylim + c(-1, 1) * bandwidth
  }

  p <- ggplot(data=x, aes(x=x, y=y))
  if (!is.null(bg_image)) {
    p <- p + eyesim_bg_layer(bg_image, xlim, ylim)
  }

  h <- rep(bandwidth, 2)
  if (density_type) {
    # Evaluate the density over the whole visible extent (plus margin).
    p <- p + expand_limits(x = xlim + c(-1, 1) * bandwidth,
                           y = ylim + c(-1, 1) * bandwidth)
  }

  if (type == "contour") {
    p <- p + stat_density_2d(aes(colour = after_stat(level)), h = h,
                             linewidth = 0.45) +
      density_scale(colours, transform, alpha_range = c(0.5, 1),
                    limits = limits, aesthetics = "colour")
  } else if (type %in% c("filled_contour", "density")) {
    # Non-overlapping density bands, each coloured at its band midpoint on
    # the continuous density scale, so the colour bar matches the pixels.
    # Shared `limits` fix the band edges, so bands agree across plots.
    n_bands <- if (type == "density") bins else 10L
    if (transform == "rank" && !is.null(limits)) {
      warning("`limits` is ignored for transform = \"rank\".", call. = FALSE)
      limits <- NULL
    }
    # The lowest band contains zero density and is left transparent.
    # The open top band (above a limit) is coloured at the midpoint of the
    # colour bar's "> limit" step, so bar and drawing agree.
    band_step <- NA_real_
    band_args <- list(
      mapping = aes(fill = after_stat(ifelse(
        level_low <= 0, NA,
        ifelse(is.finite(level_high), level_mid, level_low + band_step / 2)
      ))),
      h = h
    )
    edges <- NULL
    overflow <- FALSE
    if (!is.null(limits) && is.finite(limits[[1]]) && limits[[1]] != 0) {
      warning("Density bands start at zero; the lower limit is ignored for type = \"",
              type, "\".", call. = FALSE)
    }
    if (transform == "rank") {
      band_args$bins <- n_bands
    } else {
      # Evenly spaced edges from zero to the upper limit (or the density
      # peak), plus an open top band so density above a limit keeps the
      # top colour. Known edges give a stepped colour bar that matches
      # the bands, with the transparent zero band as its first step.
      peak <- kde_peak(x$x, x$y, h, xlim + c(-1, 1) * bandwidth, ylim + c(-1, 1) * bandwidth)
      top <- if (!is.null(limits) && is.finite(limits[[2]])) limits[[2]] else peak
      overflow <- peak > top
      edges <- seq(0, top, length.out = n_bands + 1L)
      band_step <- edges[[2]]
      band_args$breaks <- c(edges, Inf)
    }
    # Ranks have no zero: keep the lowest real band visible (the zero band
    # itself is NA and stays transparent).
    band_alpha <- if (transform == "rank") {
      c(max(alpha_range[[1]], 0.3), alpha_range[[2]])
    } else {
      alpha_range
    }
    p <- p + do.call(geom_density_2d_filled, band_args) +
      density_scale(colours, transform, alpha_range = band_alpha,
                    band_edges = edges, overflow = overflow)
    if (!is.null(bg_image) && !is.null(edges)) {
      p <- p + density_halo(h = h, breaks = edges[-1])
    }
  } else if (type == "raster") {
    p <- p + stat_density_2d(aes(fill = after_stat(density)),
                             geom = "raster", h = h, contour = FALSE, interpolate = TRUE) +
      density_scale(colours, transform, alpha_range = alpha_range, limits = limits)
    if (!is.null(bg_image)) {
      # Halo levels follow the colour scale, so shared limits share them.
      top <- if (!is.null(limits) && is.finite(limits[[2]])) {
        limits[[2]]
      } else {
        kde_peak(x$x, x$y, h, xlim + c(-1, 1) * bandwidth, ylim + c(-1, 1) * bandwidth)
      }
      p <- p + density_halo(h = h, breaks = seq(0, top, length.out = 6L)[-1])
    }
  } else if (show_path) {
    p <- p + geom_path(aes(colour = onset), linewidth = 0.5, alpha = 0.75,
                       show.legend = FALSE)
  }

  if (show_points) {
    if (density_type) {
      point_aes <- if (size_points) aes(size = duration) else aes()
      p <- p + geom_point(point_aes, shape = 21, fill = col[["ink"]], colour = "white",
                          stroke = 0.3, alpha = alpha * 0.8)
    } else {
      point_aes <- if (size_points) aes(size = duration, fill = onset) else aes(fill = onset)
      p <- p + geom_point(point_aes, shape = 21, colour = "white", stroke = 0.4,
                          alpha = alpha)
    }
    if (size_points) {
      # Point area is proportional to duration.
      p <- p + scale_size_area(
        max_size = 5,
        guide = guide_legend(
          order = 2, override.aes = list(fill = col[["muted"]], colour = "white")
        )
      ) + labs(size = "Duration")
    }
    # Ring marking the start of the scanpath, independent of its label.
    if (nrow(x) > 0L) {
      start <- x[which.min(x$onset), , drop = FALSE]
      p <- p + geom_point(data = start, aes(x = x, y = y), inherit.aes = FALSE,
                          shape = 21, size = 7, stroke = 1.8, fill = NA,
                          colour = "white") +
        geom_point(data = start, aes(x = x, y = y), inherit.aes = FALSE,
                   shape = 21, size = 7, stroke = 0.7, fill = NA,
                   colour = col[["ink"]])
    }
    if (nrow(x) < 50) {
      radius <- if (size_points) 2.5 * sqrt(x$duration / max(x$duration)) else 1.5
      lab_df <- data.frame(x = x$x, y = x$y, label = x$index, radius = radius + 0.3,
                           start = seq_len(nrow(x)) == which.min(x$onset))
      p <- p + geom_fixation_label(data = lab_df,
                                   aes(x = x, y = y, label = label, radius = radius,
                                       start = start),
                                   inherit.aes = FALSE, colour = col[["ink"]],
                                   start_r = 3.4)
    }
  }

  if (!density_type) {
    time_cols <- if (is.null(colours)) eyesim_time_colours() else colours
    p <- p + scale_colour_gradientn(colours = time_cols, aesthetics = c("colour", "fill"),
                                    guide = guide_colourbar(order = 1)) +
      labs(colour = "Onset", fill = "Onset")
  }

  p <- p + theme_eyesim_spatial()
  if (aspect == "equal") {
    p <- p + eyesim_spatial_coord(xlim, ylim, bg_image)
  } else {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim, expand = is.null(bg_image))
  }
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  p
}

#' @noRd
#' @importFrom scales trans_new
rank_trans <- scales::trans_new(name="rank",
                                transform=function(x) { rank(x, na.last = "keep") },
                                inverse=function(x) (length(x)+1) - rank(x))

#' @noRd
cuberoot_trans <- scales::trans_new(name="curoot",
                                    transform=function(x) { x^(1/3) },
                                    inverse=function(x) x^3)

#' @noRd
squareroot_trans <- scales::trans_new(name="sqroot",
                                    transform=function(x) { x^(1/2) },
                                    inverse=function(x) x^2)
