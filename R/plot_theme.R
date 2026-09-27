# Shared visual system for eyesim plots -----------------------------------

#' eyesim plot colours
#'
#' Named colour roles used by every eyesim plot. Reference (encoding) and
#' source (recall) gaze keep the same colours in every panel; time, density
#' and correspondence each have one ramp or colour.
#'
#' @param role Optional character vector of role names. `NULL` returns all.
#' @return A named character vector of hex colours.
#' @examples
#' eyesim_colours()
#' eyesim_colours(c("reference", "source"))
#' @family visualization
#' @export
eyesim_colours <- function(role = NULL) {
  roles <- c(
    reference = "#1B2A41",
    source = "#C23B4B",
    raw = "#8A8F98",
    correspondence = "#5B4FCF",
    diagnostic = "#3F4A5A",
    ink = "#1F2328",
    muted = "#5A6270",
    rule = "#D5D9DE",
    grid = "#ECEEF1"
  )
  if (is.null(role)) {
    return(roles)
  }
  unknown <- setdiff(role, names(roles))
  if (length(unknown) > 0L) {
    stop("Unknown eyesim colour role(s): ", paste(unknown, collapse = ", "))
  }
  roles[role]
}

# Sequential ramp for time / fixation order: dark start, no near-white end.
eyesim_time_colours <- function(n = 7L) {
  viridisLite::viridis(n, begin = 0.05, end = 0.82)
}

# Sequential ramp for density: light to dark, fading to transparent at zero.
eyesim_density_colours <- function(n = 9L, transparent = TRUE) {
  cols <- viridisLite::rocket(n, begin = 0.12, end = 0.97, direction = -1)
  if (!transparent) {
    return(cols)
  }
  scales::alpha(cols, seq(0, 1, length.out = n)^0.5 * 0.9)
}

#' eyesim ggplot2 theme
#'
#' `theme_eyesim()` is the house theme shared by all eyesim plots: plain
#' background, light major grid, left-aligned title block, and a small grey
#' caption for provenance notes. `theme_eyesim_spatial()` is the variant used
#' for stimulus-space plots (scanpaths, density maps): no axes or grid, and a
#' thin frame marking the stimulus extent.
#'
#' Both return ordinary ggplot2 themes, so they can be added to any plot and
#' further modified with [ggplot2::theme()].
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @return A [ggplot2::theme()] object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_eyesim()
#' @family visualization
#' @export
theme_eyesim <- function(base_size = 10, base_family = "") {
  col <- eyesim_colours()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = col[["ink"]]),
      plot.title = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(1.15), hjust = 0,
        margin = ggplot2::margin(b = 3)
      ),
      plot.subtitle = ggplot2::element_text(
        colour = col[["muted"]], size = ggplot2::rel(0.88), hjust = 0,
        lineheight = 1.1, margin = ggplot2::margin(b = 6)
      ),
      plot.caption = ggplot2::element_text(
        colour = col[["muted"]], size = ggplot2::rel(0.75), hjust = 0,
        lineheight = 1.1, margin = ggplot2::margin(t = 6)
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      axis.text = ggplot2::element_text(
        colour = col[["muted"]], size = ggplot2::rel(0.82)
      ),
      axis.title = ggplot2::element_text(
        colour = col[["muted"]], size = ggplot2::rel(0.88)
      ),
      panel.grid.major = ggplot2::element_line(
        colour = col[["grid"]], linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.82)),
      legend.text = ggplot2::element_text(
        colour = col[["muted"]], size = ggplot2::rel(0.76)
      ),
      legend.key.size = ggplot2::unit(10, "pt"),
      strip.text = ggplot2::element_text(
        face = "bold", size = ggplot2::rel(0.88), hjust = 0
      ),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

#' @rdname theme_eyesim
#' @export
theme_eyesim_spatial <- function(base_size = 10, base_family = "") {
  theme_eyesim(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        fill = NA, colour = eyesim_colours("rule"), linewidth = 0.4
      )
    )
}

#' eyesim colour scales
#'
#' `scale_fill_eyesim_density()` is the sequential density scale used by the
#' density plots. Opacity is part of the ramp and is anchored at zero
#' density, so zero is transparent (a background stimulus shows through) and
#' the colour bar shows exactly what is drawn. Set shared `limits` to compare
#' conditions on one scale. `scale_colour_eyesim_time()` is the
#' fixation-order/onset scale.
#'
#' @param ... Passed to [ggplot2::scale_fill_gradientn()] or
#'   [ggplot2::scale_colour_gradientn()].
#' @param name Legend title.
#' @param limits Scale limits. The default starts the density scale at zero.
#' @param transparent Whether low densities fade to transparent.
#' @param aesthetics Aesthetics the scale applies to.
#' @return A ggplot2 scale.
#' @examples
#' fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
#'                      onset = c(0, 200, 400), duration = c(200, 200, 200))
#' ed <- eye_density(fg, sigma = 50, xbounds = c(0, 400), ybounds = c(0, 300))
#' # Put several maps on one scale by giving them the same limits
#' plot(ed, limits = c(0, 1e-4))
#' @family visualization
#' @export
scale_fill_eyesim_density <- function(..., name = "Fixation\ndensity",
                                      limits = c(0, NA), transparent = TRUE,
                                      aesthetics = "fill") {
  ggplot2::scale_fill_gradientn(
    colours = eyesim_density_colours(transparent = transparent),
    name = name, limits = limits, aesthetics = aesthetics, ...
  )
}

#' @rdname scale_fill_eyesim_density
#' @export
scale_colour_eyesim_time <- function(..., aesthetics = "colour") {
  ggplot2::scale_colour_gradientn(
    colours = eyesim_time_colours(), aesthetics = aesthetics, ...
  )
}

# Density scale used inside the plot methods: user or eyesim colours, opacity
# ramp from alpha_range[1] (zero density) to alpha_range[2] (maximum), zero
# anchored unless the transform is a rank (whose axis has no zero).
density_scale <- function(colours = NULL, transform = "identity",
                          alpha_range = c(0, 0.9), limits = NULL,
                          aesthetics = "fill",
                          title = density_legend_title(transform),
                          band_edges = NULL, overflow = FALSE) {
  base <- if (is.null(colours)) eyesim_density_colours(transparent = FALSE) else colours
  n <- max(length(base), 9L)
  ramp <- grDevices::colorRampPalette(base)(n)
  rank <- identical(transform, "rank")
  if (rank && !is.null(limits)) {
    warning("`limits` is ignored for transform = \"rank\".", call. = FALSE)
    limits <- NULL
  }
  if (is.null(limits) && !rank) {
    limits <- c(0, NA)
  }
  if (!is.null(limits)) {
    limits[!is.finite(limits)] <- NA
  }
  trans <- resolve_transform(transform)
  position <- seq(0, 1, length.out = n)
  args <- list(
    transform = trans,
    limits = limits,
    # values above shared limits take the top colour instead of NA grey
    oob = scales::oob_squish_any,
    na.value = "transparent",
    labels = format_density,
    aesthetics = aesthetics,
    guide = ggplot2::guide_colourbar(order = 1)
  )
  if (!is.null(band_edges)) {
    # Banded densities: the colour is transparent up to the first edge (the
    # band containing zero) and ramps above it; a stepped colour bar shows
    # the band edges and colours exactly as drawn.
    shown_edges <- band_edges
    if (overflow) {
      # Density above the upper limit is drawn in the top colour; the bar
      # shows it as an extra step labelled "> limit".
      shown_edges <- c(band_edges, band_edges[[length(band_edges)]] +
                         diff(band_edges)[[1]])
    }
    tr <- scales::as.transform(trans)
    te <- tr$transform(shown_edges)
    cut <- (te[[2]] - te[[1]]) / (te[[length(te)]] - te[[1]])
    position <- c(0, cut, cut + (1 - cut) * seq(1e-6, 1, length.out = n))
    ramp <- c(ramp[[1]], ramp[[1]], ramp)
    opacity <- c(0, 0, alpha_range[[1]] + diff(alpha_range) *
                   (cut + (1 - cut) * seq(0, 1, length.out = n))^0.5)
    args$limits <- range(shown_edges)
    args$breaks <- shown_edges
    # label at most about five edges
    step <- ceiling((length(band_edges) - 1L) / 4)
    top_edge <- band_edges[[length(band_edges)]]
    args$labels <- function(x) {
      out <- format_density(x)
      idx <- seq_along(x) - 1L
      last <- length(band_edges) - 1L
      # label every `step`-th edge and the top edge, skipping a regular
      # label that would crowd the top one
      keep <- (idx %% step == 0L & (last - idx >= step | idx == last)) | idx == last
      out[!keep] <- ""
      if (overflow) {
        # one label for the limit: "> limit" at the end of the extra step
        out[idx == last] <- ""
        out[seq_along(x) == length(x)] <- paste0("> ", format_density(top_edge))
      }
      out
    }
    args$guide <- ggplot2::guide_coloursteps(order = 1, show.limits = TRUE)
  } else {
    opacity <- alpha_range[[1]] + diff(alpha_range) * position^0.5
  }
  args$colours <- scales::alpha(ramp, opacity)
  args$values <- position
  if (rank) {
    # Rank values have no meaningful tick labels; show the ramp only.
    args["breaks"] <- list(NULL)
  }
  # The legend title is set with labs() so that users can override it.
  list(
    do.call(ggplot2::scale_fill_gradientn, args),
    do.call(ggplot2::labs, stats::setNames(list(title), aesthetics[[1]]))
  )
}

# One notation per scale: scientific when the values are small.
format_density <- function(x) {
  finite <- x[is.finite(x) & x != 0]
  small <- length(finite) > 0L && max(abs(finite)) < 1e-2
  out <- if (small) {
    formatC(x, format = "e", digits = 1)
  } else {
    formatC(x, format = "g", digits = 2)
  }
  out[!is.na(x) & x == 0] <- "0"
  out
}

density_legend_title <- function(transform, base = "Fixation\ndensity") {
  switch(transform,
    identity = base,
    rank = paste0(base, "\n(rank, low to high)"),
    sqroot = paste0(base, "\n(square root)"),
    curoot = paste0(base, "\n(cube root)"),
    base
  )
}

# Wrap label text to `width` characters without splitting protected phrases
# (phrases other code and tests rely on staying intact).
eyesim_wrap <- function(text, width = 105,
                        protect = c("not a posterior", "not posterior",
                                    "posterior replay", "equal-prior episode",
                                    "gaze_info_bits")) {
  if (is.null(text) || !nzchar(text)) {
    return(text)
  }
  glue <- "\u00a0"
  for (phrase in protect) {
    text <- gsub(phrase, gsub(" ", glue, phrase, fixed = TRUE), text, fixed = TRUE)
  }
  text <- gsub(" = ", paste0(glue, "=", glue), text, fixed = TRUE)
  lines <- strwrap(text, width = width)
  gsub(glue, " ", paste(lines, collapse = "\n"), fixed = TRUE)
}
