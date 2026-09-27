# GazeWeave visual explanation ---------------------------------------------

gaze_warp_grid <- function(alignment, grid_n = 6L, line_n = 40L) {
  screen <- alignment$spec$screen
  if (is.null(screen)) {
    return(NULL)
  }
  parameters <- warp_parameters(alignment$warp)
  transform_points <- function(points) {
    sweep(points %*% t(parameters$A), 2, parameters$translation, FUN = "+")
  }
  x_values <- seq(screen$xlim[[1]], screen$xlim[[2]], length.out = grid_n)
  y_values <- seq(screen$ylim[[1]], screen$ylim[[2]], length.out = grid_n)
  vertical <- lapply(seq_along(x_values), function(i) {
    raw <- cbind(
      x = rep(x_values[[i]], line_n),
      y = seq(screen$ylim[[1]], screen$ylim[[2]], length.out = line_n)
    )
    registered <- transform_points(raw)
    data.frame(x = registered[, 1], y = registered[, 2], group = paste0("v", i))
  })
  horizontal <- lapply(seq_along(y_values), function(i) {
    raw <- cbind(
      x = seq(screen$xlim[[1]], screen$xlim[[2]], length.out = line_n),
      y = rep(y_values[[i]], line_n)
    )
    registered <- transform_points(raw)
    data.frame(x = registered[, 1], y = registered[, 2], group = paste0("h", i))
  })
  dplyr::bind_rows(c(vertical, horizontal))
}

# Shared plot vocabulary ---------------------------------------------------

# Display labels for the three gaze paths. Replay speaks of encoding/recall,
# Transport of reference/source; both use the same colour roles.
gaze_path_labels <- function(engine = c("transport", "replay")) {
  engine <- match.arg(engine)
  if (engine == "replay") {
    c(reference = "Encoding", registered = "Recall (registered)",
      raw = "Recall (raw)")
  } else {
    c(reference = "Reference", registered = "Source (registered)",
      raw = "Source (raw)")
  }
}

gaze_path_scales <- function(labels) {
  col <- eyesim_colours()
  values <- stats::setNames(
    c(col[["reference"]], col[["source"]], col[["raw"]]),
    labels[c("reference", "registered", "raw")]
  )
  linetypes <- stats::setNames(c("solid", "solid", "22"), names(values))
  shapes <- stats::setNames(c(16, 17, 1), names(values))
  list(
    ggplot2::scale_colour_manual(values = values, breaks = names(values), name = NULL),
    ggplot2::scale_linetype_manual(values = linetypes, breaks = names(values), name = NULL),
    ggplot2::scale_shape_manual(values = shapes, breaks = names(values), name = NULL)
  )
}

gaze_path_frame <- function(coords, mass, role, labels) {
  data.frame(
    x = coords[, 1], y = coords[, 2], mass = mass,
    # raw first: paths and points are drawn in level order, so the raw path
    # sits underneath; legend order is fixed by the scale breaks
    role = factor(labels[[role]], levels = labels[c("raw", "reference", "registered")])
  )
}

gaze_axis_unit <- function(spec) {
  unit <- tryCatch(transport_v3_spatial_unit(spec), error = function(e) NULL)
  if (is.null(unit) || identical(unit, "native_coordinate_units")) NULL else unit
}

# Spatial overlay of reference, raw and registered source paths plus
# correspondence arrows (opacity proportional to mass, anchored at zero).
gaze_spatial_overlay <- function(alignment, arrows, arrow_limits, labels,
                                 title, subtitle) {
  paths <- rbind(
    gaze_path_frame(alignment$source$coords, alignment$source$mass, "raw", labels),
    gaze_path_frame(alignment$reference$coords, alignment$reference$mass,
                    "reference", labels),
    gaze_path_frame(alignment$registered_source$coords,
                    alignment$registered_source$mass, "registered", labels)
  )
  grid <- gaze_warp_grid(alignment)
  unit <- gaze_axis_unit(alignment$spec)
  col <- eyesim_colours()

  plot <- ggplot2::ggplot()
  if (!is.null(grid)) {
    plot <- plot + ggplot2::geom_path(
      data = grid,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   group = .data[["group"]]),
      colour = col[["grid"]], linewidth = 0.3
    )
  }
  raw_first <- paths[order(paths$role), , drop = FALSE]
  plot +
    ggplot2::geom_path(
      data = raw_first,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   colour = .data[["role"]], linetype = .data[["role"]],
                   group = .data[["role"]]),
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = raw_first,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]], size = .data[["mass"]],
                   colour = .data[["role"]], shape = .data[["role"]]),
      stroke = 0.8
    ) +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        xend = .data[["xend"]], yend = .data[["yend"]],
        alpha = .data[["weight"]]
      ),
      colour = col[["correspondence"]], linewidth = 0.4,
      arrow = grid::arrow(length = grid::unit(0.07, "inches"), type = "closed")
    ) +
    gaze_path_scales(labels) +
    ggplot2::scale_size_area(max_size = 4.5, guide = "none") +
    ggplot2::scale_alpha_continuous(
      range = c(0.15, 0.9), limits = arrow_limits, guide = "none"
    ) +
    ggplot2::coord_equal() +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 2.2))
    ) +
    ggplot2::labs(
      title = title, subtitle = soft_wrap(subtitle),
      x = if (is.null(unit)) NULL else paste0("x (", unit, ")"),
      y = if (is.null(unit)) NULL else paste0("y (", unit, ")")
    ) +
    theme_eyesim() +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = ggplot2::unit(24, "pt"))
}

# Two-rail braid: reference on top, source below, ribbons for correspondence.
# Width and opacity are proportional to mass (zero-anchored). Both rails use
# point area for fixation mass; source points are filled in proportion to
# the share of their mass that has a correspondence (`share`, 0-1), so
# faint source points are unmatched (Transport) or background (Replay).
gaze_braid_rails <- function(ribbons, reference, source,
                             rail_labels, title, subtitle, x_label,
                             share_note) {
  col <- eyesim_colours()
  source$fill <- scales::alpha(col[["source"]], 0.1 + 0.9 * pmin(pmax(source$share, 0), 1))
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = c(0, 1), colour = col[["rule"]], linewidth = 0.5) +
    ggplot2::geom_segment(
      data = ribbons,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        xend = .data[["xend"]], yend = .data[["yend"]],
        linewidth = .data[["relative_mass"]],
        alpha = .data[["relative_mass"]]
      ),
      colour = col[["correspondence"]], lineend = "butt"
    ) +
    ggplot2::geom_point(
      data = reference,
      ggplot2::aes(x = .data[["time"]], y = .data[["rail"]],
                   size = .data[["mass"]]),
      colour = col[["reference"]], shape = 16
    ) +
    ggplot2::geom_point(
      data = source,
      ggplot2::aes(x = .data[["time"]], y = .data[["rail"]],
                   size = .data[["mass"]], fill = .data[["fill"]]),
      colour = scales::alpha(col[["source"]], 0.6), shape = 24, stroke = 0.5
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_linewidth_continuous(
      range = c(0.2, 4.5), limits = c(0, NA), guide = "none"
    ) +
    ggplot2::scale_alpha_continuous(
      range = c(0.2, 0.75), limits = c(0, NA), guide = "none"
    ) +
    ggplot2::scale_size_area(max_size = 4.5, guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = c(0, 1), labels = rail_labels, limits = c(-0.12, 1.12)
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.03)) +
    ggplot2::labs(
      title = title, subtitle = soft_wrap(subtitle),
      caption = soft_wrap(paste0("Point area = fixation mass; ", rail_labels[[1]],
                                 " fill = ", share_note, ".")),
      x = x_label, y = NULL
    ) +
    theme_eyesim() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        face = "bold", colour = col[["ink"]], size = ggplot2::rel(0.9)
      )
    )
}

# Horizontal bars for [0, 1] diagnostics with printed values.
gaze_diagnostic_bars <- function(values, title, subtitle, x_label) {
  col <- eyesim_colours()
  values$label <- ifelse(
    values$value > 0 & values$value < 0.005,
    formatC(values$value, digits = 1, format = "e"),
    ifelse(values$value < 1 & values$value >= 0.995,
           ifelse(values$value >= 0.9995, ">0.999",
                  formatC(values$value, digits = 3, format = "f")),
           formatC(values$value, digits = 2, format = "f"))
  )
  values$inside <- values$value > 0.82
  ggplot2::ggplot(
    values,
    ggplot2::aes(x = .data[["value"]], y = .data[["component"]])
  ) +
    ggplot2::geom_col(fill = col[["diagnostic"]], width = 0.62) +
    ggplot2::geom_text(
      data = values[!values$inside, , drop = FALSE],
      ggplot2::aes(label = .data[["label"]]),
      hjust = -0.2, size = 3, colour = col[["ink"]]
    ) +
    ggplot2::geom_text(
      data = values[values$inside, , drop = FALSE],
      ggplot2::aes(label = .data[["label"]]),
      hjust = 1.2, size = 3, colour = "white"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(title = title, subtitle = soft_wrap(subtitle),
                  x = x_label, y = NULL) +
    theme_eyesim() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(colour = col[["ink"]])
    )
}

# Transport panels ---------------------------------------------------------

gaze_overlay_plot <- function(alignment) {
  source_mass <- colSums(alignment$coupling)
  destinations <- t(alignment$coupling) %*% alignment$reference$coords
  destinations <- sweep(destinations, 1, source_mass, FUN = "/")
  arrows <- data.frame(
    x = alignment$registered_source$coords[, 1],
    y = alignment$registered_source$coords[, 2],
    xend = destinations[, 1],
    yend = destinations[, 2],
    weight = source_mass
  )
  gaze_spatial_overlay(
    alignment, arrows, arrow_limits = c(0, NA),
    labels = gaze_path_labels("transport"),
    title = "Registered spatial overlay",
    subtitle = paste0(
      "scale = ", signif(warp_parameters(alignment$warp)$scale, 3)
    )
  )
}

gaze_braid_plot <- function(alignment, coupling_threshold = 0.005,
                            title = "Gaze braid",
                            subtitle_prefix = "Regularized transport coupling") {
  coupling <- alignment$coupling
  relative <- coupling / sum(coupling)
  selected <- which(relative >= coupling_threshold, arr.ind = TRUE)
  if (nrow(selected) == 0L) {
    selected <- which(relative == max(relative), arr.ind = TRUE)[1, , drop = FALSE]
  }
  ribbons <- data.frame(
    x = alignment$reference$time[selected[, 1]],
    y = 1,
    xend = alignment$registered_source$time[selected[, 2]],
    yend = 0,
    mass = coupling[selected],
    relative_mass = relative[selected]
  )
  reference <- data.frame(
    time = alignment$reference$time,
    rail = 1,
    mass = alignment$reference$mass
  )
  source_mass <- alignment$registered_source$mass
  source <- data.frame(
    time = alignment$registered_source$time,
    rail = 0,
    mass = source_mass,
    share = ifelse(source_mass > 0, colSums(coupling) / source_mass, 0)
  )
  gaze_braid_rails(
    ribbons, reference, source,
    rail_labels = c("source", "reference"),
    share_note = paste0(
      "share of the fixation's mass carried by the coupling (faint = ",
      "little; MAP coverage ", signif(sum(coupling), 2), ")"
    ),
    title = title,
    subtitle = paste0(
      subtitle_prefix, "; connections shown above ", coupling_threshold,
      " of transported mass"
    ),
    x_label = "Normalized trial time"
  )
}

# Replay panels ------------------------------------------------------------

gaze_replay_overlay_plot <- function(alignment) {
  replay_probability <- rowSums(alignment$replay_posterior)
  keep <- replay_probability > 1e-8 & stats::complete.cases(alignment$barycentric)
  arrows <- data.frame(
    x = alignment$grid$coords[keep, 1],
    y = alignment$grid$coords[keep, 2],
    xend = alignment$barycentric[keep, 1],
    yend = alignment$barycentric[keep, 2],
    weight = replay_probability[keep]
  )
  gaze_spatial_overlay(
    alignment, arrows, arrow_limits = c(0, 1),
    labels = gaze_path_labels("replay"),
    title = "Registered spatial overlay",
    subtitle = paste0(
      "arrows show posterior replay destinations; scale = ",
      signif(warp_parameters(alignment$warp)$scale, 3)
    )
  )
}

gaze_replay_braid_plot <- function(alignment, coupling_threshold = 0.005) {
  posterior <- alignment$replay_posterior
  relative <- posterior / sum(posterior)
  selected <- which(relative >= coupling_threshold, arr.ind = TRUE)
  if (nrow(selected) == 0L) {
    selected <- which(relative == max(relative), arr.ind = TRUE)[1, , drop = FALSE]
  }
  ribbons <- data.frame(
    x = alignment$reference$time[selected[, 2]],
    y = 1,
    xend = alignment$grid$time[selected[, 1]],
    yend = 0,
    relative_mass = relative[selected]
  )
  reference <- data.frame(
    time = alignment$reference$time,
    rail = 1,
    mass = alignment$reference$mass
  )
  replay <- rowSums(posterior)
  n_grid <- length(alignment$grid$time)
  recall_mass <- if (length(alignment$source$mass) == n_grid) {
    alignment$source$mass
  } else {
    rep(1 / n_grid, n_grid)
  }
  source <- data.frame(
    time = alignment$grid$time,
    rail = 0,
    mass = recall_mass,
    share = replay
  )
  gaze_braid_rails(
    ribbons, reference, source,
    rail_labels = c("recall", "encoding"),
    share_note = "posterior replay probability (faint = mostly background)",
    title = "Posterior replay braid",
    subtitle = paste(
      "Ribbons are posterior correspondence mass; threshold =",
      coupling_threshold
    ),
    x_label = "Normalized gaze time"
  )
}

gaze_replay_diagnostic_plot <- function(alignment) {
  diagnostics <- data.frame(
    component = c("replay coverage", "background coverage", "restart rate"),
    value = c(
      alignment$diagnostics$replay_coverage,
      alignment$diagnostics$background_coverage,
      alignment$diagnostics$restart_rate
    )
  )
  diagnostics$component <- factor(
    diagnostics$component, levels = rev(diagnostics$component)
  )
  unit <- gaze_axis_unit(alignment$spec)
  gaze_diagnostic_bars(
    diagnostics,
    title = "Replay posterior diagnostics",
    subtitle = paste0(
      "spatial RMSE = ", signif(alignment$diagnostics$spatial_rmse, 4),
      if (is.null(unit)) "" else paste0(" ", unit),
      "; expected restarts = ",
      signif(alignment$diagnostics$expected_restarts, 3)
    ),
    x_label = "Posterior expectation"
  )
}



#' Plot a GazeWeave Replay alignment
#'
#' @param object A `gaze_replay_alignment`.
#' @param type One of `"overlay"`, `"braid"`, or `"diagnostics"`.
#' @param coupling_threshold Minimum proportion of posterior correspondence
#'   drawn in the braid.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
autoplot.gaze_replay_alignment <- function(
    object, type = c("overlay", "braid", "diagnostics"),
    coupling_threshold = 0.005, ...) {
  type <- match.arg(type)
  if (type == "overlay") return(gaze_replay_overlay_plot(object))
  if (type == "braid") {
    return(gaze_replay_braid_plot(
      object, coupling_threshold = coupling_threshold
    ))
  }
  gaze_replay_diagnostic_plot(object)
}


select_gaze_weave_result <- function(fit, row = NULL, key = NULL) {
  if (!is.null(row) && !is.null(key)) {
    stop("Use only one of row and key.")
  }
  if (!is.null(key)) {
    if (!is.list(key) || is.null(names(key)) || !all(names(key) %in% names(fit$results))) {
      stop("key must be a named list of result-column values.")
    }
    keep <- rep(TRUE, nrow(fit$results))
    for (name in names(key)) {
      keep <- keep & fit$results[[name]] %in% key[[name]]
    }
    rows <- which(keep)
    if (length(rows) != 1L) {
      stop("key must select exactly one fitted row; selected ", length(rows), ".")
    }
    return(rows)
  }
  if (is.null(row)) {
    row <- 1L
  }
  row <- as.integer(row)
  if (length(row) != 1L || is.na(row) || row < 1L || row > nrow(fit$results)) {
    stop("row must select one fitted result.")
  }
  row
}


autoplot_gaze_engine_fit <- function(object, row, key, type,
                                     coupling_threshold, engine_label) {
  selected <- select_gaze_weave_result(object, row = row, key = key)
  alignment <- object$results$alignment[[selected]]
  if (type != "combined") {
    return(ggplot2::autoplot(
      alignment, type = type, coupling_threshold = coupling_threshold
    ))
  }
  plots <- list(
    ggplot2::autoplot(alignment, type = "overlay"),
    ggplot2::autoplot(
      alignment, type = "braid", coupling_threshold = coupling_threshold
    ),
    ggplot2::autoplot(alignment, type = "diagnostics")
  )
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Install package 'patchwork' to combine GazeWeave panels; returning a list.")
    return(plots)
  }
  title <- paste0(
    engine_label, " evidence: gaze_info_bits = ",
    signif(object$results$gaze_info_bits[[selected]], 3), " bits"
  )
  gaze_combine_panels(
    plots, heights = c(2, 1.1, 0.8), title = title,
    subtitle = paste(
      "gaze_info_bits = log2(p_true / prior_true): information about the",
      "true candidate beyond its prior; 0 means none, negative means p_true",
      "fell below its prior."
    )
  )
}

# Stack panels into one figure with a shared legend and title block.
gaze_combine_panels <- function(plots, heights, title, subtitle = NULL) {
  patchwork::wrap_plots(plots, ncol = 1, heights = heights) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      subtitle = soft_wrap(subtitle),
      theme = theme_eyesim() +
        ggplot2::theme(plot.title = ggplot2::element_text(
          face = "bold", size = ggplot2::rel(1.3)
        ))
    ) &
    ggplot2::theme(legend.position = "bottom")
}

#' Plot a fitted GazeWeave Replay analysis
#'
#' @param object A `gaze_replay_fit`.
#' @param row Optional result row number.
#' @param key Optional named list selecting one result by identifiers.
#' @param type One of `"combined"`, `"overlay"`, `"braid"`, or
#'   `"diagnostics"`.
#' @param coupling_threshold Minimum posterior correspondence proportion drawn
#'   in the braid.
#' @param ... Unused.
#'
#' @return A combined patchwork when available, or an individual `ggplot`.
#' @export
autoplot.gaze_replay_fit <- function(
    object, row = NULL, key = NULL,
    type = c("combined", "overlay", "braid", "diagnostics"),
    coupling_threshold = 0.005, ...) {
  type <- match.arg(type)
  autoplot_gaze_engine_fit(
    object, row, key, type, coupling_threshold, "Replay"
  )
}
