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

gaze_overlay_plot <- function(alignment) {
  reference <- data.frame(
    x = alignment$reference$coords[, 1],
    y = alignment$reference$coords[, 2],
    mass = alignment$reference$mass,
    index = seq_len(nrow(alignment$reference$coords))
  )
  source <- data.frame(
    x = alignment$source$coords[, 1],
    y = alignment$source$coords[, 2],
    mass = alignment$source$mass,
    index = seq_len(nrow(alignment$source$coords))
  )
  registered <- data.frame(
    x = alignment$registered_source$coords[, 1],
    y = alignment$registered_source$coords[, 2],
    mass = alignment$registered_source$mass,
    index = seq_len(nrow(alignment$registered_source$coords))
  )
  source_mass <- colSums(alignment$coupling)
  destinations <- t(alignment$coupling) %*% alignment$reference$coords
  destinations <- sweep(destinations, 1, source_mass, FUN = "/")
  arrows <- data.frame(
    x = registered$x,
    y = registered$y,
    xend = destinations[, 1],
    yend = destinations[, 2],
    mass = source_mass
  )
  grid <- gaze_warp_grid(alignment)

  plot <- ggplot2::ggplot()
  if (!is.null(grid)) {
    plot <- plot + ggplot2::geom_path(
      data = grid,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   group = .data[["group"]]),
      colour = "grey88",
      linewidth = 0.3
    )
  }
  plot +
    ggplot2::geom_path(
      data = source,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = "grey65",
      linewidth = 0.7,
      linetype = 2
    ) +
    ggplot2::geom_path(
      data = reference,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = "#1B263B",
      linewidth = 0.9
    ) +
    ggplot2::geom_path(
      data = registered,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = "#D1495B",
      linewidth = 0.9
    ) +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        xend = .data[["xend"]], yend = .data[["yend"]],
        alpha = .data[["mass"]]
      ),
      colour = "#6C757D",
      linewidth = 0.35,
      arrow = grid::arrow(length = grid::unit(0.08, "inches"))
    ) +
    ggplot2::geom_point(
      data = source,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   size = .data[["mass"]]),
      colour = "grey65",
      shape = 1,
      stroke = 0.8
    ) +
    ggplot2::geom_point(
      data = reference,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   size = .data[["mass"]]),
      colour = "#1B263B"
    ) +
    ggplot2::geom_point(
      data = registered,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   size = .data[["mass"]]),
      colour = "#D1495B"
    ) +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Registered spatial overlay",
      subtitle = paste0(
        "reference (navy), raw source (grey), registered source (red); scale = ",
        signif(warp_parameters(alignment$warp)$scale, 3)
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal()
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
  source <- data.frame(
    time = alignment$registered_source$time,
    rail = 0,
    mass = alignment$registered_source$mass
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = ribbons,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        xend = .data[["xend"]], yend = .data[["yend"]],
        linewidth = .data[["relative_mass"]],
        alpha = .data[["relative_mass"]]
      ),
      colour = "#6C5CE7",
      lineend = "round"
    ) +
    ggplot2::geom_hline(yintercept = c(0, 1), colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_point(
      data = reference,
      ggplot2::aes(x = .data[["time"]], y = .data[["rail"]],
                   size = .data[["mass"]]),
      colour = "#1B263B"
    ) +
    ggplot2::geom_point(
      data = source,
      ggplot2::aes(x = .data[["time"]], y = .data[["rail"]],
                   size = .data[["mass"]]),
      colour = "#D1495B"
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 4), guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.15, 0.8), guide = "none") +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = c(0, 1),
      labels = c("source", "reference"),
      limits = c(-0.08, 1.08)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste(
        subtitle_prefix, "; connections shown above", coupling_threshold,
        "of transported mass"
      ),
      x = "Normalized trial time",
      y = NULL
    ) +
    ggplot2::theme_minimal()
}

gaze_replay_overlay_plot <- function(alignment) {
  reference <- data.frame(
    x = alignment$reference$coords[, 1],
    y = alignment$reference$coords[, 2],
    mass = alignment$reference$mass
  )
  source <- data.frame(
    x = alignment$source$coords[, 1],
    y = alignment$source$coords[, 2],
    mass = alignment$source$mass
  )
  registered <- data.frame(
    x = alignment$registered_source$coords[, 1],
    y = alignment$registered_source$coords[, 2],
    mass = alignment$registered_source$mass
  )
  replay_probability <- rowSums(alignment$replay_posterior)
  keep <- replay_probability > 1e-8 & stats::complete.cases(alignment$barycentric)
  arrows <- data.frame(
    x = alignment$grid$coords[keep, 1],
    y = alignment$grid$coords[keep, 2],
    xend = alignment$barycentric[keep, 1],
    yend = alignment$barycentric[keep, 2],
    probability = replay_probability[keep]
  )
  grid <- gaze_warp_grid(alignment)

  plot <- ggplot2::ggplot()
  if (!is.null(grid)) {
    plot <- plot + ggplot2::geom_path(
      data = grid,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   group = .data[["group"]]),
      colour = "grey88",
      linewidth = 0.3
    )
  }
  plot +
    ggplot2::geom_path(
      data = source,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = "grey65", linewidth = 0.7, linetype = 2
    ) +
    ggplot2::geom_path(
      data = reference,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = "#1B263B", linewidth = 0.9
    ) +
    ggplot2::geom_path(
      data = registered,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]]),
      colour = "#D1495B", linewidth = 0.9
    ) +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        xend = .data[["xend"]], yend = .data[["yend"]],
        alpha = .data[["probability"]]
      ),
      colour = "#6C757D", linewidth = 0.3,
      arrow = grid::arrow(length = grid::unit(0.06, "inches"))
    ) +
    ggplot2::geom_point(
      data = source,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   size = .data[["mass"]]),
      colour = "grey65", shape = 1, stroke = 0.8
    ) +
    ggplot2::geom_point(
      data = reference,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   size = .data[["mass"]]),
      colour = "#1B263B"
    ) +
    ggplot2::geom_point(
      data = registered,
      ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                   size = .data[["mass"]]),
      colour = "#D1495B"
    ) +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.08, 0.75), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Registered spatial overlay",
      subtitle = paste0(
        "arrows show posterior replay destinations; scale = ",
        signif(warp_parameters(alignment$warp)$scale, 3)
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal()
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
  source <- data.frame(
    time = alignment$grid$time,
    rail = 0,
    replay = rowSums(posterior)
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = ribbons,
      ggplot2::aes(
        x = .data[["x"]], y = .data[["y"]],
        xend = .data[["xend"]], yend = .data[["yend"]],
        linewidth = .data[["relative_mass"]],
        alpha = .data[["relative_mass"]]
      ),
      colour = "#6C5CE7", lineend = "round"
    ) +
    ggplot2::geom_hline(yintercept = c(0, 1), colour = "grey55", linewidth = 0.35) +
    ggplot2::geom_point(
      data = reference,
      ggplot2::aes(x = .data[["time"]], y = .data[["rail"]],
                   size = .data[["mass"]]),
      colour = "#1B263B"
    ) +
    ggplot2::geom_point(
      data = source,
      ggplot2::aes(x = .data[["time"]], y = .data[["rail"]],
                   alpha = .data[["replay"]]),
      colour = "#D1495B", size = 1.6
    ) +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 4), guide = "none") +
    ggplot2::scale_alpha_continuous(range = c(0.1, 0.85), guide = "none") +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = c(0, 1), labels = c("recall", "encoding"),
      limits = c(-0.08, 1.08)
    ) +
    ggplot2::labs(
      title = "Posterior replay braid",
      subtitle = paste(
        "Ribbons are posterior correspondence mass; threshold =",
        coupling_threshold
      ),
      x = "Normalized gaze time", y = NULL
    ) +
    ggplot2::theme_minimal()
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
  ggplot2::ggplot(
    diagnostics,
    ggplot2::aes(x = .data[["component"]], y = .data[["value"]])
  ) +
    ggplot2::geom_col(fill = "#3A86FF") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Replay posterior diagnostics",
      subtitle = paste0(
        "spatial RMSE = ", signif(alignment$diagnostics$spatial_rmse, 4),
        "; expected restarts = ",
        signif(alignment$diagnostics$expected_restarts, 3)
      ),
      x = NULL, y = "Posterior expectation"
    ) +
    ggplot2::theme_minimal()
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
    "gaze_info_bits = ",
    signif(object$results$gaze_info_bits[[selected]], 4),
    " (", engine_label, ")"
  )
  patchwork::wrap_plots(plots, ncol = 1, heights = c(2, 1.3, 1)) +
    patchwork::plot_annotation(title = title)
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
