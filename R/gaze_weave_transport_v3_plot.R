# Transport visual evidence ---------------------------------------------

transport_v3_raw_registered_plot <- function(alignment) {
  labels <- gaze_path_labels("transport")
  views <- c("Raw source", "Registered source")
  source <- rbind(
    cbind(gaze_path_frame(alignment$source$coords, alignment$source$mass,
                          "raw", labels), view = views[[1]]),
    cbind(gaze_path_frame(alignment$registered_source$coords,
                          alignment$registered_source$mass,
                          "registered", labels), view = views[[2]])
  )
  reference <- gaze_path_frame(alignment$reference$coords,
                               alignment$reference$mass, "reference", labels)
  reference <- rbind(cbind(reference, view = views[[1]]),
                     cbind(reference, view = views[[2]]))
  paths <- rbind(reference, source)
  paths$view <- factor(paths$view, levels = views)
  unit <- gaze_axis_unit(alignment$spec)
  ggplot2::ggplot(
    paths,
    ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                 colour = .data[["role"]], group = .data[["role"]])
  ) +
    ggplot2::geom_path(ggplot2::aes(linetype = .data[["role"]]), linewidth = 0.6) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data[["mass"]], shape = .data[["role"]]), stroke = 0.8
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["view"]])) +
    gaze_path_scales(labels) +
    ggplot2::scale_size_area(max_size = 4.5, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::coord_equal() +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(size = 2.2))
    ) +
    ggplot2::labs(
      title = "Raw versus registered source path",
      subtitle = paste0(
        "contraction scale = ", signif(alignment$diagnostics$warp_scale, 4)
      ),
      x = if (is.null(unit)) NULL else paste0("x (", unit, ")"),
      y = if (is.null(unit)) NULL else paste0("y (", unit, ")")
    ) +
    theme_eyesim() +
    ggplot2::theme(legend.position = "bottom",
                   legend.key.width = ggplot2::unit(24, "pt"))
}

transport_v3_diagnostic_plot <- function(alignment) {
  values <- data.frame(
    component = c(
      "matched coverage", "MAP coverage", "local-order preservation"
    ),
    value = c(
      alignment$diagnostics$matched_coverage,
      alignment$diagnostics$map_coverage,
      alignment$diagnostics$local_order_preservation
    ),
    stringsAsFactors = FALSE
  )
  values$component <- factor(
    values$component, levels = rev(values$component)
  )
  unit <- transport_v3_spatial_unit(alignment$spec)
  gaze_diagnostic_bars(
    values,
    title = "Transport alignment diagnostics",
    subtitle = paste0(
      "optimized correspondence, not a posterior; spatial RMSE = ",
      signif(alignment$diagnostics$spatial_rmse, 4), " ", unit,
      "; solver converged = ", isTRUE(alignment$convergence$converged)
    ),
    x_label = "Conditional diagnostic"
  )
}

transport_v3_alignment_plot <- function(alignment, type,
                                        coupling_threshold) {
  if (type == "raw_registered") {
    return(transport_v3_raw_registered_plot(alignment))
  }
  if (type == "diagnostics") {
    return(transport_v3_diagnostic_plot(alignment))
  }
  if (type == "braid") {
    return(gaze_braid_plot(
      alignment,
      coupling_threshold = coupling_threshold,
      title = "Transport gaze braid",
      subtitle_prefix = paste(
        "Ribbons are an optimized alignment correspondence, not a posterior"
      )
    ))
  }
  plot <- gaze_overlay_plot(alignment)
  plot$labels$title <- "Transport registered overlay"
  plot$labels$subtitle <- eyesim_wrap(paste0(
    gsub("\n", " ", plot$labels$subtitle, fixed = TRUE),
    "; arrows summarize optimized correspondence, not posterior probability"
  ))
  plot
}

#' Plot a Transport pair alignment
#'
#' @param object A `gaze_transport_alignment`.
#' @param type One of `"overlay"`, `"braid"`, `"raw_registered"`, or
#'   `"diagnostics"`.
#' @param coupling_threshold Minimum proportion of transported mass drawn in
#'   the braid.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
autoplot.gaze_transport_alignment <- function(
    object,
    type = c("overlay", "braid", "raw_registered", "diagnostics"),
    coupling_threshold = 0.005, ...) {
  type <- match.arg(type)
  transport_v3_alignment_plot(object, type, coupling_threshold)
}

transport_v3_select_episode <- function(record, episode = NULL) {
  true_key <- record$provenance$selected_true_candidate
  candidate <- record$candidate_evidence
  mixture <- attr(record, "candidate_result", exact = TRUE)
  if (is.null(mixture)) {
    stop("The Transport result does not retain its candidate alignment.")
  }
  episodes <- mixture$alignment$episodes
  ids <- names(episodes)
  if (is.null(episode)) {
    episode <- ids[[1L]]
  }
  episode <- as.character(episode)
  if (length(episode) != 1L || !episode %in% ids) {
    stop("episode must identify one retained study presentation: ",
         paste(ids, collapse = ", "), ".")
  }
  list(
    id = episode,
    count = length(ids),
    alignment = episodes[[episode]]$alignment,
    true_key = true_key,
    candidate = candidate
  )
}

# One score, drawn as a single labelled lollipop from zero.
transport_v3_evidence_plot <- function(record) {
  col <- eyesim_colours()
  ledger <- data.frame(
    endpoint = "gaze_info_bits",
    value = record$gaze_info_bits,
    stringsAsFactors = FALSE
  )
  ledger$label <- paste0(signif(ledger$value, 3), " bits")
  diagnostics <- record$diagnostics
  span <- max(abs(ledger$value), 1)
  ggplot2::ggplot(
    ledger,
    ggplot2::aes(x = .data[["value"]], y = .data[["endpoint"]])
  ) +
    ggplot2::geom_vline(xintercept = 0, colour = col[["muted"]], linewidth = 0.4) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data[["value"]],
                   yend = .data[["endpoint"]]),
      colour = col[["correspondence"]], linewidth = 1.2
    ) +
    ggplot2::geom_point(colour = col[["correspondence"]], size = 3.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[["label"]]),
      vjust = -1.1, size = 3.2, fontface = "bold", colour = col[["ink"]]
    ) +
    ggplot2::scale_x_continuous(limits = c(-span, span) * 1.1) +
    ggplot2::labs(
      title = "One-score evidence ledger",
      subtitle = eyesim_wrap(paste0(
        "log2(p_true / prior_true); rank ", diagnostics$template_rank,
        " of ", diagnostics$candidate_count, "; ",
        diagnostics$episode_count, " equal-prior episode(s)"
      )),
      caption = eyesim_wrap(paste(
        "Candidate-dependent calibrated evidence; alignment diagnostics",
        "do not add to this score."
      ), width = 80),
      x = "Information relative to declared candidate prior (bits)", y = NULL
    ) +
    theme_eyesim() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", colour = col[["ink"]])
    )
}

#' Plot one fitted Transport result
#'
#' `type = "evidence"` plots the single prespecified score. Alignment panels
#' display one named study episode; this selection never changes the score,
#' which always uses the fit's fixed equal-prior episode likelihood mixture.
#'
#' @param object A fitted object returned by [gaze_transport_cv()].
#' @param row Optional result row number.
#' @param key Optional named list selecting one result by identifiers.
#' @param type One of `"combined"`, `"evidence"`, `"overlay"`, `"braid"`,
#'   `"raw_registered"`, or `"diagnostics"`.
#' @param episode Optional retained episode identifier for an alignment panel.
#' @param coupling_threshold Minimum proportion of transported mass drawn in
#'   the braid.
#' @param ... Unused.
#'
#' @return A combined patchwork when available, a list of plots otherwise, or
#'   one `ggplot` for an individual panel.
#' @export
autoplot.gaze_transport_fit <- function(
    object, row = NULL, key = NULL,
    type = c(
      "combined", "evidence", "overlay", "braid", "raw_registered",
      "diagnostics"
    ),
    episode = NULL, coupling_threshold = 0.005, ...) {
  type <- match.arg(type)
  selected <- select_gaze_weave_result(object, row = row, key = key)
  record <- gaze_transport_result(object, row = selected)
  candidate <- transport_v3_true_candidate(object, selected)
  attr(record, "candidate_result") <- candidate$result
  if (type == "evidence") {
    return(transport_v3_evidence_plot(record))
  }
  selected_episode <- transport_v3_select_episode(record, episode)
  episode_note <- eyesim_wrap(paste0(
    "Alignment panel: episode ", selected_episode$id, " of ",
    selected_episode$count,
    ". gaze_info_bits retains the fixed equal-prior episode mixture."
  ), width = 110)
  if (type != "combined") {
    plot <- transport_v3_alignment_plot(
      selected_episode$alignment, type, coupling_threshold
    )
    existing <- plot$labels$caption
    plot$labels$caption <- if (is.null(existing)) {
      episode_note
    } else {
      paste(existing, episode_note, sep = "\n")
    }
    return(plot)
  }
  plots <- list(
    transport_v3_evidence_plot(record),
    transport_v3_alignment_plot(
      selected_episode$alignment, "overlay", coupling_threshold
    ),
    transport_v3_alignment_plot(
      selected_episode$alignment, "braid", coupling_threshold
    ),
    transport_v3_alignment_plot(
      selected_episode$alignment, "raw_registered", coupling_threshold
    ),
    transport_v3_alignment_plot(
      selected_episode$alignment, "diagnostics", coupling_threshold
    )
  )
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning(
      "Install package 'patchwork' to combine Transport panels; ",
      "returning a list."
    )
    return(plots)
  }
  title <- paste0(
    "Transport evidence: gaze_info_bits = ", signif(record$gaze_info_bits, 3),
    " bits",
    "; diagnostic episode ", selected_episode$id, " of ",
    selected_episode$count
  )
  gaze_combine_panels(
    plots, heights = c(0.45, 1.8, 0.9, 1.3, 0.75), title = title,
    subtitle = paste(
      "Evidence mixes all retained episodes equally; correspondences are",
      "optimized alignments, not posterior replay probabilities."
    )
  )
}
