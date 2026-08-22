# Transport visual evidence ---------------------------------------------

transport_v3_raw_registered_plot <- function(alignment) {
  raw <- data.frame(
    x = alignment$source$coords[, 1],
    y = alignment$source$coords[, 2],
    mass = alignment$source$mass,
    index = seq_len(nrow(alignment$source$coords)),
    view = "Raw source"
  )
  registered <- data.frame(
    x = alignment$registered_source$coords[, 1],
    y = alignment$registered_source$coords[, 2],
    mass = alignment$registered_source$mass,
    index = seq_len(nrow(alignment$registered_source$coords)),
    view = "Registered source"
  )
  paths <- rbind(raw, registered)
  reference <- data.frame(
    x = alignment$reference$coords[, 1],
    y = alignment$reference$coords[, 2],
    mass = alignment$reference$mass,
    index = seq_len(nrow(alignment$reference$coords))
  )
  ggplot2::ggplot(
    paths,
    ggplot2::aes(x = .data[["x"]], y = .data[["y"]])
  ) +
    ggplot2::geom_path(
      data = reference, colour = "#1B263B", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = reference,
      ggplot2::aes(size = .data[["mass"]]),
      colour = "#1B263B"
    ) +
    ggplot2::geom_path(colour = "#D1495B", linewidth = 0.8) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data[["mass"]]), colour = "#D1495B"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data[["view"]])) +
    ggplot2::scale_size_continuous(range = c(2, 6), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "Raw versus registered source path",
      subtitle = paste0(
        "navy = reference; red = source; contraction scale = ",
        signif(alignment$diagnostics$warp_scale, 4)
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal()
}

transport_v3_diagnostic_plot <- function(alignment) {
  values <- data.frame(
    diagnostic = c(
      "matched coverage", "MAP coverage", "local-order preservation"
    ),
    value = c(
      alignment$diagnostics$matched_coverage,
      alignment$diagnostics$map_coverage,
      alignment$diagnostics$local_order_preservation
    ),
    stringsAsFactors = FALSE
  )
  values$diagnostic <- factor(
    values$diagnostic, levels = rev(values$diagnostic)
  )
  unit <- transport_v3_spatial_unit(alignment$spec)
  ggplot2::ggplot(
    values,
    ggplot2::aes(x = .data[["diagnostic"]], y = .data[["value"]])
  ) +
    ggplot2::geom_col(fill = "#6C5CE7") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Transport alignment diagnostics",
      subtitle = paste0(
        "optimized correspondence, not a posterior; spatial RMSE = ",
        signif(alignment$diagnostics$spatial_rmse, 4), " ", unit,
        "; solver converged = ", isTRUE(alignment$convergence$converged)
      ),
      x = NULL, y = "Conditional diagnostic"
    ) +
    ggplot2::theme_minimal()
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
  plot$labels$subtitle <- paste0(
    plot$labels$subtitle,
    "; arrows summarize optimized correspondence, not posterior probability"
  )
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

transport_v3_evidence_plot <- function(record) {
  ledger <- data.frame(
    endpoint = "gaze_info_bits",
    value = record$gaze_info_bits,
    stringsAsFactors = FALSE
  )
  diagnostics <- record$diagnostics
  ggplot2::ggplot(
    ledger,
    ggplot2::aes(x = .data[["endpoint"]], y = .data[["value"]])
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_col(fill = "#6C5CE7", width = 0.5) +
    ggplot2::labs(
      title = "One-score evidence ledger",
      subtitle = paste0(
        "log2(p_true / prior_true); rank ", diagnostics$template_rank,
        " of ", diagnostics$candidate_count, "; ",
        diagnostics$episode_count, " equal-prior episode(s)"
      ),
      caption = paste(
        "Candidate-dependent calibrated evidence; alignment diagnostics",
        "do not add to this score."
      ),
      x = NULL, y = "Information relative to declared candidate prior (bits)"
    ) +
    ggplot2::theme_minimal()
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
  if (type != "combined") {
    plot <- transport_v3_alignment_plot(
      selected_episode$alignment, type, coupling_threshold
    )
    plot$labels$caption <- paste0(
      "Alignment panel: episode ", selected_episode$id, " of ",
      selected_episode$count,
      ". gaze_info_bits retains the fixed equal-prior episode mixture."
    )
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
    "Transport: gaze_info_bits = ", signif(record$gaze_info_bits, 4),
    "; diagnostic episode ", selected_episode$id, " of ",
    selected_episode$count
  )
  patchwork::wrap_plots(plots, ncol = 1, heights = c(1, 2, 1.3, 2, 1)) +
    patchwork::plot_annotation(
      title = title,
      subtitle = paste(
        "Evidence mixes all retained episodes equally; correspondences are",
        "optimized alignments, not posterior replay probabilities."
      )
    )
}
