# Public Transport result surface ---------------------------------------

transport_v3_spatial_unit <- function(spec) {
  if (!is.null(spec$screen)) {
    return(spec$screen$unit)
  }
  if (!is.null(spec$spatial$unit) && !identical(spec$spatial$unit, "native")) {
    return(spec$spatial$unit)
  }
  "native_coordinate_units"
}

transport_v3_true_candidate <- function(x, row) {
  candidate_table <- x$results$candidates[[row]]
  if (!is.data.frame(candidate_table) ||
      !all(c("candidate_key", "is_true") %in% names(candidate_table))) {
    stop("The selected Transport row has no candidate evidence table.")
  }
  true_key <- as.character(candidate_table$candidate_key[candidate_table$is_true])
  if (length(true_key) != 1L) {
    stop("The selected Transport row must identify one true candidate.")
  }
  candidates <- x$results$alignments[[row]]
  candidate <- candidates[[true_key]]
  if (is.null(candidate)) {
    stop("The true candidate has no retained Transport alignment.")
  }
  list(key = true_key, table = candidate_table, result = candidate)
}

transport_v3_episode_diagnostics <- function(candidate, unit) {
  episodes <- candidate$alignment$episodes
  if (!is.list(episodes) || length(episodes) == 0L) {
    stop("The selected candidate has no retained episode alignments.")
  }
  ids <- names(episodes)
  if (is.null(ids) || any(!nzchar(ids))) {
    ids <- as.character(seq_along(episodes))
  }
  weights <- candidate$diagnostics$equal_episode_weights
  weights <- weights[match(ids, names(weights))]
  value <- data.frame(
    episode = ids,
    equal_episode_weight = as.numeric(weights),
    log_score = vapply(episodes, `[[`, numeric(1), "log_score"),
    matched_coverage = vapply(episodes, function(result) {
      result$diagnostics$matched_coverage
    }, numeric(1)),
    spatial_rmse = vapply(episodes, function(result) {
      result$diagnostics$spatial_rmse
    }, numeric(1)),
    spatial_rmse_unit = rep(unit, length(episodes)),
    local_order_preservation = vapply(episodes, function(result) {
      result$diagnostics$local_order_preservation
    }, numeric(1)),
    contraction_scale = vapply(episodes, function(result) {
      result$diagnostics$warp_scale
    }, numeric(1)),
    solver_converged = vapply(episodes, function(result) {
      isTRUE(result$convergence$converged)
    }, logical(1)),
    solver_backend = vapply(episodes, function(result) {
      result$convergence$backend
    }, character(1)),
    stationarity = vapply(episodes, function(result) {
      result$diagnostics$stationarity
    }, numeric(1)),
    relative_objective_change = vapply(episodes, function(result) {
      result$diagnostics$relative_objective_change
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  if (anyNA(value$equal_episode_weight)) {
    stop("Episode weights do not match the retained episode alignments.")
  }
  value
}

#' Extract one auditable Transport result
#'
#' Transport has one primary outcome, `gaze_info_bits`. Coverage, spatial
#' error, local order, contraction, rank, and solver behavior are nested under
#' `diagnostics`; they are explanatory properties of the optimized alignment,
#' not additional evidence scores. When several study presentations are
#' available, their calibrated likelihoods receive fixed equal prior weight.
#'
#' @param x A fitted object returned by [gaze_transport_cv()].
#' @param row Optional result row number.
#' @param key Optional named list selecting one result by identifiers.
#'
#' @return A `gaze_transport_result` containing one `gaze_info_bits` value,
#'   nested diagnostics, candidate evidence, and provenance.
#' @export
gaze_transport_result <- function(x, row = NULL, key = NULL) {
  if (!inherits(x, "gaze_transport_fit")) {
    stop("x must be returned by gaze_transport_cv().")
  }
  selected <- select_gaze_weave_result(x, row = row, key = key)
  candidate <- transport_v3_true_candidate(x, selected)
  unit <- transport_v3_spatial_unit(x$spec)
  episodes <- transport_v3_episode_diagnostics(candidate$result, unit)
  result <- x$results[selected, , drop = FALSE]
  identifiers <- result[intersect(x$keys$id_columns, names(result))]
  stability <- list(
    all_converged = all(episodes$solver_converged),
    backends = unique(episodes$solver_backend),
    maximum_stationarity = max(episodes$stationarity),
    maximum_relative_objective_change = max(
      episodes$relative_objective_change
    )
  )
  structure(
    list(
      identifiers = identifiers,
      gaze_info_bits = result$gaze_info_bits[[1L]],
      diagnostics = list(
        matched_coverage = stats::weighted.mean(
          episodes$matched_coverage, episodes$equal_episode_weight
        ),
        spatial_rmse = stats::weighted.mean(
          episodes$spatial_rmse, episodes$equal_episode_weight
        ),
        spatial_rmse_unit = unit,
        local_order_preservation = stats::weighted.mean(
          episodes$local_order_preservation, episodes$equal_episode_weight
        ),
        contraction_scale = stats::weighted.mean(
          episodes$contraction_scale, episodes$equal_episode_weight
        ),
        template_rank = result$template_rank[[1L]],
        candidate_count = result$candidate_count[[1L]],
        episode_count = nrow(episodes),
        episode_normalization = "fixed_equal_prior_likelihood_mixture",
        episodes = episodes,
        solver_stability = stability
      ),
      candidate_evidence = candidate$table,
      provenance = c(
        x$provenance,
        list(
          selected_true_candidate = candidate$key,
          correspondence_semantics =
            "optimized_alignment_not_posterior_probability",
          spatial_rmse_unit = unit,
          episode_diagnostic_scope =
            "equal_weight_summary_of_true_candidate_episode_alignments"
        )
      )
    ),
    class = c("gaze_transport_result", "list")
  )
}

#' Tidy a fitted Transport analysis
#'
#' The default output contains only identifiers and the prespecified primary
#' endpoint. Use `diagnostics = TRUE` to add calibration fields and one nested
#' diagnostic record per row.
#'
#' @param x A fitted object returned by [gaze_transport_cv()].
#' @param diagnostics Include candidate/calibration columns and a nested
#'   `diagnostics` list-column.
#' @param ... Unused.
#'
#' @return A tibble. The default has identifiers and `gaze_info_bits` only.
#' @export
tidy.gaze_transport_fit <- function(x, diagnostics = FALSE, ...) {
  id_columns <- intersect(x$keys$id_columns, names(x$results))
  if (!isTRUE(diagnostics)) {
    return(tibble::as_tibble(x$results[unique(c(
      id_columns, "gaze_info_bits"
    ))]))
  }
  columns <- unique(c(
    id_columns, "gaze_info_bits", "posterior_true", "prior_true",
    "log_loss", "template_rank", "top1_credit", "candidate_count",
    "common_episode_count", "temperature", "reliability", ".cv_fold",
    "all_converged"
  ))
  columns <- intersect(columns, names(x$results))
  value <- tibble::as_tibble(x$results[columns])
  value$diagnostics <- lapply(seq_len(nrow(x$results)), function(row) {
    gaze_transport_result(x, row = row)$diagnostics
  })
  value
}

#' @export
print.gaze_transport_fit <- function(x, ...) {
  bits <- x$results$gaze_info_bits
  cat("Cross-fitted GazeWeave Transport\n")
  cat("  held-out rows:", nrow(x$results), "\n")
  cat("  primary endpoint: gaze_info_bits\n")
  cat("  finite scores:", sum(is.finite(bits)), "\n")
  if (any(is.finite(bits))) {
    cat("  mean gaze_info_bits:", signif(mean(bits, na.rm = TRUE), 5), "\n")
  }
  cat(
    "  episode normalization:",
    "fixed equal-prior likelihood mixture; one episode is the identity case\n"
  )
  invisible(x)
}
