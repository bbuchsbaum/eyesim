#' Score case-control assignments against reference maps
#' @param strata Explicit case-control locations from [sample_controls()].
#' @param refs A reference set from [reference_maps()]. Its preparation is fixed;
#'   log transforms are not applied again during scoring.
#' @param standardize Standardize each reference within subject and control draw,
#'   or retain raw reference values.
#' @param scaling Optional fitted object from [fit_score_scaling()]. Supplying it
#'   applies those parameters instead of estimating them. Explicitly combining
#'   `standardize = "none"` with a scaling object is an error.
#' @return The strata table with one column per reference and attributes
#'   `score_columns`, `reference_spec`, and `scaling` (when applied).
#' @md
#' @export
score_strata <- function(strata, refs, standardize = c("within_subject", "none"), scaling = NULL) {
  standardize <- match.arg(standardize)
  if (!is.null(scaling) && standardize == "none") {
    stop("standardize = 'none' conflicts with scaling; omit scaling for raw scores.")
  }
  if (!inherits(refs, "eyesim_references")) stop("refs must be a reference set.")
  groups <- check_strata(strata)
  check_stratum_keys(strata, groups, names(refs$keys))
  if (any(refs$refs %in% names(strata))) stop("Reference names collide with existing strata columns.")
  values <- sample_density(refs, strata)
  out <- strata
  out[refs$refs] <- values
  attr(out, "score_columns") <- refs$refs
  attr(out, "reference_spec") <- refs$spec
  if (!is.null(scaling)) return(apply_score_scaling(out, scaling))
  if (standardize == "within_subject") {
    out <- apply_score_scaling(out, fit_score_scaling(out, refs$refs,
                                                     c(refs$spec$subject, "draw")))
  }
  out
}

#' Fit and apply group-specific score standardization
#' @param scores A table with numeric score columns; missing values are allowed.
#' @param columns Score column names. Defaults to the scoring result's metadata.
#' @param by Grouping columns. Use character(0) for a single global fit.
#' @param scaling A fitted scaling object.
#' @details Fitting uses finite observations, the arithmetic mean and sample SD.
#'   A constant column maps to zero; an entirely missing column remains missing.
#'   Application rejects unseen groups and never refits. Map normalization is a
#'   separate operation. Scores that are already standardized cannot be scaled
#'   again. Fitting separately by control draw prevents additional draws from
#'   changing earlier draws' scaling.
#' @return `fit_score_scaling()` returns an `eyesim_score_scaling` list with group
#'   keys, columns, centers, scales, and counts. `apply_score_scaling()` returns
#'   a score table with the fitted object attached as `scaling`.
#' @md
#' @export
fit_score_scaling <- function(scores, columns = attr(scores, "score_columns"), by = c("subject", "draw")) {
  check_score_columns(scores, columns)
  if (!is.null(attr(scores, "scaling"))) stop("Fit scaling on raw scores, not already standardized scores.")
  if (length(by)) check_key_columns(scores, by)
  if (!nrow(scores)) stop("Cannot fit scaling on an empty table.")
  fit_scaling_values(as.matrix(scores[columns]), scaling_layout(scores, by))
}

scaling_layout <- function(scores, by) {
  groups <- if (length(by)) unique(as.data.frame(scores[by])) else data.frame(.group = 1L)
  ids <- if (length(by)) match_map_rows(scores, groups, by) else rep(1L, nrow(scores))
  list(groups = groups, ids = ids, rows = split(seq_len(nrow(scores)), ids), by = by)
}

# Reuse group membership for replay draws instead of reconstructing key tables.
fit_scaling_values <- function(values, layout, keep = NULL) {
  columns <- colnames(values)
  center <- scale <- counts <- matrix(NA_real_, nrow(layout$groups), length(columns), dimnames = list(NULL, columns))
  rows <- layout$rows
  if (!is.null(keep)) rows <- lapply(rows, function(r) r[keep[r]])
  for (g in seq_len(nrow(layout$groups))) for (j in seq_along(columns)) {
    vals <- values[rows[[as.character(g)]], j]
    vals <- vals[!is.na(vals)]
    counts[g, j] <- length(vals)
    if (length(vals)) {
      center[g, j] <- mean(vals)
      sd <- if (length(vals) > 1L) stats::sd(vals) else 0
      scale[g, j] <- if (sd > 0) sd else 1
    }
  }
  structure(list(groups = layout$groups, by = layout$by, columns = columns, center = center,
                 scale = scale, n = counts), class = "eyesim_score_scaling")
}

#' @rdname fit_score_scaling
#' @md
#' @export
apply_score_scaling <- function(scores, scaling) {
  if (!inherits(scaling, "eyesim_score_scaling")) stop("scaling must be a fitted score scaling object.")
  if (!is.null(attr(scores, "scaling"))) stop("Scores are already standardized.")
  check_score_columns(scores, scaling$columns)
  ids <- if (length(scaling$by)) match_map_rows(scores, scaling$groups, scaling$by, strict = TRUE) else rep(1L, nrow(scores))
  if (anyNA(ids)) stop("Scaling has no fitted parameters for some groups.")
  out <- scores
  for (j in seq_along(scaling$columns)) {
    nm <- scaling$columns[j]
    out[[nm]] <- (out[[nm]] - scaling$center[ids, j]) / scaling$scale[ids, j]
  }
  attr(out, "scaling") <- scaling
  out
}

#' Summarize case minus mean-control reference scores
#' @param scored A table returned by [score_strata()].
#' @param by Return one contrast per fixation/control draw or an equally weighted
#'   mean of fixation contrasts per subject/trial/control draw.
#' @param subject,trial Subject and trial column names. Subject defaults to the
#'   reference specification when available.
#' @param columns Score columns, defaulting to the score metadata.
#' @return A data frame with one contrast per reference. Control draws remain
#'   separate. Missing scores propagate rather than silently changing the pool.
#' @md
#' @export
case_control_contrast <- function(scored, by = c("fixation", "trial"),
                                  subject = NULL, trial = "trial",
                                  columns = attr(scored, "score_columns")) {
  by <- match.arg(by)
  check_score_columns(scored, columns)
  groups <- check_strata(scored)
  layout <- contrast_layout(scored, groups, by, subject, trial)
  contrast_values(as.matrix(scored[columns]), layout)
}

# Cache row membership once; exclusion changes counts, never case identities.
contrast_layout <- function(scored, groups, by, subject, trial) {
  case_rows <- vapply(groups, function(g) g[scored$is_case[g]], integer(1))
  ids <- integer(nrow(scored))
  ids[unlist(groups, use.names = FALSE)] <- rep(seq_along(groups), lengths(groups))
  out <- as.data.frame(scored[case_rows, , drop = FALSE])
  trial_ids <- trial_counts <- NULL
  output_rows <- case_rows
  if (by == "trial") {
    if (is.null(subject)) subject <- attr(scored, "reference_spec")$subject %||% "subject"
    keys <- unique(c(subject, trial, "draw"))
    check_key_columns(out, keys)
    trial_ids <- key_group_ids(out, keys)
    trial_counts <- tabulate(trial_ids)
    first <- match(seq_along(trial_counts), trial_ids)
    out <- out[first, keys, drop = FALSE]
    out$n_fixations <- trial_counts
    output_rows <- case_rows[first]
  }
  rownames(out) <- NULL
  # Contrasts no longer carry row-level fits or donor-source tables.
  attr(out, "scaling") <- NULL
  attr(out, "control_sources") <- NULL
  list(case_rows = case_rows, ids = ids, controls = which(!scored$is_case),
       trial_ids = trial_ids, trial_counts = trial_counts,
       output_rows = output_rows, template = out)
}

# The public contrast and repeated-null paths share the same reduction.
contrast_values <- function(values, layout, keep = NULL, scaling = NULL, scaling_ids = NULL) {
  controls <- layout$controls
  if (!is.null(keep)) controls <- controls[keep[controls]]
  counts <- tabulate(layout$ids[controls], length(layout$case_rows))
  if (any(counts == 0L)) {
    stop("Replacement-item exclusion leaves a stratum without controls. Increase n in sample_controls() or explicitly choose control_policy = 'fixed' for a different baseline.")
  }
  # Subtract before summing: a large common offset must not erase a tiny
  # contrast that becomes appreciable after dividing by a small score SD.
  differences <- values[layout$case_rows[layout$ids[controls]], , drop = FALSE] -
    values[controls, , drop = FALSE]
  result <- rowsum(differences, layout$ids[controls], reorder = TRUE) / counts
  out <- layout$template
  if (!is.null(layout$trial_ids)) {
    result <- rowsum(result, layout$trial_ids, reorder = TRUE) / layout$trial_counts
  } else {
    out$n_controls <- counts
  }
  # Centering cancels in case-minus-control differences. Subject/draw scaling
  # is constant within a trial, so divide the reduced values once.
  if (!is.null(scaling)) result <- result / scaling$scale[scaling_ids, , drop = FALSE]
  out[colnames(values)] <- as.data.frame(result)
  attr(out, "score_columns") <- colnames(values)
  out
}

#' Evaluate replayable mismatched-template draws
#' @inheritParams mismatch_assignments
#' @param assignments Saved output of [mismatch_assignments()]. If NULL, create
#'   `n` draws using `seed`; supplied assignments are evaluated without RNG use.
#' @param replace_ref Replace only self_item or the entire reference tuple.
#' @param restandardize Refit scaling in every null draw. FALSE applies the
#'   original observed fit. Each paired matched/mismatched comparison uses the
#'   same fit (estimated on retained mismatched scores when TRUE). This choice
#'   alone does not establish an exact randomization test: calibration also
#'   requires an appropriate exchangeability assumption.
#' @param standardize Score standardization, as in [score_strata()].
#' @param by Contrast resolution, as in [case_control_contrast()].
#' @param control_policy By default, `"exclude_replacement"` removes observed
#'   controls donated by any trial of the replacement item, using the saved
#'   donor metadata. Cases are always retained. `"fixed"` retains all controls:
#'   this item-discrimination baseline can be negative under a genuine matched
#'   effect because it includes the replacement item's own test fixations.
#' @return A list with `observed` contrasts on the original controls,
#'   `paired_observed` matched contrasts on each draw's retained controls and
#'   scaling, `contrasts` for mismatched templates, assignments, original
#'   observed scaling, and settings. Both paired tables carry `null_draw` and
#'   identical row keys. Compare `paired_observed` with `contrasts` to use the
#'   same controls and units. Evaluation proceeds one draw at a time.
#' @details Exclusion uses donor items, including repeated trials of that item,
#'   rather than assuming trial labels equal item labels. It does not redraw
#'   controls. A stratum with no remaining controls is an error. Generated
#'   uniform/map controls have no donor item and are retained. For older saved
#'   strata without donor metadata, complete case rows can supply the donor
#'   keys; otherwise resample and save the strata or choose the fixed policy.
#'   No automatic p-values or calibrated effect decomposition are computed.
#' @md
#' @export
mismatch_draws <- function(strata, refs, n = 100, trial = "trial", seed = 1,
                           assignments = NULL, replace_ref = c("self_item", "all"),
                           restandardize = TRUE, standardize = c("within_subject", "none"),
                           by = c("fixation", "trial"),
                           control_policy = c("exclude_replacement", "fixed")) {
  replace_ref <- match.arg(replace_ref)
  standardize <- match.arg(standardize)
  by <- match.arg(by)
  control_policy <- match.arg(control_policy)
  if (!is.logical(restandardize) || length(restandardize) != 1L || is.na(restandardize)) stop("restandardize must be TRUE or FALSE.")
  if (replace_ref == "self_item" && !"self_item" %in% refs$refs) stop("self_item is not in this reference set.")
  subject <- refs$spec$subject
  item <- refs$spec$item
  raw <- score_strata(strata, refs, standardize = "none")
  layout <- contrast_layout(raw, check_strata(strata), by, subject, trial)
  donor_item <- if (control_policy == "exclude_replacement") control_item_ids(strata, refs) else NULL
  idx <- density_grid_index(strata, refs$x, refs$y)
  cell <- idx[, 1] + (idx[, 2] - 1L) * length(refs$x)
  values <- as.matrix(raw[refs$refs])
  scale_layout <- if (standardize == "none") NULL else scaling_layout(raw, c(subject, "draw"))
  scaling <- if (is.null(scale_layout)) NULL else fit_scaling_values(values, scale_layout)
  scaling_ids <- scale_layout$ids[layout$output_rows]
  observed <- contrast_values(values, layout, scaling = scaling, scaling_ids = scaling_ids)
  if (is.null(assignments)) assignments <- mismatch_assignments(strata, refs, n, trial, seed)
  keys <- unique(c(subject, trial))
  trials <- unique(as.data.frame(strata[unique(c(keys, item))]))
  trial_row <- match_map_rows(strata, trials, keys, strict = TRUE)
  check_key_columns(assignments, unique(c(keys, item, "null_draw", "replacement_item")))
  draw_ids <- unique(assignments$null_draw)
  if (!length(draw_ids)) stop("assignments must contain at least one null draw.")
  out <- lapply(draw_ids, function(d) {
    assn <- assignments[assignments$null_draw == d, , drop = FALSE]
    ai <- match_map_rows(trials, assn, keys, strict = TRUE)
    if (anyNA(ai) || length(unique(ai)) != nrow(assn)) stop("Assignments must cover exactly the scored trials.")
    if (any(!equal_key_values(assn[[item]][ai], trials[[item]])) ||
        any(equal_key_values(assn$replacement_item, assn[[item]]))) {
      stop("Assignments must preserve original items and replace them with different items.")
    }
    replacement_keys <- trials[names(refs$keys)]
    replacement_keys[[item]] <- assn$replacement_item[ai]
    target <- match_map_rows(replacement_keys, refs$keys, names(refs$keys), strict = TRUE)
    if (anyNA(target) || any(!refs$valid[target])) stop("Replacement items must have valid maps for the same subject.")
    target <- target[trial_row]
    keep <- if (is.null(donor_item)) NULL else strata$is_case | donor_item != refs$ii[target]
    scored <- values
    columns <- if (replace_ref == "all") refs$refs else "self_item"
    scored[, columns] <- as.matrix(evaluate_references(refs, target, cell, columns))
    fit <- if (!is.null(scaling) && restandardize) fit_scaling_values(scored, scale_layout, keep) else scaling
    contrasts <- contrast_values(scored, layout, keep, fit, scaling_ids)
    matched <- contrast_values(values, layout, keep, fit, scaling_ids)
    contrasts$null_draw <- d
    matched$null_draw <- d
    list(matched = matched, contrasts = contrasts)
  })
  list(observed = observed,
       paired_observed = dplyr::bind_rows(lapply(out, `[[`, "matched")),
       contrasts = dplyr::bind_rows(lapply(out, `[[`, "contrasts")), assignments = assignments, scaling = scaling,
       settings = list(replace_ref = replace_ref, restandardize = restandardize,
                       standardize = standardize, by = by, control_policy = control_policy))
}

# Resolve donor identity independently of the case metadata carried by every
# assignment row. Stored source rows survive case subsetting and RDS replay.
control_item_ids <- function(strata, refs) {
  if (!"control_id" %in% names(strata)) stop("Replacement-item exclusion requires control_id donor provenance.")
  donors <- strata$control_id
  controls <- !strata$is_case & !is.na(donors)
  if (!is.numeric(donors) || any(!is.finite(donors[controls])) ||
      any(donors[controls] != trunc(donors[controls]))) stop("control_id must index saved donor rows.")
  if (identical(attr(strata, "sampling")$method, "observed") && anyNA(donors[!strata$is_case])) {
    stop("Observed controls require donor provenance.")
  }
  out <- integer(nrow(strata))
  if (!any(controls)) {
    return(out)
  }
  keys <- names(refs$keys)
  sources <- attr(strata, "control_sources")
  if (is.null(sources)) {
    sources <- unique(as.data.frame(strata[strata$is_case, c("case_id", keys), drop = FALSE]))
    if (anyDuplicated(sources$case_id)) stop("Donor case_id has inconsistent reference keys.")
    rows <- match(donors[controls], sources$case_id)
  } else {
    rows <- donors[controls]
    rows[rows < 1L | rows > nrow(sources)] <- NA_integer_
  }
  if (anyNA(rows)) stop("Missing donor provenance for replacement-item exclusion; retain control_sources or all donor case rows.")
  check_key_columns(sources, keys)
  if (any(!equal_key_values(sources[[refs$spec$subject]][rows], strata[[refs$spec$subject]][controls]))) {
    stop("Observed donor subjects must match case subjects.")
  }
  out[controls] <- match(sources[[refs$spec$item]][rows], unique(refs$keys[[refs$spec$item]]), nomatch = 0L)
  out
}

check_score_columns <- function(scores, columns) {
  if (!is.data.frame(scores) || !is.character(columns) || !length(columns) ||
      anyDuplicated(columns) || !all(columns %in% names(scores)) ||
      !all(vapply(scores[columns], is.numeric, logical(1))) ||
      any(vapply(scores[columns], function(x) any(is.infinite(x)), logical(1)))) {
    stop("columns must identify distinct numeric score columns without infinite values.")
  }
}

check_strata <- function(strata) {
  check_key_columns(strata, c("case_id", "stratum", "draw"))
  if (!nrow(strata) || !is.logical(strata$is_case) || anyNA(strata$is_case)) {
    stop("Strata require rows and a non-missing logical is_case column.")
  }
  groups <- split(seq_len(nrow(strata)), key_group_ids(strata, c("draw", "stratum")))
  if (any(vapply(groups, function(g) length(g) < 2L || sum(strata$is_case[g]) != 1L ||
                 length(unique(strata$case_id[g])) != 1L, logical(1)))) {
    stop("Each stratum requires exactly one case, at least one control, and one case_id.")
  }
  if (anyDuplicated(strata[strata$is_case, c("draw", "case_id")])) stop("Each case_id must have one stratum per draw.")
  groups
}

check_stratum_keys <- function(strata, groups, keys) {
  check_key_columns(strata, keys)
  ids <- key_group_ids(strata, keys)
  if (any(vapply(groups, function(g) length(unique(ids[g])) != 1L, logical(1)))) {
    stop("All locations in a stratum must carry the case's reference keys.")
  }
}

#' @export
#' @noRd
print.eyesim_score_scaling <- function(x, ...) {
  cat("Score scaling:", nrow(x$groups), "groups;", length(x$columns), "reference columns\n")
  cat("Groups:", if (length(x$by)) paste(x$by, collapse = ", ") else "global", "\n")
  cat("References:", paste(x$columns, collapse = ", "), "\n")
  invisible(x)
}
