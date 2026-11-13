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

# Groups are kept in order of first appearance, as unique() would return them.
scaling_layout <- function(scores, by) {
  if (!length(by)) return(list(groups = data.frame(.group = 1L), ids = rep(1L, nrow(scores)), by = by))
  key <- key_group_ids(scores, by)
  first <- which(!duplicated(key))
  list(groups = as.data.frame(scores[by])[first, , drop = FALSE],
       ids = match(key, key[first]), by = by)
}

# Group counts, sums and sums of squares of values shifted by `shift` (one row
# per group). Shifting by a value near the group mean keeps the one-pass
# variance accurate, and the sums can be updated by removing rows.
shifted_moments <- function(x, ids, n, shift) {
  d <- x - shift[ids, , drop = FALSE]
  ok <- !is.na(d)
  d[!ok] <- 0
  group_sums(cbind(ok + 0, d, d^2), ids, n)
}

# The scaling estimate of fit_scaling_values() from shifted moments.
moment_scale <- function(moments, columns) {
  k <- length(columns)
  counts <- moments[, seq_len(k), drop = FALSE]
  s1 <- moments[, k + seq_len(k), drop = FALSE]
  s2 <- moments[, 2L * k + seq_len(k), drop = FALSE]
  sd <- sqrt(pmax(s2 - s1^2 / counts, 0) / (counts - 1))
  scale <- sd
  scale[!(counts > 1 & is.finite(sd) & sd > 0)] <- 1
  scale[counts == 0] <- NA_real_
  dimnames(scale) <- list(NULL, columns)
  scale
}

# Reuse group membership for replay draws instead of reconstructing key tables.
# Vectorized over groups and columns: two-pass mean (with one correction pass,
# as mean() does) and sample SD of the non-missing values in each group.
fit_scaling_values <- function(values, layout, keep = NULL) {
  ng <- nrow(layout$groups)
  rows <- if (is.null(keep)) seq_len(nrow(values)) else which(keep)
  ids <- layout$ids[rows]
  x <- values[rows, , drop = FALSE]
  ok <- !is.na(x)
  x[!ok] <- 0
  counts <- group_sums(ok + 0, ids, ng)
  center <- group_sums(x, ids, ng) / counts
  deviation <- function() {
    d <- x - center[ids, , drop = FALSE]
    d[!ok] <- 0
    d
  }
  center <- center + group_sums(deviation(), ids, ng) / counts
  sd <- sqrt(group_sums(deviation()^2, ids, ng) / (counts - 1))
  scale <- sd
  scale[!(counts > 1 & is.finite(sd) & sd > 0)] <- 1
  center[counts == 0] <- NA_real_
  scale[counts == 0] <- NA_real_
  structure(list(groups = layout$groups, by = layout$by, columns = colnames(values),
                 center = center, scale = scale, n = counts), class = "eyesim_score_scaling")
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
  ids <- layout$ids[controls]
  n <- length(layout$case_rows)
  # Subtract before summing: a large common offset must not erase a tiny
  # contrast that becomes appreciable after dividing by a small score SD.
  differences <- values[layout$case_rows[ids], , drop = FALSE] - values[controls, , drop = FALSE]
  finish_contrast(group_sums(differences, ids, n), tabulate(ids, n), layout, scaling, scaling_ids)
}

check_control_counts <- function(counts) {
  if (any(counts == 0L)) {
    stop("Item exclusion leaves a stratum without controls. Donors must include items other than both the case and replacement items. Increase n only if such donors exist, or explicitly choose control_policy = 'fixed' for a different baseline.")
  }
}

# Column sums by group, with a zero row for every group absent from `ids`.
group_sums <- function(x, ids, n) {
  out <- matrix(0, n, ncol(x), dimnames = list(NULL, colnames(x)))
  if (length(ids)) {
    s <- rowsum(x, ids, reorder = TRUE)
    out[as.integer(rownames(s)), ] <- s
  }
  out
}

# Stratum sums of case-minus-control differences become (scaled) contrasts.
finish_contrast <- function(sums, counts, layout, scaling = NULL, scaling_ids = NULL) {
  check_control_counts(counts)
  result <- sums / counts
  out <- layout$template
  if (!is.null(layout$trial_ids)) {
    result <- group_sums(result, layout$trial_ids, length(layout$trial_counts)) / layout$trial_counts
  } else {
    out$n_controls <- counts
  }
  # Centering cancels in case-minus-control differences. Subject/draw scaling
  # is constant within a trial, so divide the reduced values once.
  if (!is.null(scaling)) result <- result / scaling$scale[scaling_ids, , drop = FALSE]
  out[colnames(sums)] <- as.data.frame(result)
  attr(out, "score_columns") <- colnames(sums)
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
#'   controls donated by any trial of either the case's item or the replacement
#'   item, using the saved donor metadata. Cases are always retained.
#'   `"fixed"` retains all controls:
#'   this item-discrimination baseline can be negative under a genuine matched
#'   effect because it includes the replacement item's own test fixations.
#' @return A list with `observed` contrasts on the original controls,
#'   `paired_observed` matched contrasts on each draw's retained controls and
#'   scaling, `contrasts` for mismatched templates, assignments, original
#'   observed scaling, and settings. Both paired tables carry `null_draw` and
#'   identical row keys. Compare `paired_observed` with `contrasts` to use the
#'   same controls and units. Evaluation proceeds one draw at a time.
#' @details Exclusion uses donor items, including repeated trials of either item,
#'   rather than assuming trial labels equal item labels. It does not redraw
#'   controls. A stratum with no remaining controls is an error: a donor pool
#'   containing only the case and replacement items cannot support the paired
#'   comparison, regardless of the number of controls drawn. In particular,
#'   two test trials fail when they contain those two items. Generated
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
  keep_own <- if (is.null(donor_item)) NULL else
    strata$is_case | donor_item != match(strata[[item]], unique(refs$keys[[item]]))
  idx <- density_grid_index(strata, refs$x, refs$y)
  cell <- idx[, 1] + (idx[, 2] - 1L) * length(refs$x)
  values <- as.matrix(raw[refs$refs])
  scale_layout <- if (standardize == "none") NULL else scaling_layout(raw, c(subject, "draw"))
  scaling <- if (is.null(scale_layout)) NULL else fit_scaling_values(values, scale_layout)
  scaling_ids <- scale_layout$ids[layout$output_rows]
  observed <- contrast_values(values, layout, scaling = scaling, scaling_ids = scaling_ids)
  if (is.null(assignments)) assignments <- mismatch_assignments(strata, refs, n, trial, seed)
  keys <- unique(c(subject, trial))
  # Every row carries its case's keys, so case rows list the trials in order.
  trials <- unique(as.data.frame(strata[strata$is_case, unique(c(keys, item)), drop = FALSE]))
  trial_row <- match_map_rows(strata, trials, keys, strict = TRUE)
  check_key_columns(assignments, unique(c(keys, item, "null_draw", "replacement_item")))
  draw_ids <- unique(assignments$null_draw)
  if (!length(draw_ids)) stop("assignments must contain at least one null draw.")
  # Validate and resolve every draw's replacement maps in one pass.
  draw_of <- match(assignments$null_draw, draw_ids)
  trial_of <- match_map_rows(assignments, trials, keys)
  if (anyNA(trial_of) || anyDuplicated(data.frame(draw_of, trial_of)) ||
      any(tabulate(draw_of, length(draw_ids)) != nrow(trials))) {
    stop("Assignments must cover exactly the scored trials.")
  }
  if (any(!equal_key_values(assignments[[item]], trials[[item]][trial_of])) ||
      any(equal_key_values(assignments$replacement_item, assignments[[item]]))) {
    stop("Assignments must preserve original items and replace them with different items.")
  }
  replacement_keys <- trials[trial_of, names(refs$keys), drop = FALSE]
  replacement_keys[[item]] <- assignments$replacement_item
  replacement <- match_map_rows(replacement_keys, refs$keys, names(refs$keys), strict = TRUE)
  if (anyNA(replacement) || any(!refs$valid[replacement])) stop("Replacement items must have valid maps for the same subject.")
  draw_rows <- split(seq_len(nrow(assignments)), draw_of)
  columns <- if (replace_ref == "all") refs$refs else "self_item"
  # Matched differences do not change across draws; only the retained controls
  # do. Reduce them once, then remove each draw's few excluded donors.
  controls <- layout$controls
  if (!is.null(keep_own)) controls <- controls[keep_own[controls]]
  stratum <- layout$ids[controls]
  n_strata <- length(layout$case_rows)
  matched_diff <- values[layout$case_rows[stratum], , drop = FALSE] - values[controls, , drop = FALSE]
  matched_total <- group_sums(matched_diff, stratum, n_strata)
  count_total <- tabulate(stratum, n_strata)
  refit <- !is.null(scaling) && restandardize
  if (refit) {
    # Unchanged columns: moments over cases and base controls, from which each
    # draw removes its dropped donors. Replaced columns are refitted per draw.
    fixed_columns <- setdiff(refs$refs, columns)
    shift <- scaling$center
    shift[is.na(shift)] <- 0
    case_idx <- which(strata$is_case)
    case_ids <- scale_layout$ids[case_idx]
    stratum_group <- scale_layout$ids[layout$case_rows]
    base_rows <- c(case_idx, controls)
    fixed_moments <- shifted_moments(values[base_rows, fixed_columns, drop = FALSE],
                                     scale_layout$ids[base_rows], nrow(shift),
                                     shift[, fixed_columns, drop = FALSE])
  }
  out <- lapply(seq_along(draw_ids), function(k) {
    rows <- draw_rows[[k]]
    by_trial <- integer(nrow(trials))
    by_trial[trial_of[rows]] <- replacement[rows]
    target <- by_trial[trial_row]
    dropped <- if (is.null(donor_item)) integer(0) else
      which(donor_item[controls] == refs$ii[target[controls]])
    counts <- count_total - tabulate(stratum[dropped], n_strata)
    check_control_counts(counts)
    matched_sums <- matched_total -
      group_sums(matched_diff[dropped, , drop = FALSE], stratum[dropped], n_strata)
    kept <- if (length(dropped)) controls[-dropped] else controls
    kept_stratum <- if (length(dropped)) stratum[-dropped] else stratum
    rescored <- matrix(vapply(columns, function(ref) reference_values(refs, target, cell, ref),
                              numeric(length(cell))), ncol = length(columns),
                       dimnames = list(NULL, columns))
    kept_diff <- rescored[layout$case_rows[kept_stratum], , drop = FALSE] - rescored[kept, , drop = FALSE]
    mismatch_sums <- matched_sums
    fit <- scaling
    if (!refit) {
      mismatch_sums[, columns] <- group_sums(kept_diff, kept_stratum, n_strata)
    } else {
      scale <- matrix(NA_real_, nrow(shift), length(refs$refs), dimnames = list(NULL, refs$refs))
      if (length(fixed_columns)) {
        lost <- controls[dropped]
        moments <- fixed_moments - shifted_moments(
          values[lost, fixed_columns, drop = FALSE], scale_layout$ids[lost], nrow(shift),
          shift[, fixed_columns, drop = FALSE])
        scale[, fixed_columns] <- moment_scale(moments, fixed_columns)
      }
      # Shift replaced columns by each group's mean rescored case score, which
      # lies near the group mean; the observed centre can be missing or far off,
      # and the one-pass variance needs a shift close to the mean.
      case_values <- rescored[case_idx, , drop = FALSE]
      case_present <- !is.na(case_values)
      case_values[!case_present] <- 0
      replaced_shift <- group_sums(case_values, case_ids, nrow(shift)) /
        group_sums(case_present + 0, case_ids, nrow(shift))
      replaced_shift[!is.finite(replaced_shift)] <- 0
      # One pass over the kept controls yields both the contrast sums and the
      # replaced columns' moments by stratum; strata nest in scaling groups.
      shifted <- rescored[kept, , drop = FALSE] - replaced_shift[scale_layout$ids[kept], , drop = FALSE]
      present <- !is.na(shifted)
      shifted[!present] <- 0
      nc <- length(columns)
      sums <- group_sums(cbind(kept_diff, present + 0, shifted, shifted^2), kept_stratum, n_strata)
      mismatch_sums[, columns] <- sums[, seq_len(nc), drop = FALSE]
      moments <- group_sums(sums[, -seq_len(nc), drop = FALSE], stratum_group, nrow(shift)) +
        shifted_moments(rescored[case_idx, , drop = FALSE], case_ids, nrow(shift), replaced_shift)
      scale[, columns] <- moment_scale(moments, columns)
      fit <- list(scale = scale)
    }
    contrasts <- finish_contrast(mismatch_sums, counts, layout, fit, scaling_ids)
    matched <- finish_contrast(matched_sums, counts, layout, fit, scaling_ids)
    contrasts$null_draw <- draw_ids[k]
    matched$null_draw <- draw_ids[k]
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
