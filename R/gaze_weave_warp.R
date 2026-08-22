# Cross-fitted GazeWeave registration ---------------------------------------

gaze_key <- function(tab, columns, label = "key") {
  if (!is.character(columns) || length(columns) == 0L ||
      !all(columns %in% names(tab))) {
    stop(label, " must name one or more columns in the table.")
  }
  if (anyNA(tab[columns])) {
    stop(label, " columns cannot contain missing values.")
  }
  encoded <- lapply(tab[columns], function(column) {
    value <- enc2utf8(as.character(column))
    paste0(nchar(value, type = "bytes"), ":", value)
  })
  do.call(paste, c(encoded, sep = "|"))
}

identity_gaze_warp_model <- function() {
  structure(
    list(
      type = "none",
      fit_by = NULL,
      group_models = list(all = list(
        A = diag(2L),
        translation = c(0, 0),
        center = c(0, 0),
        shift = c(0, 0),
        scale = 1
      )),
      info = list(type = "none", groups = "all")
    ),
    class = c("gaze_warp_model", "list")
  )
}

resolve_gaze_warp_center <- function(warp_spec, screen) {
  if (is.character(warp_spec$center)) {
    if (is.null(screen)) {
      stop("A screen-centred warp requires screen geometry.")
    }
    return(screen$center)
  }
  as.numeric(warp_spec$center)
}

fit_gaze_warp_model <- function(ref_tab, source_tab, match_on, refvar, sourcevar,
                                spec) {
  warp_spec <- spec$warp
  if (identical(warp_spec$type, "none")) {
    return(identity_gaze_warp_model())
  }
  if (!identical(warp_spec$type, "contraction")) {
    stop("Unsupported GazeWeave warp: ", warp_spec$type, ".")
  }
  if (!refvar %in% names(ref_tab) || !sourcevar %in% names(source_tab)) {
    stop("refvar and sourcevar must identify fixation-group columns.")
  }

  ref_match_key <- gaze_key(ref_tab, match_on, "match_on")
  source_match_key <- gaze_key(source_tab, match_on, "match_on")
  if (anyDuplicated(ref_match_key)) {
    stop("Reference match_on keys must be unique for GazeWeave.")
  }
  ref_index <- match(source_match_key, ref_match_key)
  if (anyNA(ref_index)) {
    stop("Every source path used for warp fitting must have a reference match.")
  }

  fit_by <- warp_spec$fit_by
  if (is.null(fit_by)) {
    source_group <- rep("all", nrow(source_tab))
  } else {
    if (!all(fit_by %in% names(ref_tab)) || !all(fit_by %in% names(source_tab))) {
      stop("warp fit_by columns must exist in both reference and source tables.")
    }
    source_group <- gaze_key(source_tab, fit_by, "warp fit_by")
    expected_ref_group <- gaze_key(ref_tab[ref_index, , drop = FALSE], fit_by, "warp fit_by")
    if (any(source_group != expected_ref_group)) {
      stop("Matched reference and source paths disagree on warp fit_by values.")
    }
  }

  center <- resolve_gaze_warp_center(warp_spec, spec$screen)
  group_levels <- sort(unique(source_group))
  group_models <- list()
  group_info <- vector("list", length(group_levels))

  for (g in seq_along(group_levels)) {
    group_key <- group_levels[[g]]
    source_rows <- which(source_group == group_key)
    paired_ref_rows <- ref_index[source_rows]
    ref_measures <- lapply(ref_tab[[refvar]][paired_ref_rows], as_gaze_measure,
                           chronology = spec$chronology)
    source_measures <- lapply(source_tab[[sourcevar]][source_rows], as_gaze_measure,
                              chronology = spec$chronology)
    ref_moments <- gaze_measure_moments(ref_measures)
    source_moments <- gaze_measure_moments(source_measures)

    ref_radius <- sum(diag(ref_moments$cov))
    source_radius <- sum(diag(source_moments$cov))
    scale <- sqrt((ref_radius + warp_spec$shrink) /
                    (source_radius + warp_spec$shrink))
    A <- diag(scale, 2L)

    if (warp_spec$translation) {
      translation <- as.numeric(ref_moments$mean - A %*% source_moments$mean)
      shift <- as.numeric(ref_moments$mean - center -
                            scale * (source_moments$mean - center))
    } else {
      shift <- c(0, 0)
      translation <- as.numeric(center - A %*% center)
    }

    group_models[[group_key]] <- list(
      A = A,
      translation = translation,
      center = center,
      shift = shift,
      scale = scale,
      matched_pairs = length(source_rows)
    )
    group_info[[g]] <- list(
      group = group_key,
      matched_pairs = length(source_rows),
      scale = scale,
      translation = translation,
      shift = shift,
      reference_radius = ref_radius,
      source_radius = source_radius
    )
  }

  structure(
    list(
      type = "contraction",
      fit_by = fit_by,
      group_models = group_models,
      info = list(
        type = "contraction",
        fit_by = fit_by,
        center = center,
        translation = warp_spec$translation,
        shrink = warp_spec$shrink,
        groups = group_info
      )
    ),
    class = c("gaze_warp_model", "list")
  )
}

subset_gaze_warp_model <- function(model, group_key = NULL) {
  if (identical(model$type, "none")) {
    return(model)
  }
  if (is.null(group_key)) {
    if (length(model$group_models) != 1L) {
      stop("A grouped warp model requires a group key.")
    }
    group_key <- names(model$group_models)[[1]]
  }
  group_model <- model$group_models[[group_key]]
  if (is.null(group_model)) {
    stop("No fitted GazeWeave warp is available for group '", group_key, "'.")
  }
  structure(
    list(
      type = model$type,
      fit_by = model$fit_by,
      group_key = group_key,
      group_models = stats::setNames(list(group_model), group_key),
      info = model$info
    ),
    class = c("gaze_warp_model", "list")
  )
}

warp_parameters <- function(model) {
  group_model <- model$group_models[[1]]
  list(
    A = group_model$A,
    translation = group_model$translation,
    center = group_model$center,
    shift = group_model$shift,
    scale = group_model$scale
  )
}

apply_gaze_warp_model <- function(measure, model) {
  if (!inherits(model, "gaze_warp_model")) {
    stop("warp_model must be a fitted gaze_warp_model.")
  }
  if (length(model$group_models) != 1L) {
    stop("apply_gaze_warp_model() requires a single-group warp model.")
  }
  parameters <- warp_parameters(model)
  transform_gaze_measure(
    measure,
    A = parameters$A,
    translation = parameters$translation
  )
}
