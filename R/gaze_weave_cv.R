# Shared cross-validation utilities ----------------------------------------

resolve_gaze_weave_filter <- function(tab, filter_spec, label) {
  if (is.null(filter_spec)) {
    return(rep(TRUE, nrow(tab)))
  }
  values <- if (is.function(filter_spec)) filter_spec(tab) else filter_spec
  if (!is.logical(values) || length(values) != nrow(tab) || anyNA(values)) {
    stop(
      label,
      " must be NULL, a complete logical vector, or a function returning one."
    )
  }
  values
}

# Seed GazeWeave's random draws reproducibly: the generator kinds are fixed
# (R's defaults) whatever the caller has set, and when the calling function
# exits the caller's generator kinds and random number stream are restored.
local_gaze_weave_seed <- function(seed, envir = parent.frame()) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("seed must be a single number.")
  }
  restore_stream <- snapshot_session_rng()
  kinds <- RNGkind()
  restore <- function() {
    suppressWarnings(RNGkind(kinds[[1L]], kinds[[2L]], kinds[[3L]]))
    restore_stream()
  }
  do.call(on.exit, list(as.call(list(restore)), add = TRUE), envir = envir)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  invisible(NULL)
}

make_gaze_weave_folds <- function(source_tab, split_on, contrast_on = NULL,
                                  n_folds = NULL, seed = 1) {
  split_key <- gaze_key(source_tab, split_on, "split_on")
  stratum_key <- if (is.null(contrast_on)) {
    rep("all", nrow(source_tab))
  } else {
    gaze_key(source_tab, contrast_on, "contrast_on")
  }

  group_map <- unique(data.frame(
    split_key = split_key,
    stratum_key = stratum_key,
    stringsAsFactors = FALSE
  ))
  strata_per_group <- table(group_map$split_key)
  if (any(strata_per_group > 1L)) {
    stop("Each split_on group must belong to exactly one contrast_on stratum.")
  }
  groups_per_stratum <- table(group_map$stratum_key)
  minimum_groups <- min(groups_per_stratum)
  if (minimum_groups < 2L) {
    stop(
      "Cross-fitted GazeWeave requires at least two split groups per ",
      "contrast stratum."
    )
  }
  if (is.null(n_folds)) {
    n_folds <- min(5L, minimum_groups)
  }
  n_folds <- as.integer(n_folds)
  if (length(n_folds) != 1L || is.na(n_folds) || n_folds < 2L ||
      n_folds > minimum_groups) {
    stop(
      "n_folds must be between two and the smallest number of split groups ",
      "in a contrast stratum."
    )
  }

  # Fold assignment is a random draw. Run it under `seed` with a fixed
  # generator and locale-independent ordering, so the same seed gives the same
  # folds everywhere, and then restore the caller's random number stream.
  local_gaze_weave_seed(seed)
  group_map$fold <- NA_integer_
  for (stratum in sort(unique(group_map$stratum_key), method = "radix")) {
    rows <- which(group_map$stratum_key == stratum)
    group_names <- sort(group_map$split_key[rows], method = "radix")
    shuffled <- sample(group_names, length(group_names))
    fold_for_group <- stats::setNames(
      rep(seq_len(n_folds), length.out = length(shuffled)),
      shuffled
    )
    group_map$fold[rows] <- unname(
      fold_for_group[group_map$split_key[rows]]
    )
  }
  fold_lookup <- stats::setNames(group_map$fold, group_map$split_key)

  list(
    fold_id = unname(fold_lookup[split_key]),
    split_key = split_key,
    stratum_key = stratum_key,
    n_folds = n_folds,
    group_map = group_map
  )
}

# Folds for methods that score each held-out row only against the candidates
# in its own fold (Replay and the baseline court): every fold x contrast
# stratum holding evaluation rows needs the true candidate and at least one
# nonmatch, i.e. two distinct match_on keys. With `n_folds = NULL`, the
# default is lowered until that holds; an explicit `n_folds` that violates it
# is an error naming the problem.
make_gaze_weave_candidate_folds <- function(source_tab, split_on, contrast_on,
                                            n_folds, seed, match_on,
                                            eval_mask) {
  match_key <- gaze_key(source_tab, match_on, "match_on")[eval_mask]
  stratum_key <- if (is.null(contrast_on)) {
    rep("all", sum(eval_mask))
  } else {
    gaze_key(source_tab, contrast_on, "contrast_on")[eval_mask]
  }
  enough_candidates <- function(folds) {
    if (length(match_key) == 0L) return(TRUE)
    cell <- paste(folds$fold_id[eval_mask], stratum_key, sep = "\r")
    counts <- tapply(match_key, cell, function(key) length(unique(key)))
    all(counts >= 2L)
  }
  folds <- make_gaze_weave_folds(source_tab, split_on, contrast_on, n_folds, seed)
  if (is.null(n_folds)) {
    while (!enough_candidates(folds) && folds$n_folds > 2L) {
      folds <- make_gaze_weave_folds(
        source_tab, split_on, contrast_on, folds$n_folds - 1L, seed
      )
    }
  }
  if (!enough_candidates(folds)) {
    stop(
      "With ", folds$n_folds, " folds, some fold has fewer than two ",
      "candidate items (the true match and a nonmatch) in a contrast ",
      "stratum. Use fewer folds or larger contrast strata."
    )
  }
  folds
}

#' Cross-fitted GazeWeave analysis
#'
#' `gaze_weave_cv()` is the common entry point for the two GazeWeave engines.
#' Use Transport for symmetric explanatory alignment and Replay for a
#' directional encoding-to-recall model. There is no universal engine default;
#' the supplied engine-specific specification selects the scientific model.
#'
#' Both engines return held-out candidate probabilities and
#' `gaze_info_bits = log2(p_true / prior_true)`. Replay correspondence is a
#' posterior probability. Transport correspondence is an optimized alignment
#' and must not be interpreted as statistical uncertainty.
#'
#' @param ref_tab,source_tab Reference and source tables.
#' @param match_on One or more columns defining the true template match.
#' @param contrast_on Optional columns restricting nonmatching candidates.
#' @param refvar,sourcevar List columns containing fixation groups.
#' @param spec A [gaze_transport_spec()] or [gaze_replay_spec()].
#' @param engine Either `"transport"` or `"replay"`. When omitted, it is
#'   inferred from `spec`.
#' @param split_on Columns defining the held-out unit. Defaults to `match_on`.
#' @param n_folds Number of cross-fitting folds. `NULL` uses up to five; for
#'   Replay, which scores each held-out row against the candidates in its own
#'   fold, the default is lowered when needed so that every fold keeps a true
#'   candidate and a nonmatch in each contrast stratum.
#'   Calibration needs more: each training fold must hold at least
#'   `2 * calibration folds` items (4 by default) per contrast stratum, or that
#'   fold falls back to an uncalibrated score and a warning of class
#'   `gaze_calibration_fallback` is raised.
#' @param seed Fold-assignment seed. `NULL` uses the engine-specific default.
#'   Folds depend on `seed` alone: the draw uses R's default generator and
#'   locale-independent ordering whatever `RNGkind()` and locale are in effect,
#'   and the caller's generator and random number stream are left unchanged.
#' @param fit_source_filter,eval_source_filter Optional logical vectors or
#'   functions selecting fitting and evaluation source rows.
#' @param episode_on Optional reference columns identifying separate study
#'   presentations for Transport. Their likelihoods receive equal prior weight.
#' @param priorvar Optional Transport reference column containing one positive
#'   design-prior weight per candidate.
#' @param template_on Optional reference columns identifying separate encoding
#'   templates of one candidate for Replay (the counterpart of `episode_on`);
#'   see [gaze_replay_cv()].
#'
#' @return An engine-specific fitted object. Use `broom::tidy()` for the sole
#'   primary endpoint, `gaze_info_bits`.
#' @section Experimental:
#' GazeWeave is experimental. Its interface, defaults and result columns may
#' change in future releases; record the package version with any results.
#' @export
gaze_weave_cv <- function(
    ref_tab, source_tab, match_on, contrast_on = NULL,
    refvar = "fixgroup", sourcevar = "fixgroup",
    spec = NULL, engine = NULL, split_on = match_on,
    n_folds = NULL, seed = NULL,
    fit_source_filter = NULL, eval_source_filter = NULL,
    episode_on = NULL, priorvar = NULL, template_on = NULL) {
  inferred <- if (inherits(spec, "gaze_replay_spec")) {
    "replay"
  } else if (inherits(spec, "gaze_transport_spec")) {
    "transport"
  } else {
    NULL
  }
  if (is.null(engine)) {
    if (is.null(inferred)) {
      stop(
        "No GazeWeave engine is selected. Supply gaze_transport_spec() for ",
        "symmetric alignment or gaze_replay_spec() for directional replay."
      )
    }
    engine <- inferred
  }
  engine <- match.arg(engine, c("transport", "replay"))
  if (!is.null(inferred) && !identical(engine, inferred)) {
    stop(
      "engine = '", engine, "' is incompatible with the supplied ",
      class(spec)[[1L]], "."
    )
  }

  if (identical(engine, "replay")) {
    if (is.null(seed)) seed <- 1L
    if (is.null(spec)) spec <- gaze_replay_spec()
    if (!is.null(episode_on) || !is.null(priorvar)) {
      stop("episode_on and priorvar are available only for Transport; ",
           "Replay takes template_on.")
    }
    return(gaze_replay_cv(
      ref_tab = ref_tab,
      source_tab = source_tab,
      match_on = match_on,
      contrast_on = contrast_on,
      refvar = refvar,
      sourcevar = sourcevar,
      spec = spec,
      split_on = split_on,
      n_folds = n_folds,
      seed = seed,
      fit_source_filter = fit_source_filter,
      eval_source_filter = eval_source_filter,
      template_on = template_on
    ))
  }

  if (!is.null(template_on)) {
    stop("template_on is available only for Replay; Transport takes episode_on.")
  }
  if (is.null(seed)) seed <- 20260822L
  if (is.null(spec)) spec <- gaze_transport_spec()
  gaze_transport_cv(
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = match_on,
    contrast_on = contrast_on,
    refvar = refvar,
    sourcevar = sourcevar,
    spec = spec,
    split_on = split_on,
    n_folds = n_folds,
    seed = seed,
    fit_source_filter = fit_source_filter,
    eval_source_filter = eval_source_filter,
    episode_on = episode_on,
    priorvar = priorvar
  )
}
