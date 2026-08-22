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

  set.seed(seed)
  group_map$fold <- NA_integer_
  for (stratum in sort(unique(group_map$stratum_key))) {
    rows <- which(group_map$stratum_key == stratum)
    group_names <- sort(group_map$split_key[rows])
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
#' @param n_folds Number of cross-fitting folds.
#' @param seed Fold-assignment seed. `NULL` uses the engine-specific default.
#' @param fit_source_filter,eval_source_filter Optional logical vectors or
#'   functions selecting fitting and evaluation source rows.
#' @param episode_on Optional reference columns identifying separate study
#'   presentations for Transport. Their likelihoods receive equal prior weight.
#' @param priorvar Optional Transport reference column containing one positive
#'   design-prior weight per candidate.
#'
#' @return An engine-specific fitted object. Use `broom::tidy()` for the sole
#'   primary endpoint, `gaze_info_bits`.
#' @export
gaze_weave_cv <- function(
    ref_tab, source_tab, match_on, contrast_on = NULL,
    refvar = "fixgroup", sourcevar = "fixgroup",
    spec = NULL, engine = NULL, split_on = match_on,
    n_folds = NULL, seed = NULL,
    fit_source_filter = NULL, eval_source_filter = NULL,
    episode_on = NULL, priorvar = NULL) {
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
      stop("episode_on and priorvar are available only for Transport.")
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
      eval_source_filter = eval_source_filter
    ))
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
