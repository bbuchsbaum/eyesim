# Nested calibration and exhaustive-candidate Transport -----------------

transport_v3_path_quality <- function(path) {
  measure <- as_gaze_measure(path, gaze_local_order())
  fixations <- measure$fixations
  coalesced <- coalesce_adjacent_gaze_fixations(
    fixations, distance = 0, kept_index = seq_len(nrow(fixations))
  )$fixations
  mass <- coalesced$duration / sum(coalesced$duration)
  entropy_bits <- -sum(mass * log2(mass))
  normalized_entropy <- if (length(mass) == 1L) {
    0
  } else {
    entropy_bits / log2(length(mass))
  }
  center <- c(sum(coalesced$x * mass), sum(coalesced$y * mass))
  radius_sq <- (coalesced$x - center[[1L]])^2 +
    (coalesced$y - center[[2L]])^2
  concentration <- sum(mass^2)
  list(
    raw_fixation_count = as.integer(nrow(fixations)),
    coalesced_fixation_count = as.integer(nrow(coalesced)),
    effective_fixations = as.numeric(1 / concentration),
    total_duration = as.numeric(sum(coalesced$duration)),
    duration_concentration = as.numeric(concentration),
    duration_entropy_bits = as.numeric(entropy_bits),
    normalized_duration_entropy = as.numeric(normalized_entropy),
    spatial_dispersion = as.numeric(sqrt(sum(mass * radius_sq))),
    generic_gaze_quality = as.numeric(normalized_entropy)
  )
}

transport_v3_reliability <- function(effective_fixations, kappa) {
  gaze_replay_reliability(effective_fixations, kappa)
}

transport_v3_candidate_score <- function(episode_scores, temperature) {
  if (!is.numeric(episode_scores) || length(episode_scores) == 0L ||
      anyNA(episode_scores) || any(!is.finite(episode_scores))) {
    stop("episode_scores must contain finite values.")
  }
  transport_v3_log_mean_exp(episode_scores / temperature)
}

transport_v3_profile_scores <- function(profile, temperature) {
  if (!is.list(profile) || length(profile) < 2L) {
    stop("profile must contain episode scores for at least two candidates.")
  }
  vapply(profile, transport_v3_candidate_score, numeric(1),
         temperature = temperature)
}

transport_v3_calibration_loss <- function(
    temperature, kappa, profile_sets, true_index, prior_sets,
    effective_fixations) {
  reliability <- transport_v3_reliability(effective_fixations, kappa)
  losses <- vapply(seq_along(profile_sets), function(index) {
    score <- transport_v3_profile_scores(
      profile_sets[[index]], temperature
    )
    base <- gaze_candidate_probabilities(score, prior_sets[[index]], 1)
    posterior <- reliability[[index]] * base +
      (1 - reliability[[index]]) * prior_sets[[index]]
    -log(posterior[[true_index[[index]]]])
  }, numeric(1))
  mean(losses)
}

fit_transport_v3_calibration <- function(
    profile_sets, true_index, effective_fixations, prior = NULL,
    spec = gaze_transport_spec()) {
  if (!is.list(profile_sets) || length(profile_sets) < 2L) {
    stop("profile_sets must contain at least two inner out-of-fold rows.")
  }
  score_sets <- lapply(profile_sets, transport_v3_profile_scores,
                       temperature = 1)
  true_index <- validate_gaze_truth(true_index, score_sets)
  prior_sets <- as_gaze_prior_sets(prior, score_sets)
  if (!is.numeric(effective_fixations) ||
      length(effective_fixations) != length(profile_sets) ||
      any(!is.finite(effective_fixations)) || any(effective_fixations <= 0)) {
    stop("effective_fixations must provide positive support for every row.")
  }
  bounds <- spec$temperature_bounds
  temperature_penalty <- if (is.infinite(
    spec$calibration$log_temperature_prior_sd
  )) 0 else {
    1 / (
      2 * length(profile_sets) *
        spec$calibration$log_temperature_prior_sd^2
    )
  }
  temperature_fit <- stats::optimize(
    function(log_temperature) {
      transport_v3_calibration_loss(
        exp(log_temperature), 0, profile_sets, true_index, prior_sets,
        effective_fixations
      ) + temperature_penalty * log_temperature^2
    },
    interval = log(bounds)
  )
  temperature_only <- exp(temperature_fit$minimum)
  temperature_only_loss <- transport_v3_calibration_loss(
    temperature_only, 0, profile_sets, true_index, prior_sets,
    effective_fixations
  )

  temperature <- temperature_only
  kappa <- 0
  selected_loss <- temperature_only_loss
  selected_boundary <- TRUE
  optimizer <- NULL
  if (identical(spec$reliability, "effective_fixations")) {
    kappa_bounds <- spec$reliability_kappa_bounds
    lower <- c(log(bounds[[1L]]), 0)
    upper <- c(log(bounds[[2L]]), log1p(kappa_bounds[[2L]]))
    kappa_penalty <- if (is.infinite(
      spec$calibration$log1p_kappa_prior_sd
    )) 0 else {
      1 / (
        2 * length(profile_sets) *
          spec$calibration$log1p_kappa_prior_sd^2
      )
    }
    objective <- function(parameter) {
      candidate_temperature <- exp(parameter[[1L]])
      candidate_kappa <- expm1(parameter[[2L]])
      transport_v3_calibration_loss(
        candidate_temperature, candidate_kappa, profile_sets, true_index,
        prior_sets, effective_fixations
      ) + temperature_penalty * parameter[[1L]]^2 +
        kappa_penalty * parameter[[2L]]^2
    }
    starts <- rbind(
      c(log(temperature_only), 0),
      c(log(temperature_only), log1p(kappa_bounds[[2L]]) / 2),
      c(0, log1p(kappa_bounds[[2L]]) / 4)
    )
    fits <- lapply(seq_len(nrow(starts)), function(index) {
      tryCatch(
        stats::optim(
          starts[index, ], objective, method = "L-BFGS-B",
          lower = lower, upper = upper
        ),
        error = function(condition) list(
          par = starts[index, ], value = objective(starts[index, ]),
          convergence = 999L, message = conditionMessage(condition)
        )
      )
    })
    values <- vapply(fits, function(fit) {
      if (is.null(fit$value) || !is.finite(fit$value)) Inf else fit$value
    }, numeric(1))
    optimizer <- fits[[which.min(values)]]
    proposed_temperature <- exp(optimizer$par[[1L]])
    proposed_kappa <- expm1(optimizer$par[[2L]])
    proposed_loss <- transport_v3_calibration_loss(
      proposed_temperature, proposed_kappa, profile_sets, true_index,
      prior_sets, effective_fixations
    )
    if (proposed_loss <= temperature_only_loss + 1e-10) {
      temperature <- proposed_temperature
      kappa <- proposed_kappa
      selected_loss <- proposed_loss
      selected_boundary <- kappa <= 1e-10
    }
  }

  reliability <- transport_v3_reliability(effective_fixations, kappa)
  evidence <- lapply(seq_along(profile_sets), function(index) {
    score_gaze_candidates(
      transport_v3_profile_scores(profile_sets[[index]], temperature),
      true_index = true_index[[index]],
      candidate_key = names(profile_sets[[index]]),
      prior = prior_sets[[index]],
      reliability = reliability[[index]],
      candidate_pool_id = "inner-out-of-fold"
    )
  })
  structure(
    list(
      temperature = as.numeric(temperature),
      kappa = as.numeric(kappa),
      log_loss = as.numeric(selected_loss),
      temperature_only_temperature = as.numeric(temperature_only),
      temperature_only_log_loss = as.numeric(temperature_only_loss),
      reliability_worsening = as.numeric(
        selected_loss - temperature_only_loss
      ),
      selected_temperature_only_boundary = selected_boundary,
      training_n = length(profile_sets),
      bounds = bounds,
      kappa_bounds = spec$reliability_kappa_bounds,
      reliability_policy = spec$reliability,
      scheme = if (identical(spec$reliability, "effective_fixations")) {
        paste(
          "inner_oof_episode_scale_temperature",
          "response_blind_effective_fixation_shrinkage",
          sep = "_"
        )
      } else {
        "inner_oof_episode_scale_temperature_only"
      },
      evidence = evidence,
      reliability_curve = gaze_reliability_table(evidence),
      optimizer = optimizer
    ),
    class = c("gaze_transport_calibration", "list")
  )
}

transport_v3_reference_bank <- function(
    ref_tab, source_row, match_on, contrast_on, refvar, episode_on,
    priorvar) {
  ref_key <- gaze_key(ref_tab, match_on, "match_on")
  source_key <- gaze_key(source_row, match_on, "match_on")
  source_contrast <- if (is.null(contrast_on)) {
    "all"
  } else {
    gaze_key(source_row, contrast_on, "contrast_on")
  }
  ref_contrast <- if (is.null(contrast_on)) {
    rep("all", nrow(ref_tab))
  } else {
    gaze_key(ref_tab, contrast_on, "contrast_on")
  }
  candidate_rows <- which(ref_contrast == source_contrast)
  candidate_keys <- sort(unique(ref_key[candidate_rows]))
  if (length(candidate_keys) < 2L || !source_key %in% candidate_keys) {
    stop("Every row requires a true candidate and at least one permitted nonmatch.")
  }
  candidates <- lapply(candidate_keys, function(candidate_key) {
    rows <- candidate_rows[ref_key[candidate_rows] == candidate_key]
    episode_id <- if (is.null(episode_on)) {
      rep("episode_1", length(rows))
    } else {
      gaze_key(ref_tab[rows, , drop = FALSE], episode_on, "episode_on")
    }
    if (anyDuplicated(episode_id)) {
      stop("episode_on must uniquely identify presentations within a candidate.")
    }
    stats::setNames(ref_tab[[refvar]][rows], episode_id)
  })
  names(candidates) <- candidate_keys
  prior <- if (is.null(priorvar)) {
    rep(1 / length(candidate_keys), length(candidate_keys))
  } else {
    vapply(candidate_keys, function(candidate_key) {
      rows <- candidate_rows[ref_key[candidate_rows] == candidate_key]
      value <- unique(ref_tab[[priorvar]][rows])
      if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
          value <= 0) {
        stop("priorvar must provide one finite positive design weight per candidate.")
      }
      value
    }, numeric(1))
  }
  prior <- normalize_gaze_prior(prior, length(candidate_keys))
  list(
    candidates = candidates,
    candidate_keys = candidate_keys,
    true_key = source_key,
    true_index = match(source_key, candidate_keys),
    prior = prior,
    contrast_key = source_contrast,
    pool_id = paste0("exhaustive:", source_contrast)
  )
}

transport_v3_warp_reference_table <- function(
    ref_tab, match_on, refvar, episode_on, spec) {
  ref_key <- gaze_key(ref_tab, match_on, "match_on")
  keys <- sort(unique(ref_key))
  rows <- lapply(keys, function(key) {
    indices <- which(ref_key == key)
    if (is.null(episode_on) && length(indices) != 1L) {
      stop("Duplicate reference keys require episode_on.")
    }
    measures <- lapply(
      ref_tab[[refvar]][indices], as_gaze_measure,
      chronology = spec$chronology
    )
    episode_weight <- 1 / length(measures)
    pooled <- do.call(rbind, lapply(measures, function(measure) {
      fixation <- measure$fixations
      fixation$duration <- episode_weight *
        fixation$duration / sum(fixation$duration)
      fixation
    }))
    result <- ref_tab[indices[[1L]], unique(c(match_on, spec$warp$fit_by)),
                      drop = FALSE]
    result[[refvar]] <- list(fixation_group(
      x = pooled$x, y = pooled$y, duration = pooled$duration,
      onset = seq(0, length.out = nrow(pooled))
    ))
    result
  })
  dplyr::bind_rows(rows)
}

fit_transport_v3_warp <- function(
    ref_tab, source_tab, match_on, refvar, sourcevar, episode_on, spec) {
  pooled_reference <- transport_v3_warp_reference_table(
    ref_tab, match_on, refvar, episode_on, spec
  )
  fit_gaze_warp_model(
    pooled_reference, source_tab, match_on, refvar, sourcevar, spec
  )
}

transport_v3_row_warp <- function(warp, source_row, spec) {
  if (identical(warp$type, "none")) return(warp)
  group <- if (is.null(spec$warp$fit_by)) {
    names(warp$group_models)[[1L]]
  } else {
    gaze_key(source_row, spec$warp$fit_by, "warp fit_by")
  }
  subset_gaze_warp_model(warp, group)
}

score_transport_v3_cv_row <- function(
    source_row, ref_tab, match_on, contrast_on, refvar, sourcevar,
    episode_on, priorvar, spec, warp, temperature = 1, kappa = 0) {
  bank <- transport_v3_reference_bank(
    ref_tab, source_row, match_on, contrast_on, refvar, episode_on, priorvar
  )
  row_warp <- transport_v3_row_warp(warp, source_row, spec)
  aligner <- function(reference, source, spec, candidate_key) {
    gaze_transport_align(
      reference, source, spec, warp_model = row_warp,
      candidate_key = candidate_key
    )
  }
  quality <- transport_v3_path_quality(source_row[[sourcevar]][[1L]])
  reliability <- if (identical(spec$reliability, "effective_fixations")) {
    transport_v3_reliability(quality$effective_fixations, kappa)
  } else {
    1
  }
  scored <- score_transport_v3_episode_candidates(
    source = source_row[[sourcevar]][[1L]],
    reference_candidates = bank$candidates,
    chronology = spec$chronology,
    spec = spec,
    true_key = bank$true_key,
    prior = bank$prior,
    temperature = temperature,
    reliability = reliability,
    candidate_pool_id = bank$pool_id,
    aligner = aligner
  )
  profile <- lapply(scored$candidates, function(candidate) {
    candidate$diagnostics$episode_scores
  })
  list(
    evidence = scored$evidence,
    profile = profile,
    quality = quality,
    reliability = reliability,
    candidates = scored$candidates,
    episode_evidence = scored$episode_evidence,
    common_episode_ids = scored$common_episode_ids,
    omitted_by_candidate = scored$omitted_by_candidate,
    prior = bank$prior,
    true_index = bank$true_index,
    candidate_keys = bank$candidate_keys,
    pool_id = bank$pool_id,
    warp = row_warp,
    all_converged = all(vapply(scored$candidates, function(candidate) {
      candidate$convergence$converged
    }, logical(1)))
  )
}

transport_v3_identity_calibration <- function(spec, reason) {
  structure(
    list(
      temperature = 1,
      kappa = 0,
      log_loss = NA_real_,
      temperature_only_temperature = 1,
      temperature_only_log_loss = NA_real_,
      reliability_worsening = 0,
      selected_temperature_only_boundary = TRUE,
      training_n = 0L,
      bounds = spec$temperature_bounds,
      kappa_bounds = spec$reliability_kappa_bounds,
      reliability_policy = spec$reliability,
      scheme = "identity_boundary_fallback",
      reason = reason,
      evidence = list(),
      reliability_curve = data.frame()
    ),
    class = c("gaze_transport_calibration", "list")
  )
}

fit_transport_v3_inner_calibration <- function(
    ref_tab, source_tab, match_on, contrast_on, refvar, sourcevar,
    episode_on, priorvar, spec, fold_contrast_on = contrast_on) {
  source_key <- gaze_key(source_tab, match_on, "match_on")
  inner <- tryCatch(
    make_gaze_weave_folds(
      source_tab, split_on = match_on, contrast_on = fold_contrast_on,
      n_folds = spec$calibration$folds, seed = spec$calibration$seed
    ),
    error = function(condition) condition
  )
  if (inherits(inner, "error")) {
    return(list(
      calibration = transport_v3_identity_calibration(
        spec, conditionMessage(inner)
      ),
      receipts = list()
    ))
  }
  profiles <- vector("list", nrow(source_tab))
  true_index <- integer(nrow(source_tab))
  priors <- vector("list", nrow(source_tab))
  effective_fixations <- numeric(nrow(source_tab))
  receipts <- vector("list", inner$n_folds)
  for (fold in seq_len(inner$n_folds)) {
    eval_rows <- which(inner$fold_id == fold)
    train_rows <- which(inner$fold_id != fold)
    train_keys <- unique(source_key[train_rows])
    eval_keys <- unique(source_key[eval_rows])
    ref_key <- gaze_key(ref_tab, match_on, "match_on")
    warp_reference <- ref_tab[ref_key %in% train_keys, , drop = FALSE]
    warp_key <- gaze_key(
      warp_reference, c(match_on, episode_on), "warp reference"
    )
    warp_reference <- warp_reference[!duplicated(warp_key), , drop = FALSE]
    warp <- fit_transport_v3_warp(
      warp_reference,
      source_tab[train_rows, , drop = FALSE],
      match_on, refvar, sourcevar, episode_on, spec
    )
    for (row in eval_rows) {
      scored <- score_transport_v3_cv_row(
        source_tab[row, , drop = FALSE], ref_tab, match_on, contrast_on,
        refvar, sourcevar, episode_on, priorvar, spec, warp
      )
      profiles[[row]] <- scored$profile
      true_index[[row]] <- scored$true_index
      priors[[row]] <- scored$prior
      effective_fixations[[row]] <- scored$quality$effective_fixations
    }
    receipts[[fold]] <- list(
      fold = fold,
      train_rows = source_tab$..gaze_row_id[train_rows],
      eval_rows = source_tab$..gaze_row_id[eval_rows],
      train_match_keys = train_keys,
      eval_match_keys = eval_keys,
      overlap_match_n = length(intersect(train_keys, eval_keys)),
      warp = warp$info,
      candidate_policy = "exhaustive_permitted_with_actual_prior",
      episode_policy = "equal_prior_common_valid_intersection"
    )
  }
  calibration <- fit_transport_v3_calibration(
    profiles, true_index, effective_fixations, priors, spec
  )
  calibration$folds <- receipts
  calibration$seed <- spec$calibration$seed
  list(calibration = calibration, receipts = receipts)
}

transport_v3_expected_calibration_error <- function(reliability_table) {
  if (!is.data.frame(reliability_table) || nrow(reliability_table) == 0L) {
    return(NA_real_)
  }
  sum(
    reliability_table$n * abs(
      reliability_table$mean_probability -
        reliability_table$observed_frequency
    )
  ) / sum(reliability_table$n)
}

#' Nested-calibrated exhaustive-candidate Transport
#'
#' Outer folds hold out participant-item cells. Within each outer-training set,
#' warp fitting and episode-scale temperature/reliability calibration are
#' repeated on item-grouped inner folds. Every held-out row is scored against
#' every permitted candidate in its contrast stratum with the actual declared
#' candidate prior and one common solver, warp, and equal-episode policy.
#'
#' @inheritParams gaze_weave_cv
#' @param spec A [gaze_transport_spec()].
#' @param episode_on Optional columns identifying the separate study
#'   presentations within a candidate. Their likelihoods receive equal prior
#'   weight and are never concatenated.
#' @param priorvar Optional reference-table column containing one positive
#'   design-prior weight per candidate. `NULL` declares a uniform design.
#'
#' @return A `gaze_transport_fit` with held-out candidate probabilities,
#'   fold receipts, response-blind quality diagnostics, and calibration checks.
#' @export
gaze_transport_cv <- function(
    ref_tab, source_tab, match_on, contrast_on = NULL,
    refvar = "fixgroup", sourcevar = "fixgroup",
    spec = gaze_transport_spec(), split_on = match_on,
    n_folds = NULL, seed = 20260822L,
    fit_source_filter = NULL, eval_source_filter = NULL,
    episode_on = NULL, priorvar = NULL) {
  if (!inherits(spec, "gaze_transport_spec")) {
    stop("spec must be created by gaze_transport_spec().")
  }
  required_ref <- unique(c(
    match_on, contrast_on, episode_on, priorvar, spec$warp$fit_by, refvar
  ))
  required_source <- unique(c(
    match_on, contrast_on, split_on, spec$warp$fit_by, sourcevar
  ))
  if (!all(required_ref %in% names(ref_tab)) ||
      !all(required_source %in% names(source_tab))) {
    stop("Transport matching, episode, prior, warp, and path columns must exist.")
  }
  if (is.null(episode_on) && anyDuplicated(gaze_key(
    ref_tab, match_on, "match_on"
  ))) {
    stop("Duplicate reference candidates require episode_on.")
  }
  source_tab <- dplyr::ungroup(source_tab)
  source_tab$..gaze_row_id <- seq_len(nrow(source_tab))
  source_key <- gaze_key(source_tab, match_on, "match_on")
  ref_key <- gaze_key(ref_tab, match_on, "match_on")
  if (any(!source_key %in% ref_key)) {
    stop("Every Transport source key must exist in the reference table.")
  }
  fit_mask <- resolve_gaze_weave_filter(
    source_tab, fit_source_filter, "fit_source_filter"
  )
  eval_mask <- resolve_gaze_weave_filter(
    source_tab, eval_source_filter, "eval_source_filter"
  )
  folds <- make_gaze_weave_folds(
    source_tab, split_on, contrast_on, n_folds, seed
  )
  fold_results <- vector("list", folds$n_folds)
  fold_info <- vector("list", folds$n_folds)
  all_evidence <- list()
  for (fold in seq_len(folds$n_folds)) {
    eval_rows <- which(folds$fold_id == fold & eval_mask)
    if (length(eval_rows) == 0L) next
    eval_keys <- unique(source_key[eval_rows])
    train_candidates <- which(folds$fold_id != fold & fit_mask)
    train_rows <- train_candidates[
      !source_key[train_candidates] %in% eval_keys
    ]
    train_keys <- unique(source_key[train_rows])
    if (length(train_rows) < 2L) {
      stop("Fold ", fold, " has fewer than two Transport training rows.")
    }
    warp <- fit_transport_v3_warp(
      ref_tab[ref_key %in% train_keys, , drop = FALSE],
      source_tab[train_rows, , drop = FALSE],
      match_on, refvar, sourcevar, episode_on, spec
    )
    inner <- fit_transport_v3_inner_calibration(
      ref_tab, source_tab[train_rows, , drop = FALSE], match_on,
      contrast_on, refvar, sourcevar, episode_on, priorvar, spec
    )
    calibration <- inner$calibration
    scored <- lapply(eval_rows, function(row) {
      score_transport_v3_cv_row(
        source_tab[row, , drop = FALSE], ref_tab, match_on, contrast_on,
        refvar, sourcevar, episode_on, priorvar, spec, warp,
        temperature = calibration$temperature,
        kappa = calibration$kappa
      )
    })
    result_fold <- source_tab[eval_rows, , drop = FALSE]
    result_fold$.cv_fold <- fold
    for (field in c(
      "gaze_info_bits", "base_gaze_info_bits", "odds_bits",
      "posterior_true", "prior_true", "log_loss", "base_log_loss",
      "brier_score", "template_rank", "top1_credit", "candidate_count",
      "temperature", "reliability"
    )) {
      result_fold[[field]] <- vapply(scored, function(result) {
        result$evidence[[field]]
      }, numeric(1))
    }
    result_fold$candidate_count <- as.integer(result_fold$candidate_count)
    quality_fields <- names(scored[[1L]]$quality)
    for (field in quality_fields) {
      prototype <- scored[[1L]]$quality[[field]]
      result_fold[[field]] <- vapply(
        scored, function(result) result$quality[[field]],
        if (is.integer(prototype)) integer(1) else numeric(1)
      )
    }
    result_fold$common_episode_count <- vapply(
      scored, function(result) length(result$common_episode_ids), integer(1)
    )
    result_fold$all_converged <- vapply(
      scored, `[[`, logical(1), "all_converged"
    )
    result_fold$candidates <- lapply(scored, function(result) {
      result$evidence$candidates
    })
    result_fold$episode_evidence <- lapply(scored, `[[`, "episode_evidence")
    result_fold$alignments <- lapply(scored, `[[`, "candidates")
    fold_results[[fold]] <- result_fold
    all_evidence <- c(all_evidence, lapply(scored, `[[`, "evidence"))
    fold_info[[fold]] <- list(
      fold = fold,
      train_rows = source_tab$..gaze_row_id[train_rows],
      eval_rows = source_tab$..gaze_row_id[eval_rows],
      train_match_keys = train_keys,
      eval_match_keys = eval_keys,
      overlap_match_n = length(intersect(train_keys, eval_keys)),
      calibration = calibration,
      inner_receipts = inner$receipts,
      warp = warp$info,
      candidate_policy = "exhaustive_permitted_with_actual_prior",
      episode_policy = "equal_prior_common_valid_intersection",
      heldout_cell_contributed_to_fit = FALSE
    )
  }
  results <- dplyr::bind_rows(fold_results)
  if (nrow(results) == 0L) {
    stop("No held-out Transport rows were scored.")
  }
  results <- results[order(results$..gaze_row_id), , drop = FALSE]
  results$..gaze_row_id <- NULL
  reliability_curve <- gaze_reliability_table(all_evidence)
  pool_sensitivity <- stats::aggregate(
    cbind(gaze_info_bits, log_loss) ~ candidate_count,
    data = results, FUN = mean
  )
  structure(
    list(
      results = results,
      spec = spec,
      folds = fold_info,
      calibration = list(
        heldout_log_loss = mean(results$log_loss),
        heldout_reliability_curve = reliability_curve,
        heldout_ece = transport_v3_expected_calibration_error(
          reliability_curve
        ),
        candidate_pool_size_sensitivity = pool_sensitivity
      ),
      keys = list(
        match_on = match_on,
        contrast_on = contrast_on,
        split_on = split_on,
        episode_on = episode_on,
        priorvar = priorvar,
        id_columns = unique(c(match_on, contrast_on, split_on))
      ),
      provenance = list(
        engine = "transport",
        primary_score = "gaze_info_bits=log2(p_true/prior_true)",
        candidate_pool = "exhaustive_permitted_within_contrast",
        candidate_prior = if (is.null(priorvar)) {
          "declared_uniform_design"
        } else {
          paste0("actual_design_prior:", priorvar)
        },
        calibration = "inner_oof_episode_scale_temperature_and_reliability",
        reliability = spec$reliability,
        generic_gaze_quality = "separate_response_blind_diagnostic_channel",
        seed = seed,
        n_folds = folds$n_folds
      )
    ),
    class = c("gaze_transport_fit", "list")
  )
}
