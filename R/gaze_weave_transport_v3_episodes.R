# Multi-presentation Transport -------------------------------------------

prepare_transport_v3_episodes <- function(episodes, chronology,
                                          episode_ids = NULL) {
  if (!inherits(chronology, "gaze_chronology_spec")) {
    stop("chronology must be a gaze chronology specification.")
  }
  if (inherits(episodes, "fixation_group") || inherits(episodes, "gaze_measure")) {
    episodes <- list(episodes)
  }
  if (!is.list(episodes) || length(episodes) == 0L) {
    stop("episodes must contain at least one fixation path.")
  }
  if (is.null(episode_ids)) {
    candidate_names <- names(episodes)
    episode_ids <- if (!is.null(candidate_names) &&
        all(!is.na(candidate_names)) && all(nzchar(candidate_names))) {
      candidate_names
    } else {
      paste0("episode_", seq_along(episodes))
    }
  }
  if (length(episode_ids) != length(episodes) || anyNA(episode_ids) ||
      any(!nzchar(as.character(episode_ids))) ||
      anyDuplicated(as.character(episode_ids))) {
    stop("episode_ids must contain one complete unique value per episode.")
  }
  episode_ids <- as.character(episode_ids)

  prepared <- vector("list", length(episodes))
  errors <- rep(NA_character_, length(episodes))
  for (index in seq_along(episodes)) {
    result <- tryCatch(
      {
        if (inherits(episodes[[index]], "gaze_measure")) {
          episodes[[index]]
        } else {
          as_gaze_measure(episodes[[index]], chronology)
        }
      },
      error = function(condition) condition
    )
    if (inherits(result, "error")) {
      errors[[index]] <- conditionMessage(result)
    } else {
      prepared[[index]] <- result
    }
  }
  names(prepared) <- episode_ids
  valid <- !vapply(prepared, is.null, logical(1))

  structure(
    list(
      episodes = prepared[valid],
      episode_ids = episode_ids,
      valid_ids = episode_ids[valid],
      invalid_ids = episode_ids[!valid],
      invalid_reasons = stats::setNames(errors[!valid], episode_ids[!valid]),
      requested_count = length(episodes),
      valid_count = sum(valid),
      provenance = list(
        preparation = "independent_per_episode",
        concatenated = FALSE,
        cross_episode_edges = FALSE,
        missing_rule = paste(
          "equal weights over common valid presentation indices",
          "across the candidate pool"
        )
      )
    ),
    class = c("gaze_transport_episodes", "list")
  )
}

transport_v3_common_episode_ids <- function(prepared_candidates) {
  if (!is.list(prepared_candidates) || length(prepared_candidates) == 0L ||
      !all(vapply(
        prepared_candidates,
        inherits,
        logical(1),
        "gaze_transport_episodes"
      ))) {
    stop("prepared_candidates must contain prepared Transport episodes.")
  }
  common <- prepared_candidates[[1]]$valid_ids
  if (length(prepared_candidates) > 1L) {
    for (index in 2:length(prepared_candidates)) {
      common <- common[common %in% prepared_candidates[[index]]$valid_ids]
    }
  }
  common
}

transport_v3_log_mean_exp <- function(values) {
  if (!is.numeric(values) || length(values) == 0L || anyNA(values) ||
      any(values == Inf) || !any(is.finite(values))) {
    stop("episode scores must contain a finite value and no NA or Inf.")
  }
  gaze_log_sum_exp(values) - log(length(values))
}

score_transport_v3_episode_candidate <- function(
    source, prepared_reference, spec, candidate_key,
    common_episode_ids = prepared_reference$valid_ids,
    temperature = 1, aligner) {
  if (!inherits(prepared_reference, "gaze_transport_episodes")) {
    stop("prepared_reference must be prepared Transport episodes.")
  }
  if (!is.function(aligner)) stop("aligner must be a function.")
  if (!is.numeric(temperature) || length(temperature) != 1L ||
      !is.finite(temperature) || temperature <= 0) {
    stop("temperature must be one finite positive value.")
  }
  common_episode_ids <- as.character(common_episode_ids)
  if (length(common_episode_ids) == 0L) {
    stop("No episode index is valid across the candidate pool.")
  }
  if (any(!common_episode_ids %in% prepared_reference$valid_ids)) {
    stop("common_episode_ids must be valid for this candidate.")
  }
  source_measure <- if (inherits(source, "gaze_measure")) {
    source
  } else {
    as_gaze_measure(source, spec$chronology)
  }
  episode_results <- lapply(common_episode_ids, function(episode_id) {
    aligner(
      prepared_reference$episodes[[episode_id]],
      source_measure,
      spec,
      candidate_key = candidate_key
    )
  })
  if (!all(vapply(episode_results, inherits, logical(1), "gaze_engine_result"))) {
    stop("aligner must return one gaze_engine_result per episode.")
  }
  episode_scores <- vapply(episode_results, `[[`, numeric(1), "log_score")
  if (anyNA(episode_scores) || any(!is.finite(episode_scores))) {
    stop("Every common episode must have a finite alignment score.")
  }
  converged <- vapply(episode_results, function(result) {
    isTRUE(result$convergence$converged)
  }, logical(1))
  scaled_scores <- episode_scores / temperature
  mixture_score <- transport_v3_log_mean_exp(scaled_scores)
  omitted_ids <- setdiff(prepared_reference$episode_ids, common_episode_ids)

  new_gaze_engine_result(
    engine = "transport",
    candidate_key = candidate_key,
    log_score = mixture_score,
    diagnostics = list(
      episode_count = length(common_episode_ids),
      requested_episode_count = prepared_reference$requested_count,
      common_episode_ids = common_episode_ids,
      omitted_episode_ids = omitted_ids,
      episode_scores = stats::setNames(episode_scores, common_episode_ids),
      scaled_episode_scores = stats::setNames(scaled_scores, common_episode_ids),
      equal_episode_weights = stats::setNames(
        rep(1 / length(common_episode_ids), length(common_episode_ids)),
        common_episode_ids
      ),
      temperature = temperature
    ),
    alignment = list(
      episodes = stats::setNames(episode_results, common_episode_ids),
      combination = "equal_prior_likelihood_mixture"
    ),
    convergence = list(
      converged = all(converged),
      episode_converged = stats::setNames(converged, common_episode_ids)
    ),
    provenance = list(
      engine_version = 3L,
      directionality = "symmetric",
      score_semantics = "equal_prior_episode_likelihood_mixture",
      duration_semantics = "unit_duration_mass_within_separate_episodes",
      candidate_invariant = TRUE,
      concatenated = FALSE,
      episode_weighting = "fixed_equal_prior",
      missing_rule = paste(
        "renormalize over intersection of valid presentation indices",
        "across candidate pool"
      )
    )
  )
}

score_transport_v3_episode_candidates <- function(
    source, reference_candidates, chronology, spec, true_key,
    prior = NULL, temperature = 1, reliability = 1,
    candidate_pool_id = "declared", aligner) {
  if (!is.list(reference_candidates) || length(reference_candidates) < 2L) {
    stop("reference_candidates must contain at least two candidates.")
  }
  candidate_keys <- names(reference_candidates)
  if (is.null(candidate_keys) || any(!nzchar(candidate_keys)) ||
      anyDuplicated(candidate_keys)) {
    stop("reference_candidates must have complete unique candidate names.")
  }
  prepared <- lapply(reference_candidates, function(episodes) {
    if (inherits(episodes, "gaze_transport_episodes")) {
      episodes
    } else {
      prepare_transport_v3_episodes(episodes, chronology)
    }
  })
  common_episode_ids <- transport_v3_common_episode_ids(prepared)
  if (length(common_episode_ids) == 0L) {
    stop("No episode index is valid across the candidate pool.")
  }
  results <- lapply(seq_along(prepared), function(index) {
    score_transport_v3_episode_candidate(
      source = source,
      prepared_reference = prepared[[index]],
      spec = spec,
      candidate_key = candidate_keys[[index]],
      common_episode_ids = common_episode_ids,
      temperature = temperature,
      aligner = aligner
    )
  })
  evidence <- score_gaze_engine_results(
    results,
    true_key = true_key,
    prior = prior,
    temperature = 1,
    reliability = reliability,
    candidate_pool_id = candidate_pool_id
  )
  episode_evidence <- lapply(common_episode_ids, function(episode_id) {
    episode_results <- lapply(results, function(result) {
      result$alignment$episodes[[episode_id]]
    })
    score_gaze_engine_results(
      episode_results,
      true_key = true_key,
      prior = prior,
      temperature = temperature,
      reliability = reliability,
      candidate_pool_id = paste0(candidate_pool_id, ":", episode_id)
    )
  })
  names(episode_evidence) <- common_episode_ids

  list(
    evidence = evidence,
    candidates = stats::setNames(results, candidate_keys),
    episode_evidence = episode_evidence,
    common_episode_ids = common_episode_ids,
    omitted_by_candidate = lapply(prepared, `[[`, "invalid_ids"),
    provenance = list(
      candidate_pool_id = candidate_pool_id,
      episode_weighting = "fixed_equal_prior",
      common_valid_rule = "intersection_then_equal_renormalization"
    )
  )
}
