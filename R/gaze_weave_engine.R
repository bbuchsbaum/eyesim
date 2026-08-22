# Shared GazeWeave engine contract ------------------------------------------

gaze_engine_contract <- function() {
  list(
    version = 2L,
    primary_endpoint = "gaze_info_bits",
    score_scale = "natural_log",
    score_direction = "larger_is_better",
    estimands = c("coverage", "selection", "correspondence", "evidence"),
    engines = list(
      transport = list(directionality = "symmetric"),
      replay = list(directionality = "encoding_to_recall")
    ),
    required_provenance = c(
      "engine_version",
      "directionality",
      "score_semantics",
      "duration_semantics",
      "candidate_invariant"
    )
  )
}

new_gaze_engine_result <- function(engine, candidate_key, log_score,
                                   diagnostics = list(), alignment = NULL,
                                   convergence = list(converged = TRUE),
                                   provenance) {
  contract <- gaze_engine_contract()
  if (!is.character(engine) || length(engine) != 1L ||
      !engine %in% names(contract$engines)) {
    stop("engine must be one of: ", paste(names(contract$engines), collapse = ", "), ".")
  }
  if (length(candidate_key) != 1L || is.na(candidate_key)) {
    stop("candidate_key must be a complete scalar value.")
  }
  if (!is.numeric(log_score) || length(log_score) != 1L || is.na(log_score) ||
      log_score == Inf) {
    stop("log_score must be a numeric scalar that is finite or -Inf.")
  }
  if (!is.list(diagnostics) || is.null(names(diagnostics)) ||
      any(!nzchar(names(diagnostics)))) {
    stop("diagnostics must be a fully named list.")
  }
  if (!is.list(convergence) || is.null(convergence$converged) ||
      !is.logical(convergence$converged) || length(convergence$converged) != 1L ||
      is.na(convergence$converged)) {
    stop("convergence must contain one complete logical value named converged.")
  }
  if (!is.list(provenance)) {
    stop("provenance must be a named list.")
  }
  missing_provenance <- setdiff(contract$required_provenance, names(provenance))
  if (length(missing_provenance) > 0L) {
    stop("provenance is missing: ", paste(missing_provenance, collapse = ", "), ".")
  }
  expected_direction <- contract$engines[[engine]]$directionality
  if (!identical(provenance$directionality, expected_direction)) {
    stop("provenance directionality does not match the engine contract.")
  }
  if (!identical(provenance$candidate_invariant, TRUE)) {
    stop("candidate-invariant preprocessing must be confirmed in provenance.")
  }

  structure(
    list(
      engine = engine,
      candidate_key = candidate_key,
      log_score = as.numeric(log_score),
      diagnostics = diagnostics,
      alignment = alignment,
      convergence = convergence,
      provenance = provenance
    ),
    class = c("gaze_engine_result", "list")
  )
}

validate_gaze_engine_results <- function(results) {
  if (!is.list(results) || length(results) < 2L ||
      !all(vapply(results, inherits, logical(1), "gaze_engine_result"))) {
    stop("results must contain at least two gaze_engine_result objects.")
  }
  keys <- vapply(results, function(result) as.character(result$candidate_key), character(1))
  if (anyDuplicated(keys)) {
    stop("candidate keys must be unique within an evaluation set.")
  }
  engines <- vapply(results, `[[`, character(1), "engine")
  if (length(unique(engines)) != 1L) {
    stop("all candidate results must come from the same engine.")
  }
  scores <- vapply(results, `[[`, numeric(1), "log_score")
  if (!any(is.finite(scores))) {
    stop("at least one candidate log_score must be finite.")
  }
  invisible(TRUE)
}

# Candidate evidence --------------------------------------------------------

gaze_log_sum_exp <- function(x) {
  finite <- is.finite(x)
  if (!any(finite)) {
    return(-Inf)
  }
  maximum <- max(x[finite])
  maximum + log(sum(exp(x[finite] - maximum)))
}

normalize_gaze_prior <- function(prior, n_candidates) {
  if (is.null(prior)) {
    return(rep(1 / n_candidates, n_candidates))
  }
  if (!is.numeric(prior) || length(prior) != n_candidates ||
      any(!is.finite(prior)) || any(prior <= 0) || sum(prior) <= 0) {
    stop("prior must contain one finite positive value per candidate.")
  }
  as.numeric(prior / sum(prior))
}

gaze_candidate_log_probabilities <- function(log_score, prior = NULL,
                                             temperature = 1) {
  if (!is.numeric(log_score) || length(log_score) < 2L || anyNA(log_score) ||
      any(log_score == Inf) || !any(is.finite(log_score))) {
    stop("log_score must contain at least two values, at least one finite, and no NA or Inf.")
  }
  if (!is.numeric(temperature) || length(temperature) != 1L ||
      !is.finite(temperature) || temperature <= 0) {
    stop("temperature must be a finite positive scalar.")
  }
  prior <- normalize_gaze_prior(prior, length(log_score))
  log_weight <- log(prior) + log_score / temperature
  normalizer <- gaze_log_sum_exp(log_weight)
  log_posterior <- log_weight - normalizer
  log_posterior[!is.finite(log_weight)] <- -Inf
  log_posterior
}

gaze_candidate_probabilities <- function(log_score, prior = NULL,
                                         temperature = 1) {
  log_posterior <- gaze_candidate_log_probabilities(
    log_score, prior = prior, temperature = temperature
  )
  posterior <- exp(log_posterior)
  posterior / sum(posterior)
}

gaze_rank_summary <- function(posterior, true_index,
                              tolerance = sqrt(.Machine$double.eps)) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a finite non-negative scalar.")
  }
  true_probability <- posterior[[true_index]]
  threshold <- tolerance * max(1, max(abs(posterior)))
  better <- sum(posterior > true_probability + threshold)
  tied <- sum(abs(posterior - true_probability) <= threshold)
  c(
    template_rank = 1 + better + (tied - 1) / 2,
    top1_credit = if (better == 0L) 1 / tied else 0,
    tied_candidates = tied
  )
}

gaze_logit <- function(probability) {
  log(probability) - log1p(-probability)
}

score_gaze_candidates <- function(log_score, true_index, candidate_key = NULL,
                                  prior = NULL, temperature = 1,
                                  reliability = 1,
                                  candidate_pool_id = "declared",
                                  tie_tolerance = sqrt(.Machine$double.eps)) {
  n_candidates <- length(log_score)
  if (is.null(candidate_key)) {
    candidate_key <- seq_len(n_candidates)
  }
  if (length(candidate_key) != n_candidates || anyNA(candidate_key) ||
      anyDuplicated(as.character(candidate_key))) {
    stop("candidate_key must contain one complete unique value per candidate.")
  }
  true_index <- as.integer(true_index)
  if (length(true_index) != 1L || is.na(true_index) || true_index < 1L ||
      true_index > n_candidates) {
    stop("true_index must identify exactly one candidate.")
  }
  if (!is.character(candidate_pool_id) || length(candidate_pool_id) != 1L ||
      is.na(candidate_pool_id) || !nzchar(candidate_pool_id)) {
    stop("candidate_pool_id must be a non-empty string.")
  }
  if (!is.numeric(reliability) || length(reliability) != 1L ||
      !is.finite(reliability) || reliability < 0 || reliability > 1) {
    stop("reliability must be one finite value between zero and one.")
  }

  prior <- normalize_gaze_prior(prior, n_candidates)
  base_log_posterior <- gaze_candidate_log_probabilities(
    log_score, prior, temperature
  )
  base_posterior <- exp(base_log_posterior)
  base_posterior <- base_posterior / sum(base_posterior)
  log_posterior <- if (reliability == 1) {
    base_log_posterior
  } else if (reliability == 0) {
    log(prior)
  } else {
    vapply(seq_len(n_candidates), function(i) {
      gaze_log_sum_exp(c(
        log(reliability) + base_log_posterior[[i]],
        log1p(-reliability) + log(prior[[i]])
      ))
    }, numeric(1))
  }
  posterior <- exp(log_posterior)
  true_probability <- posterior[[true_index]]
  true_log_probability <- log_posterior[[true_index]]
  true_prior <- prior[[true_index]]
  rank <- gaze_rank_summary(posterior, true_index, tie_tolerance)
  truth <- rep(0, n_candidates)
  truth[[true_index]] <- 1

  candidate_table <- data.frame(
    candidate_key = as.character(candidate_key),
    log_score = as.numeric(log_score),
    prior = prior,
    base_log_posterior = base_log_posterior,
    base_posterior = base_posterior,
    log_posterior = log_posterior,
    posterior = posterior,
    is_true = seq_len(n_candidates) == true_index,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      gaze_info_bits = (true_log_probability - log(true_prior)) / log(2),
      base_gaze_info_bits = (
        base_log_posterior[[true_index]] - log(true_prior)
      ) / log(2),
      odds_bits = (
        true_log_probability - gaze_log_sum_exp(log_posterior[-true_index]) -
          gaze_logit(true_prior)
      ) / log(2),
      posterior_true = true_probability,
      prior_true = true_prior,
      log_loss = -true_log_probability,
      base_log_loss = -base_log_posterior[[true_index]],
      brier_score = sum((posterior - truth)^2),
      template_rank = unname(rank[["template_rank"]]),
      top1_credit = unname(rank[["top1_credit"]]),
      tied_candidates = as.integer(rank[["tied_candidates"]]),
      candidate_count = n_candidates,
      temperature = as.numeric(temperature),
      reliability = as.numeric(reliability),
      candidate_pool_id = candidate_pool_id,
      candidates = candidate_table
    ),
    class = c("gaze_candidate_evidence", "list")
  )
}

score_gaze_engine_results <- function(results, true_key, prior = NULL,
                                      temperature = 1,
                                      reliability = 1,
                                      candidate_pool_id = "declared",
                                      tie_tolerance = sqrt(.Machine$double.eps)) {
  validate_gaze_engine_results(results)
  keys <- vapply(results, function(result) as.character(result$candidate_key), character(1))
  true_index <- match(as.character(true_key), keys)
  if (is.na(true_index)) {
    stop("true_key is not present in the candidate results.")
  }
  score_gaze_candidates(
    log_score = vapply(results, `[[`, numeric(1), "log_score"),
    true_index = true_index,
    candidate_key = keys,
    prior = prior,
    temperature = temperature,
    reliability = reliability,
    candidate_pool_id = candidate_pool_id,
    tie_tolerance = tie_tolerance
  )
}

as_gaze_score_sets <- function(log_scores) {
  score_sets <- if (is.matrix(log_scores) || is.data.frame(log_scores)) {
    lapply(seq_len(nrow(log_scores)), function(i) as.numeric(log_scores[i, ]))
  } else {
    log_scores
  }
  if (!is.list(score_sets) || length(score_sets) < 2L ||
      any(vapply(score_sets, length, integer(1)) < 2L)) {
    stop("log_scores must contain at least two candidate sets with at least two candidates each.")
  }
  valid <- vapply(score_sets, function(scores) {
    is.numeric(scores) && !anyNA(scores) && !any(scores == Inf) && any(is.finite(scores))
  }, logical(1))
  if (!all(valid)) {
    stop("each candidate score set must contain at least one finite score and no NA or Inf.")
  }
  score_sets
}

as_gaze_prior_sets <- function(prior, score_sets) {
  if (is.null(prior)) {
    return(lapply(score_sets, function(scores) normalize_gaze_prior(NULL, length(scores))))
  }
  if (is.matrix(prior) || is.data.frame(prior)) {
    if (nrow(prior) != length(score_sets)) {
      stop("prior rows must match candidate score sets.")
    }
    prior <- lapply(seq_len(nrow(prior)), function(i) as.numeric(prior[i, ]))
  } else if (is.numeric(prior)) {
    prior <- rep(list(prior), length(score_sets))
  }
  if (!is.list(prior) || length(prior) != length(score_sets)) {
    stop("prior must be NULL, a reusable vector, or one vector per score set.")
  }
  Map(normalize_gaze_prior, prior, vapply(score_sets, length, integer(1)))
}

validate_gaze_truth <- function(true_index, score_sets) {
  true_index <- as.integer(true_index)
  if (length(true_index) != length(score_sets) || anyNA(true_index) ||
      any(true_index < 1L) ||
      any(true_index > vapply(score_sets, length, integer(1)))) {
    stop("true_index must identify one candidate in every score set.")
  }
  true_index
}

gaze_temperature_loss <- function(temperature, score_sets, true_index,
                                  prior_sets) {
  losses <- vapply(seq_along(score_sets), function(i) {
    log_posterior <- gaze_candidate_log_probabilities(
      score_sets[[i]], prior_sets[[i]], temperature
    )
    -log_posterior[[true_index[[i]]]]
  }, numeric(1))
  mean(losses)
}

fit_gaze_temperature <- function(log_scores, true_index, prior = NULL,
                                 bounds = c(0.05, 20),
                                 log_temperature_prior_sd = 1) {
  score_sets <- as_gaze_score_sets(log_scores)
  true_index <- validate_gaze_truth(true_index, score_sets)
  prior_sets <- as_gaze_prior_sets(prior, score_sets)
  if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds)) ||
      bounds[[1]] <= 0 || bounds[[2]] <= bounds[[1]]) {
    stop("bounds must contain two ordered finite positive temperatures.")
  }
  if (!is.numeric(log_temperature_prior_sd) ||
      length(log_temperature_prior_sd) != 1L ||
      is.na(log_temperature_prior_sd) || log_temperature_prior_sd <= 0) {
    stop("log_temperature_prior_sd must be one positive scalar.")
  }
  penalty_scale <- if (is.infinite(log_temperature_prior_sd)) {
    0
  } else {
    1 / (2 * length(score_sets) * log_temperature_prior_sd^2)
  }

  fit <- stats::optimize(
    function(log_temperature) {
      gaze_temperature_loss(
        exp(log_temperature), score_sets, true_index, prior_sets
      ) + penalty_scale * log_temperature^2
    },
    interval = log(bounds)
  )
  temperature <- exp(fit$minimum)
  boundary_tolerance <- 1e-4 * diff(log(bounds))

  structure(
    list(
      temperature = temperature,
      log_loss = gaze_temperature_loss(
        temperature, score_sets, true_index, prior_sets
      ),
      penalized_objective = fit$objective,
      uncalibrated_log_loss = gaze_temperature_loss(
        1, score_sets, true_index, prior_sets
      ),
      training_n = length(score_sets),
      bounds = as.numeric(bounds),
      log_temperature_prior_sd = log_temperature_prior_sd,
      at_boundary = any(abs(fit$minimum - log(bounds)) <= boundary_tolerance)
    ),
    class = c("gaze_temperature_fit", "list")
  )
}

crossfit_gaze_temperature <- function(log_scores, true_index, fold_id,
                                      prior = NULL, bounds = c(0.05, 20),
                                      candidate_pool_id = "declared") {
  score_sets <- as_gaze_score_sets(log_scores)
  true_index <- validate_gaze_truth(true_index, score_sets)
  prior_sets <- as_gaze_prior_sets(prior, score_sets)
  if (length(fold_id) != length(score_sets) || anyNA(fold_id) ||
      length(unique(fold_id)) < 2L) {
    stop("fold_id must assign every score set to at least two complete folds.")
  }

  evidence <- vector("list", length(score_sets))
  temperature <- numeric(length(score_sets))
  fold_info <- list()
  fold_levels <- sort(unique(fold_id))
  for (fold in fold_levels) {
    eval_rows <- which(fold_id == fold)
    train_rows <- which(fold_id != fold)
    fit <- fit_gaze_temperature(
      score_sets[train_rows],
      true_index[train_rows],
      prior = prior_sets[train_rows],
      bounds = bounds
    )
    temperature[eval_rows] <- fit$temperature
    for (i in eval_rows) {
      evidence[[i]] <- score_gaze_candidates(
        score_sets[[i]],
        true_index = true_index[[i]],
        prior = prior_sets[[i]],
        temperature = fit$temperature,
        candidate_pool_id = candidate_pool_id
      )
    }
    fold_info[[as.character(fold)]] <- list(
      fold = fold,
      train_rows = train_rows,
      eval_rows = eval_rows,
      overlap_n = length(intersect(train_rows, eval_rows)),
      fit = fit
    )
  }

  structure(
    list(
      evidence = evidence,
      temperature = temperature,
      folds = fold_info,
      candidate_pool_id = candidate_pool_id
    ),
    class = c("gaze_temperature_cv", "list")
  )
}

gaze_reliability_table <- function(evidence, bins = 10L) {
  if (!is.list(evidence) || length(evidence) == 0L ||
      !all(vapply(evidence, inherits, logical(1), "gaze_candidate_evidence"))) {
    stop("evidence must be a non-empty list of gaze_candidate_evidence objects.")
  }
  bins <- as.integer(bins)
  if (length(bins) != 1L || is.na(bins) || bins < 2L) {
    stop("bins must be an integer of at least two.")
  }
  candidates <- do.call(rbind, lapply(seq_along(evidence), function(i) {
    table <- evidence[[i]]$candidates
    table$observation <- i
    table
  }))
  candidates$bin <- cut(
    candidates$posterior,
    breaks = seq(0, 1, length.out = bins + 1L),
    include.lowest = TRUE,
    right = TRUE
  )
  split_rows <- split(candidates, candidates$bin, drop = TRUE)
  result <- do.call(rbind, lapply(names(split_rows), function(bin) {
    rows <- split_rows[[bin]]
    data.frame(
      bin = bin,
      n = nrow(rows),
      mean_probability = mean(rows$posterior),
      observed_frequency = mean(rows$is_true),
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}
