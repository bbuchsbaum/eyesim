# Shared candidate calibration, revision 2026.10 (plan A3) -----------------
#
# Both engines map a held-out row's candidate scores to probabilities with
#
#   p_ik  proportional to  prior_ik * exp(beta_i * s_ik),
#   beta_i = beta * (n_i / n_ref)^gamma,        T_i = 1 / beta_i,
#
# where s_ik is the candidate score after the optional typicality offset,
# n_i the row's evidence count and n_ref the geometric mean evidence of the
# calibration rows. For Transport, whose candidate score is a log mean over
# study episodes, beta_i multiplies each episode score inside the log mean,
# exactly as the temperature did before. gamma = 0 is one global temperature.
# (beta, gamma) are fitted on inner out-of-fold rows only. Revision 2026.08
# never reaches this file.

#' Calibration controls for GazeWeave revision 2026.10
#'
#' Controls how Transport and Replay turn held-out candidate scores into
#' candidate probabilities under revision `"2026.10"`. Pass the result as
#' `calibration_control` to [gaze_transport_spec()] or [gaze_replay_spec()].
#'
#' **Evidence-scaled temperature.** Row \eqn{i} is scored at temperature
#' \eqn{T_i = T (\bar n / n_i)^\gamma}, where \eqn{n_i} is the row's evidence
#' count (Transport: duration-effective fixations of the recall; Replay: the
#' number of recall fixations the HMM observes) and \eqn{\bar n} is the
#' geometric mean over the calibration rows. \eqn{T} and \eqn{\gamma} are
#' fitted jointly by minimising the candidate log loss of inner out-of-fold
#' rows; no held-out row enters the fit. \eqn{\gamma = 0} is one global
#' temperature.
#'
#' **Ranking.** Calibration changes probabilities and bits only. Under
#' revision `"2026.10"`, `template_rank` and `top1_credit` (and any AUC
#' computed from ranks) always use one fixed reference ranking: the log
#' posterior at inverse temperature one, \eqn{\log \pi_k + r_k}, where
#' \eqn{r_k} is the engine's native candidate score after any typicality
#' offset (Transport: the log mean of the episode scores; Replay: the trial
#' log likelihood) and \eqn{\pi_k} the declared prior (uniform unless
#' `priorvar` is given). This holds for both methods and every fitted
#' \eqn{T}, \eqn{\gamma} or Stein factor, including a calibration that
#' returns the prior. The ranking cannot be read off the calibrated
#' probabilities: for multi-episode Transport the inverse temperature acts
#' inside the log mean over episodes, which ranks by the arithmetic episode
#' mean as \eqn{1/T \to 0} and by the log mean at \eqn{T = 1}. Transport
#' per-episode ranks use the same rule. The reference score uses no label,
#' so relabelling a row permutes it but never changes it. Ties use a
#' relative tolerance of \eqn{\sqrt{\epsilon}} on this log scale.
#' Replay's trial log likelihoods are often near-ties (see `typicality`), so
#' a Replay top-1 or AUC should be computed from the candidate tables'
#' `ranking_score` with a tie tolerance declared before the analysis (for
#' example 0.01 nats), splitting credit among tied candidates.
#'
#' The fit penalises the standardised inverse temperature
#' \eqn{\tilde\beta = \sigma / T} with a Gaussian prior centred at zero, where
#' \eqn{\sigma} is the median within-row standard deviation of the calibration
#' rows' candidate scores (label-free). The prior therefore only ever pulls
#' toward the declared candidate prior (less confidence), never toward a
#' fixed temperature. The earlier prior on \eqn{\log T}, centred at
#' \eqn{T = 1}, pulled toward overconfidence whenever the scores' natural
#' scale exceeded one nat (Replay's total trial log likelihood). The
#' temperature has no upper bound (\eqn{T = \infty} returns the prior);
#' the engine's lower temperature bound still caps every row's confidence.
#'
#' A calibration fitted on a few dozen inner rows is noisy, and under a null
#' every positive inverse temperature is overconfidence on held-out rows.
#' The fitted \eqn{\tilde\beta} is therefore shrunk by the positive-part
#' Stein factor \eqn{\max(0, 1 - p / LR)}, where \eqn{LR = 2N(L_0 - L)} is
#' the inner likelihood-ratio statistic against the declared prior
#' (\eqn{L_0}: its log loss) and \eqn{p} the number of free parameters (two
#' when \eqn{\gamma} is fitted). Under a null \eqn{LR} is about
#' \eqn{\chi^2_p} and the factor is small or zero; with real evidence it is
#' close to one. A hard selection test (keep the fit only if it beats the
#' prior by AIC) was rejected: conditioning on passing inflates the selected
#' fit, which made weak signals overconfident.
#'
#' The effective-fixation reliability shrink (\eqn{\kappa}) is not used by
#' `"evidence_scaled"`: it could only make sparse rows less confident, never
#' rich rows more confident, and with \eqn{\gamma} it is a second,
#' weakly identified parameter on the same evidence axis.
#'
#' **Typicality offset.** With `typicality = TRUE`, each candidate's score
#' has its typicality \eqn{b_k} subtracted before the softmax. \eqn{b_k} is
#' the mean score of candidate \eqn{k} against a fixed, seeded subsample of at
#' most `typicality_sources` recall paths from the fold's training rows whose
#' item (the `typicality_item_on` key) differs from candidate \eqn{k}'s item.
#' It is computed from the candidate's encodings and training recalls only,
#' so it never depends on the held-out recall or on which candidate is
#' labelled true, and relabelling a held-out row leaves every candidate score
#' unchanged. The subsample is drawn round-robin over items. Scores are
#' centred within each source before averaging; the removed per-source
#' constant is common to all candidates and cancels in the softmax, but for
#' Replay (total log likelihood) it would otherwise dominate the sampling
#' noise, because it grows with the recall's fixation count. Each mean is
#' shrunk toward the common mean by the empirical-Bayes factor
#' \eqn{\tau^2 / (\tau^2 + s_k^2)}, where \eqn{s_k^2} is the sampling variance
#' of candidate \eqn{k}'s mean and \eqn{\tau^2} the between-candidate variance
#' beyond sampling noise, so that offsets made mostly of noise do not perturb
#' the ranking. When a candidate has fewer than two other-item sources, no
#' offset is applied in that fold and the status is recorded. Inner
#' calibration rows use offsets estimated from their own inner training rows.
#' Offset-adjusted scores are ranking scores, not normalised likelihoods.
#' Two approximations: Transport offsets and scales are estimated over each
#' candidate's valid study episodes, while a scored row uses the episodes
#' valid for its whole pool; and Replay offsets always use the training
#' background protocol, also under `background_support = "held_out"`. Subtracting a mean equalises the candidates'
#' average scores under generic gaze but not their score variances; the fit
#' records the per-candidate standard deviation over the typicality sources.
#'
#' @param method `"evidence_scaled"` (default) fits \eqn{T} and \eqn{\gamma}
#'   as described. `"global"` reproduces the pre-A3 revision `"2026.10"`
#'   calibration: one temperature with a log-normal prior centred at
#'   one and, when the engine's `reliability` is `"effective_fixations"`,
#'   the \eqn{\kappa} shrink. The pre-change Transport default used that
#'   shrink, so reproducing pre-A3 `"2026.10"` requires
#'   `reliability = "effective_fixations"` passed explicitly; with the
#'   (new) default `reliability = "none"` the global path fits a different
#'   model. Probabilities and bits then match exactly; the specification
#'   additionally carries this control. Rank and top-1 use the reference
#'   ranking described under Ranking, so they differ from the pre-change
#'   values only in rows where the fitted temperature or the \eqn{\kappa}
#'   mixture had reordered the candidates (multi-episode Transport, or a
#'   non-uniform prior).
#' @param gamma_bounds Bounds for \eqn{\gamma}; equal values fix it (for
#'   example `c(0, 0)` for one global temperature with the new prior).
#' @param inverse_temperature_prior_sd Standard deviation of the zero-centred
#'   Gaussian prior on the standardised inverse temperature \eqn{\tilde\beta}.
#'   `Inf` removes the penalty.
#' @param stein_shrinkage Logical. When `TRUE` (default), the fitted inverse
#'   temperature is multiplied by \eqn{\max(0, 1 - p / LR)}, where \eqn{LR}
#'   is the inner likelihood-ratio statistic of the fit against the declared
#'   prior and \eqn{p} its number of free parameters (see Details).
#' @param typicality `"none"`, `"mean"` (subtract the offset) or
#'   `"standardized"` (subtract the offset from the row-centred score and
#'   divide by the candidate's shrunk score standard deviation over the
#'   sources). `NULL` (the default) uses the engine default: `"standardized"`
#'   for Transport and `"none"` for Replay. On a centre-biased simulated null
#'   with six candidates per pool, `"standardized"` brought Transport's
#'   central-candidate argmax share per candidate from 0.40 to 0.16 (chance
#'   0.167, MC error 0.020) without lowering top-1 or AUC on signal data.
#'   For Replay it lowered the share only from 0.30 to about 0.25-0.27
#'   (MC error 0.015), and it is therefore opt-in for Replay. The residual
#'   is not a heavy-tail effect. Under these nulls the HMM background state
#'   absorbs most recall fixations, so the candidates' likelihoods nearly
#'   coincide: in 57\% of held-out rows every candidate lay within 0.01
#'   nats of the others. The very large per-candidate SD ratio across rows
#'   comes from these near-degenerate pools, and the central bias from the
#'   argmax among near-ties: central candidates won 0.54 of the tie rows
#'   against a 0.37 share of the candidates (few rows, so imprecise). An
#'   offset cannot remove it; a pre-declared tie tolerance for top-1 and
#'   AUC (see Ranking) is the appropriate treatment.
#' @param typicality_sources Maximum number of training recalls used per
#'   fold to estimate every candidate's offset.
#' @param typicality_item_on Columns defining an "item" for the
#'   other-item rule. `NULL` uses `setdiff(match_on, contrast_on)`, or
#'   `match_on` when that is empty.
#'
#' @return A `gaze_calibration_control`.
#' @export
gaze_calibration_control <- function(
    method = c("evidence_scaled", "global"),
    gamma_bounds = c(0, 1),
    inverse_temperature_prior_sd = 3,
    stein_shrinkage = TRUE,
    typicality = NULL,
    typicality_sources = 24L,
    typicality_item_on = NULL) {
  if (!is.null(typicality)) {
    typicality <- match.arg(typicality, c("none", "mean", "standardized"))
  }
  method <- match.arg(method)
  if (!is.logical(stein_shrinkage) || length(stein_shrinkage) != 1L ||
      is.na(stein_shrinkage)) {
    stop("stein_shrinkage must be TRUE or FALSE.")
  }
  if (!is.numeric(gamma_bounds) || length(gamma_bounds) != 2L ||
      any(!is.finite(gamma_bounds)) || gamma_bounds[[1L]] > gamma_bounds[[2L]]) {
    stop("gamma_bounds must be two finite non-decreasing values.")
  }
  if (!is.numeric(inverse_temperature_prior_sd) ||
      length(inverse_temperature_prior_sd) != 1L ||
      is.na(inverse_temperature_prior_sd) ||
      inverse_temperature_prior_sd <= 0) {
    stop("inverse_temperature_prior_sd must be one positive value.")
  }
  typicality_sources <- as.integer(typicality_sources)
  if (length(typicality_sources) != 1L || is.na(typicality_sources) ||
      typicality_sources < 2L) {
    stop("typicality_sources must be an integer of at least two.")
  }
  if (!is.null(typicality_item_on) &&
      (!is.character(typicality_item_on) || !length(typicality_item_on) ||
       anyNA(typicality_item_on))) {
    stop("typicality_item_on must be NULL or column names.")
  }
  if (identical(method, "global")) {
    if (!is.null(typicality) && !identical(typicality, "none")) {
      stop("The typicality offset requires method = \"evidence_scaled\".")
    }
    typicality <- "none"
  }
  structure(
    list(
      method = method,
      gamma_bounds = as.numeric(gamma_bounds),
      inverse_temperature_prior_sd = as.numeric(inverse_temperature_prior_sd),
      stein_shrinkage = stein_shrinkage,
      typicality = typicality,
      typicality_sources = typicality_sources,
      typicality_item_on = typicality_item_on
    ),
    class = c("gaze_calibration_control", "list")
  )
}

#' @export
print.gaze_calibration_control <- function(x, ...) {
  cat("GazeWeave calibration control (revision 2026.10)\n")
  cat("  method:", x$method, "\n")
  if (identical(x$method, "evidence_scaled")) {
    cat("  gamma bounds:", paste(format(x$gamma_bounds), collapse = " to "),
        "\n")
    cat("  inverse-temperature prior sd:", format(x$inverse_temperature_prior_sd),
        "\n")
  }
  cat("  typicality offset:", if (is.null(x$typicality)) {
    "engine default"
  } else if (!identical(x$typicality, "none")) {
    paste0(x$typicality, " (", x$typicality_sources, " training recalls)")
  } else "none", "\n")
  invisible(x)
}

# Engine defaults for the typicality offset (see gaze_calibration_control()).
gaze_default_typicality <- c(transport = "standardized", replay = "none")

# Resolve a spec's calibration control. Revision 2026.08 specs have none.
resolve_gaze_calibration_control <- function(calibration_control, revision,
                                             reliability, engine) {
  if (identical(revision, "2026.08")) {
    if (!is.null(calibration_control)) {
      stop("calibration_control requires revision \"2026.10\".")
    }
    return(NULL)
  }
  if (is.null(calibration_control)) {
    calibration_control <- gaze_calibration_control()
  }
  if (!inherits(calibration_control, "gaze_calibration_control")) {
    stop("calibration_control must be created by gaze_calibration_control().")
  }
  if (is.null(calibration_control$typicality)) {
    calibration_control["typicality"] <- list(
      gaze_default_typicality[[engine]]
    )
  }
  if (identical(calibration_control$method, "evidence_scaled") &&
      identical(reliability, "effective_fixations")) {
    stop(
      "reliability = \"effective_fixations\" (the kappa shrink) is replaced ",
      "by the evidence-scaled temperature under revision \"2026.10\". Use ",
      "reliability = \"none\", or gaze_calibration_control(method = ",
      "\"global\") to keep the kappa shrink."
    )
  }
  calibration_control
}

gaze_calibration_evidence_scaled <- function(spec) {
  control <- spec$calibration$control
  !is.null(control) && identical(control$method, "evidence_scaled")
}

gaze_calibration_typicality <- function(spec) {
  gaze_calibration_evidence_scaled(spec) &&
    !identical(spec$calibration$control$typicality, "none")
}

gaze_calibration_standardized <- function(spec) {
  gaze_calibration_evidence_scaled(spec) &&
    identical(spec$calibration$control$typicality, "standardized")
}

# Ranking scores of one row (a list of per-candidate score vectors, e.g.
# Transport episode scores; one value per candidate for Replay) given a
# typicality object (offsets: named per-candidate vectors matching the score
# names, or scalars; scale: named per-candidate factors or NULL).
#   "mean":         s - b
#   "standardized": (s - row centre - (b - grand)) * a
# The row centre is the mean over the row's candidates (per episode), which
# uses no label; it is needed only because the scale factors differ across
# candidates, so a row-level constant would no longer cancel.
gaze_typicality_adjust <- function(profile, typicality) {
  if (is.null(typicality)) return(profile)
  keys <- names(profile)
  offsets <- typicality$offsets[keys]
  adjusted <- stats::setNames(lapply(keys, function(key) {
    scores <- profile[[key]]
    offset <- offsets[[key]]
    if (!is.null(names(scores)) && !is.null(names(offset))) {
      offset <- offset[names(scores)]
    }
    scores - offset
  }), keys)
  if (is.null(typicality$scale)) return(adjusted)
  # Row centre per element (episode) across the row's candidates; every
  # candidate's profile has the same element names in the same order.
  centre <- Reduce(`+`, lapply(keys, function(key) unname(profile[[key]]))) /
    length(keys)
  stats::setNames(lapply(keys, function(key) {
    (adjusted[[key]] - centre + typicality$grand) * typicality$scale[[key]]
  }), keys)
}

gaze_log_mean_exp <- function(values) {
  gaze_log_sum_exp(values) - log(length(values))
}

# Candidate logits of one row at inverse temperature beta. `profile` holds,
# per candidate, its (offset-adjusted) episode scores; a single-score engine
# passes one value per candidate.
gaze_calibrated_logits <- function(profile, beta) {
  if (beta == 0) return(rep(0, length(profile)))
  vapply(profile, function(scores) gaze_log_mean_exp(beta * scores),
         numeric(1))
}

gaze_calibration_row_beta <- function(calibration, evidence) {
  if (!is.numeric(evidence) || any(!is.finite(evidence)) ||
      any(evidence <= 0)) {
    stop("evidence counts must be finite and positive.")
  }
  beta <- calibration$inverse_temperature *
    (evidence / calibration$evidence_reference)^calibration$gamma
  # The engine's lower temperature bound caps every row's confidence, also
  # rows whose evidence exceeds the reference.
  bounds <- calibration$temperature_bounds
  if (!is.null(bounds)) beta <- pmin(beta, 1 / bounds[[1L]])
  beta
}

gaze_calibration_row_loss <- function(profile, beta, prior, true_index) {
  log_posterior <- gaze_candidate_log_probabilities(
    gaze_calibrated_logits(profile, beta), prior, 1
  )
  -log_posterior[[true_index]]
}

# Fit (beta, gamma) on inner out-of-fold rows.
#
# profiles: list over rows of lists over candidates of numeric scores
#   (already offset-adjusted); evidence: positive evidence counts per row.
fit_gaze_evidence_calibration <- function(profiles, true_index, evidence,
                                          prior = NULL,
                                          control = gaze_calibration_control(),
                                          temperature_bounds = c(0.05, 100)) {
  if (!is.list(profiles) || length(profiles) < 2L ||
      any(lengths(profiles) < 2L)) {
    stop("profiles must contain at least two rows of at least two candidates.")
  }
  score_sets <- lapply(profiles, function(profile) {
    vapply(profile, gaze_log_mean_exp, numeric(1))
  })
  true_index <- validate_gaze_truth(true_index, score_sets)
  prior_sets <- as_gaze_prior_sets(prior, score_sets)
  if (!is.numeric(evidence) || length(evidence) != length(profiles) ||
      any(!is.finite(evidence)) || any(evidence <= 0)) {
    stop("evidence must give one finite positive count per row.")
  }
  n <- length(profiles)
  spread <- vapply(score_sets, stats::sd, numeric(1))
  score_scale <- stats::median(spread[is.finite(spread)])
  if (!is.finite(score_scale) || score_scale <= 0) score_scale <- 1
  evidence_reference <- exp(mean(log(evidence)))
  relative <- evidence / evidence_reference
  beta_max <- score_scale / temperature_bounds[[1L]]
  prior_sd <- control$inverse_temperature_prior_sd
  penalty_scale <- if (is.infinite(prior_sd)) 0 else 1 / (2 * n * prior_sd^2)

  loss <- function(standard_beta, gamma) {
    beta <- pmin(standard_beta / score_scale * relative^gamma,
                 1 / temperature_bounds[[1L]])
    mean(vapply(seq_len(n), function(i) {
      gaze_calibration_row_loss(profiles[[i]], beta[[i]], prior_sets[[i]],
                                true_index[[i]])
    }, numeric(1)))
  }
  objective <- function(standard_beta, gamma) {
    loss(standard_beta, gamma) + penalty_scale * standard_beta^2
  }
  # Profile out beta for a fixed gamma; the zero boundary is checked exactly.
  profile_beta <- function(gamma) {
    fit <- stats::optimize(function(b) objective(b, gamma),
                           interval = c(0, beta_max))
    candidates <- c(0, fit$minimum, beta_max)
    values <- vapply(candidates, objective, numeric(1), gamma = gamma)
    best <- which.min(values)
    list(standard_beta = candidates[[best]], objective = values[[best]])
  }
  bounds <- control$gamma_bounds
  if (bounds[[1L]] == bounds[[2L]]) {
    gamma <- bounds[[1L]]
    fitted <- profile_beta(gamma)
  } else {
    grid <- seq(bounds[[1L]], bounds[[2L]], length.out = 11L)
    grid_fit <- lapply(grid, profile_beta)
    grid_value <- vapply(grid_fit, `[[`, numeric(1), "objective")
    best <- which.min(grid_value)
    step <- diff(grid[1:2])
    refine <- stats::optimize(
      function(g) profile_beta(g)$objective,
      interval = c(max(bounds[[1L]], grid[[best]] - step),
                   min(bounds[[2L]], grid[[best]] + step))
    )
    if (refine$objective < grid_value[[best]]) {
      gamma <- refine$minimum
      fitted <- profile_beta(gamma)
    } else {
      gamma <- grid[[best]]
      fitted <- grid_fit[[best]]
    }
  }
  # Positive-part Stein shrinkage of the inverse temperature by the
  # likelihood-ratio statistic against the declared prior (beta = 0):
  # beta <- beta * max(0, 1 - p / LR), with p free parameters. Under a null,
  # LR is about chi-square(p) and the fitted beta is mostly noise that turns
  # into held-out overconfidence; with real evidence LR >> p and the factor
  # is near one. Unlike a hard selection test it is continuous in the data,
  # so it has no winner's-curse jump for weak signals.
  free_parameters <- 1L + as.integer(bounds[[1L]] < bounds[[2L]])
  unshrunk_beta <- fitted$standard_beta
  lr_statistic <- 2 * n * (loss(0, 0) - loss(unshrunk_beta, gamma))
  shrinkage_factor <- if (isTRUE(control$stein_shrinkage)) {
    if (is.finite(lr_statistic) && lr_statistic > 0) {
      max(0, 1 - free_parameters / lr_statistic)
    } else 0
  } else {
    1
  }
  if (shrinkage_factor < 1) {
    standard_beta <- unshrunk_beta * shrinkage_factor
    fitted <- list(standard_beta = standard_beta,
                   objective = objective(standard_beta, gamma))
  }
  # A fit whose inverse temperature is zero has no evidence scaling.
  if (fitted$standard_beta == 0) gamma <- 0
  global <- if (bounds[[1L]] <= 0 && bounds[[2L]] >= 0) {
    profile_beta(0)
  } else {
    NULL
  }
  inverse_temperature <- fitted$standard_beta / score_scale
  structure(
    list(
      inverse_temperature = inverse_temperature,
      temperature = if (inverse_temperature > 0) 1 / inverse_temperature else Inf,
      gamma = gamma,
      evidence_reference = evidence_reference,
      score_scale = score_scale,
      standardized_inverse_temperature = fitted$standard_beta,
      log_loss = loss(fitted$standard_beta, gamma),
      penalized_objective = fitted$objective,
      global_temperature = if (is.null(global)) NA_real_ else
        if (global$standard_beta > 0) score_scale / global$standard_beta else Inf,
      global_log_loss = if (is.null(global)) NA_real_ else
        loss(global$standard_beta, 0),
      uniform_log_loss = loss(0, 0),
      training_n = n,
      gamma_bounds = bounds,
      temperature_bounds = as.numeric(temperature_bounds),
      inverse_temperature_prior_sd = prior_sd,
      at_confidence_bound = fitted$standard_beta >= beta_max * (1 - 1e-6),
      at_prior = fitted$standard_beta == 0,
      unshrunk_standardized_inverse_temperature = unshrunk_beta,
      likelihood_ratio = lr_statistic,
      free_parameters = free_parameters,
      stein_factor = shrinkage_factor,
      kappa = 0,
      scheme = "inner_oof_evidence_scaled_temperature"
    ),
    class = c("gaze_evidence_calibration", "list")
  )
}

# Calibration used when too few inner rows exist: the declared prior is
# returned (no confidence), which can never be overconfident.
gaze_evidence_calibration_fallback <- function(control, temperature_bounds,
                                               reason) {
  structure(
    list(
      inverse_temperature = 0, temperature = Inf, gamma = 0,
      evidence_reference = 1, score_scale = NA_real_,
      standardized_inverse_temperature = 0, log_loss = NA_real_,
      penalized_objective = NA_real_, global_temperature = NA_real_,
      global_log_loss = NA_real_, uniform_log_loss = NA_real_,
      training_n = 0L, gamma_bounds = control$gamma_bounds,
      temperature_bounds = as.numeric(temperature_bounds),
      inverse_temperature_prior_sd = control$inverse_temperature_prior_sd,
      at_confidence_bound = FALSE, at_prior = TRUE, kappa = 0,
      scheme = "declared_prior_fallback", reason = reason
    ),
    class = c("gaze_evidence_calibration", "list")
  )
}

# Held-out evidence for one row under a fitted evidence calibration. The
# calibration sets the probabilities and bits; rank and top-1 always come
# from the fixed reference ranking (see gaze_reference_rank()).
score_gaze_calibrated_row <- function(profile, true_index, evidence,
                                      calibration, candidate_key = NULL,
                                      prior = NULL,
                                      candidate_pool_id = "declared") {
  beta <- gaze_calibration_row_beta(calibration, evidence)
  result <- score_gaze_candidates(
    gaze_calibrated_logits(profile, beta),
    true_index = true_index,
    candidate_key = candidate_key,
    prior = prior,
    temperature = 1,
    candidate_pool_id = candidate_pool_id
  )
  result$temperature <- if (beta > 0) 1 / beta else Inf
  result$inverse_temperature <- beta
  # Untempered, offset-adjusted candidate score (a ranking score, not a
  # normalised likelihood when an offset was subtracted).
  ranking <- vapply(profile, gaze_log_mean_exp, numeric(1))
  result$candidates$ranking_score <- ranking
  gaze_apply_reference_rank(result, ranking, true_index)
}

# The calibration-independent ranking of one row: the log posterior at the
# reference inverse temperature one, log prior_k + r_k, where r_k is the
# engine's native candidate score after any typicality offset (Transport:
# the log mean of the episode scores; Replay: the trial log likelihood).
# A fitted inverse temperature below one would move a multi-episode
# Transport ranking toward the arithmetic episode mean, and a zero one would
# tie every candidate, so rank and top-1 never use the fitted calibration.
# Ties use a relative tolerance on this log scale, declared here once. The
# score uses no label, so relabelling a row permutes, never changes, it.
gaze_reference_rank <- function(ranking_score, prior, true_index,
                                tolerance = sqrt(.Machine$double.eps)) {
  value <- log(normalize_gaze_prior(prior, length(ranking_score))) +
    ranking_score
  target <- value[[true_index]]
  if (target == -Inf) {
    better <- sum(value > -Inf)
    tied <- sum(value == -Inf)
  } else {
    threshold <- tolerance * max(1, abs(target))
    better <- sum(value > target + threshold)
    tied <- sum(abs(value - target) <= threshold)
  }
  c(
    template_rank = 1 + better + (tied - 1) / 2,
    top1_credit = if (better == 0) 1 / tied else 0,
    tied_candidates = tied
  )
}

# Replace the rank fields of a scored row with the reference ranking.
gaze_apply_reference_rank <- function(evidence, ranking_score, true_index) {
  rank <- gaze_reference_rank(
    ranking_score, evidence$candidates$prior, true_index
  )
  evidence$template_rank <- unname(rank[["template_rank"]])
  evidence$top1_credit <- unname(rank[["top1_credit"]])
  evidence$tied_candidates <- as.integer(rank[["tied_candidates"]])
  evidence
}

# Cross-fitted calibration summary shared by both engines' CV results.
gaze_calibration_fit_summary <- function(fold_info, results, typicality) {
  folds <- Filter(Negate(is.null), fold_info)
  fold_table <- do.call(rbind, lapply(folds, function(fold) {
    calibration <- fold$calibration
    if (is.null(calibration)) {
      # Fitted without a contrastive calibration set: temperature one.
      calibration <- gaze_replay_identity_calibration()
    }
    field <- function(name) {
      value <- calibration[[name]]
      if (is.null(value)) NA_real_ else value
    }
    data.frame(
      fold = fold$fold,
      temperature = calibration$temperature,
      gamma = calibration$gamma,
      inverse_temperature = calibration$inverse_temperature,
      evidence_reference = calibration$evidence_reference,
      at_prior = isTRUE(calibration$at_prior),
      inner_log_loss = field("log_loss"),
      inner_global_log_loss = field("global_log_loss"),
      inner_uniform_log_loss = field("uniform_log_loss"),
      likelihood_ratio = field("likelihood_ratio"),
      stein_factor = field("stein_factor"),
      training_n = field("training_n"),
      scheme = calibration$scheme,
      stringsAsFactors = FALSE
    )
  }))
  out <- list(fold_calibration = fold_table)
  # Scores are centred within each held-out row (a row constant is common to
  # every candidate and cancels), then summarised per candidate: the spread
  # of candidate means is what the offset removes; the per-candidate
  # standard deviations are what a mean offset cannot equalise.
  candidates <- do.call(rbind, lapply(results$candidates, function(table) {
    data.frame(
      candidate_key = table$candidate_key,
      raw_score = table$raw_score - mean(table$raw_score),
      ranking_score = table$ranking_score - mean(table$ranking_score),
      stringsAsFactors = FALSE
    )
  }))
  if (!is.null(candidates) && nrow(candidates)) {
    by_key <- split(candidates, candidates$candidate_key)
    dispersion <- do.call(rbind, lapply(names(by_key), function(key) {
      rows <- by_key[[key]]
      data.frame(
        candidate_key = key, n = nrow(rows),
        raw_mean = mean(rows$raw_score),
        raw_sd = if (nrow(rows) > 1L) stats::sd(rows$raw_score) else NA_real_,
        ranking_mean = mean(rows$ranking_score),
        ranking_sd = if (nrow(rows) > 1L) {
          stats::sd(rows$ranking_score)
        } else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
    out$heldout_score_dispersion <- dispersion
  }
  if (typicality) {
    out$typicality <- do.call(rbind, lapply(folds, `[[`, "typicality"))
  }
  out
}

# Typicality sources --------------------------------------------------------

gaze_typicality_item_on <- function(control, match_on, contrast_on) {
  if (!is.null(control$typicality_item_on)) return(control$typicality_item_on)
  item_on <- setdiff(match_on, contrast_on)
  if (length(item_on)) item_on else match_on
}

# A fixed, seeded subsample of training rows used as typicality sources for a
# whole fold. It depends on the training rows only. Rows are drawn round-robin
# over items (in a seeded order), so that every candidate keeps sources of
# several other items.
gaze_typicality_source_rows <- function(source_item, cap, seed) {
  n_train <- length(source_item)
  if (n_train <= cap) return(seq_len(n_train))
  restore_session_rng <- snapshot_session_rng()
  on.exit(restore_session_rng(), add = TRUE)
  set.seed(seed)
  items <- sort(unique(source_item))
  items <- items[sample.int(length(items))]
  queues <- lapply(items, function(item) {
    rows <- which(source_item == item)
    rows[sample.int(length(rows))]
  })
  depth <- max(lengths(queues))
  order <- unlist(lapply(seq_len(depth), function(level) {
    unlist(lapply(queues, function(rows) if (length(rows) >= level) rows[[level]]))
  }))
  sort(order[seq_len(cap)])
}

# Summarise a candidate-by-source score matrix (NA where a source shares the
# candidate's item) into offsets and the diagnostics the plan asks for.
gaze_typicality_source_centre <- function(score_matrix) {
  apply(score_matrix, 2, function(column) {
    if (any(!is.na(column))) mean(column, na.rm = TRUE) else NA_real_
  })
}

gaze_typicality_centre_sources <- function(score_matrix) {
  sweep(score_matrix, 2, gaze_typicality_source_centre(score_matrix))
}

# When any candidate has fewer than two other-item sources, no offset is
# applied in that fold (all offsets zero, status recorded): offsetting only
# some candidates would itself favour candidates.
#
# Each candidate's mean is shrunk toward the grand mean by the empirical-Bayes
# factor lambda_k = tau^2 / (tau^2 + se_k^2), where se_k^2 is the sampling
# variance of its mean over sources and tau^2 the between-candidate variance
# left after removing sampling noise (method of moments, floored at zero).
# Offsets that are mostly sampling noise would only add noise to the ranking;
# a real typicality difference (tau^2 >> se^2) is removed almost entirely.
# The grand mean is common to every candidate and cancels in the softmax.
gaze_typicality_summary <- function(score_matrix, candidate_key) {
  # Scores are first centred within each source (column): a source-level
  # constant (for Replay, mostly the recall's fixation count) is common to
  # all candidates, cancels in the softmax, and would otherwise dominate the
  # sampling variance of the offsets.
  # A non-finite score (e.g. a -Inf likelihood) carries no usable mean.
  score_matrix[!is.finite(score_matrix)] <- NA_real_
  score_matrix <- gaze_typicality_centre_sources(score_matrix)
  used <- rowSums(!is.na(score_matrix))
  enough <- all(used >= 2L)
  source_sd <- vapply(seq_len(nrow(score_matrix)), function(index) {
    if (used[[index]] >= 2L) stats::sd(score_matrix[index, ], na.rm = TRUE)
    else NA_real_
  }, numeric(1))
  source_mean <- vapply(seq_len(nrow(score_matrix)), function(index) {
    if (used[[index]] >= 1L) mean(score_matrix[index, ], na.rm = TRUE)
    else NA_real_
  }, numeric(1))
  shrinkage <- rep(0, length(candidate_key))
  grand_mean <- 0
  if (enough) {
    sampling_variance <- source_sd^2 / used
    grand_mean <- mean(source_mean)
    between <- if (length(source_mean) > 1L) {
      max(0, stats::var(source_mean) - mean(sampling_variance))
    } else 0
    shrinkage <- if (isTRUE(between > 0)) {
      between / (between + sampling_variance)
    } else {
      rep(0, length(candidate_key))
    }
  }
  offset <- if (enough) {
    grand_mean + shrinkage * (source_mean - grand_mean)
  } else {
    rep(0, length(candidate_key))
  }
  # Scale factors (used by typicality = "standardized"): each candidate's log
  # standard deviation over the sources, shrunk toward their mean by the same
  # empirical-Bayes rule (sampling variance of a log sd ~ 1 / (2 (n - 1))).
  # A candidate whose scores vary more under generic gaze is scaled down.
  scale <- rep(1, length(candidate_key))
  if (enough && isTRUE(all(source_sd > 0))) {
    log_sd <- log(source_sd)
    log_sd_variance <- 1 / (2 * (used - 1))
    centre <- mean(log_sd)
    between <- if (length(log_sd) > 1L) {
      max(0, stats::var(log_sd) - mean(log_sd_variance))
    } else 0
    weight <- if (isTRUE(between > 0)) between / (between + log_sd_variance) else 0
    scale <- exp(-(weight * (log_sd - centre)))
  }
  data.frame(
    candidate_key = candidate_key,
    offset = as.numeric(offset),
    scale = as.numeric(scale),
    source_mean = source_mean,
    source_sd = source_sd,
    n_sources = as.integer(used),
    shrinkage = as.numeric(shrinkage),
    grand_mean = grand_mean,
    status = if (enough) "applied" else "insufficient_sources_not_applied",
    stringsAsFactors = FALSE
  )
}
