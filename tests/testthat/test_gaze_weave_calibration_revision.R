# Shared calibration, revision 2026.10 (plan A3): evidence-scaled temperature
# and the typicality offset, for Transport and Replay. All data are
# simulated. Heavy engine-level checks run only with EYESIM_SLOW_TESTS=true.

# Evaluate `code` under `seed`, then restore the caller's random stream.
calibration_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) {
    assign(".Random.seed", old, envir = globalenv())
  } else {
    rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  code
}

# Score-level simulation ------------------------------------------------------

# Candidate scores on a per-mass scale, as Transport produces them: the
# signal d does not grow with the recall's fixation count n, the noise sd
# shrinks as 1 / sqrt(n). The ideal inverse temperature is then proportional
# to n (gamma = 1). The true candidate is candidate 1.
calibration_score_rows <- function(n_rows, n_candidates, signal,
                                   evidence = c(2, 10), noise = 1) {
  lapply(seq_len(n_rows), function(i) {
    n <- evidence[[sample.int(length(evidence), 1)]]
    score <- stats::rnorm(n_candidates, 0, noise / sqrt(n))
    score[[1]] <- score[[1]] + signal
    list(profile = as.list(score), evidence = n, truth = 1L)
  })
}

fit_score_rows <- function(rows, ...) {
  eyesim:::fit_gaze_evidence_calibration(
    lapply(rows, `[[`, "profile"), vapply(rows, `[[`, integer(1), "truth"),
    vapply(rows, `[[`, numeric(1), "evidence"),
    control = gaze_calibration_control(...)
  )
}

score_rows <- function(rows, calibration) {
  evidence <- lapply(rows, function(row) {
    eyesim:::score_gaze_calibrated_row(
      row$profile, row$truth, row$evidence, calibration
    )
  })
  data.frame(
    evidence = vapply(rows, `[[`, numeric(1), "evidence"),
    bits = vapply(evidence, `[[`, numeric(1), "gaze_info_bits"),
    log_loss = vapply(evidence, `[[`, numeric(1), "log_loss"),
    top1 = vapply(evidence, `[[`, numeric(1), "top1_credit"),
    max_p = vapply(evidence, function(e) max(e$candidates$posterior),
                   numeric(1)),
    correct = vapply(evidence, function(e) {
      which.max(e$candidates$posterior) == 1L
    }, logical(1))
  )
}

calibration_cache <- local({
  cache <- list()
  function(name, value) {
    if (is.null(cache[[name]])) cache[[name]] <<- value()
    cache[[name]]
  }
})

sparse_rich_fixture <- function() {
  calibration_cache("sparse_rich", function() {
    calibration_with_seed(11, {
      train <- calibration_score_rows(600, 4, signal = 0.6)
      test <- calibration_score_rows(3000, 4, signal = 0.6)
      list(
        scaled = fit_score_rows(train),
        global = fit_score_rows(train, gamma_bounds = c(0, 0)),
        test = test
      )
    })
  })
}

test_that("specifications carry a calibration control only under 2026.10", {
  expect_null(gaze_transport_spec(revision = "2026.08")$calibration$control)
  expect_null(gaze_replay_spec(revision = "2026.08")$calibration$control)
  expect_error(
    gaze_transport_spec(revision = "2026.08",
                        calibration_control = gaze_calibration_control()),
    "requires revision"
  )
  expect_error(
    gaze_replay_spec(revision = "2026.08",
                     calibration_control = gaze_calibration_control()),
    "requires revision"
  )
  control <- gaze_transport_spec()$calibration$control
  expect_s3_class(control, "gaze_calibration_control")
  expect_identical(control$method, "evidence_scaled")
  expect_identical(control$typicality, "standardized")
  expect_identical(
    suppressMessages(gaze_replay_spec())$calibration$control$typicality,
    "none"
  )
  expect_null(gaze_calibration_control()$typicality)
  expect_identical(
    gaze_transport_spec(
      calibration_control = gaze_calibration_control(typicality = "none")
    )$calibration$control$typicality,
    "none"
  )
  expect_identical(gaze_transport_spec()$reliability, "none")
  expect_identical(gaze_transport_spec(revision = "2026.08")$reliability,
                   "effective_fixations")
  # The kappa shrink is superseded by the evidence-scaled temperature.
  expect_error(gaze_transport_spec(reliability = "effective_fixations"),
               "replaced")
  expect_error(gaze_replay_spec(reliability = "effective_fixations"),
               "replaced")
  expect_s3_class(
    gaze_replay_spec(
      reliability = "effective_fixations",
      calibration_control = gaze_calibration_control(method = "global")
    ),
    "gaze_replay_spec"
  )
  expect_error(gaze_calibration_control(method = "global", typicality = "mean"),
               "requires")
})

test_that("gamma = 0 gives one temperature for every evidence count", {
  fixture <- sparse_rich_fixture()
  global <- fixture$global
  beta <- eyesim:::gaze_calibration_row_beta(global, c(1, 2, 10, 50))
  expect_identical(global$gamma, 0)
  expect_true(all(beta == beta[[1]]))
  expect_equal(global$temperature, 1 / beta[[1]])
})

test_that("evidence scaling improves held-out log loss in sparse and rich strata", {
  fixture <- sparse_rich_fixture()
  scaled <- score_rows(fixture$test, fixture$scaled)
  global <- score_rows(fixture$test, fixture$global)

  expect_gt(fixture$scaled$gamma, 0.5)
  for (n in c(2, 10)) {
    stratum <- scaled$evidence == n
    gain <- global$log_loss[stratum] - scaled$log_loss[stratum]
    expect_gt(mean(gain), 2 * stats::sd(gain) / sqrt(sum(stratum)))
  }
  # Calibration only: within a row the ranking never changes.
  expect_identical(scaled$top1, global$top1)
})

test_that("evidence-scaled probabilities are calibrated by stratum and confidence", {
  fixture <- sparse_rich_fixture()
  scaled <- score_rows(fixture$test, fixture$scaled)
  # Tolerance 0.05: about three binomial standard errors of accuracy for
  # the smallest reported group (n >= 250).
  for (n in c(2, 10)) {
    stratum <- scaled[scaled$evidence == n, ]
    expect_lt(abs(mean(stratum$max_p) - mean(stratum$correct)), 0.05)
  }
  bins <- cut(scaled$max_p, c(0.25, 0.4, 0.55, 0.7, 0.85, 1),
              include.lowest = TRUE)
  for (bin in levels(bins)) {
    rows <- scaled[bins == bin, ]
    if (nrow(rows) < 250) next
    expect_lt(abs(mean(rows$max_p) - mean(rows$correct)), 0.05)
  }
  # One global temperature is miscalibrated in opposite directions.
  global <- score_rows(fixture$test, fixture$global)
  gap <- vapply(c(2, 10), function(n) {
    stratum <- global[global$evidence == n, ]
    mean(stratum$max_p) - mean(stratum$correct)
  }, numeric(1))
  expect_gt(gap[[1]], 0.05)
  expect_lt(gap[[2]], -0.05)
})

test_that("the inverse-temperature prior only pulls toward the declared prior", {
  calibration_with_seed(3, {
    rows <- calibration_score_rows(30, 3, signal = 0.8)
  })
  free <- fit_score_rows(rows, inverse_temperature_prior_sd = Inf,
                         gamma_bounds = c(0, 0))
  tight <- fit_score_rows(rows, inverse_temperature_prior_sd = 0.25,
                          gamma_bounds = c(0, 0))
  expect_lt(tight$inverse_temperature, free$inverse_temperature)
  expect_gt(tight$temperature, free$temperature)
  # Scores on a large scale (Replay's total log likelihood) do not change
  # the fitted probabilities: the prior is on the standardised scale.
  scaled_rows <- lapply(rows, function(row) {
    row$profile <- lapply(row$profile, `*`, 40)
    row
  })
  big <- fit_score_rows(scaled_rows, gamma_bounds = c(0, 0))
  small <- fit_score_rows(rows, gamma_bounds = c(0, 0))
  expect_equal(big$temperature, 40 * small$temperature, tolerance = 1e-3)
})

test_that("score-level nulls are neither above chance nor overconfident", {
  # Inner calibration sets as small as a Replay fold's (16 rows, pools of
  # two); held-out pools of four.
  outcome <- calibration_with_seed(5, {
    lapply(seq_len(80), function(rep) {
      train <- calibration_score_rows(16, 2, signal = 0)
      test <- calibration_score_rows(40, 4, signal = 0)
      list(
        stein = score_rows(test, fit_score_rows(train)),
        plain = score_rows(test, fit_score_rows(train,
                                                stein_shrinkage = FALSE))
      )
    })
  })
  stein <- do.call(rbind, lapply(outcome, `[[`, "stein"))
  plain <- do.call(rbind, lapply(outcome, `[[`, "plain"))
  mc_error <- sqrt(0.25 * 0.75 / nrow(stein))
  expect_lt(abs(mean(stein$top1) - 0.25), 3 * mc_error)
  expect_gte(mean(stein$bits), -0.03)
  # Without the Stein factor the noisy inner fit is overconfident.
  expect_lt(mean(plain$bits), mean(stein$bits))
})

test_that("the Stein factor leaves strong calibrations nearly untouched", {
  rows <- calibration_with_seed(8, calibration_score_rows(40, 4, signal = 1))
  fit <- fit_score_rows(rows)
  expect_gt(fit$stein_factor, 0.9)
  expect_equal(fit$standardized_inverse_temperature,
               fit$stein_factor *
                 fit$unshrunk_standardized_inverse_temperature)
  null <- calibration_with_seed(9, calibration_score_rows(16, 2, signal = 0))
  expect_lt(fit_score_rows(null)$stein_factor, 0.5)
})

test_that("rank and top-1 follow the ranking score when the prior is returned", {
  calibration <- eyesim:::gaze_evidence_calibration_fallback(
    gaze_calibration_control(), c(0.05, 100), "test"
  )
  evidence <- eyesim:::score_gaze_calibrated_row(
    list(0.1, 2, -1), true_index = 2L, evidence = 5, calibration = calibration
  )
  expect_equal(evidence$candidates$posterior, rep(1 / 3, 3))
  expect_identical(evidence$top1_credit, 1)
  expect_identical(evidence$template_rank, 1)
  expect_equal(evidence$gaze_info_bits, 0)
})

test_that("rank and top-1 do not depend on the inverse temperature", {
  # Two-episode candidates. As beta -> 0 the tempered log mean ranks by the
  # arithmetic episode mean (a: -5 < b: -3), at beta = 1 by the log mean
  # (a: log((1 + e^-10) / 2) ~ -0.69 > b: -3). Rank and top-1 must use one
  # fixed score whatever the fitted beta; only the probabilities change.
  profile <- list(a = c(0, -10), b = c(-3, -3))
  calibration_at <- function(beta) {
    list(inverse_temperature = beta, gamma = 0, evidence_reference = 1,
         temperature_bounds = NULL)
  }
  rows <- lapply(c(0, 1e-6, 0.5, 1), function(beta) {
    eyesim:::score_gaze_calibrated_row(profile, true_index = 2L, evidence = 3,
                                       calibration = calibration_at(beta))
  })
  for (row in rows) {
    expect_identical(row$template_rank, 2)
    expect_identical(row$top1_credit, 0)
    expect_identical(row$tied_candidates, 1L)
  }
  expect_equal(rows[[1]]$candidates$posterior, c(0.5, 0.5))
  expect_gt(rows[[4]]$candidates$posterior[[1]], 0.5)
  # Label independence: relabelling a as true gives the complementary rank.
  other <- eyesim:::score_gaze_calibrated_row(
    profile, true_index = 1L, evidence = 3, calibration = calibration_at(1e-6)
  )
  expect_identical(other$template_rank, 1)
  expect_identical(other$top1_credit, 1)
})

# Engine fixtures ------------------------------------------------------------

calibration_screen_w <- 30
calibration_screen_h <- 20

calibration_path <- function(x, y, duration = rep(250, length(x))) {
  make_gaze_fixations(cbind(x, y), duration = duration,
                      onset = cumsum(c(0, duration[-length(duration)] + 20)))
}

calibration_uniform_path <- function(n) {
  calibration_path(stats::runif(n, 1, calibration_screen_w - 1),
                   stats::runif(n, 1, calibration_screen_h - 1))
}

calibration_central_path <- function(n, sd) {
  calibration_path(
    pmin(pmax(stats::rnorm(n, calibration_screen_w / 2, sd), 0),
         calibration_screen_w),
    pmin(pmax(stats::rnorm(n, calibration_screen_h / 2, sd), 0),
         calibration_screen_h)
  )
}

calibration_peripheral_path <- function(n) {
  side <- sample.int(4, n, TRUE)
  x <- ifelse(side == 1, stats::runif(n, 0.5, 3),
              ifelse(side == 2, stats::runif(n, 27, 29.5),
                     stats::runif(n, 1, 29)))
  y <- ifelse(side == 3, stats::runif(n, 0.5, 3),
              ifelse(side == 4, stats::runif(n, 17, 19.5),
                     stats::runif(n, 1, 19)))
  calibration_path(x, y)
}

calibration_transport_spec <- function(...) {
  gaze_transport_spec(
    coverage_nodes = 2, entropy_schedule = 0.03, maxit = 40,
    tolerance = 1e-3, projection_maxit = 300, projection_tolerance = 1e-7,
    backend = "optimized", ...
  )
}

# Participants x items. Item k of every participant has its own encoding.
# `recall` is "signal" (a noisy subset of the encoding, `keep` fixations) or
# "centre" (independent centre-biased gaze). `layout` "mixed" makes items 1-2
# central, 3-4 uniform and 5-6 peripheral encodings.
simulate_calibration_transport <- function(n_participants, n_items,
                                           recall = c("signal", "centre"),
                                           layout = c("uniform", "mixed"),
                                           keep = c(2, 8), n_encoding = 8,
                                           noise = 2.5, seed = 1) {
  recall <- match.arg(recall)
  layout <- match.arg(layout)
  calibration_with_seed(seed, {
    rows <- expand.grid(item = seq_len(n_items),
                        participant = paste0("p", seq_len(n_participants)),
                        stringsAsFactors = FALSE)
    reference <- lapply(seq_len(nrow(rows)), function(r) {
      item <- rows$item[[r]]
      if (identical(layout, "mixed") && item <= 2) {
        calibration_central_path(n_encoding, 3)
      } else if (identical(layout, "mixed") && item >= 5) {
        calibration_peripheral_path(n_encoding)
      } else {
        calibration_uniform_path(n_encoding)
      }
    })
    source <- lapply(seq_len(nrow(rows)), function(r) {
      if (identical(recall, "centre")) {
        return(calibration_central_path(8, 5))
      }
      encoding <- reference[[r]]
      n <- keep[[sample.int(length(keep), 1)]]
      index <- sort(sample.int(nrow(encoding), n))
      calibration_path(encoding$x[index] + stats::rnorm(n, 0, noise),
                       encoding$y[index] + stats::rnorm(n, 0, noise))
    })
  })
  list(
    ref = tibble::tibble(participant = rows$participant, item = rows$item,
                         fixgroup = reference),
    src = tibble::tibble(participant = rows$participant, item = rows$item,
                         fixgroup = source)
  )
}

run_calibration_transport <- function(data, ..., n_folds = 2, seed = 1) {
  gaze_transport_cv(
    data$ref, data$src, match_on = c("participant", "item"),
    contrast_on = "participant", n_folds = n_folds, seed = seed,
    spec = calibration_transport_spec(...)
  )
}

transport_small <- function() {
  calibration_cache("transport_small", function() {
    simulate_calibration_transport(4, 4, seed = 2)
  })
}

test_that("Transport fits the calibration on inner rows and scales by evidence", {
  data <- transport_small()
  fit <- run_calibration_transport(data)
  results <- fit$results

  expect_true(all(c("temperature", "inverse_temperature", "evidence_count")
                  %in% names(results)))
  expect_equal(results$evidence_count, results$effective_fixations)
  for (fold in fit$folds) {
    calibration <- fold$calibration
    expect_s3_class(calibration, "gaze_evidence_calibration")
    # Inner rows are the outer training rows only.
    expect_identical(calibration$training_n, length(fold$train_rows))
    rows <- results$.cv_fold == fold$fold
    expected <- 1 / eyesim:::gaze_calibration_row_beta(
      calibration, results$evidence_count[rows]
    )
    expect_equal(results$temperature[rows], expected)
  }
  expect_true(is.data.frame(fit$calibration$fold_calibration))
  expect_match(fit$provenance$calibration, "evidence_scaled")
})

test_that("held-out labels never enter the Transport calibration", {
  data <- transport_small()
  fit <- run_calibration_transport(data)
  # Swap the recalls of two held-out rows of the same participant and fold:
  # their labels change, the training rows do not.
  fold <- fit$folds[[1]]
  held <- fold$eval_rows
  pair <- held[data$src$participant[held] == data$src$participant[held[[1]]]]
  expect_gte(length(pair), 2L)
  swapped <- data
  swapped$src$fixgroup[pair[1:2]] <- data$src$fixgroup[rev(pair[1:2])]
  refit <- run_calibration_transport(swapped)

  # Fold 1's calibration and typicality offsets use only its training rows.
  expect_identical(refit$folds[[1]]$calibration$inverse_temperature,
                   fold$calibration$inverse_temperature)
  expect_identical(refit$folds[[1]]$calibration$gamma, fold$calibration$gamma)
  expect_identical(refit$folds[[1]]$typicality, fold$typicality)
  expect_false(is.null(fold$typicality))
  # The swapped rows are training rows of fold 2, whose fit may change.
  rows <- match(pair[1:2], seq_len(nrow(data$src)))
  expect_identical(
    lapply(refit$results$candidates[rows], `[[`, "ranking_score"),
    rev(lapply(fit$results$candidates[rows], `[[`, "ranking_score"))
  )
})

test_that("a pre-A3 2026.10 specification reproduces its calibration exactly", {
  data <- transport_small()
  global <- calibration_transport_spec(
    reliability = "effective_fixations",
    calibration_control = gaze_calibration_control(method = "global")
  )
  saved <- global
  saved$calibration$control <- NULL
  strip <- function(fit) {
    fit$spec <- NULL
    fit$results$alignments <- NULL
    fit
  }
  run <- function(spec) {
    gaze_transport_cv(data$ref, data$src, match_on = c("participant", "item"),
                      contrast_on = "participant", n_folds = 2, spec = spec)
  }
  expect_identical(strip(run(saved)), strip(run(global)))
})

# Two study presentations per candidate with unrelated encodings, and
# random recalls: episode scores disagree, so a tempered log mean would rank
# differently at different temperatures.
simulate_two_episode_transport <- function(seed, n_participants = 3,
                                           n_items = 4) {
  calibration_with_seed(seed, {
    ref <- expand.grid(item = seq_len(n_items),
                       participant = paste0("p", seq_len(n_participants)),
                       presentation = 1:2, stringsAsFactors = FALSE)
    ref$fixgroup <- lapply(seq_len(nrow(ref)), function(i) {
      calibration_uniform_path(6)
    })
    src <- expand.grid(item = seq_len(n_items),
                       participant = paste0("p", seq_len(n_participants)),
                       stringsAsFactors = FALSE)
    src$fixgroup <- lapply(seq_len(nrow(src)), function(i) {
      calibration_uniform_path(sample(c(2, 8), 1))
    })
    list(ref = tibble::as_tibble(ref), src = tibble::as_tibble(src))
  })
}

test_that("multi-episode Transport top-1 and AUC do not depend on the calibration", {
  data <- simulate_two_episode_transport(2)
  run <- function(control, reliability = "none") {
    gaze_transport_cv(
      data$ref, data$src, match_on = c("participant", "item"),
      contrast_on = "participant", n_folds = 2, seed = 1,
      episode_on = "presentation",
      spec = calibration_transport_spec(reliability = reliability,
                                        calibration_control = control)
    )
  }
  none <- function(...) gaze_calibration_control(typicality = "none", ...)
  fits <- list(
    scaled = run(none()),
    unshrunk = run(none(stein_shrinkage = FALSE)),
    global = run(gaze_calibration_control(method = "global"),
                 reliability = "effective_fixations")
  )
  # A Stein factor of zero in every fold returns the declared prior.
  fit_original <- eyesim:::fit_gaze_evidence_calibration
  local_mocked_bindings(
    fit_gaze_evidence_calibration = function(...) {
      fit <- fit_original(...)
      fit$inverse_temperature <- 0
      fit$standardized_inverse_temperature <- 0
      fit$gamma <- 0
      fit$stein_factor <- 0
      fit
    },
    .package = "eyesim"
  )
  fits$prior <- run(none())
  expect_true(all(fits$prior$results$gaze_info_bits == 0))
  expect_true(any(fits$unshrunk$results$inverse_temperature > 0))
  # The probabilities differ across calibrations ...
  expect_false(isTRUE(all.equal(fits$unshrunk$results$posterior_true,
                                fits$prior$results$posterior_true)))
  # ... but rank, top-1 and the per-row AUC are identical.
  auc <- function(fit) {
    (fit$results$candidate_count - fit$results$template_rank) /
      (fit$results$candidate_count - 1)
  }
  for (fit in fits[-1]) {
    expect_identical(fit$results$template_rank, fits$scaled$results$template_rank)
    expect_identical(fit$results$top1_credit, fits$scaled$results$top1_credit)
    expect_identical(auc(fit), auc(fits$scaled))
  }
  expect_true(all(fits$scaled$results$common_episode_count == 2L))
})

test_that("typicality offsets use training recalls of other items only", {
  data <- simulate_calibration_transport(4, 6, recall = "centre",
                                         layout = "mixed", seed = 3)
  fit <- run_calibration_transport(
    data, calibration_control = gaze_calibration_control(
      typicality = "standardized", typicality_sources = 8
    )
  )
  table <- fit$calibration$typicality
  expect_true(all(table$status == "applied"))
  expect_true(all(table$n_sources >= 2L))
  expect_true(all(is.finite(table$source_sd)))
  # Standardized ranking score: (raw - row mean - (offset - grand)) * scale.
  # One study episode per candidate here, so the candidate score is linear.
  candidates <- fit$results$candidates[[1]]
  grand <- table$grand_mean[table$fold == fit$results$.cv_fold[[1]]][[1]]
  expected <- (candidates$raw_score - mean(candidates$raw_score) -
                 (candidates$typicality_offset - grand)) *
    candidates$typicality_scale
  expect_equal(candidates$ranking_score, expected)
  expect_true(all(candidates$typicality_scale > 0))
  expect_match(fit$provenance$calibration, "typicality")

  mean_fit <- run_calibration_transport(
    data, calibration_control = gaze_calibration_control(
      typicality = "mean", typicality_sources = 8
    )
  )
  candidates <- mean_fit$results$candidates[[1]]
  expect_equal(candidates$ranking_score,
               candidates$raw_score - candidates$typicality_offset)
  expect_true(all(candidates$typicality_scale == 1))
})

test_that("relabelling a held-out row leaves Transport candidate scores unchanged", {
  data <- simulate_calibration_transport(4, 6, recall = "centre",
                                         layout = "mixed", seed = 4)
  spec <- calibration_transport_spec(
    calibration_control = gaze_calibration_control(
      typicality = "standardized", typicality_sources = 8
    )
  )
  train <- data$src[data$src$participant != "p1", ]
  warp <- eyesim:::fit_transport_v3_warp(
    data$ref, train, c("participant", "item"), "fixgroup", "fixgroup",
    NULL, spec
  )
  keys <- gaze_key(data$ref[data$ref$participant == "p1", ],
                   c("participant", "item"))
  offsets <- eyesim:::transport_v3_typicality(
    keys, data$ref, train, c("participant", "item"), "participant",
    "fixgroup", "fixgroup", NULL, spec, warp
  )
  calibration <- structure(
    list(inverse_temperature = 3, gamma = 0.5, evidence_reference = 4),
    class = c("gaze_evidence_calibration", "list")
  )
  score_as <- function(item) {
    row <- data$src[data$src$participant == "p1" & data$src$item == 1, ]
    row$item <- item
    scored <- eyesim:::score_transport_v3_cv_row(
      row, data$ref, c("participant", "item"), "participant", "fixgroup",
      "fixgroup", NULL, NULL, spec, warp
    )
    eyesim:::transport_v3_calibrated_evidence(
      scored, calibration, offsets, "p1"
    )$evidence$candidates
  }
  first <- score_as(1L)
  second <- score_as(4L)
  expect_identical(first$ranking_score, second$ranking_score)
  expect_identical(first$log_score, second$log_score)
  expect_identical(which(first$is_true), 1L)
  expect_identical(which(second$is_true), 4L)
})

# Replay -------------------------------------------------------------------

calibration_replay_w <- 1024
calibration_replay_h <- 768

calibration_replay_path <- function(xy, duration) {
  make_gaze_fixations(xy, duration = duration,
                      onset = cumsum(c(0, duration[-length(duration)])))
}

calibration_replay_point <- function(center, scale) {
  repeat {
    point <- center + stats::rnorm(2) * scale
    if (point[[1]] >= 0 && point[[1]] <= calibration_replay_w &&
        point[[2]] >= 0 && point[[2]] <= calibration_replay_h) {
      return(point)
    }
  }
}

# Replay data. "signal": recalls revisit the item's encoding fixations with
# probability `rho`, otherwise participant background. "centre": recalls are
# independent centre-biased gaze. Layout "mixed": items with k %% 3 == 1 are
# central, k %% 3 == 0 peripheral, the rest uniform.
simulate_calibration_replay <- function(n_participants, n_items,
                                        recall = c("signal", "centre"),
                                        layout = c("uniform", "mixed"),
                                        n_fixations = c(8, 24), rho = 0.5,
                                        n_reference = 6, seed = 1) {
  recall <- match.arg(recall)
  layout <- match.arg(layout)
  centre <- c(calibration_replay_w, calibration_replay_h) / 2
  calibration_with_seed(seed, {
    rows <- expand.grid(item = seq_len(n_items),
                        participant = paste0("p", seq_len(n_participants)),
                        stringsAsFactors = FALSE)
    reference <- source <- vector("list", nrow(rows))
    for (r in seq_len(nrow(rows))) {
      kind <- if (identical(layout, "uniform")) 2L else rows$item[[r]] %% 3L
      xy <- switch(
        as.character(kind),
        "1" = t(replicate(n_reference, calibration_replay_point(centre, 90))),
        "0" = cbind(
          ifelse(stats::runif(n_reference) < 0.5,
                 stats::runif(n_reference, 20, 120),
                 stats::runif(n_reference, 904, 1004)),
          stats::runif(n_reference, 20, 748)
        ),
        cbind(stats::runif(n_reference, 60, 964),
              stats::runif(n_reference, 60, 708))
      )
      reference[[r]] <- calibration_replay_path(
        xy, stats::rlnorm(n_reference, log(250), 0.4)
      )
      n <- n_fixations[[sample.int(length(n_fixations), 1)]]
      out <- t(vapply(seq_len(n), function(t) {
        if (identical(recall, "signal") && stats::runif(1) < rho) {
          calibration_replay_point(xy[sample.int(n_reference, 1), ], 40)
        } else {
          calibration_replay_point(centre, 180)
        }
      }, numeric(2)))
      source[[r]] <- calibration_replay_path(
        out, stats::rlnorm(n, log(230), 0.4)
      )
    }
  })
  list(
    ref = tibble::tibble(participant = rows$participant, image_id = rows$item,
                         fixgroup = reference),
    src = tibble::tibble(participant = rows$participant, image_id = rows$item,
                         fixgroup = source)
  )
}

calibration_replay_spec <- function(...) {
  suppressMessages(gaze_replay_spec(
    max_skip = 2L, student_df = 4, scale_floor = 6,
    transition_grid = list(
      background = c(0.03, 0.10), restart = c(0.02, 0.10),
      advance = c(0.25, 0.55), background_stay = c(0.85, 0.95)
    ),
    screen = gaze_screen(calibration_replay_w, calibration_replay_h),
    background_by = "participant", ...
  ))
}

run_calibration_replay <- function(data, ..., seed = 1) {
  suppressMessages(gaze_replay_cv(
    data$ref, data$src, match_on = c("participant", "image_id"),
    contrast_on = "participant", n_folds = 2, seed = seed,
    spec = calibration_replay_spec(...)
  ))
}

test_that("Replay fits an evidence-scaled calibration on recall fixation counts", {
  data <- simulate_calibration_replay(3, 8, seed = 6)
  fit <- run_calibration_replay(data)
  results <- fit$results

  expect_equal(results$evidence_count, as.numeric(results$raw_fixation_count))
  for (fold in fit$folds) {
    expect_s3_class(fold$calibration, "gaze_evidence_calibration")
    expect_identical(fold$calibration$evidence_measure, "recall_fixations")
    rows <- results$.cv_fold == fold$fold
    expected <- 1 / eyesim:::gaze_calibration_row_beta(
      fold$calibration, results$evidence_count[rows]
    )
    expect_equal(results$temperature[rows], expected)
  }
  expect_match(fit$provenance$temperature, "evidence-scaled")
})

test_that("relabelling a held-out row leaves Replay candidate scores unchanged", {
  data <- simulate_calibration_replay(3, 8, recall = "centre",
                                      layout = "mixed", seed = 7)
  spec <- calibration_replay_spec(
    calibration_control = gaze_calibration_control(typicality = "standardized")
  )
  held <- data$src$participant == "p1" & data$src$image_id <= 4
  train <- data$src[!held, ]
  ref_train <- data$ref[!(data$ref$participant == "p1" &
                            data$ref$image_id <= 4), ]
  ref_eval <- data$ref[data$ref$participant == "p1" &
                         data$ref$image_id <= 4, ]
  model <- suppressMessages(fit_gaze_replay_model(
    ref_train, train, c("participant", "image_id"), "participant",
    spec = spec
  ))
  offsets <- eyesim:::replay_typicality(
    model, ref_eval, train, c("participant", "image_id"), "participant",
    "fixgroup", "fixgroup"
  )
  expect_true(all(offsets$table$status == "applied"))
  row <- data$src[data$src$participant == "p1" & data$src$image_id == 1, ]
  score_as <- function(item) {
    relabelled <- row
    relabelled$image_id <- item
    eyesim:::score_gaze_replay_row(
      relabelled, ref_eval, c("participant", "image_id"), "participant",
      "fixgroup", "fixgroup", model, typicality = offsets
    )$candidates
  }
  first <- score_as(1L)
  second <- score_as(3L)
  expect_identical(first$ranking_score, second$ranking_score)
  expect_identical(first$calibrated_logit, second$calibrated_logit)
  expect_identical(which(first$is_true), 1L)
  expect_identical(which(second$is_true), 3L)
})

# Pre-change reference fits ---------------------------------------------------
#
# fixtures/calibration_base_0d9fbcf.rds was produced by
# calibration_base_fixture() below with the package at master commit 0d9fbcf
# (before A3), passing that commit's default specifications:
#   transport_08 = calibration_base_transport_spec(revision = "2026.08")
#   transport_10 = calibration_base_transport_spec()
#   replay_10    = calibration_replay_spec()
# Transport uses the pure-R reference backend so that the stored numbers do
# not depend on how the native code was compiled.

calibration_base_transport_spec <- function(...) {
  gaze_transport_spec(
    coverage_nodes = 2, entropy_schedule = 0.03, maxit = 40,
    tolerance = 1e-3, projection_maxit = 300, projection_tolerance = 1e-7,
    backend = "reference", ...
  )
}

calibration_base_fixture <- function(transport_08, transport_10, replay_10) {
  transport_data <- simulate_calibration_transport(2, 3, n_encoding = 4, keep = c(2, 4),
                                                   seed = 81)
  replay_data <- simulate_calibration_replay(3, 8, n_fixations = c(6, 10),
                                             seed = 82)
  run_transport <- function(spec) {
    fit <- gaze_transport_cv(
      transport_data$ref, transport_data$src,
      match_on = c("participant", "item"), contrast_on = "participant",
      n_folds = 2, seed = 3, spec = spec
    )
    list(
      results = fit$results[, setdiff(names(fit$results), "alignments")],
      calibration = lapply(fit$folds, function(fold) {
        fold$calibration[c("temperature", "kappa", "log_loss")]
      }),
      summary = fit$calibration,
      provenance = fit$provenance
    )
  }
  replay_fit <- suppressMessages(gaze_replay_cv(
    replay_data$ref, replay_data$src,
    match_on = c("participant", "image_id"), contrast_on = "participant",
    n_folds = 2, seed = 4, spec = replay_10
  ))
  list(
    transport_08 = run_transport(transport_08),
    transport_10 = run_transport(transport_10),
    replay_10 = list(
      results = replay_fit$results[, setdiff(names(replay_fit$results),
                                             "alignment")],
      calibration = lapply(replay_fit$folds, function(fold) {
        fold$calibration[c("temperature", "log_loss")]
      }),
      provenance = replay_fit$provenance
    )
  )
}

test_that("revision 2026.08 and the global method reproduce pre-A3 fits exactly", {
  golden <- readRDS(test_path("fixtures", "calibration_base_0d9fbcf.rds"))
  global <- gaze_calibration_control(method = "global")
  current <- calibration_base_fixture(
    transport_08 = calibration_base_transport_spec(revision = "2026.08"),
    transport_10 = calibration_base_transport_spec(
      reliability = "effective_fixations", calibration_control = global
    ),
    replay_10 = calibration_replay_spec(calibration_control = global)
  )

  expect_identical(current$transport_08, golden$transport_08)
  expect_identical(current$transport_10, golden$transport_10)
  expect_identical(current$replay_10, golden$replay_10)
  expect_null(gaze_transport_spec(revision = "2026.08")$calibration$control)
  expect_null(gaze_replay_spec(revision = "2026.08")$calibration$control)
})

# Slow checks (set EYESIM_SLOW_TESTS=true) ------------------------------------

skip_unless_slow_calibration <- function() {
  skip_unless_slow_tests("slow calibration check")
}

# Pool held-out rows over seeds. `top1` credit against chance 1/K.
pool_null <- function(fits) {
  results <- do.call(rbind, lapply(fits, function(fit) {
    fit$results[, c("top1_credit", "gaze_info_bits", "candidate_count")]
  }))
  chance <- 1 / results$candidate_count
  list(
    top1 = mean(results$top1_credit), chance = mean(chance),
    mc_error = sqrt(sum(chance * (1 - chance))) / nrow(results),
    bits = mean(results$gaze_info_bits),
    bits_se = stats::sd(results$gaze_info_bits) / sqrt(nrow(results)),
    n = nrow(results)
  )
}

# Share of held-out rows whose argmax candidate has a central encoding,
# per central candidate (chance 1 / K).
# The argmax uses the (offset-adjusted) ranking score, which is also what
# rank and top-1 use under every calibration.
central_share <- function(fits, is_central) {
  wins <- unlist(lapply(fits, function(fit) {
    vapply(fit$results$candidates, function(candidates) {
      score <- if (is.null(candidates$ranking_score)) {
        candidates$log_score
      } else {
        candidates$ranking_score
      }
      top <- which.max(score)
      is_central(candidates[top, , drop = FALSE])
    }, logical(1))
  }))
  n_central <- unlist(lapply(fits, function(fit) {
    vapply(fit$results$candidates, function(candidates) {
      sum(vapply(seq_len(nrow(candidates)), function(i) {
        is_central(candidates[i, , drop = FALSE])
      }, logical(1)))
    }, integer(1))
  }))
  k <- unlist(lapply(fits, function(fit) fit$results$candidate_count))
  p0 <- n_central / k
  list(share_per_candidate = mean(wins) / mean(n_central),
       chance_per_candidate = mean(1 / k),
       mc_error = sqrt(sum(p0 * (1 - p0))) / length(wins) / mean(n_central),
       wins = mean(wins), expected = mean(p0))
}

test_that("Transport nulls are at chance and not overconfident", {
  skip_unless_slow_calibration()
  fits <- lapply(1:4, function(seed) {
    run_calibration_transport(
      simulate_calibration_transport(8, 6, recall = "centre", seed = 20 + seed),
      seed = seed
    )
  })
  null <- pool_null(fits)
  expect_lt(abs(null$top1 - null$chance), 3 * null$mc_error)
  expect_gte(null$bits, -0.03)
})

test_that("Replay nulls are at chance and not overconfident", {
  skip_unless_slow_calibration()
  fits <- lapply(1:4, function(seed) {
    run_calibration_replay(
      simulate_calibration_replay(8, 8, recall = "centre", seed = 30 + seed),
      seed = seed
    )
  })
  null <- pool_null(fits)
  expect_lt(abs(null$top1 - null$chance), 3 * null$mc_error)
  expect_gte(null$bits, -0.03)
})

test_that("Transport evidence scaling improves held-out log loss on sparse and rich recalls", {
  skip_unless_slow_calibration()
  data <- simulate_calibration_transport(10, 6, keep = c(2, 8), seed = 41)
  scaled <- run_calibration_transport(data)
  global <- run_calibration_transport(
    data, calibration_control = gaze_calibration_control(gamma_bounds = c(0, 0))
  )
  gain <- global$results$log_loss - scaled$results$log_loss
  expect_gt(mean(gain), 0)
  expect_identical(scaled$results$top1_credit, global$results$top1_credit)
})

test_that("the typicality offset removes the central-candidate advantage", {
  skip_unless_slow_calibration()
  # Transport candidate keys encode the item last ("...|<len>:<item>").
  central <- function(candidate) {
    as.integer(sub(".*:", "", candidate$candidate_key)) <= 2L
  }
  run <- function(typicality) {
    lapply(1:3, function(seed) {
      run_calibration_transport(
        simulate_calibration_transport(8, 6, recall = "centre",
                                       layout = "mixed", seed = 50 + seed),
        seed = seed,
        calibration_control = gaze_calibration_control(typicality = typicality)
      )
    })
  }
  without <- central_share(run("none"), central)
  standardized <- central_share(run("standardized"), central)
  # Without an offset central candidates win far above chance.
  expect_gt(without$share_per_candidate - without$chance_per_candidate,
            3 * without$mc_error)
  expect_lt(abs(standardized$share_per_candidate -
                  standardized$chance_per_candidate),
            3 * standardized$mc_error)
})

# Known limitation, pinned so that a fix is noticed: for Replay the offset
# lowers the central share (0.30 -> 0.27 in the A3 study) but does not bring
# it within MC error of 1 / K. Under these nulls the HMM background absorbs
# most recalls and most pools are near-ties (57% of rows within 0.01 nats);
# the residual bias is the argmax among near-ties, which no offset removes.
# The offset is therefore opt-in for Replay.
test_that("the Replay typicality offset reduces but does not remove central bias", {
  skip_unless_slow_calibration()
  central <- function(candidate) candidate$image_id %% 3 == 1
  run <- function(typicality) {
    lapply(1:3, function(seed) {
      run_calibration_replay(
        simulate_calibration_replay(6, 12, recall = "centre",
                                    layout = "mixed", seed = 60 + seed),
        seed = seed,
        calibration_control = gaze_calibration_control(typicality = typicality)
      )
    })
  }
  without <- central_share(run("none"), central)
  standardized <- central_share(run("standardized"), central)
  expect_lt(standardized$share_per_candidate, without$share_per_candidate)
  expect_gt(standardized$share_per_candidate -
              standardized$chance_per_candidate,
            3 * standardized$mc_error)
})

test_that("the typicality offset does not reduce discrimination on signal data", {
  skip_unless_slow_calibration()
  data <- simulate_calibration_transport(8, 6, keep = c(2, 8), seed = 71)
  off <- run_calibration_transport(
    data, calibration_control = gaze_calibration_control(typicality = "none")
  )
  on <- run_calibration_transport(data)  # Transport default: "standardized"
  expect_identical(on$spec$calibration$control$typicality, "standardized")
  difference <- on$results$top1_credit - off$results$top1_credit
  expect_gt(mean(difference),
            -2 * stats::sd(difference) / sqrt(length(difference)))

  replay <- simulate_calibration_replay(6, 12, seed = 72)
  replay_off <- run_calibration_replay(replay)
  replay_on <- run_calibration_replay(
    replay, calibration_control = gaze_calibration_control(typicality = "standardized")
  )
  replay_difference <- replay_on$results$top1_credit -
    replay_off$results$top1_credit
  expect_gt(mean(replay_difference),
            -2 * stats::sd(replay_difference) /
              sqrt(length(replay_difference)))
})
