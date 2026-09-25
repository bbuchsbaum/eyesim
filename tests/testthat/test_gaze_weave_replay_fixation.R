# Replay revision 2026.10: fixation-level observation model (plan A2b).
# One hidden-state visit per recall fixation; each fixation emits a position
# and a duration. All data are simulated.

fixation_screen_w <- 1024
fixation_screen_h <- 768

fixation_rtrunc_t <- function(center, scale, degrees = 4) {
  repeat {
    point <- center + stats::rnorm(2) * scale /
      sqrt(stats::rgamma(1, degrees / 2, degrees / 2))
    if (point[[1]] >= 0 && point[[1]] <= fixation_screen_w &&
        point[[2]] >= 0 && point[[2]] <= fixation_screen_h) {
      return(point)
    }
  }
}

fixation_rtrunc_normal <- function(center, sd) {
  repeat {
    point <- center + stats::rnorm(2) * sd
    if (point[[1]] >= 0 && point[[1]] <= fixation_screen_w &&
        point[[2]] >= 0 && point[[2]] <= fixation_screen_h) {
      return(point)
    }
  }
}

fixation_path <- function(xy, duration) {
  make_gaze_fixations(xy, duration = duration,
                      onset = cumsum(c(0, duration[-length(duration)])))
}

fixation_true_duration <- list(
  replay_intercept = 1.2, replay_slope = 0.8, replay_sd = 0.2,
  background_mean = 5.3, background_sd = 0.45
)

# Exact generative model of the fixation-level Replay HMM: restart by
# encoding duration mass, screen-truncated Student replay emissions,
# log-normal durations with mean a + b log(encoding duration).
simulate_fixation_replay <- function(n_participants = 4L, n_items = 10L,
                                     n_fixations = 24L, n_reference = 6L,
                                     scale = 40,
                                     parameters = list(
                                       background = 0.15, restart = 0.05,
                                       advance = 0.5, background_stay = 0.7
                                     ),
                                     duration = fixation_true_duration,
                                     mode = c("signal", "uniform", "group"),
                                     seed = 1L) {
  mode <- match.arg(mode)
  set.seed(seed)
  centre <- c(fixation_screen_w, fixation_screen_h) / 2
  background_center <- lapply(seq_len(n_participants), function(p) {
    centre + c(stats::runif(1, -120, 120), stats::runif(1, -80, 80))
  })
  layout <- lapply(seq_len(n_items), function(k) {
    cbind(stats::runif(n_reference, 60, fixation_screen_w - 60),
          stats::runif(n_reference, 60, fixation_screen_h - 60))
  })
  group_points <- do.call(rbind, layout)
  rows <- expand.grid(item = seq_len(n_items),
                      participant = seq_len(n_participants))
  reference <- source <- vector("list", nrow(rows))
  background_share <- numeric(nrow(rows))
  for (r in seq_len(nrow(rows))) {
    xy <- layout[[rows$item[[r]]]] +
      matrix(stats::rnorm(2 * n_reference, 0, 5), n_reference)
    encoding_duration <- stats::rlnorm(n_reference, log(250), 0.5)
    reference[[r]] <- fixation_path(xy, encoding_duration)
    out <- matrix(NA_real_, n_fixations, 2)
    recall_duration <- numeric(n_fixations)
    if (mode == "signal") {
      mass <- encoding_duration / sum(encoding_duration)
      transition <- gaze_replay_transition(mass, parameters, 2L)
      state <- sample.int(n_reference + 1L, 1, prob = transition$initial)
      visited <- integer(n_fixations)
      for (t in seq_len(n_fixations)) {
        if (t > 1L) {
          state <- sample.int(n_reference + 1L, 1,
                              prob = transition$transition[state, ])
        }
        visited[[t]] <- state
        if (state == 1L) {
          out[t, ] <- fixation_rtrunc_normal(
            background_center[[rows$participant[[r]]]], c(150, 120)
          )
          recall_duration[[t]] <- stats::rlnorm(
            1, duration$background_mean, duration$background_sd
          )
        } else {
          out[t, ] <- fixation_rtrunc_t(xy[state - 1L, ], scale)
          recall_duration[[t]] <- stats::rlnorm(
            1, duration$replay_intercept + duration$replay_slope *
              log(encoding_duration[[state - 1L]]),
            duration$replay_sd
          )
        }
      }
      background_share[[r]] <- mean(visited == 1L)
    } else {
      for (t in seq_len(n_fixations)) {
        out[t, ] <- if (mode == "uniform") {
          c(stats::runif(1, 0, fixation_screen_w),
            stats::runif(1, 0, fixation_screen_h))
        } else {
          # Group density: gaze follows the pooled layout of every item,
          # independently of which item is being recalled.
          fixation_rtrunc_t(
            group_points[sample.int(nrow(group_points), 1), ], scale
          )
        }
      }
      recall_duration <- stats::rlnorm(n_fixations, log(230), 0.45)
      background_share[[r]] <- 1
    }
    source[[r]] <- fixation_path(out, recall_duration)
  }
  participant <- paste0("p", rows$participant)
  list(
    ref = tibble::tibble(participant = participant, image_id = rows$item,
                         fixgroup = reference),
    src = tibble::tibble(participant = participant, image_id = rows$item,
                         fixgroup = source),
    layout = layout,
    true_background = mean(background_share),
    true_scale = scale
  )
}

fixation_spec <- function(...) {
  suppressMessages(gaze_replay_spec(
    max_skip = 2L, student_df = 4, scale_floor = 6,
    transition_grid = list(
      background = c(0.03, 0.10), restart = c(0.02, 0.10),
      advance = c(0.25, 0.55), background_stay = c(0.85, 0.95)
    ),
    screen = gaze_screen(fixation_screen_w, fixation_screen_h),
    background_by = "participant", ...
  ))
}

fit_fixation_model <- function(data, spec = fixation_spec()) {
  suppressMessages(fit_gaze_replay_model(
    data$ref, data$src, match_on = c("participant", "image_id"), spec = spec
  ))
}

fixation_cache <- local({
  cache <- list()
  function(name, value) {
    if (is.null(cache[[name]])) cache[[name]] <<- value()
    cache[[name]]
  }
})

fixation_signal <- function() {
  fixation_cache("signal", function() simulate_fixation_replay(seed = 1L))
}
fixation_signal_model <- function() {
  fixation_cache("signal_model", function() {
    fit_fixation_model(fixation_signal())
  })
}
fixation_uniform <- function() {
  fixation_cache("uniform", function() {
    simulate_fixation_replay(n_items = 8L, n_fixations = 25L,
                             mode = "uniform", seed = 2L)
  })
}
fixation_uniform_model <- function() {
  fixation_cache("uniform_model", function() {
    fit_fixation_model(fixation_uniform())
  })
}

fixation_key <- function(value, column) {
  tab <- data.frame(value, stringsAsFactors = FALSE)
  names(tab) <- column
  gaze_key(tab, column)
}

test_that("each recall fixation is one observation with position and duration", {
  model <- fixation_signal_model()
  data <- fixation_signal()
  result <- gaze_replay_align(data$ref$fixgroup[[3]], data$src$fixgroup[[3]],
                              model, background_key = "p1")
  n_fixations <- nrow(data$src$fixgroup[[3]])

  expect_identical(nrow(result$alignment$posterior), n_fixations)
  expect_equal(result$alignment$grid$duration,
               data$src$fixgroup[[3]]$duration, tolerance = 0)
  # The score is the total trial log likelihood: no division by a grid.
  expect_equal(result$log_score, result$alignment$log_likelihood,
               tolerance = 0)
  expect_identical(result$provenance$score_semantics,
                   "total_trial_log_likelihood")
  expect_null(model$training$em$repeated_grid_rows)
})

test_that("spatial and duration emissions are normalised densities", {
  model <- fixation_signal_model()
  emission <- model$emission_models[[1]]
  nx <- 512L
  ny <- 384L
  gx <- (seq_len(nx) - 0.5) * fixation_screen_w / nx
  gy <- (seq_len(ny) - 0.5) * fixation_screen_h / ny
  points <- as.matrix(expand.grid(gx, gy))
  cell <- (fixation_screen_w / nx) * (fixation_screen_h / ny)
  centers <- rbind(
    c(fixation_screen_w / 2, fixation_screen_h / 2),
    c(1, fixation_screen_h / 2), c(1, 1),
    c(fixation_screen_w + 5, fixation_screen_h - 2)
  )
  spatial <- gaze_replay_truncated_t_log_density(
    points, centers, emission$replay_scale, model$spec$student_df,
    emission$screen
  )
  background <- gaze_replay_background_log_density(
    points, model$background, fixation_key("p1", "participant"),
    fixation_key(3L, "image_id")
  )$log_density
  expect_equal(unname(colSums(exp(cbind(background, spatial))) * cell),
               rep(1, 5), tolerance = 1e-3)

  encoding <- c(80, 250, 900)
  duration_integral <- vapply(seq_len(length(encoding) + 1L), function(j) {
    stats::integrate(function(d) {
      exp(gaze_replay_duration_log_density(
        log(d), log(encoding), emission$duration
      )[, j])
    }, 0, Inf, rel.tol = 1e-10)$value
  }, numeric(1))
  expect_equal(duration_integral, rep(1, length(encoding) + 1L),
               tolerance = 1e-6)
})

test_that("fixation forward-backward and Baum-Welch counts match a dense HMM", {
  data <- fixation_signal()
  model <- fixation_signal_model()
  reference <- as_gaze_measure(data$ref$fixgroup[[5]], gaze_local_order())
  prepared <- prepare_gaze_replay_source(
    data$src$fixgroup[[5]], model, NULL, fixation_key("p1", "participant"),
    fixation_key(5L, "image_id")
  )
  transition <- gaze_replay_transition(
    reference$mass, model$parameters, model$spec$max_skip
  )
  log_emission <- gaze_replay_emissions(
    reference, prepared$grid, model$emission_models[[1]],
    model$spec$student_df
  )
  fast <- gaze_replay_forward_backward(log_emission, transition)
  dense <- gaze_hmm_forward_backward(
    log_emission, transition$initial, transition$transition
  )

  expect_identical(nrow(log_emission), nrow(data$src$fixgroup[[5]]))
  expect_equal(fast$log_likelihood, dense$log_likelihood, tolerance = 1e-10)
  expect_lt(max(abs(fast$posterior - dense$posterior)), 1e-10)

  counts <- gaze_replay_expected_counts(log_emission, transition, fast)
  xi <- apply(dense$transition_posterior, c(2, 3), sum)
  total <- transition$transition
  share <- function(component) {
    out <- matrix(0, nrow(total), ncol(total))
    out[total > 0] <- component[total > 0] / total[total > 0]
    sum(xi * out)
  }
  expect_equal(counts$counts[["bb"]], xi[1, 1], tolerance = 1e-10)
  expect_equal(counts$counts[["replay_b"]], sum(xi[-1, 1]), tolerance = 1e-10)
  expect_equal(counts$counts[["replay_restart"]],
               share(transition$restart_component), tolerance = 1e-10)
  expect_equal(counts$replay_local, share(transition$local_component),
               tolerance = 1e-10)
})

test_that("training EM recovers spatial, occupancy and duration parameters", {
  data <- fixation_signal()
  model <- fixation_signal_model()
  duration <- model$emission_models[[1]]$duration
  trace <- model$training$em$log_likelihood_trace
  # The replay duration mean at the typical encoding duration is better
  # determined than the intercept alone.
  at_typical <- function(value) {
    value$replay_intercept + value$replay_slope * log(250)
  }

  expect_true(model$training$em$converged)
  expect_true(all(diff(trace) > -1e-8 * abs(trace[-1])))
  expect_equal(model$emission_models[[1]]$replay_scale, data$true_scale,
               tolerance = 0.15)
  expect_lt(abs(model$training$em$background_fraction -
                  data$true_background), 0.1)
  expect_lt(abs(duration$replay_slope -
                  fixation_true_duration$replay_slope), 0.15)
  expect_lt(abs(at_typical(duration) - at_typical(fixation_true_duration)),
            0.05)
  expect_equal(duration$replay_sd, fixation_true_duration$replay_sd,
               tolerance = 0.2)
  expect_lt(abs(duration$background_mean -
                  fixation_true_duration$background_mean), 0.1)
  expect_equal(duration$background_sd, fixation_true_duration$background_sd,
               tolerance = 0.15)
})

test_that("a uniform null at typical recall lengths is fitted as background", {
  model <- fixation_uniform_model()
  signal_model <- fixation_signal_model()
  signal <- fixation_signal()
  set.seed(31)
  candidates <- lapply(1:4, function(k) {
    fixation_path(cbind(stats::runif(6, 60, fixation_screen_w - 60),
                        stats::runif(6, 60, fixation_screen_h - 60)),
                  stats::rlnorm(6, log(250), 0.5))
  })
  null_spread <- vapply(1:10, function(i) {
    n <- sample(20:30, 1)
    recall <- fixation_path(
      cbind(stats::runif(n, 0, fixation_screen_w),
            stats::runif(n, 0, fixation_screen_h)),
      stats::rlnorm(n, log(230), 0.45)
    )
    score <- vapply(candidates, function(candidate) {
      gaze_replay_align(candidate, recall, model)$log_score
    }, numeric(1))
    diff(range(score))
  }, numeric(1))
  # Signal: true-candidate margin over the other items of the same
  # participant under the model fitted to replaying recalls.
  signal_margin <- vapply(1:10, function(r) {
    score <- vapply(1:10, function(k) {
      gaze_replay_align(signal$ref$fixgroup[[k]], signal$src$fixgroup[[r]],
                        signal_model)$log_score
    }, numeric(1))
    score[[r]] - max(score[-r])
  }, numeric(1))

  # A single uniform fixation can land on an encoding location, so a null
  # spread is not exactly zero; it must stay negligible against the signal.
  expect_gt(model$training$em$background_fraction, 0.8)
  expect_lt(stats::median(null_spread), 0.02 * stats::median(signal_margin))
  expect_lt(max(null_spread), 0.1 * stats::median(signal_margin))
})

test_that("grid_size is ignored under revision 2026.10, with a one-time message", {
  data <- simulate_fixation_replay(n_participants = 2L, n_items = 4L,
                                   n_fixations = 12L, seed = 4L)
  reset_gaze_replay_notices()
  expect_message(
    coarse <- gaze_replay_spec(
      grid_size = 16L, scale_floor = 6, background_by = "participant",
      screen = gaze_screen(fixation_screen_w, fixation_screen_h)
    ),
    "grid_size is ignored"
  )
  expect_silent(
    fine <- gaze_replay_spec(
      grid_size = 128L, scale_floor = 6, background_by = "participant",
      screen = gaze_screen(fixation_screen_w, fixation_screen_h)
    )
  )
  expect_silent(gaze_replay_spec(scale_floor = 6))
  messages <- character()
  first <- withCallingHandlers(
    fit_gaze_replay_model(data$ref, data$src, c("participant", "image_id"),
                          spec = coarse),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  second <- suppressMessages(fit_gaze_replay_model(
    data$ref, data$src, c("participant", "image_id"), spec = fine
  ))

  expect_false(any(grepl("repeat", messages)))
  expect_identical(first$parameters, second$parameters)
  expect_identical(first$emission_models, second$emission_models)
  fit <- suppressMessages(gaze_replay_cv(
    data$ref, data$src, match_on = c("participant", "image_id"),
    contrast_on = "participant", n_folds = 2, seed = 2, spec = coarse
  ))
  expect_null(fit$provenance$repeated_grid_rows)
  expect_identical(fit$provenance$score_semantics,
                   "total_trial_log_likelihood")
})

relabel_fixture <- function() {
  fixation_cache("relabel", function() {
    set.seed(11)
    w <- fixation_screen_w
    h <- fixation_screen_h
    layout <- lapply(1:2, function(k) {
      cbind(stats::runif(5, 80, w - 80), stats::runif(5, 80, h - 80))
    })
    rows <- expand.grid(image_id = 1:2, participant = paste0("p", 1:40),
                        stringsAsFactors = FALSE)
    rows <- do.call(rbind, lapply(split(rows, rows$participant), function(d) {
      d[sample(nrow(d), 2), ]
    }))
    clamp <- function(xy) pmin(pmax(xy, 1), w - 1)
    source <- lapply(rows$image_id, function(k) {
      xy <- clamp(layout[[k]][sample(5, 24, TRUE), ] +
                    matrix(stats::rnorm(48, 0, 30), 24))
      noise <- sample(24, 10)
      xy[noise, ] <- cbind(stats::runif(10, 1, w - 1),
                           stats::runif(10, 1, h - 1))
      fixation_path(xy, stats::rlnorm(24, log(230), 0.4))
    })
    reference <- lapply(rows$image_id, function(k) {
      fixation_path(layout[[k]] + stats::rnorm(10, 0, 5),
                    stats::rlnorm(5, log(250), 0.4))
    })
    model <- fit_fixation_model(list(
      ref = tibble::tibble(participant = rows$participant,
                           image_id = rows$image_id, fixgroup = reference),
      src = tibble::tibble(participant = rows$participant,
                           image_id = rows$image_id, fixgroup = source)
    ))
    ref_eval <- tibble::tibble(
      participant = "pnew", image_id = 1:2,
      fixgroup = lapply(layout, function(xy) {
        fixation_path(xy + stats::rnorm(10, 0, 5),
                      stats::rlnorm(5, log(250), 0.4))
      })
    )
    list(model = model, layout = layout, ref_eval = ref_eval)
  })
}

test_that("candidate scores do not depend on which candidate is labelled true", {
  fixture <- relabel_fixture()
  all_layout <- do.call(rbind, fixture$layout)
  score_as <- function(recall, target) {
    row <- tibble::tibble(participant = "pnew", image_id = target,
                          fixgroup = list(recall))
    score_gaze_replay_row(row, fixture$ref_eval, c("participant", "image_id"),
                          "participant", "fixgroup", "fixgroup", fixture$model)
  }
  set.seed(5)
  outcomes <- t(replicate(30, {
    xy <- all_layout[sample(nrow(all_layout), 24, TRUE), ] +
      matrix(stats::rnorm(48, 0, 30), 24)
    recall <- fixation_path(pmin(pmax(xy, 1), fixation_screen_w - 1),
                            stats::rlnorm(24, log(230), 0.4))
    first <- score_as(recall, 1L)
    second <- score_as(recall, 2L)
    c(
      invariant = max(abs(first$candidates$log_score -
                            second$candidates$log_score)),
      top1 = first$evidence$top1_credit
    )
  }))

  expect_lt(max(outcomes[, "invariant"]), 1e-12)
  expect_lt(abs(mean(outcomes[, "top1"]) - 0.5), 3 * sqrt(0.25 / 30))
})

# Legacy revision 2026.08 ------------------------------------------------------

# A full cross-fitted 2026.08 analysis. The golden fixture was produced by
# this function at commit 1934252 (before the fixation-level model); the
# frozen revision must reproduce it bit for bit.
legacy_replay_fixture_fit <- function() {
  set.seed(21)
  w <- 1024
  h <- 768
  fg <- function(xy, dur) {
    fixation_group(x = xy[, 1], y = xy[, 2], duration = dur,
                   onset = c(0, cumsum(dur)[-length(dur)]))
  }
  rows <- expand.grid(image_id = 1:6, participant = c("p1", "p2"),
                      stringsAsFactors = FALSE)
  ref <- lapply(seq_len(nrow(rows)), function(i) {
    fg(cbind(stats::runif(6, 50, w - 50), stats::runif(6, 50, h - 50)),
       stats::rgamma(6, 4, 4 / 250))
  })
  src <- lapply(ref, function(r) {
    fg(cbind(r$x[c(2, 3, 5)] + stats::rnorm(3, 0, 40),
             r$y[c(2, 3, 5)] + stats::rnorm(3, 0, 40)),
       stats::rgamma(3, 4, 4 / 250))
  })
  spec <- gaze_replay_spec(
    grid_size = 16L, max_skip = 2L, student_df = 4, scale_floor = 6,
    transition_grid = list(background = c(0.03, 0.10), restart = 0.05,
                           advance = c(0.25, 0.55), background_stay = 0.9),
    screen = gaze_screen(w, h), reliability = "effective_fixations",
    warp = gaze_warp_contraction(center = "screen", translation = TRUE),
    revision = "2026.08"
  )
  fit <- gaze_replay_cv(
    tibble::tibble(participant = rows$participant, image_id = rows$image_id,
                   fixgroup = ref),
    tibble::tibble(participant = rows$participant, image_id = rows$image_id,
                   fixgroup = src),
    match_on = c("participant", "image_id"), contrast_on = "participant",
    n_folds = 2, seed = 3, spec = spec
  )
  list(
    results = fit$results[, c("gaze_info_bits", "posterior_true", "log_loss",
                              "reliability", "top1_credit")],
    candidate_scores = lapply(fit$results$candidates, `[[`, "log_score"),
    posteriors = lapply(fit$results$alignment, `[[`, "posterior"),
    folds = lapply(fit$folds, function(fold) {
      fold[c("transition_parameters", "emission_models")]
    }),
    temperature = lapply(fit$folds, function(fold) {
      fold$calibration$temperature
    }),
    provenance = fit$provenance
  )
}

test_that("revision 2026.08 reproduces its frozen cross-fitted fit exactly", {
  golden <- readRDS(test_path("fixtures", "replay_legacy_2026_08.rds"))
  current <- legacy_replay_fixture_fit()

  expect_identical(current, golden)
})

# Slow checks (set EYESIM_SLOW_TESTS=true) ------------------------------------

test_that("a group-density null gives no candidate preference beyond MC error", {
  skip_unless_slow_tests("slow Replay null check")
  top1 <- bits <- chance <- numeric()
  for (seed in 1:4) {
    data <- simulate_fixation_replay(n_participants = 4L, n_items = 8L,
                                     n_fixations = 24L, mode = "group",
                                     seed = 40L + seed)
    fit <- suppressMessages(gaze_replay_cv(
      data$ref, data$src, match_on = c("participant", "image_id"),
      contrast_on = "participant", n_folds = 2, seed = seed,
      spec = fixation_spec()
    ))
    top1 <- c(top1, fit$results$top1_credit)
    bits <- c(bits, fit$results$gaze_info_bits)
    chance <- c(chance, 1 / fit$results$candidate_count)
  }
  # Binomial MC error of the top-1 rate around chance (1 / K per row).
  mc_error <- sqrt(sum(chance * (1 - chance))) / length(top1)

  expect_lt(abs(mean(top1) - mean(chance)), 3 * mc_error)
  # No overconfidence beyond MC error. When the Stein factor returns the
  # declared prior in every fold, every row has exactly zero bits and the
  # MC error is zero, so the boundary itself (equality) is the ideal outcome.
  expect_lte(mean(bits), 3 * stats::sd(bits) / sqrt(length(bits)))
})
