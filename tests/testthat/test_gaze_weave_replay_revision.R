# Replay revision 2026.10: screen normalisation, background, and EM fitting.
# All data are simulated from the Replay generative model itself. The
# fixation-level observation model (plan A2b) has its own file,
# test_gaze_weave_replay_fixation.R.

replay_revision_under_test <- "2026.10"
revision_screen_w <- 1024
revision_screen_h <- 768

revision_rtrunc_t <- function(center, scale, degrees = 4) {
  repeat {
    point <- center + stats::rnorm(2) * scale /
      sqrt(stats::rgamma(1, degrees / 2, degrees / 2))
    if (point[[1]] >= 0 && point[[1]] <= revision_screen_w &&
        point[[2]] >= 0 && point[[2]] <= revision_screen_h) {
      return(point)
    }
  }
}

revision_rtrunc_background <- function(center, sd) {
  repeat {
    point <- center + stats::rnorm(2) * sd
    if (point[[1]] >= 0 && point[[1]] <= revision_screen_w &&
        point[[2]] >= 0 && point[[2]] <= revision_screen_h) {
      return(point)
    }
  }
}

revision_path <- function(xy) {
  make_gaze_fixations(xy, duration = rep(1, nrow(xy)),
                      onset = seq_len(nrow(xy)) - 1)
}

# One recall fixation per hidden-state visit (equal durations), so the
# fixation-level Replay HMM is the exact generative model of the recalls.
simulate_replay_hmm <- function(n_participants = 4L, n_items = 12L,
                                n_bins = 32L, n_reference = 6L, scale = 40,
                                parameters = list(
                                  background = 0.15, restart = 0.05,
                                  advance = 0.5, background_stay = 0.7
                                ),
                                null = FALSE, seed = 1L) {
  set.seed(seed)
  centre <- c(revision_screen_w, revision_screen_h) / 2
  background_center <- lapply(seq_len(n_participants), function(p) {
    centre + c(stats::runif(1, -120, 120), stats::runif(1, -80, 80))
  })
  rows <- expand.grid(item = seq_len(n_items),
                      participant = seq_len(n_participants))
  reference <- source <- vector("list", nrow(rows))
  background_share <- numeric(nrow(rows))
  mass <- rep(1 / n_reference, n_reference)
  transition <- gaze_replay_transition(mass, parameters, 2L)
  for (r in seq_len(nrow(rows))) {
    xy <- cbind(
      stats::runif(n_reference, 60, revision_screen_w - 60),
      stats::runif(n_reference, 60, revision_screen_h - 60)
    )
    reference[[r]] <- revision_path(xy)
    state <- if (null) 1L else
      sample.int(n_reference + 1L, 1, prob = transition$initial)
    out <- matrix(NA_real_, n_bins, 2)
    visited <- integer(n_bins)
    for (t in seq_len(n_bins)) {
      if (t > 1L && !null) {
        state <- sample.int(n_reference + 1L, 1,
                            prob = transition$transition[state, ])
      }
      visited[[t]] <- state
      out[t, ] <- if (state == 1L) {
        revision_rtrunc_background(
          background_center[[rows$participant[[r]]]], c(150, 120)
        )
      } else {
        revision_rtrunc_t(xy[state - 1L, ], scale)
      }
    }
    background_share[[r]] <- mean(visited == 1L)
    source[[r]] <- revision_path(out)
  }
  participant <- paste0("p", rows$participant)
  list(
    ref = tibble::tibble(participant = participant, image_id = rows$item,
                         fixgroup = reference),
    src = tibble::tibble(participant = participant, image_id = rows$item,
                         fixgroup = source),
    true_background = mean(background_share),
    true_scale = scale
  )
}

make_revision_spec <- function(revision = replay_revision_under_test,
                               n_bins = 32L) {
  args <- list(
    grid_size = n_bins, max_skip = 2L, student_df = 4, scale_floor = 6,
    transition_grid = list(
      background = c(0.03, 0.10), restart = c(0.02, 0.10),
      advance = c(0.25, 0.55), background_stay = c(0.85, 0.95)
    ),
    screen = gaze_screen(revision_screen_w, revision_screen_h),
    revision = revision
  )
  if (identical(revision, "2026.10")) args$background_by <- "participant"
  do.call(gaze_replay_spec, args)
}

fit_revision_model <- function(data, revision = replay_revision_under_test) {
  fit_gaze_replay_model(
    data$ref, data$src, match_on = c("participant", "image_id"),
    spec = make_revision_spec(revision)
  )
}

# Posterior share of observations assigned to the background state on the
# training rows. Revision 2026.10 records the out-of-fold EM value; the frozen
# revision is summarised with its own alignment posteriors.
fitted_background_fraction <- function(model, data) {
  if (!is.null(model$training$em)) {
    return(model$training$em$background_fraction)
  }
  mean(vapply(seq_len(nrow(data$src)), function(i) {
    gaze_replay_align(data$ref$fixgroup[[i]], data$src$fixgroup[[i]],
                      model)$diagnostics$background_coverage
  }, numeric(1)))
}

revision_key <- function(value, column) {
  tab <- data.frame(value, stringsAsFactors = FALSE)
  names(tab) <- column
  gaze_key(tab, column)
}

revision_emission_grid <- function(model, points) {
  grid <- list(coords = points, log_duration = rep(0, nrow(points)))
  if (identical(model$revision, "2026.10")) {
    grid <- gaze_replay_prepare_background(
      grid, model, revision_key("p1", "participant"),
      revision_key(3L, "image_id")
    )
  }
  grid
}

revision_cache <- local({
  cache <- list()
  function(name, value) {
    if (is.null(cache[[name]])) cache[[name]] <<- value()
    cache[[name]]
  }
})

signal_data <- function() {
  revision_cache("signal", function() simulate_replay_hmm(seed = 1L))
}
signal_model <- function() {
  revision_cache("signal_model", function() fit_revision_model(signal_data()))
}

test_that("every Replay emission state integrates to one over the screen", {
  skip_on_cran()
  model <- signal_model()
  nx <- 512L
  ny <- 384L
  gx <- (seq_len(nx) - 0.5) * revision_screen_w / nx
  gy <- (seq_len(ny) - 0.5) * revision_screen_h / ny
  points <- as.matrix(expand.grid(gx, gy))
  cell <- (revision_screen_w / nx) * (revision_screen_h / ny)
  # Centre, edge midpoint, corner, and a fixation recorded just off screen.
  reference <- list(
    coords = rbind(
      c(revision_screen_w / 2, revision_screen_h / 2),
      c(1, revision_screen_h / 2),
      c(1, 1),
      c(revision_screen_w + 5, revision_screen_h - 2)
    ),
    fixations = data.frame(duration = c(1, 2, 3, 4))
  )
  grid <- revision_emission_grid(model, points)
  log_emission <- gaze_replay_emissions(
    reference, grid, model$emission_models[[1]], model$spec$student_df
  )
  # Every point has the same duration, so dividing out each state's duration
  # density leaves its spatial density, which must integrate to one.
  duration <- gaze_replay_duration_log_density(
    grid$log_duration[[1]], log(reference$fixations$duration),
    model$emission_models[[1]]$duration
  )
  integral <- unname(colSums(exp(sweep(log_emission, 2, duration))) * cell)

  expect_equal(integral, rep(1, ncol(log_emission)), tolerance = 1e-3)
})

test_that("edge-heavy true candidates are not penalised against central ones", {
  model <- signal_model()
  shape <- rbind(c(0, 0), c(70, 15), c(25, 80), c(95, 95))
  placements <- list(
    corner = sweep(shape, 2, c(2, 2), "+"),
    edge = sweep(shape, 2, c(2, revision_screen_h / 2 - 45), "+"),
    centre = sweep(shape, 2, c(revision_screen_w / 2 - 45,
                               revision_screen_h / 2 - 45), "+")
  )
  rect <- list(xlim = c(0, revision_screen_w),
               ylim = c(0, revision_screen_h),
               area = revision_screen_w * revision_screen_h)
  state <- rep(seq_len(nrow(shape)), each = 8L)
  # Per-fixation shortfall of the model score (a total log likelihood) from
  # the exact screen-truncated generating log density of the replayed
  # recall. Identical relative geometry must give the same shortfall
  # wherever the candidate sits.
  shortfall <- function(candidate) {
    recall <- t(vapply(state, function(s) {
      revision_rtrunc_t(candidate[s, ], 40)
    }, numeric(2)))
    generating <- gaze_replay_truncated_t_log_density(
      recall, candidate, 40, 4, rect
    )[cbind(seq_along(state), state)]
    (sum(generating) - gaze_replay_align(
      revision_path(candidate), revision_path(recall), model
    )$log_score) / length(state)
  }
  set.seed(9)
  values <- t(replicate(20, vapply(placements, shortfall, numeric(1))))
  mean_shortfall <- colMeans(values)

  expect_lt(abs(mean_shortfall[["corner"]] - mean_shortfall[["centre"]]), 0.1)
  expect_lt(abs(mean_shortfall[["edge"]] - mean_shortfall[["centre"]]), 0.1)
})

test_that("training fits recover the simulated replay scale and background share", {
  data <- signal_data()
  model <- signal_model()

  expect_equal(model$emission_models[[1]]$replay_scale, data$true_scale,
               tolerance = 0.15)
  expect_lt(abs(fitted_background_fraction(model, data) -
                  data$true_background), 0.1)
})

null_data <- function() {
  revision_cache("null", function() {
    simulate_replay_hmm(null = TRUE, n_items = 8L, seed = 2L)
  })
}
null_model <- function() {
  revision_cache("null_model", function() fit_revision_model(null_data()))
}

test_that("a full null is fitted as background", {
  skip_on_cran()
  data <- null_data()
  model <- null_model()
  # Long-run background share implied by the fitted two-level chain
  # (background <-> replay); "no replay" must be representable.
  stationary <- with(model$parameters, background /
                       (background + 1 - background_stay))

  expect_gt(fitted_background_fraction(model, data), 0.8)
  expect_gt(stationary, 0.8)
})

test_that("forward-backward and Baum-Welch counts match a dense HMM", {
  data <- signal_data()
  model <- signal_model()
  reference <- as_gaze_measure(data$ref$fixgroup[[5]], gaze_local_order())
  source <- as_gaze_measure(data$src$fixgroup[[5]], gaze_local_order())
  grid <- revision_emission_grid(model, source$coords)
  transition <- gaze_replay_transition(
    reference$mass, model$parameters, model$spec$max_skip
  )
  log_emission <- gaze_replay_emissions(
    reference, grid, model$emission_models[[1]], model$spec$student_df
  )
  fast <- gaze_replay_forward_backward(log_emission, transition)
  dense <- gaze_hmm_forward_backward(
    log_emission, transition$initial, transition$transition
  )

  expect_equal(fast$log_likelihood, dense$log_likelihood, tolerance = 1e-10)
  expect_lt(max(abs(fast$posterior - dense$posterior)), 1e-10)

  counts <- gaze_replay_expected_counts(log_emission, transition, fast)
  xi <- apply(dense$transition_posterior, c(2, 3), sum)
  local <- transition$local_component
  restart <- transition$restart_component
  total <- transition$transition
  share <- function(component) {
    out <- matrix(0, nrow(total), ncol(total))
    out[total > 0] <- component[total > 0] / total[total > 0]
    sum(xi * out)
  }
  expect_equal(counts$counts[["bb"]], xi[1, 1], tolerance = 1e-10)
  expect_equal(counts$counts[["b_replay"]], sum(xi[1, -1]), tolerance = 1e-10)
  expect_equal(counts$counts[["replay_b"]], sum(xi[-1, 1]), tolerance = 1e-10)
  expect_equal(counts$counts[["replay_restart"]], share(restart),
               tolerance = 1e-10)
  expect_equal(counts$replay_local, share(local), tolerance = 1e-10)
})

test_that("revision 2026.10 EM increases the training likelihood monotonically", {
  model <- signal_model()
  trace <- model$training$em$log_likelihood_trace

  expect_true(model$training$em$converged)
  expect_true(all(diff(trace) > -1e-8 * abs(trace[-1])))
  expect_identical(model$revision, "2026.10")
  expect_identical(model$version, 5L)
  expect_null(model$training$transition_candidates)
})

test_that("backgrounds are participant-level and exclude the target item", {
  model <- signal_model()
  points <- rbind(c(300, 300), c(500, 400), c(900, 700))
  pooled <- gaze_replay_background_log_density(points, model$background)
  participant <- revision_key("p1", "participant")
  own <- gaze_replay_background_log_density(
    points, model$background, participant
  )
  held_out <- gaze_replay_background_log_density(
    points, model$background, participant, revision_key(3L, "image_id")
  )
  source <- revision_path(rbind(c(300, 300), c(500, 400)))
  keyed <- gaze_replay_align(
    revision_path(rbind(c(310, 305), c(505, 395))), source, model,
    background_key = "p1"
  )

  expect_identical(pooled$level, "pooled")
  expect_identical(own$level, "participant")
  expect_identical(held_out$trials, own$trials - 1L)
  expect_false(isTRUE(all.equal(own$log_density, held_out$log_density)))
  expect_setequal(
    unique(model$background$table$exclude_key),
    revision_key(seq_len(12), "image_id")
  )
  expect_identical(keyed$alignment$grid$background_level, "participant")
})

test_that("revision 2026.10 cross-fitting scores every held-out row", {
  data <- simulate_replay_hmm(n_participants = 2L, n_items = 8L,
                              n_bins = 16L, seed = 3L)
  spec <- make_revision_spec(n_bins = 16L)
  fit <- gaze_replay_cv(
    data$ref, data$src, match_on = c("participant", "image_id"),
    contrast_on = "participant", n_folds = 2, seed = 5, spec = spec
  )

  expect_true(all(is.finite(fit$results$gaze_info_bits)))
  expect_true(all(fit$results$all_converged))
  expect_identical(fit$provenance$revision, "2026.10")
  expect_gt(mean(fit$results$top1_credit), 0.5)
})

test_that("revision 2026.08 keeps the frozen grid fit and rejects background_by", {
  expect_identical(gaze_replay_spec()$revision, "2026.10")
  expect_error(
    gaze_replay_spec(revision = "2026.08", background_by = "participant"),
    "background_by requires"
  )
  spec <- make_revision_spec("2026.08", n_bins = 16L)
  data <- simulate_replay_hmm(n_participants = 1L, n_items = 3L,
                              n_bins = 16L, seed = 4L)
  model <- fit_gaze_replay_model(data$ref, data$src,
                                 c("participant", "image_id"), spec = spec)

  expect_identical(model$revision, "2026.08")
  expect_null(model$training$em)
  expect_equal(nrow(model$training$transition_candidates), 16L)
  expect_true(all(unlist(model$parameters) %in%
                    unlist(spec$transition_grid)))
})

# Review regressions ---------------------------------------------------------

review_path <- function(xy, duration = rep(1, nrow(xy))) {
  make_gaze_fixations(xy, duration = duration,
                      onset = cumsum(c(0, duration[-length(duration)])))
}

# Training where every participant has two trials, so participant-level
# backgrounds fall back to the pooled density at scoring time.
ambiguous_null_fixture <- function(n_items = 2L, seed = 11L) {
  set.seed(seed)
  w <- revision_screen_w
  h <- revision_screen_h
  layout <- lapply(seq_len(n_items), function(k) {
    cbind(stats::runif(5, 80, w - 80), stats::runif(5, 80, h - 80))
  })
  rows <- expand.grid(image_id = seq_len(n_items),
                      participant = paste0("p", 1:40),
                      stringsAsFactors = FALSE)
  rows <- do.call(rbind, lapply(split(rows, rows$participant), function(d) {
    d[sample(nrow(d), 2), ]
  }))
  clamp <- function(xy) pmin(pmax(xy, 1), w - 1)
  source <- lapply(rows$image_id, function(k) {
    xy <- clamp(layout[[k]][sample(5, 24, TRUE), ] +
                  matrix(stats::rnorm(48, 0, 30), 24))
    noise <- sample(24, 10)
    xy[noise, ] <- cbind(stats::runif(10, 1, w - 1), stats::runif(10, 1, h - 1))
    review_path(xy)
  })
  reference <- lapply(rows$image_id, function(k) {
    review_path(layout[[k]] + stats::rnorm(10, 0, 5))
  })
  spec <- gaze_replay_spec(
    grid_size = 24L, max_skip = 2L, student_df = 4, scale_floor = 6,
    screen = gaze_screen(w, h), background_by = "participant"
  )
  model <- suppressMessages(fit_gaze_replay_model(
    tibble::tibble(participant = rows$participant, image_id = rows$image_id,
                   fixgroup = reference),
    tibble::tibble(participant = rows$participant, image_id = rows$image_id,
                   fixgroup = source),
    c("participant", "image_id"), spec = spec
  ))
  ref_eval <- tibble::tibble(
    participant = "pnew", image_id = seq_len(n_items),
    fixgroup = lapply(layout, function(xy) {
      review_path(xy + stats::rnorm(10, 0, 5))
    })
  )
  list(model = model, layout = layout, ref_eval = ref_eval)
}

test_that("scoring backgrounds never depend on which candidate is true", {
  skip_on_cran()
  fixture <- ambiguous_null_fixture()
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
    recall <- review_path(pmin(pmax(xy, 1), revision_screen_w - 1))
    first <- score_as(recall, 1L)
    second <- score_as(recall, 2L)
    c(
      invariant = max(abs(first$candidates$log_score -
                            second$candidates$log_score)),
      top1 = first$evidence$top1_credit,
      # Per recall fixation: the scale of the former one-bin-per-fixation
      # score.
      margin = (first$candidates$log_score[[1]] -
        first$candidates$log_score[[2]]) / 24
    )
  }))

  expect_lt(max(outcomes[, "invariant"]), 1e-12)
  expect_lt(abs(mean(outcomes[, "top1"]) - 0.5), 0.25)
  expect_lt(abs(mean(outcomes[, "margin"])), 0.2)
})

offscreen_fixture <- function(seed, offscreen, background_share) {
  set.seed(seed)
  w <- revision_screen_w
  h <- revision_screen_h
  path <- function(xy) {
    make_gaze_fixations(
      xy, duration = stats::runif(nrow(xy), 100, 400),
      onset = cumsum(c(0, stats::runif(nrow(xy) - 1, 200, 400)))
    )
  }
  rows <- expand.grid(image_id = 1:8, participant = paste0("p", 1:4),
                      stringsAsFactors = FALSE)
  layout <- lapply(1:8, function(k) {
    cbind(stats::runif(6, -offscreen, w + offscreen),
          stats::runif(6, -offscreen, h + offscreen))
  })
  reference <- lapply(rows$image_id, function(k) {
    path(layout[[k]] + stats::rnorm(12, 0, 10))
  })
  source <- lapply(rows$image_id, function(k) {
    xy <- layout[[k]][sample(6, 12, TRUE), ] + stats::rnorm(24, 0, 35)
    noise <- stats::runif(12) < background_share
    xy[noise, ] <- cbind(stats::runif(sum(noise), 0, w),
                         stats::runif(sum(noise), 0, h))
    path(pmin(pmax(xy, 0), w))
  })
  spec <- gaze_replay_spec(
    grid_size = 24L, max_skip = 2L, student_df = 4, scale_floor = 6,
    screen = gaze_screen(w, h), background_by = "participant"
  )
  suppressMessages(fit_gaze_replay_model(
    tibble::tibble(participant = rows$participant, image_id = rows$image_id,
                   fixgroup = reference),
    tibble::tibble(participant = rows$participant, image_id = rows$image_id,
                   fixgroup = source),
    c("participant", "image_id"), spec = spec
  ))
}

test_that("EM stays monotone when encoding fixations lie off screen", {
  model <- offscreen_fixture(3L, offscreen = 150, background_share = 0.9)
  trace <- model$training$em$log_likelihood_trace

  expect_true(all(diff(trace) > -1e-8 * abs(trace[-1])))
  expect_true(model$training$em$converged)
})

uniform_null_model <- function(grid_size, seed = 6L) {
  set.seed(seed)
  w <- revision_screen_w
  h <- revision_screen_h
  path <- function(xy) {
    make_gaze_fixations(xy, duration = rep(250, nrow(xy)),
                        onset = 300 * (seq_len(nrow(xy)) - 1))
  }
  rows <- expand.grid(image_id = 1:8, participant = paste0("p", 1:4),
                      stringsAsFactors = FALSE)
  layout <- lapply(1:8, function(k) {
    cbind(stats::runif(6, 0, w), stats::runif(6, 0, h))
  })
  reference <- lapply(rows$image_id, function(k) {
    path(layout[[k]] + stats::rnorm(12, 0, 10))
  })
  source <- lapply(rows$image_id, function(k) {
    path(cbind(stats::runif(12, 0, w), stats::runif(12, 0, h)))
  })
  spec <- gaze_replay_spec(
    grid_size = grid_size, max_skip = 2L, student_df = 4, scale_floor = 6,
    screen = gaze_screen(w, h), background_by = "participant"
  )
  suppressMessages(fit_gaze_replay_model(
    tibble::tibble(participant = rows$participant, image_id = rows$image_id,
                   fixgroup = reference),
    tibble::tibble(participant = rows$participant, image_id = rows$image_id,
                   fixgroup = source),
    c("participant", "image_id"), spec = spec
  ))
}

test_that("a uniform on-screen null is fitted as background", {
  skip_on_cran()
  model <- uniform_null_model(grid_size = 12L)
  bound <- min(revision_screen_w, revision_screen_h) / 6

  expect_gt(model$training$em$background_fraction, 0.8)
  expect_lte(model$emission_models[[1]]$replay_scale, bound + 1e-8)
})

test_that("a grid longer than the recalls no longer unidentifies the null", {
  skip_on_cran()
  # Under the former duration-bin observation model a grid of 24 bins for
  # 12-fixation recalls repeated every fixation, and the uniform null was
  # fitted with a background share below 0.8. The fixation-level model
  # observes each fixation once, so grid_size cannot matter.
  model <- uniform_null_model(grid_size = 24L)

  expect_gt(model$training$em$background_fraction, 0.8)
  expect_identical(model$parameters,
                   uniform_null_model(grid_size = 12L)$parameters)
})

test_that("an all-background null leaves the candidates exchangeable", {
  model <- null_model()
  set.seed(12)
  candidates <- lapply(1:4, function(k) {
    revision_path(cbind(stats::runif(6, 60, revision_screen_w - 60),
                        stats::runif(6, 60, revision_screen_h - 60)))
  })
  spread <- vapply(1:10, function(i) {
    recall <- simulate_replay_hmm(n_participants = 1L, n_items = 1L,
                                  null = TRUE, seed = 100L + i)$src$fixgroup[[1]]
    score <- vapply(candidates, function(candidate) {
      gaze_replay_align(candidate, recall, model)$log_score
    }, numeric(1))
    posterior <- exp(score * 32 - max(score * 32))
    posterior <- posterior / sum(posterior)
    c(range = diff(range(score)), posterior = max(abs(posterior - 0.25)))
  }, numeric(2))

  expect_lt(max(spread["range", ]), 0.01)
  expect_lt(max(spread["posterior", ]), 0.05)
})

held_out_protocol_data <- function() {
  data <- simulate_replay_hmm(n_participants = 3L, n_items = 16L,
                              n_bins = 16L, seed = 9L)
  data$ref$block <- data$ref$image_id %% 2L
  data$src$block <- data$src$image_id %% 2L
  data
}

test_that("participant holdout uses the population background and says so", {
  skip_on_cran()
  data <- held_out_protocol_data()
  expect_message(
    fit <- gaze_replay_cv(
      data$ref, data$src, match_on = c("participant", "image_id"),
      contrast_on = NULL, split_on = "participant",
      n_folds = 3, seed = 1, spec = make_revision_spec(n_bins = 16L)
    ),
    "absent from their training fold"
  )
  expect_identical(fit$provenance$background_fallback$unseen_level_rows,
                   nrow(data$src))
  # Every item is a candidate here, so excluding the candidate pool would
  # leave no population data; the unexcluded population density is used.
  expect_true(all(fit$results$background_level == "pooled_all_items"))
})

test_that("held-out support uses only non-candidate recalls of the participant", {
  skip_on_cran()
  data <- held_out_protocol_data()
  spec <- gaze_replay_spec(
    grid_size = 16L, max_skip = 2L, student_df = 4, scale_floor = 6,
    screen = gaze_screen(revision_screen_w, revision_screen_h),
    background_by = "participant", background_support = "held_out"
  )
  fit <- suppressMessages(gaze_replay_cv(
    data$ref, data$src, match_on = c("participant", "image_id"),
    contrast_on = c("participant", "block"),
    split_on = c("participant", "image_id"), n_folds = 2, seed = 1,
    spec = spec
  ))
  expect_true(all(fit$results$background_level == "held_out_support"))
  expect_true(all(is.finite(fit$results$gaze_info_bits)))

  model <- suppressMessages(fit_gaze_replay_model(
    data$ref[data$ref$participant != "p1", ],
    data$src[data$src$participant != "p1", ],
    c("participant", "image_id"), spec = spec
  ))
  eval_source <- data$src[data$src$participant == "p1", ]
  ref_eval <- data$ref[data$ref$participant == "p1", ]
  row <- eval_source[eval_source$image_id == 2L, ]
  relabelled <- row
  relabelled$image_id <- 4L
  first <- score_gaze_replay_row(
    row, ref_eval, c("participant", "image_id"), c("participant", "block"),
    "fixgroup", "fixgroup", model, support_tab = eval_source
  )
  second <- score_gaze_replay_row(
    relabelled, ref_eval, c("participant", "image_id"),
    c("participant", "block"), "fixgroup", "fixgroup", model,
    support_tab = eval_source
  )
  expect_identical(first$alignment$grid$background_level, "held_out_support")
  expect_equal(first$candidates$log_score, second$candidates$log_score,
               tolerance = 1e-12)
  expect_error(
    gaze_replay_spec(background_support = "held_out"),
    "requires revision"
  )
})

test_that("saved objects without a revision score as frozen 2026.08", {
  data <- simulate_replay_hmm(n_participants = 1L, n_items = 3L,
                              n_bins = 16L, seed = 4L)
  spec <- make_revision_spec("2026.08", n_bins = 16L)
  model <- fit_gaze_replay_model(data$ref, data$src,
                                 c("participant", "image_id"), spec = spec)
  old_spec <- unclass(spec)
  old_spec[c("revision", "background_by", "background_support")] <- NULL
  old_spec <- structure(old_spec, class = class(spec))
  old_model <- unclass(model)
  old_model[c("revision", "background", "match_on")] <- NULL
  old_model$spec <- old_spec
  old_model <- structure(old_model, class = class(model))
  old_model <- unserialize(serialize(old_model, NULL))

  reference <- data$ref$fixgroup[[1]]
  source <- data$src$fixgroup[[2]]
  expect_identical(
    gaze_replay_align(reference, source, old_model)$log_score,
    gaze_replay_align(reference, source, model)$log_score
  )
  refit <- fit_gaze_replay_model(data$ref, data$src,
                                 c("participant", "image_id"),
                                 spec = old_spec)
  expect_identical(refit$parameters, model$parameters)
  expect_identical(refit$revision, "2026.08")

  future <- model
  future$revision <- "2027.01"
  expect_error(gaze_replay_align(reference, source, future),
               "Unknown Replay revision")
  future_spec <- spec
  future_spec$revision <- "2027.01"
  expect_error(
    fit_gaze_replay_model(data$ref, data$src, c("participant", "image_id"),
                          spec = future_spec),
    "Unknown Replay revision"
  )
})

test_that("undeclared screens and pooled background fallbacks are reported", {
  data <- simulate_replay_hmm(n_participants = 1L, n_items = 4L,
                              n_bins = 16L, seed = 7L)
  plane <- gaze_replay_spec(grid_size = 16L, scale_floor = 6)
  reset_gaze_replay_notices()
  expect_message(
    fit_gaze_replay_model(data$ref, data$src, c("participant", "image_id"),
                          spec = plane),
    "whole plane"
  )
  expect_silent(
    fit_gaze_replay_model(data$ref, data$src, c("participant", "image_id"),
                          spec = plane)
  )

  sparse <- simulate_replay_hmm(n_participants = 4L, n_items = 4L,
                                n_bins = 16L, seed = 8L)
  two_each <- list(ref = sparse$ref[sparse$ref$image_id <= 2L, ],
                   src = sparse$src[sparse$src$image_id <= 2L, ])
  expect_message(fit_revision_model(two_each), "pooled background")
  fit <- suppressMessages(gaze_replay_cv(
    sparse$ref, sparse$src, match_on = c("participant", "image_id"),
    contrast_on = "participant", n_folds = 2, seed = 2,
    spec = make_revision_spec(n_bins = 16L)
  ))
  expect_gt(fit$provenance$background_fallback$scored_rows, 0L)
  expect_true("pooled" %in% fit$results$background_level)
})

test_that("grids longer than the recall's fixation count change nothing", {
  # The repeated-grid report of the duration-bin model is retired: each
  # recall fixation is observed exactly once.
  data <- simulate_replay_hmm(n_participants = 2L, n_items = 4L,
                              n_bins = 16L, seed = 10L)
  exact <- make_revision_spec(n_bins = 16L)
  repeated <- make_revision_spec(n_bins = 32L)

  expect_no_message(
    exact_model <- fit_gaze_replay_model(
      data$ref, data$src, c("participant", "image_id"), spec = exact
    ),
    message = "repeat"
  )
  expect_no_message(
    repeated_model <- fit_gaze_replay_model(
      data$ref, data$src, c("participant", "image_id"), spec = repeated
    ),
    message = "repeat"
  )
  expect_identical(repeated_model$parameters, exact_model$parameters)
  expect_identical(repeated_model$emission_models,
                   exact_model$emission_models)
  fit <- suppressMessages(gaze_replay_cv(
    data$ref, data$src, match_on = c("participant", "image_id"),
    contrast_on = "participant", n_folds = 2, seed = 2, spec = repeated
  ))
  expect_false("repeated_grid_rows" %in% names(fit$provenance))
})

test_that("the participant-holdout message does not recommend held-out support", {
  skip_on_cran()
  data <- held_out_protocol_data()
  messages <- character()
  withCallingHandlers(
    gaze_replay_cv(
      data$ref, data$src, match_on = c("participant", "image_id"),
      contrast_on = NULL, split_on = "participant",
      n_folds = 3, seed = 1, spec = make_revision_spec(n_bins = 16L)
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  holdout <- grep("absent from their training fold", messages, value = TRUE)

  expect_length(holdout, 1L)
  expect_match(holdout, "item-level split", fixed = TRUE)
})
