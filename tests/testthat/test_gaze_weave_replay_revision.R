# Replay revision 2026.10: screen normalisation, background, and EM fitting.
# All data are simulated from the Replay generative model itself.

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

# One fixation per duration bin (grid_size = n_bins), so the Replay HMM is
# the exact generative model of the simulated recalls.
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

# Posterior share of duration bins assigned to the background state on the
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
  grid <- list(coords = points)
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
  model <- signal_model()
  nx <- 512L
  ny <- 384L
  gx <- (seq_len(nx) - 0.5) * revision_screen_w / nx
  gy <- (seq_len(ny) - 0.5) * revision_screen_h / ny
  points <- as.matrix(expand.grid(gx, gy))
  cell <- (revision_screen_w / nx) * (revision_screen_h / ny)
  # Centre, edge midpoint, corner, and a fixation recorded just off screen.
  reference <- list(coords = rbind(
    c(revision_screen_w / 2, revision_screen_h / 2),
    c(1, revision_screen_h / 2),
    c(1, 1),
    c(revision_screen_w + 5, revision_screen_h - 2)
  ))
  log_emission <- gaze_replay_emissions(
    reference, revision_emission_grid(model, points),
    model$emission_models[[1]], model$spec$student_df
  )
  integral <- colSums(exp(log_emission)) * cell

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
  # Per-bin shortfall of the model score from the exact screen-truncated
  # generating log density of the replayed recall. Identical relative
  # geometry must give the same shortfall wherever the candidate sits.
  shortfall <- function(candidate) {
    recall <- t(vapply(state, function(s) {
      revision_rtrunc_t(candidate[s, ], 40)
    }, numeric(2)))
    generating <- gaze_replay_truncated_t_log_density(
      recall, candidate, 40, 4, rect
    )[cbind(seq_along(state), state)]
    sum(generating) / length(state) - gaze_replay_align(
      revision_path(candidate), revision_path(recall), model
    )$log_score
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

test_that("a full null is fitted as background", {
  data <- simulate_replay_hmm(null = TRUE, n_items = 8L, seed = 2L)
  model <- fit_revision_model(data)
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
  grid <- gaze_duration_grid(source, model$spec$grid_size)
  grid <- revision_emission_grid(model, grid$coords)
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
