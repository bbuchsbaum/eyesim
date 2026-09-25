# GazeWeave Replay revision 2026.10 ----------------------------------------
#
# Screen-truncated emissions, participant-level out-of-fold background
# densities, and Baum-Welch (EM) parameter fitting. Revision "2026.08" keeps
# the frozen behaviour in gaze_weave_replay.R.

# Objects saved before revisions existed carry no revision and are frozen
# revision 2026.08; any other unknown revision is refused.
gaze_replay_revision_of <- function(x) {
  revision <- x$revision
  if (is.null(revision)) return("2026.08")
  if (!is.character(revision) || length(revision) != 1L ||
      !revision %in% c("2026.08", "2026.10")) {
    stop("Unknown Replay revision '", paste(revision, collapse = ", "),
         "'; this eyesim supports \"2026.08\" and \"2026.10\".")
  }
  revision
}

gaze_replay_is_legacy <- function(spec) {
  identical(gaze_replay_revision_of(spec), "2026.08")
}

gaze_replay_em_control <- function() {
  list(
    max_iterations = 200L,
    tolerance = 1e-6,
    probability_bounds = c(1e-6, 1 - 1e-6),
    background_floor = 0.05,
    background_floor_bounds = c(1e-3, 1 - 1e-6),
    background_min_trials = 3L,
    quadrature_nodes = 40L,
    max_scale_fraction = 1 / 6
  )
}

# One-time notices ------------------------------------------------------------

gaze_replay_notices <- new.env(parent = emptyenv())

reset_gaze_replay_notices <- function() {
  rm(list = ls(gaze_replay_notices), envir = gaze_replay_notices)
  invisible(NULL)
}

note_gaze_replay_plane <- function(spec) {
  if (gaze_replay_is_legacy(spec) || !is.null(spec$screen) ||
      isTRUE(gaze_replay_notices$plane)) {
    return(invisible(FALSE))
  }
  assign("plane", TRUE, envir = gaze_replay_notices)
  message(
    "Replay revision 2026.10: no screen declared, so emission densities are ",
    "normalised over the whole plane. Declare gaze_screen() for ",
    "screen-truncated emissions. (Shown once per session.)"
  )
  invisible(TRUE)
}

gaze_replay_repeated_grid_message <- function(count, total, what,
                                              grid_size) {
  paste0(
    "Replay: ", count, " of ", total, " ", what, " have fewer fixations ",
    "than grid_size = ", grid_size, ", so fixations repeat across ",
    "duration bins. Replay and background are then not identified: the ",
    "fitted background share is not interpretable, and candidate scores can ",
    "differ under a pure background null."
  )
}

# Screen geometry ------------------------------------------------------------

gaze_replay_screen_rect <- function(screen) {
  list(
    xlim = as.numeric(screen$xlim),
    ylim = as.numeric(screen$ylim),
    area = diff(screen$xlim) * diff(screen$ylim)
  )
}

# Without a declared screen the support is the whole plane: emissions are then
# ordinary (untruncated) densities. A rectangle inferred from training
# coordinates is deliberately not used, because held-out paths can lie outside
# it and clamping would distort their geometry.
resolve_gaze_replay_screen <- function(spec) {
  if (!is.null(spec$screen)) {
    rect <- gaze_replay_screen_rect(spec$screen)
    rect$source <- "spec"
    return(rect)
  }
  list(xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), area = Inf, source = "plane")
}

gaze_replay_clamp_to_screen <- function(coords, rect) {
  coords[, 1] <- pmin(pmax(coords[, 1], rect$xlim[[1]]), rect$xlim[[2]])
  coords[, 2] <- pmin(pmax(coords[, 2], rect$ylim[[1]]), rect$ylim[[2]])
  coords
}

# Truncated bivariate Student emissions --------------------------------------

gaze_replay_gamma_quadrature <- local({
  cache <- list()
  function(degrees, n_nodes) {
    key <- paste(format(degrees, digits = 17), n_nodes)
    if (!is.null(cache[[key]])) return(cache[[key]])
    # Generalized Gauss-Laguerre (Golub-Welsch) for x^alpha exp(-x).
    alpha <- degrees / 2 - 1
    k <- seq_len(n_nodes) - 1L
    diagonal <- 2 * k + alpha + 1
    off <- sqrt(seq_len(n_nodes - 1L) * (seq_len(n_nodes - 1L) + alpha))
    jacobi <- diag(diagonal, n_nodes)
    jacobi[cbind(seq_len(n_nodes - 1L), seq_len(n_nodes - 1L) + 1L)] <- off
    jacobi[cbind(seq_len(n_nodes - 1L) + 1L, seq_len(n_nodes - 1L))] <- off
    decomposition <- eigen(jacobi, symmetric = TRUE)
    weights <- decomposition$vectors[1, ]^2
    value <- list(
      # W ~ Gamma(shape = degrees / 2, rate = degrees / 2).
      precision = decomposition$values / (degrees / 2),
      weights = weights / sum(weights)
    )
    cache[[key]] <<- value
    value
  }
})

# Probability that an isotropic bivariate Student variable centred at each row
# of `center` with common `scale` falls inside the screen rectangle.
gaze_replay_t_screen_mass <- function(center, scale, degrees, rect,
                                      n_nodes = gaze_replay_em_control()$quadrature_nodes) {
  quadrature <- gaze_replay_gamma_quadrature(degrees, n_nodes)
  root <- sqrt(quadrature$precision) / scale
  mass_axis <- function(values, limits) {
    upper <- outer(limits[[2]] - values, root)
    lower <- outer(limits[[1]] - values, root)
    stats::pnorm(upper) - stats::pnorm(lower)
  }
  mass <- mass_axis(center[, 1], rect$xlim) * mass_axis(center[, 2], rect$ylim)
  pmax(as.numeric(mass %*% quadrature$weights), .Machine$double.xmin)
}

gaze_replay_truncated_t_log_density <- function(observed, centers, scale,
                                                degrees, rect) {
  # Encoding fixations recorded off screen are projected onto it, so that no
  # state has vanishing on-screen mass.
  centers <- gaze_replay_clamp_to_screen(centers, rect)
  log_mass <- log(gaze_replay_t_screen_mass(centers, scale, degrees, rect))
  out <- matrix(NA_real_, nrow(observed), nrow(centers))
  for (j in seq_len(nrow(centers))) {
    out[, j] <- gaze_bivariate_t_log_density(
      observed, centers[j, ], scale, degrees
    ) - log_mass[[j]]
  }
  out
}

# Participant-level screen-truncated background ------------------------------

gaze_replay_background_table <- function(measures, background_key,
                                         exclude_key) {
  rows <- lapply(seq_along(measures), function(i) {
    measure <- measures[[i]]
    data.frame(
      x = measure$coords[, 1],
      y = measure$coords[, 2],
      mass = measure$mass,
      trial = i,
      background_key = background_key[[i]],
      exclude_key = exclude_key[[i]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gaze_replay_kde_bandwidth <- function(points, weights, n_points, floor) {
  weights <- weights / sum(weights)
  center <- colSums(points * weights)
  variance <- colSums(sweep(points, 2, center)^2 * weights)
  spread <- sqrt(mean(variance))
  max(floor, spread * max(n_points, 1)^(-1 / 6))
}

gaze_replay_background_subset <- function(background, background_key = NULL,
                                          exclude_key = NULL, min_trials,
                                          support = NULL) {
  if (!is.null(support) && nrow(support)) {
    support <- support[!support$exclude_key %in% exclude_key, , drop = FALSE]
    if (length(unique(support$trial)) >= min_trials) {
      return(list(table = support, level = "held_out_support"))
    }
  }
  table <- background$table
  keep <- rep(TRUE, nrow(table))
  if (length(exclude_key)) keep <- keep & !table$exclude_key %in% exclude_key
  level <- "pooled"
  if (length(exclude_key) &&
      length(unique(table$trial[keep])) < min_trials) {
    # Excluding the whole candidate pool left too little population data
    # (e.g. every item is a candidate). Excluding nothing is equally
    # independent of which candidate is true.
    keep <- rep(TRUE, nrow(table))
    level <- "pooled_all_items"
  }
  if (!is.null(background_key)) {
    own <- keep & table$background_key == background_key
    if (length(unique(table$trial[own])) >= min_trials) {
      keep <- own
      level <- "participant"
    }
  }
  list(table = table[keep, , drop = FALSE], level = level)
}

gaze_replay_background_log_density <- function(observed, background,
                                               background_key = NULL,
                                               exclude_key = NULL,
                                               support = NULL) {
  control <- background$control
  rect <- background$screen
  floor_log_density <- if (is.finite(rect$area)) {
    rep(-log(rect$area), nrow(observed))
  } else {
    # Plane support: a broad Gaussian over the pooled training recalls.
    rowSums(stats::dnorm(
      observed, rep(background$floor_center, each = nrow(observed)),
      background$floor_sd, log = TRUE
    ))
  }
  subset <- gaze_replay_background_subset(
    background, background_key, exclude_key, control$background_min_trials,
    support
  )
  table <- subset$table
  floor_weight <- background$floor_weight
  if (is.null(floor_weight)) floor_weight <- control$background_floor
  if (!nrow(table)) {
    return(list(
      log_density = floor_log_density,
      log_kde = NULL,
      log_floor = floor_log_density,
      level = "floor", bandwidth = NA_real_, trials = 0L
    ))
  }
  trials <- unique(table$trial)
  # Every trial carries equal total weight; duration mass within trial.
  weight <- table$mass / stats::ave(table$mass, table$trial, FUN = sum) /
    length(trials)
  points <- cbind(table$x, table$y)
  bandwidth <- gaze_replay_kde_bandwidth(
    points, weight, nrow(points), background$bandwidth_floor
  )
  mass_axis <- function(values, limits) {
    stats::pnorm((limits[[2]] - values) / bandwidth) -
      stats::pnorm((limits[[1]] - values) / bandwidth)
  }
  kernel_mass <- pmax(
    mass_axis(points[, 1], rect$xlim) * mass_axis(points[, 2], rect$ylim),
    .Machine$double.xmin
  )
  scaled <- weight / kernel_mass / (2 * pi * bandwidth^2)
  dx <- outer(observed[, 1], points[, 1], FUN = "-")
  dy <- outer(observed[, 2], points[, 2], FUN = "-")
  kde <- as.numeric(exp(-(dx^2 + dy^2) / (2 * bandwidth^2)) %*% scaled)
  log_kde <- log(pmax(kde, .Machine$double.xmin))
  list(
    log_density = gaze_replay_mix_background(
      log_kde, floor_log_density, floor_weight
    ),
    log_kde = log_kde,
    log_floor = floor_log_density,
    level = subset$level,
    bandwidth = bandwidth,
    trials = length(trials)
  )
}

# Background = (1 - floor_weight) * KDE + floor_weight * floor (uniform on the
# screen, or a broad Gaussian on the plane). floor_weight is fitted by EM.
gaze_replay_mix_background <- function(log_kde, log_floor, floor_weight) {
  if (is.null(log_kde)) return(log_floor)
  a <- log1p(-floor_weight) + log_kde
  b <- log(floor_weight) + log_floor
  top <- pmax(a, b)
  top + log(exp(a - top) + exp(b - top))
}

# Revision emissions ----------------------------------------------------------

gaze_replay_emissions_revision <- function(reference, grid, emission_model,
                                           degrees) {
  if (is.null(grid$log_background)) {
    stop("Replay revision 2026.10 emissions need a prepared background density.")
  }
  cbind(
    grid$log_background,
    gaze_replay_truncated_t_log_density(
      grid$coords, reference$coords, emission_model$replay_scale, degrees,
      emission_model$screen
    )
  )
}

# Baum-Welch expected counts ---------------------------------------------------

gaze_replay_expected_counts <- function(log_emission, transition, fit) {
  n_time <- nrow(log_emission)
  n_states <- ncol(log_emission)
  mass <- transition$reference_mass
  emission_max <- apply(log_emission, 1, max)
  emission_scaled <- exp(log_emission - emission_max)
  parameters <- transition$parameters
  outer_replay <- matrix(0, n_states - 1L, n_states - 1L)
  counts <- c(bb = 0, b_replay = 0, replay_b = 0, replay_restart = 0)
  if (n_time > 1L) {
    for (time in seq_len(n_time - 1L)) {
      w <- emission_scaled[time + 1L, ] * fit$backward[time + 1L, ] /
        fit$scaling[[time + 1L]]
      alpha <- fit$forward[time, ]
      replay_alpha <- alpha[-1]
      replay_total <- sum(replay_alpha)
      restart_target <- sum(mass * w[-1])
      counts[["bb"]] <- counts[["bb"]] +
        alpha[[1]] * parameters$background_stay * w[[1]]
      counts[["b_replay"]] <- counts[["b_replay"]] +
        alpha[[1]] * (1 - parameters$background_stay) * restart_target
      counts[["replay_b"]] <- counts[["replay_b"]] +
        replay_total * parameters$background * w[[1]]
      counts[["replay_restart"]] <- counts[["replay_restart"]] +
        replay_total * parameters$restart * restart_target
      outer_replay <- outer_replay + outer(replay_alpha, w[-1])
    }
  }
  local <- transition$local_component[-1, -1, drop = FALSE] * outer_replay
  n_reference <- nrow(local)
  has_successor <- seq_len(n_reference) < n_reference
  stay <- diag(local)
  list(
    counts = counts,
    replay_local = sum(local),
    advance = sum(local) - sum(stay),
    stay_eligible = sum(stay[has_successor]),
    initial_background = fit$posterior[1, 1],
    posterior = fit$posterior
  )
}

gaze_replay_m_step_transitions <- function(stats, previous, bounds) {
  clamp <- function(value, fallback) {
    if (!is.finite(value)) return(fallback)
    min(max(value, bounds[[1]]), bounds[[2]])
  }
  replay_exit <- stats$replay_b + stats$initial_background
  replay_total <- stats$replay_b + stats$replay_restart +
    stats$replay_local + stats$initial_n
  background <- clamp(replay_exit / replay_total, previous$background)
  restart_share <- clamp(
    stats$replay_restart / (stats$replay_restart + stats$replay_local),
    previous$restart / (1 - previous$background)
  )
  advance <- clamp(
    stats$advance / (stats$advance + stats$stay_eligible),
    previous$advance
  )
  background_stay <- clamp(
    stats$bb / (stats$bb + stats$b_replay),
    previous$background_stay
  )
  list(
    background = background,
    restart = (1 - background) * restart_share,
    advance = advance,
    background_stay = background_stay
  )
}

gaze_replay_m_step_scale <- function(items, degrees, rect, floor, current,
                                     upper) {
  objective <- function(log_scale) {
    scale <- exp(log_scale)
    total <- 0
    for (item in items) {
      replay_weight <- item$weight * item$posterior[, -1, drop = FALSE]
      if (!any(replay_weight > 0)) next
      log_density <- lgamma((degrees + 2) / 2) - lgamma(degrees / 2) -
        log(degrees * pi) - 2 * log(scale) -
        ((degrees + 2) / 2) * log1p(item$distance_sq / (degrees * scale^2))
      log_mass <- log(gaze_replay_t_screen_mass(
        item$centers, scale, degrees, rect
      ))
      total <- total + sum(replay_weight * log_density) -
        sum(colSums(replay_weight) * log_mass)
    }
    total
  }
  lower <- log(floor)
  if (lower >= log(upper)) return(floor)
  best <- stats::optimize(
    objective, c(lower, log(upper)), maximum = TRUE, tol = 1e-6
  )
  candidates <- c(best$maximum, log(current))
  values <- vapply(candidates, objective, numeric(1))
  exp(candidates[[which.max(values)]])
}

# EM fit on training pairs ------------------------------------------------------

gaze_replay_initial_parameters <- function(spec) {
  values <- lapply(spec$transition_grid, stats::median)
  if (values$background + values$restart >= 1) {
    values$restart <- (1 - values$background) / 2
  }
  values
}

fit_gaze_replay_em <- function(pairs, pair_episode, groups, spec, rect,
                               scale_upper) {
  control <- gaze_replay_em_control()
  degrees <- spec$student_df
  bounds <- control$probability_bounds
  parameters <- gaze_replay_initial_parameters(spec)
  scale <- lapply(groups, function(group) {
    member <- vapply(pairs, `[[`, character(1), "group") == group
    nearest <- unlist(lapply(pairs[member], function(pair) {
      sqrt(apply(pair$distance_sq, 1, min))
    }), use.names = FALSE)
    min(max(stats::median(nearest), spec$scale_floor), scale_upper)
  })
  names(scale) <- groups
  floor_weight <- control$background_floor
  trace <- numeric()
  converged <- FALSE
  decreased <- FALSE
  n_bins <- sum(vapply(split(seq_along(pairs), pair_episode), function(rows) {
    nrow(pairs[[rows[[1]]]]$grid)
  }, numeric(1)))

  run_e_step <- function(parameters, scale, floor_weight) {
    lapply(pairs, function(pair) {
      transition <- gaze_replay_transition(
        pair$reference_mass, parameters, spec$max_skip
      )
      log_background <- gaze_replay_mix_background(
        pair$log_kde, pair$log_floor, floor_weight
      )
      log_emission <- cbind(
        log_background,
        gaze_replay_truncated_t_log_density(
          pair$grid, pair$centers, scale[[pair$group]], degrees, rect
        )
      )
      fit <- gaze_replay_forward_backward(log_emission, transition)
      list(
        log_likelihood = fit$log_likelihood,
        stats = gaze_replay_expected_counts(log_emission, transition, fit)
      )
    })
  }

  e_step <- run_e_step(parameters, scale, floor_weight)
  for (iteration in seq_len(control$max_iterations)) {
    pair_ll <- vapply(e_step, `[[`, numeric(1), "log_likelihood")
    episode_rows <- split(seq_along(pairs), pair_episode)
    weights <- numeric(length(pairs))
    total <- 0
    for (rows in episode_rows) {
      mixture <- gaze_log_sum_exp(pair_ll[rows] - log(length(rows)))
      weights[rows] <- exp(pair_ll[rows] - log(length(rows)) - mixture)
      total <- total + mixture
    }
    trace <- c(trace, total)
    if (length(trace) > 1L) {
      change <- trace[[length(trace)]] - trace[[length(trace) - 1L]]
      # EM must not decrease the likelihood; a decrease beyond tolerance
      # signals an inconsistent update and is never reported as convergence.
      if (change < -control$tolerance * n_bins) decreased <- TRUE
      if (abs(change) <= control$tolerance * n_bins) {
        converged <- !decreased
        break
      }
    }

    stats <- list(bb = 0, b_replay = 0, replay_b = 0, replay_restart = 0,
                  replay_local = 0, advance = 0, stay_eligible = 0,
                  initial_background = 0, initial_n = 0)
    for (i in seq_along(pairs)) {
      value <- e_step[[i]]$stats
      w <- weights[[i]]
      for (name in names(value$counts)) {
        stats[[name]] <- stats[[name]] + w * value$counts[[name]]
      }
      stats$replay_local <- stats$replay_local + w * value$replay_local
      stats$advance <- stats$advance + w * value$advance
      stats$stay_eligible <- stats$stay_eligible + w * value$stay_eligible
      stats$initial_background <- stats$initial_background +
        w * value$initial_background
      stats$initial_n <- stats$initial_n + w
    }
    parameters <- gaze_replay_m_step_transitions(stats, parameters, bounds)
    # Exact M-step for the floor weight: maximise the expected background
    # emission log density sum_t gamma_t(bg) log((1 - w) kde_t + w floor_t),
    # which is concave in w.
    mixable <- which(!vapply(pairs, function(pair) is.null(pair$log_kde),
                             logical(1)))
    if (length(mixable)) {
      gamma <- unlist(lapply(mixable, function(i) {
        weights[[i]] * e_step[[i]]$stats$posterior[, 1]
      }))
      log_kde <- unlist(lapply(mixable, function(i) pairs[[i]]$log_kde))
      log_floor <- unlist(lapply(mixable, function(i) pairs[[i]]$log_floor))
      floor_objective <- function(w) {
        sum(gamma * gaze_replay_mix_background(log_kde, log_floor, w))
      }
      if (sum(gamma) > 0) {
        candidate <- stats::optimize(
          floor_objective, control$background_floor_bounds, maximum = TRUE,
          tol = 1e-8
        )$maximum
        if (floor_objective(candidate) >= floor_objective(floor_weight)) {
          floor_weight <- candidate
        }
      }
    }
    for (group in groups) {
      member <- which(vapply(pairs, `[[`, character(1), "group") == group)
      items <- lapply(member, function(i) list(
        weight = weights[[i]],
        posterior = e_step[[i]]$stats$posterior,
        distance_sq = pairs[[i]]$distance_sq,
        centers = pairs[[i]]$centers
      ))
      scale[[group]] <- gaze_replay_m_step_scale(
        items, degrees, rect, spec$scale_floor, scale[[group]], scale_upper
      )
    }
    e_step <- run_e_step(parameters, scale, floor_weight)
  }

  posterior_background <- vapply(e_step, function(value) {
    mean(value$stats$posterior[, 1])
  }, numeric(1))
  pair_ll <- vapply(e_step, `[[`, numeric(1), "log_likelihood")
  episode_rows <- split(seq_along(pairs), pair_episode)
  weights <- numeric(length(pairs))
  for (rows in episode_rows) {
    mixture <- gaze_log_sum_exp(pair_ll[rows] - log(length(rows)))
    weights[rows] <- exp(pair_ll[rows] - log(length(rows)) - mixture)
  }
  list(
    parameters = parameters,
    replay_scale = scale,
    log_likelihood_trace = trace,
    iterations = length(trace),
    converged = converged,
    monotone = !decreased,
    floor_weight = floor_weight,
    background_fraction = sum(weights * posterior_background) /
      length(episode_rows),
    pair_background_coverage = posterior_background
  )
}

fit_gaze_replay_revision <- function(episodes, source_tab, match_on, spec,
                                     group_key) {
  control <- gaze_replay_em_control()
  note_gaze_replay_plane(spec)
  background_by <- spec$background_by
  if (!is.null(background_by) && !all(background_by %in% names(source_tab))) {
    stop("Replay background_by columns must exist in the source table.")
  }
  exclude_on <- setdiff(match_on, background_by)
  if (!length(exclude_on)) exclude_on <- match_on
  background_key <- if (is.null(background_by)) {
    rep("all", nrow(source_tab))
  } else {
    gaze_key(source_tab, background_by, "background_by")
  }
  exclude_key <- gaze_key(source_tab, exclude_on, "background exclusion")

  all_coords <- do.call(rbind, c(
    lapply(episodes, function(episode) episode$source$coords),
    unlist(lapply(episodes, function(episode) {
      lapply(episode$references, `[[`, "coords")
    }), recursive = FALSE)
  ))
  rect <- resolve_gaze_replay_screen(spec)
  all_coords <- gaze_replay_clamp_to_screen(all_coords, rect)
  # A replay state is local: its scale may not exceed a fixed fraction of the
  # shorter side of the screen (or of the training extent without a screen).
  # Otherwise a near-uniform replay state is indistinguishable from the
  # background and "no replay" is not identifiable.
  extent <- if (is.finite(rect$area)) {
    c(diff(rect$xlim), diff(rect$ylim))
  } else {
    apply(all_coords, 2, function(v) diff(range(v)))
  }
  scale_upper <- max(
    control$max_scale_fraction * min(extent), 2 * spec$scale_floor
  )
  source_coords <- gaze_replay_clamp_to_screen(do.call(rbind, lapply(
    episodes, function(episode) episode$source$coords
  )), rect)
  background <- list(
    table = gaze_replay_background_table(
      lapply(episodes, function(episode) {
        measure <- episode$source
        measure$coords <- gaze_replay_clamp_to_screen(measure$coords, rect)
        measure
      }),
      background_key, exclude_key
    ),
    screen = rect,
    control = control,
    bandwidth_floor = spec$scale_floor,
    floor_center = colMeans(source_coords),
    floor_sd = max(
      4 * sqrt(mean(apply(source_coords, 2, stats::var))), spec$scale_floor,
      na.rm = TRUE
    ),
    background_by = background_by,
    exclude_on = exclude_on,
    support = spec$background_support
  )

  pairs <- list()
  pair_episode <- integer()
  background_level <- character(length(episodes))
  for (i in seq_along(episodes)) {
    episode <- episodes[[i]]
    grid <- gaze_duration_grid(episode$source, spec$grid_size)
    coords <- gaze_replay_clamp_to_screen(grid$coords, rect)
    bg <- gaze_replay_background_log_density(
      coords, background, background_key[[i]], exclude_key[[i]]
    )
    background_level[[i]] <- bg$level
    for (reference in episode$references) {
      # The same on-screen centres used by the emissions (see
      # gaze_replay_truncated_t_log_density) enter the scale M-step.
      centers <- gaze_replay_clamp_to_screen(reference$coords, rect)
      dx <- outer(coords[, 1], centers[, 1], FUN = "-")
      dy <- outer(coords[, 2], centers[, 2], FUN = "-")
      pairs[[length(pairs) + 1L]] <- list(
        grid = coords,
        centers = centers,
        reference_mass = reference$mass,
        distance_sq = dx^2 + dy^2,
        log_kde = bg$log_kde,
        log_floor = bg$log_floor,
        group = episode$group
      )
      pair_episode <- c(pair_episode, i)
    }
  }
  groups <- sort(unique(group_key))
  em <- fit_gaze_replay_em(pairs, pair_episode, groups, spec, rect, scale_upper)
  emission_models <- lapply(groups, function(group) {
    list(
      revision = "2026.10",
      replay_scale = em$replay_scale[[group]],
      screen = rect,
      training_grid_points = sum(vapply(pairs, function(pair) {
        if (identical(pair$group, group)) nrow(pair$grid) else 0L
      }, integer(1)))
    )
  })
  names(emission_models) <- groups
  background$floor_weight <- em$floor_weight
  fixation_count <- vapply(episodes, function(episode) {
    nrow(episode$source$coords)
  }, integer(1))
  repeated_n <- sum(fixation_count < spec$grid_size)
  if (repeated_n > 0L) {
    message(gaze_replay_repeated_grid_message(
      repeated_n, length(fixation_count), "training recalls", spec$grid_size
    ))
  }
  fallback_n <- sum(background_level != "participant")
  if (!is.null(background_by) && fallback_n > 0L) {
    message(
      "Replay: ", fallback_n, " of ", length(background_level),
      " training rows used the pooled background because their ",
      "background_by level had fewer than ", control$background_min_trials,
      " other training trials."
    )
  }
  list(
    emission_models = emission_models,
    parameters = em$parameters,
    background = background,
    screen = rect,
    em = c(em[setdiff(names(em), c("parameters", "replay_scale"))],
           list(repeated_grid_rows = repeated_n)),
    background_level = background_level
  )
}

gaze_replay_prepare_background <- function(grid, model, background_key = NULL,
                                           exclude_key = NULL,
                                           support = NULL) {
  rect <- model$background$screen
  grid$raw_coords <- grid$coords
  grid$coords <- gaze_replay_clamp_to_screen(grid$coords, rect)
  bg <- gaze_replay_background_log_density(
    grid$coords, model$background, background_key, exclude_key, support
  )
  grid$log_background <- bg$log_density
  grid$background_level <- bg$level
  grid$background_bandwidth <- bg$bandwidth
  grid
}

# Encode user-facing background level values (in `background_by` column order)
# with the same key encoding used during fitting.
encode_gaze_replay_background_key <- function(model, value) {
  if (is.null(value)) return(NULL)
  if (!identical(model$revision, "2026.10")) {
    stop("background_key applies only to revision 2026.10 Replay models.")
  }
  by <- model$background$background_by
  if (is.null(by)) {
    stop("This Replay model was fitted without background_by.")
  }
  if (length(value) != length(by) || anyNA(value)) {
    stop("background_key must give one value per background_by column.")
  }
  tab <- as.data.frame(
    stats::setNames(lapply(value, identity), by),
    stringsAsFactors = FALSE
  )
  gaze_key(tab, by, "background_key")
}

# Held-out support rows (background_support = "held_out"): recalls of the
# scored row's own background_by level that are available at evaluation time,
# restricted later to items outside the candidate pool. They are registered
# with the fitted warp exactly like the scored recall.
gaze_replay_support_table <- function(model, support_tab, sourcevar,
                                      background_key) {
  background <- model$background
  if (!identical(background$support, "held_out") || is.null(support_tab) ||
      is.null(background_key) || !nrow(support_tab)) {
    return(NULL)
  }
  keys <- gaze_key(support_tab, background$background_by, "background_by")
  rows <- which(keys == background_key)
  if (!length(rows)) return(NULL)
  exclude <- gaze_key(
    support_tab[rows, , drop = FALSE], background$exclude_on,
    "background exclusion"
  )
  measures <- lapply(rows, function(row) {
    group <- if (is.null(model$spec$warp$fit_by)) {
      names(model$emission_models)[[1]]
    } else {
      gaze_key(support_tab[row, , drop = FALSE], model$spec$warp$fit_by,
               "warp fit_by")
    }
    prepared <- prepare_gaze_replay_source(
      support_tab[[sourcevar]][[row]], model, group
    )
    measure <- prepared$registered_source
    measure$coords <- gaze_replay_clamp_to_screen(
      measure$coords, background$screen
    )
    measure
  })
  gaze_replay_background_table(measures, rep(background_key, length(rows)),
                               exclude)
}
