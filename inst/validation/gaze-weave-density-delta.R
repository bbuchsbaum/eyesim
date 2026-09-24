# Density Delta court: held-out predictive gain of the participant's own
# four-presentation encoding density over group and background densities.
#
# Protocol: inst/validation/GAZEWEAVE-DENSITY-DELTA.md (frozen section above
# the FROZEN-PROTOCOL-END marker). Run from the eyesim source root after
# devtools::load_all():
#
#   source("inst/validation/gaze-weave-density-delta.R")
#   density_delta_freeze()                 # once, before any real scoring
#   run_density_delta_simulation()         # synthetic sweep (no real data)
#   run_density_delta_null1()              # semi-synthetic null (1)
#   run_density_delta_court()              # primary, nulls 2-3, sensitivity
#
# Participant-linked outputs stay in the Git-ignored results directory.

# Amendment 1 (1.1.0): explicit observation model, fold-disjoint background
# support, fixation-sum aggregation sensitivity, and an FPR-based null rule.
# Version 1.0.0 outputs remain at the top of the results directory.
density_delta_protocol <- "density-delta/1.1.0"
density_delta_version <- "1.1.0"
density_delta_seed <- 20260924L
density_delta_result_dir <- file.path(
  "inst", "validation", "gaze-weave-density-delta-results",
  paste0("v", density_delta_version)
)
density_delta_freeze_dir <- file.path(
  "inst", "validation", "gaze-weave-density-delta-freeze"
)
density_delta_simulation_dir <- file.path(
  "inst", "validation", "gaze-weave-density-delta-simulation"
)
density_delta_protocol_file <- file.path(
  "inst", "validation", "GAZEWEAVE-DENSITY-DELTA.md"
)
density_delta_script_file <- file.path(
  "inst", "validation", "gaze-weave-density-delta.R"
)
# The 1.1.0 frozen text runs through the amendment marker, so it contains the
# unchanged 1.0.0 section plus Amendment 1.
density_delta_marker <- "<!-- AMENDMENT-1-END -->"

# Reuse the own-group court: importer, cohort, folds, candidate plan, donor
# matching, and the crossed participant x item bootstrap helpers.
density_delta_dependency <- system.file(
  "validation", "gaze-weave-own-group-decomposition.R", package = "eyesim"
)
if (!nzchar(density_delta_dependency)) {
  density_delta_dependency <- file.path(
    "inst", "validation", "gaze-weave-own-group-decomposition.R"
  )
}
if (!file.exists(density_delta_dependency)) {
  stop("The own-group decomposition helpers are unavailable.")
}
source(density_delta_dependency, local = TRUE)

density_delta_config <- function() {
  primary_grid <- 20 * sqrt(2)^(0:6)
  list(
    protocol = density_delta_protocol,
    seed = density_delta_seed,
    screen = c(width = 800, height = 600),
    # Group and background bandwidths are selected on this grid; the own
    # bandwidth is selected on the same values, and the extended grid adds
    # two geometric steps on each side for the bandwidth sensitivity.
    bandwidth_grid = primary_grid,
    own_bandwidth_grid = 20 * sqrt(2)^(-2:8),
    own_primary_index = 3:9,
    uniform_floor = 0.01,
    em_max_iter = 5000L,
    em_tol = 1e-10,
    windows = list(combined = c(0, 3000), delay = c(500, 3000)),
    minimum_study_fixations = 3L,
    minimum_retrieval_fixations = 3L,
    alpha_one_sided = 0.025,
    bootstrap_draws = 2000L,
    simulation_bootstrap_draws = 499L,
    null_delta_tolerance = 0.01,
    # Trial aggregation: "duration_mean" (primary) or "fixation_sum"
    # (sensitivity). Fitting and scoring always use the same aggregation.
    aggregation = "duration_mean",
    # Background support: the participant's retrieval trials on items of the
    # other item fold only (never an evaluated item or candidate in the fold).
    background_support = "other_item_fold",
    null1_replicates = 200L,
    # FPR tolerance = alpha + 2 Monte-Carlo SE at alpha for the replicate count.
    null1_fpr_tolerance = 0.025 + 2 * sqrt(0.025 * 0.975 / 200),
    bandwidth_shifts = c(half = -2L, double = 2L),
    simulation_replicates = 100L,
    simulation_cells = density_delta_simulation_cells()
  )
}

density_delta_simulation_cells <- function() {
  strengths <- c(0, 0.05, 0.1, 0.2, 0.3)
  base <- list(participants = 36L, items = 36L, overlap = "low",
               heterogeneity = 0)
  variants <- list(
    base = base,
    overlap_high = utils::modifyList(base, list(overlap = "high")),
    heterogeneous = utils::modifyList(base, list(heterogeneity = 1)),
    n_24x36 = utils::modifyList(base, list(participants = 24L)),
    n_12x24 = utils::modifyList(base, list(participants = 12L, items = 24L))
  )
  cells <- lapply(names(variants), function(name) {
    v <- variants[[name]]
    data.frame(
      variant = name, strength = strengths,
      participants = v$participants, items = v$items,
      overlap = v$overlap, heterogeneity = v$heterogeneity,
      stringsAsFactors = FALSE
    )
  })
  cells <- do.call(rbind, cells)
  cells$cell <- seq_len(nrow(cells))
  cells
}

# ---------------------------------------------------------------------------
# Screen-normalised kernels and templates
# ---------------------------------------------------------------------------

density_delta_kernel_mass <- function(mx, my, h, screen) {
  (stats::pnorm((screen[["width"]] - mx) / h) - stats::pnorm(-mx / h)) *
    (stats::pnorm((screen[["height"]] - my) / h) - stats::pnorm(-my / h))
}

# Density of a weighted template at query points, one column per bandwidth.
# Every isotropic Gaussian kernel is divided by its on-screen mass, so each
# kernel, and therefore the template, integrates to one over the screen.
density_delta_template_density <- function(template, qx, qy, bandwidths,
                                           screen) {
  n <- length(qx)
  if (!n) return(matrix(numeric(0), 0L, length(bandwidths)))
  if (!length(template$x)) stop("A template has no fixations.")
  d2 <- outer(qx, template$x, "-")^2 + outer(qy, template$y, "-")^2
  out <- matrix(0, n, length(bandwidths))
  for (k in seq_along(bandwidths)) {
    h <- bandwidths[[k]]
    mass <- density_delta_kernel_mass(template$x, template$y, h, screen)
    kernel <- exp(-d2 / (2 * h^2)) / (2 * pi * h^2)
    out[, k] <- drop(kernel %*% (template$w / mass))
  }
  out
}

# Each episode is duration-normalised to unit mass before averaging, so a long
# presentation does not outweigh a short one.
# Templates are lists of equal-length x, y, w vectors with sum(w) == 1.
density_delta_episode_template <- function(episodes) {
  episodes <- Filter(function(e) length(e$x) > 0L && sum(e$duration) > 0,
                     episodes)
  if (!length(episodes)) stop("A template has no usable episode.")
  n <- length(episodes)
  list(
    x = unlist(lapply(episodes, `[[`, "x"), use.names = FALSE),
    y = unlist(lapply(episodes, `[[`, "y"), use.names = FALSE),
    w = unlist(lapply(episodes, function(e) e$duration / sum(e$duration) / n),
               use.names = FALSE)
  )
}

# Average already-normalised templates with equal weight (e.g. one per donor).
density_delta_pool_templates <- function(templates) {
  if (!length(templates)) stop("No templates to pool.")
  n <- length(templates)
  list(
    x = unlist(lapply(templates, `[[`, "x"), use.names = FALSE),
    y = unlist(lapply(templates, `[[`, "y"), use.names = FALSE),
    w = unlist(lapply(templates, function(t) t$w / sum(t$w) / n),
               use.names = FALSE)
  )
}

density_delta_sample_template <- function(template, n, h, screen) {
  if (!n) return(data.frame(x = numeric(0), y = numeric(0)))
  centre <- sample.int(length(template$x), n, replace = TRUE, prob = template$w)
  density_delta_truncated_normal(
    template$x[centre], template$y[centre], h, screen
  )
}

density_delta_truncated_normal <- function(mx, my, sd, screen) {
  n <- length(mx)
  sd <- rep_len(sd, n)
  x <- stats::rnorm(n, mx, sd)
  y <- stats::rnorm(n, my, sd)
  bad <- x < 0 | x > screen[["width"]] | y < 0 | y > screen[["height"]]
  while (any(bad)) {
    x[bad] <- stats::rnorm(sum(bad), mx[bad], sd[bad])
    y[bad] <- stats::rnorm(sum(bad), my[bad], sd[bad])
    bad <- x < 0 | x > screen[["width"]] | y < 0 | y > screen[["height"]]
  }
  data.frame(x = x, y = y)
}

# ---------------------------------------------------------------------------
# Design: per-fixation component densities for every trial and bandwidth
# ---------------------------------------------------------------------------

density_delta_features <- function(y, templates, bandwidths, screen) {
  do.call(rbind, lapply(seq_along(y), function(t) {
    density_delta_template_density(
      templates[[t]], y[[t]]$x, y[[t]]$y, bandwidths, screen
    )
  }))
}

density_delta_design <- function(trials, y, own, group, background, config) {
  n <- nrow(trials)
  if (length(y) != n || length(own) != n || length(group) != n ||
      length(background) != n) {
    stop("Design inputs must align with the trial table.")
  }
  if (any(vapply(y, nrow, integer(1)) < 1L)) {
    stop("Every scored trial needs at least one fixation.")
  }
  fix <- data.frame(
    trial = rep(seq_len(n), vapply(y, nrow, integer(1))),
    x = unlist(lapply(y, `[[`, "x"), use.names = FALSE),
    y = unlist(lapply(y, `[[`, "y"), use.names = FALSE),
    # Observation weight of each fixation within its trial. duration_mean:
    # d_j / D (trial score is the duration-weighted mean log density).
    # fixation_sum: 1 (trial score is the summed fixation log density).
    a = unlist(lapply(y, function(e) {
      if (identical(config$aggregation, "fixation_sum")) {
        rep(1, length(e$duration))
      } else {
        e$duration / sum(e$duration)
      }
    }), use.names = FALSE)
  )
  if (!config$aggregation %in% c("duration_mean", "fixation_sum")) {
    stop("Unknown trial aggregation.")
  }
  screen <- config$screen
  if (any(fix$x < 0 | fix$x > screen[["width"]] |
          fix$y < 0 | fix$y > screen[["height"]])) {
    stop("Retrieval fixations must lie on the screen rectangle.")
  }
  list(
    trials = trials, y = y, fix = fix,
    templates = list(own = own, group = group, background = background),
    g = density_delta_features(y, group, config$bandwidth_grid, screen),
    b = density_delta_features(y, background, config$bandwidth_grid, screen),
    o = density_delta_features(y, own, config$own_bandwidth_grid, screen),
    config = config
  )
}

# Replace the retrieval fixations (e.g. by simulated gaze) and recompute every
# feature against the unchanged templates.
density_delta_redesign <- function(design, y) {
  t <- design$templates
  density_delta_design(design$trials, y, t$own, t$group, t$background,
                       design$config)
}

# ---------------------------------------------------------------------------
# Population-level mixture fit (training fold only) and held-out scoring
# ---------------------------------------------------------------------------

# Maximise sum_j a_j log(floor * u + (1 - floor) * F_j w) over the simplex.
# The a_j sum to one within each trial, so the objective is the sum of
# per-trial duration-weighted mean log densities.
density_delta_em <- function(F, a, u, floor, max_iter, tol) {
  K <- ncol(F)
  w <- rep(1 / K, K)
  base <- floor * u
  previous <- -Inf
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    mix <- drop(F %*% w)
    p <- base + (1 - floor) * mix
    objective <- sum(a * log(p))
    if (objective - previous <= tol * (abs(objective) + 1)) {
      converged <- TRUE
      break
    }
    previous <- objective
    w <- w * colSums(F * (a / p))
    w <- w / sum(w)
  }
  list(weights = w, objective = objective, iterations = iter,
       converged = converged)
}

density_delta_fit <- function(design, train, own_shift = 0L) {
  cfg <- design$config
  rows <- design$fix$trial %in% train
  if (!any(rows)) stop("The training fold is empty.")
  a <- design$fix$a[rows]
  u <- 1 / (cfg$screen[["width"]] * cfg$screen[["height"]])
  G <- design$g[rows, , drop = FALSE]
  B <- design$b[rows, , drop = FALSE]
  O <- design$o[rows, , drop = FALSE]
  em <- function(F) {
    density_delta_em(F, a, u, cfg$uniform_floor, cfg$em_max_iter, cfg$em_tol)
  }
  grid <- expand.grid(ig = seq_along(cfg$bandwidth_grid),
                      ib = seq_along(cfg$bandwidth_grid))
  fits0 <- lapply(seq_len(nrow(grid)), function(k) {
    em(cbind(G[, grid$ig[[k]]], B[, grid$ib[[k]]], u))
  })
  best0 <- which.max(vapply(fits0, `[[`, numeric(1), "objective"))
  ig <- grid$ig[[best0]]
  ib <- grid$ib[[best0]]
  own_candidates <- cfg$own_primary_index
  fits1 <- lapply(own_candidates, function(io) em(cbind(G[, ig], B[, ib], O[, io], u)))
  best1 <- which.max(vapply(fits1, `[[`, numeric(1), "objective"))
  io <- own_candidates[[best1]]
  fit1 <- fits1[[best1]]
  if (own_shift != 0L) {
    io <- io + as.integer(own_shift)
    if (io < 1L || io > length(cfg$own_bandwidth_grid)) {
      stop("The own-bandwidth shift leaves the extended grid.")
    }
    fit1 <- em(cbind(G[, ig], B[, ib], O[, io], u))
  }
  fit0 <- fits0[[best0]]
  list(
    ig = ig, ib = ib, io = io,
    h_group = cfg$bandwidth_grid[[ig]],
    h_background = cfg$bandwidth_grid[[ib]],
    h_own = cfg$own_bandwidth_grid[[io]],
    weights0 = stats::setNames(fit0$weights, c("group", "background", "uniform")),
    weights1 = stats::setNames(
      fit1$weights, c("group", "background", "own", "uniform")
    ),
    train_trials = length(unique(design$fix$trial[rows])),
    train_objective0 = fit0$objective,
    train_objective1 = fit1$objective,
    converged = fit0$converged && fit1$converged
  )
}

density_delta_mixture <- function(design, fit, rows, model = c("base", "full")) {
  model <- match.arg(model)
  cfg <- design$config
  u <- 1 / (cfg$screen[["width"]] * cfg$screen[["height"]])
  g <- design$g[rows, fit$ig]
  b <- design$b[rows, fit$ib]
  if (model == "base") {
    w <- fit$weights0
    mix <- w[["group"]] * g + w[["background"]] * b + w[["uniform"]] * u
  } else {
    w <- fit$weights1
    mix <- w[["group"]] * g + w[["background"]] * b +
      w[["own"]] * design$o[rows, fit$io] + w[["uniform"]] * u
  }
  cfg$uniform_floor * u + (1 - cfg$uniform_floor) * mix
}

density_delta_score <- function(design, fit, eval) {
  rows <- which(design$fix$trial %in% eval)
  trial <- design$fix$trial[rows]
  a <- design$fix$a[rows]
  p0 <- density_delta_mixture(design, fit, rows, "base")
  p1 <- density_delta_mixture(design, fit, rows, "full")
  ll0 <- rowsum(a * log(p0), trial)
  ll1 <- rowsum(a * log(p1), trial)
  index <- as.integer(rownames(ll0))
  u <- 1 / prod(design$config$screen)
  data.frame(
    trial = index,
    loglik_base = drop(ll0), loglik_full = drop(ll1),
    delta = drop(ll1 - ll0), base_gain_over_uniform = drop(ll0) - log(u)
  )
}

density_delta_crossfit <- function(design, folds, own_shift = 0L) {
  trials <- design$trials
  scores <- list()
  fits <- list()
  for (fold in folds) {
    in_p <- as.character(trials$participant) %in% as.character(fold$eval_participants)
    in_i <- trials$item %in% fold$eval_items
    eval <- which(in_p & in_i)
    train <- which(!in_p & !in_i)
    if (!length(eval)) next
    fit <- density_delta_fit(design, train, own_shift)
    s <- density_delta_score(design, fit, eval)
    s$outer_fold <- fold$id
    scores[[length(scores) + 1L]] <- s
    fits[[length(fits) + 1L]] <- data.frame(
      outer_fold = fold$id, train_trials = fit$train_trials,
      eval_trials = length(eval),
      h_group = fit$h_group, h_background = fit$h_background,
      h_own = fit$h_own,
      w0_group = fit$weights0[["group"]],
      w0_background = fit$weights0[["background"]],
      w0_uniform = fit$weights0[["uniform"]],
      w1_group = fit$weights1[["group"]],
      w1_background = fit$weights1[["background"]],
      w1_own = fit$weights1[["own"]],
      w1_uniform = fit$weights1[["uniform"]],
      converged = fit$converged
    )
  }
  scores <- do.call(rbind, scores)
  if (anyDuplicated(scores$trial)) stop("A trial was evaluated twice.")
  scores <- cbind(trials[scores$trial, c("participant", "item"), drop = FALSE],
                  scores)
  rownames(scores) <- NULL
  list(scores = scores, fits = do.call(rbind, fits))
}

# ---------------------------------------------------------------------------
# Crossed participant x item inference
# ---------------------------------------------------------------------------

density_delta_preserve_rng <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  force(expr)
}

density_delta_boot_weights <- function(tab, draws, seed) {
  density_delta_preserve_rng({
    plan <- transport_exploratory_bootstrap_plan(tab, draws, seed)
    transport_exploratory_align_weights(tab, plan)
  })
}

density_delta_weighted_means <- function(value, weights) {
  total <- colSums(weights)
  out <- colSums(weights * value) / total
  out[total <= 0] <- NA_real_
  out
}

density_delta_crossed_test <- function(tab, value, label, draws, seed,
                                       alpha = 0.025) {
  if (length(value) != nrow(tab)) stop("Values must align with the trials.")
  weights <- density_delta_boot_weights(tab, draws, seed)
  boot <- density_delta_weighted_means(value, weights)
  interval <- transport_exploratory_interval(boot, 1 - 2 * alpha)
  data.frame(
    contrast = label, trials = nrow(tab),
    participants = length(unique(tab$participant)),
    items = length(unique(tab$item)),
    estimate = mean(value), lower_95 = interval[["lower"]],
    upper_95 = interval[["upper"]],
    probability_le_zero = mean(boot <= 0, na.rm = TRUE),
    reject_one_sided = interval[["lower"]] > 0,
    stringsAsFactors = FALSE
  )
}

# Unpaired difference of two trial sets (different participant-item pairs)
# under one crossed bootstrap over the union of participants and items.
density_delta_contrast_test <- function(first, second, label, draws, seed,
                                        alpha = 0.025) {
  combined <- rbind(
    data.frame(participant = first$participant, item = first$item,
               value = first$delta, set = 1L),
    data.frame(participant = second$participant, item = second$item,
               value = second$delta, set = 2L)
  )
  weights <- density_delta_boot_weights(combined, draws, seed)
  one <- combined$set == 1L
  boot <- density_delta_weighted_means(combined$value[one], weights[one, , drop = FALSE]) -
    density_delta_weighted_means(combined$value[!one], weights[!one, , drop = FALSE])
  interval <- transport_exploratory_interval(boot, 1 - 2 * alpha)
  data.frame(
    contrast = label, trials = nrow(combined),
    participants = length(unique(combined$participant)),
    items = length(unique(combined$item)),
    estimate = mean(first$delta) - mean(second$delta),
    lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
    probability_le_zero = mean(boot <= 0, na.rm = TRUE),
    reject_one_sided = interval[["lower"]] > 0,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Synthetic generator (tests and simulation sweep; no participant data)
# ---------------------------------------------------------------------------

density_delta_simulate <- function(
    participants = 36L, items = 36L, strength = 0, overlap = c("low", "high"),
    heterogeneity = 0, seed = 1L, screen = c(width = 800, height = 600),
    study_own_share = 0.5, group_donors = 9L, hotspot_sd = 35,
    background_sd = 130) {
  overlap <- match.arg(overlap)
  set.seed(seed)
  W <- screen[["width"]]
  H <- screen[["height"]]
  hot <- lapply(seq_len(items), function(i) {
    list(x = stats::runif(4, 0.1 * W, 0.9 * W),
         y = stats::runif(4, 0.1 * H, 0.9 * H),
         w = stats::rgamma(4, 2))
  })
  bg <- lapply(seq_len(participants), function(s) {
    c(x = W / 2 + stats::rnorm(1, 0, 60), y = H / 2 + stats::rnorm(1, 0, 45))
  })
  multiplier <- exp(heterogeneity * stats::rnorm(participants) -
                      heterogeneity^2 / 2)
  trials <- expand.grid(item = seq_len(items), participant = seq_len(participants))
  trials <- trials[c("participant", "item")]
  trials$participant <- sprintf("s%02d", trials$participant)
  own <- lapply(seq_len(nrow(trials)), function(t) {
    i <- trials$item[[t]]
    if (overlap == "low") {
      list(x = stats::runif(2, 0.1 * W, 0.9 * W),
           y = stats::runif(2, 0.1 * H, 0.9 * H), w = c(1, 1))
    } else {
      pref <- stats::rgamma(4, 0.3) + 1e-6
      list(x = hot[[i]]$x, y = hot[[i]]$y, w = pref)
    }
  })
  draw_spots <- function(spots, n) {
    k <- sample.int(length(spots$x), n, replace = TRUE, prob = spots$w)
    density_delta_truncated_normal(spots$x[k], spots$y[k], hotspot_sd, screen)
  }
  draw_bg <- function(s, n) {
    density_delta_truncated_normal(rep(bg[[s]][["x"]], n),
                                   rep(bg[[s]][["y"]], n), background_sd, screen)
  }
  draw_mixture <- function(n, probs, s, t) {
    source <- sample.int(3L, n, replace = TRUE, prob = probs)
    xy <- data.frame(x = numeric(n), y = numeric(n))
    if (any(source == 1L)) xy[source == 1L, ] <- draw_spots(hot[[trials$item[[t]]]], sum(source == 1L))
    if (any(source == 2L)) xy[source == 2L, ] <- draw_bg(s, sum(source == 2L))
    if (any(source == 3L)) xy[source == 3L, ] <- draw_spots(own[[t]], sum(source == 3L))
    xy$duration <- stats::rgamma(n, shape = 4, scale = 75) + 80
    xy
  }
  s_index <- match(trials$participant, sort(unique(trials$participant)))
  study <- lapply(seq_len(nrow(trials)), function(t) {
    lapply(1:4, function(p) {
      n <- 1L + stats::rpois(1, 4)
      draw_mixture(n, c((1 - study_own_share) * 0.9, 0.1, study_own_share),
                   s_index[[t]], t)
    })
  })
  own_strength <- pmin(strength * multiplier[s_index], 0.95)
  y <- lapply(seq_len(nrow(trials)), function(t) {
    n <- 3L + stats::rpois(1, 2.5)
    rest <- 1 - own_strength[[t]]
    draw_mixture(n, c(0.6 * rest, 0.4 * rest, own_strength[[t]]),
                 s_index[[t]], t)
  })
  list(trials = trials, study = study, y = y, own_strength = own_strength,
       group_donors = group_donors, screen = screen)
}

# own_source = "pseudo" replaces each own template by one other participant's
# four presentations of the same item and removes that donor from the group
# template, mirroring the real-data pseudo-own control.
density_delta_synthetic_design <- function(sim, folds,
                                           config = density_delta_config(),
                                           own_source = c("own", "pseudo")) {
  own_source <- match.arg(own_source)
  trials <- sim$trials
  study_templates <- lapply(sim$study, density_delta_episode_template)
  by_item <- split(seq_len(nrow(trials)), trials$item)
  donors <- lapply(seq_len(nrow(trials)), function(t) {
    others <- by_item[[as.character(trials$item[[t]])]]
    others <- others[trials$participant[others] != trials$participant[[t]]]
    # Deterministic donor order mirroring the per-version donor rotation.
    offset <- match(trials$participant[[t]], sort(unique(trials$participant)))
    others <- own_group_rotate(others, offset)
    others[seq_len(min(length(others), sim$group_donors + 1L))]
  })
  own <- if (own_source == "own") {
    study_templates
  } else {
    lapply(donors, function(d) study_templates[[d[[1L]]]])
  }
  group <- lapply(donors, function(d) {
    if (own_source == "pseudo") d <- d[-1L]
    d <- d[seq_len(min(length(d), sim$group_donors))]
    density_delta_pool_templates(study_templates[d])
  })
  item_fold <- density_delta_item_fold_map(folds)
  by_participant <- split(seq_len(nrow(trials)), trials$participant)
  background <- lapply(seq_len(nrow(trials)), function(t) {
    others <- by_participant[[trials$participant[[t]]]]
    support <- density_delta_background_support(
      trials$item[others], trials$item[[t]], item_fold
    )
    density_delta_episode_template(sim$y[others[support]])
  })
  density_delta_design(trials, sim$y, own, group, background, config)
}

# Item -> item-fold label, derived from the outer folds (distinct eval-item
# sets receive distinct labels).
density_delta_item_fold_map <- function(folds) {
  sets <- lapply(folds, function(f) sort(as.integer(f$eval_items)))
  keys <- vapply(sets, paste, character(1), collapse = ",")
  labels <- match(keys, unique(keys))
  map <- integer(0)
  for (k in seq_along(sets)) {
    map[as.character(sets[[k]])] <- labels[[k]]
  }
  map
}

# Amendment 1: background support is the participant's retrieval trials on
# items of the *other* item fold. No support trial is ever evaluated, or used
# as a candidate or target, in the fold where the target is evaluated, and the
# same rule applies to training, evaluation, null and simulated trials.
density_delta_background_support <- function(items, target, item_fold) {
  target_fold <- item_fold[as.character(target)]
  fold <- item_fold[as.character(items)]
  if (is.na(target_fold)) stop("The target item has no item fold.")
  unname(!is.na(fold) & fold != target_fold)
}

density_delta_synthetic_folds <- function(trials, seed) {
  density_delta_preserve_rng({
    set.seed(seed)
    participants <- sample(sort(unique(trials$participant)))
    items <- sample(sort(unique(trials$item)))
    pf <- rep(1:2, length.out = length(participants))
    itf <- rep(1:2, length.out = length(items))
    folds <- list()
    for (p in 1:2) for (i in 1:2) {
      folds[[length(folds) + 1L]] <- list(
        id = length(folds) + 1L,
        eval_participants = participants[pf == p],
        eval_items = items[itf == i]
      )
    }
    folds
  })
}

density_delta_simulation_replicate <- function(cell, replicate, config) {
  seed <- config$seed + 1000L * cell$cell + replicate
  sim <- density_delta_simulate(
    participants = cell$participants, items = cell$items,
    strength = cell$strength, overlap = cell$overlap,
    heterogeneity = cell$heterogeneity, seed = seed, screen = config$screen
  )
  folds <- density_delta_synthetic_folds(sim$trials, seed)
  own <- density_delta_crossfit(
    density_delta_synthetic_design(sim, folds, config), folds
  )
  pseudo <- density_delta_crossfit(
    density_delta_synthetic_design(sim, folds, config, "pseudo"), folds
  )
  if (!identical(own$scores$trial, pseudo$scores$trial)) {
    stop("Own and pseudo-own scores do not align.")
  }
  draws <- config$simulation_bootstrap_draws
  alpha <- config$alpha_one_sided
  tests <- list(
    delta = density_delta_crossed_test(
      own$scores, own$scores$delta, "delta", draws, seed, alpha),
    pseudo = density_delta_crossed_test(
      pseudo$scores, pseudo$scores$delta, "pseudo", draws, seed, alpha),
    own_minus_pseudo = density_delta_crossed_test(
      own$scores, own$scores$delta - pseudo$scores$delta, "own_minus_pseudo",
      draws, seed, alpha)
  )
  do.call(rbind, lapply(names(tests), function(test) {
    data.frame(
      cell = cell$cell, replicate = replicate, test = test,
      estimate = tests[[test]]$estimate, lower_95 = tests[[test]]$lower_95,
      upper_95 = tests[[test]]$upper_95,
      reject = tests[[test]]$reject_one_sided,
      mean_own_weight = mean(own$fits$w1_own),
      mean_pseudo_weight = mean(pseudo$fits$w1_own),
      stringsAsFactors = FALSE
    )
  }))
}

density_delta_rate_interval <- function(successes, n) {
  ci <- stats::binom.test(successes, n)$conf.int
  c(lower = ci[[1L]], upper = ci[[2L]])
}

density_delta_summarise_simulation <- function(reps, cells) {
  groups <- split(reps, list(reps$cell, reps$test), drop = TRUE)
  rows <- lapply(groups, function(part) {
    ci <- density_delta_rate_interval(sum(part$reject), nrow(part))
    data.frame(
      cell = part$cell[[1L]], test = part$test[[1L]], replicates = nrow(part),
      rejection_rate = mean(part$reject),
      rate_lower_95 = ci[["lower"]], rate_upper_95 = ci[["upper"]],
      mean_delta = mean(part$estimate), sd_delta = stats::sd(part$estimate),
      mean_own_weight = mean(part$mean_own_weight)
    )
  })
  out <- merge(cells, do.call(rbind, rows), by = "cell")
  out$quantity <- ifelse(
    out$strength == 0 | out$test == "pseudo", "false_positive_rate", "power"
  )
  out[order(out$cell, out$test), ]
}

run_density_delta_simulation <- function(
    config = density_delta_config(), cores = 6L,
    output_dir = density_delta_simulation_dir,
    replicates = config$simulation_replicates) {
  cells <- config$simulation_cells
  jobs <- expand.grid(replicate = seq_len(replicates), cell = cells$cell)
  started <- proc.time()[["elapsed"]]
  reps <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
    cell <- cells[cells$cell == jobs$cell[[k]], , drop = FALSE]
    density_delta_simulation_replicate(cell, jobs$replicate[[k]], config)
  }, mc.cores = cores, mc.preschedule = FALSE)
  failed <- vapply(reps, inherits, logical(1), "try-error")
  if (any(failed)) stop("Simulation replicates failed: ", sum(failed))
  reps <- do.call(rbind, reps)
  summary <- density_delta_summarise_simulation(reps, cells)
  summary$elapsed_seconds_total <- proc.time()[["elapsed"]] - started
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, file.path(output_dir, "simulation-summary.csv"),
                   row.names = FALSE)
  utils::write.csv(reps, file.path(output_dir, "simulation-replicates.csv"),
                   row.names = FALSE)
  invisible(summary)
}

# ---------------------------------------------------------------------------
# Real data: importer, cohort, folds (reused), and templates
# ---------------------------------------------------------------------------

density_delta_data_dir <- function() {
  Sys.getenv("EYESIM_WYNN_DATA_DIR",
              file.path("test_data", "wynn_probe_delay"))
}

density_delta_read_inputs <- function(data_dir = density_delta_data_dir()) {
  probe_delay_verify_inputs(data_dir)
  files <- probe_delay_data_files(data_dir)
  study_raw <- probe_delay_read_csv(files[["study"]], "study")
  retrieval_raw <- probe_delay_read_csv(files[["retrieval"]], "retrieval")
  list(
    study = probe_delay_study_paths(study_raw),
    retrieval = probe_delay_retrieval_paths(retrieval_raw, "combined"),
    retrieval_delay = probe_delay_retrieval_paths(retrieval_raw, "delay"),
    trials = probe_delay_trial_tables(study_raw, retrieval_raw)
  )
}

density_delta_context <- function(data_dir = density_delta_data_dir(),
                                  cache_dir = density_delta_result_dir) {
  cache <- file.path(cache_dir, "context-cache.rds")
  if (file.exists(cache)) {
    context <- readRDS(cache)
    if (identical(context$protocol, density_delta_protocol)) return(context)
  }
  raw <- density_delta_read_inputs(data_dir)
  cohort <- full_recognition_select_cohort(raw)
  fold_plan <- full_recognition_fold_plan(cohort)
  candidate_plan <- full_recognition_candidate_plan(cohort)
  full_recognition_validate_design(cohort, candidate_plan, fold_plan)
  context <- list(
    protocol = density_delta_protocol, raw = raw, cohort = cohort,
    fold_plan = fold_plan, candidate_plan = candidate_plan,
    study = own_group_complete_study(raw$study)
  )
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(context, cache, version = 3)
  context
}

density_delta_path_df <- function(fixgroup) {
  data.frame(x = fixgroup$x, y = fixgroup$y, duration = fixgroup$duration)
}

density_delta_study_index <- function(study) {
  key <- paste(study$participant, study$item, study$study_image_version, sep = ":")
  split(seq_len(nrow(study)), key)
}

density_delta_own_template <- function(context, index, participant, item, version) {
  rows <- index[[paste(participant, item, version, sep = ":")]]
  if (length(rows) != 4L) stop("Own template requires four presentations.")
  part <- context$study[rows, , drop = FALSE]
  part <- part[order(part$presentation), , drop = FALSE]
  density_delta_episode_template(lapply(part$fixgroup, density_delta_path_df))
}

density_delta_donors <- function(context, item, version) {
  keep <- context$study$item == item &
    as.character(context$study$study_image_version) == as.character(version)
  sort(unique(as.character(context$study$participant[keep])))
}

density_delta_group_template <- function(context, index, item, version, exclude) {
  donors <- setdiff(density_delta_donors(context, item, version), exclude)
  if (!length(donors)) stop("No group donor is available.")
  list(
    template = density_delta_pool_templates(lapply(donors, function(d) {
      density_delta_own_template(context, index, d, item, version)
    })),
    donors = length(donors)
  )
}

density_delta_retrieval_table <- function(context, window) {
  if (window == "combined") context$raw$retrieval else context$raw$retrieval_delay
}

density_delta_background_templates <- function(context, window) {
  retrieval <- density_delta_retrieval_table(context, window)
  by_participant <- split(seq_len(nrow(retrieval)), retrieval$participant)
  lapply(by_participant, function(rows) {
    list(items = retrieval$item[rows],
         paths = lapply(retrieval$fixgroup[rows], density_delta_path_df))
  })
}

density_delta_background <- function(bg, participant, item, item_fold) {
  part <- bg[[as.character(participant)]]
  keep <- density_delta_background_support(part$items, item, item_fold)
  density_delta_episode_template(part$paths[keep])
}

density_delta_trials <- function(context, family, window) {
  retrieval <- density_delta_retrieval_table(context, window)
  cohort <- context$cohort
  if (family == "old_lure") {
    trials <- cohort$pairs[c("participant", "item", "probe_type")]
  } else {
    comb <- context$raw$retrieval
    keep <- comb$probe_type == "newtest" &
      comb$nfix >= context$config_min_retrieval &
      as.character(comb$participant) %in% as.character(cohort$participants) &
      comb$item %in% cohort$item_map$item
    trials <- comb[keep, c("participant", "item", "probe_type",
                           "retrieval_image_version")]
  }
  key <- paste(trials$participant, trials$item, sep = ":")
  rkey <- paste(retrieval$participant, retrieval$item, sep = ":")
  if (anyDuplicated(rkey[rkey %in% key])) stop("Retrieval trials are not unique.")
  position <- match(key, rkey)
  trials <- trials[!is.na(position), , drop = FALSE]
  trials$row <- position[!is.na(position)]
  trials <- trials[order(trials$participant, trials$item), , drop = FALSE]
  rownames(trials) <- NULL
  trials
}

# variant: primary | wrong_item | pseudo_own | donor_support
density_delta_real_design <- function(
    context, family = c("old_lure", "newtest"),
    variant = c("primary", "wrong_item", "pseudo_own", "donor_support"),
    window = c("combined", "delay"), config = density_delta_config()) {
  family <- match.arg(family)
  variant <- match.arg(variant)
  window <- match.arg(window)
  if (family == "newtest" && variant != "pseudo_own") {
    stop("Newtest trials support only pseudo-own templates.")
  }
  context$config_min_retrieval <- config$minimum_retrieval_fixations
  trials <- density_delta_trials(context, family, window)
  retrieval <- density_delta_retrieval_table(context, window)
  index <- density_delta_study_index(context$study)
  bg <- density_delta_background_templates(context, window)
  item_fold <- stats::setNames(
    as.integer(context$cohort$item_map$item_fold),
    as.character(context$cohort$item_map$item)
  )
  plan <- context$candidate_plan
  n <- nrow(trials)
  own <- group <- background <- vector("list", n)
  y <- vector("list", n)
  trials$group_donors <- NA_integer_
  trials$template_item <- trials$item
  for (t in seq_len(n)) {
    p <- as.character(trials$participant[[t]])
    item <- as.integer(trials$item[[t]])
    y[[t]] <- density_delta_path_df(retrieval$fixgroup[[trials$row[[t]]]])
    background[[t]] <- density_delta_background(bg, p, item, item_fold)
    version <- if (family == "old_lure") {
      own_group_version_lookup(context$study, p)(item)
    } else {
      own_group_newtest_version_lookup(
        data.frame(retrieval_image_version = trials$retrieval_image_version[[t]])
      )(item)
    }
    exclude <- p
    if (variant == "pseudo_own") {
      pseudo <- own_group_select_donors(
        context$study, p, item, version
      )$donor_participant[[1L]]
      own[[t]] <- density_delta_own_template(context, index, pseudo, item, version)
      exclude <- c(p, pseudo)
    } else if (variant == "wrong_item") {
      rows <- plan[plan$candidate_set_id == own_group_trial_key(p, item) &
                     !plan$is_true, , drop = FALSE]
      wrong <- rows$candidate_item[which.min(rows$candidate_position)]
      trials$template_item[[t]] <- wrong
      wrong_version <- own_group_version_lookup(context$study, p)(wrong)
      own[[t]] <- density_delta_own_template(context, index, p, wrong, wrong_version)
    } else {
      own[[t]] <- density_delta_own_template(context, index, p, item, version)
    }
    if (variant == "donor_support") {
      donors <- own_group_select_donors(context$study, p, item, version)
      group[[t]] <- density_delta_episode_template(
        lapply(donors$fixgroup, density_delta_path_df)
      )
      trials$group_donors[[t]] <- length(unique(donors$donor_participant))
    } else {
      g <- density_delta_group_template(context, index, item, version, exclude)
      group[[t]] <- g$template
      trials$group_donors[[t]] <- g$donors
    }
  }
  design <- density_delta_design(trials, y, own, group, background, config)
  design$family <- family
  design$variant <- variant
  design$window <- window
  design
}

density_delta_real_folds <- function(context) context$fold_plan$folds

# ---------------------------------------------------------------------------
# Freeze: protocol text, script, and configuration hashed before scoring
# ---------------------------------------------------------------------------

density_delta_protocol_text <- function(path = density_delta_protocol_file) {
  lines <- readLines(path, warn = FALSE)
  end <- match(density_delta_marker, lines)
  if (is.na(end)) stop("The protocol lacks the frozen-section marker.")
  lines[seq_len(end)]
}

density_delta_hash_lines <- function(lines) {
  path <- tempfile(fileext = ".txt")
  on.exit(unlink(path))
  writeLines(lines, path, useBytes = TRUE)
  probe_delay_sha256(path)
}

density_delta_hashes <- function(config = density_delta_config()) {
  c(
    protocol_section = density_delta_hash_lines(density_delta_protocol_text()),
    script = probe_delay_sha256(density_delta_script_file),
    configuration = density_delta_hash_lines(
      utils::capture.output(dput(config, control = c("keepNA", "keepInteger")))
    )
  )
}

# Version 1.0.0 was frozen to freeze-record.csv; later versions are suffixed.
density_delta_freeze_file <- function(version = density_delta_version) {
  if (identical(version, "1.0.0")) "freeze-record.csv" else
    paste0("freeze-record-", version, ".csv")
}

density_delta_freeze <- function(output_dir = density_delta_freeze_dir,
                                 config = density_delta_config(),
                                 data_dir = density_delta_data_dir()) {
  path <- file.path(output_dir, density_delta_freeze_file())
  if (file.exists(path)) stop("The density-delta court is already frozen.")
  hashes <- density_delta_hashes(config)
  manifest <- utils::read.csv(file.path(data_dir, "manifest.csv"),
                              stringsAsFactors = FALSE)
  record <- data.frame(
    component = c(names(hashes), paste0("input:", manifest$file)),
    sha256 = c(unname(hashes), manifest$sha256),
    stringsAsFactors = FALSE
  )
  record$protocol <- density_delta_protocol
  record$frozen_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(record, path, row.names = FALSE)
  writeLines(
    utils::capture.output(dput(config, control = c("keepNA", "keepInteger"))),
    file.path(output_dir, paste0("configuration-", density_delta_version, ".txt"))
  )
  invisible(record)
}

density_delta_verify_freeze <- function(output_dir = density_delta_freeze_dir,
                                        config = density_delta_config()) {
  path <- file.path(output_dir, density_delta_freeze_file())
  if (!file.exists(path)) stop("Freeze the protocol before any real scoring.")
  record <- utils::read.csv(path, stringsAsFactors = FALSE)
  hashes <- density_delta_hashes(config)
  frozen <- stats::setNames(record$sha256, record$component)[names(hashes)]
  if (!identical(unname(frozen), unname(hashes))) {
    stop("Protocol, script, or configuration changed after the freeze: ",
         paste(names(hashes)[frozen != hashes], collapse = ", "))
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Null (1): simulated group + background gaze, no own contribution
# ---------------------------------------------------------------------------

density_delta_null1_draw <- function(design, fits, seed) {
  set.seed(seed)
  cfg <- design$config
  fold_of <- design$fold_of
  lapply(seq_len(nrow(design$trials)), function(t) {
    fit <- fits[[as.character(fold_of[[t]])]]
    old <- design$y[[t]]
    n <- nrow(old)
    w <- c(cfg$uniform_floor,
           (1 - cfg$uniform_floor) * unname(fit$weights0))
    source <- sample.int(4L, n, replace = TRUE, prob = w)
    xy <- data.frame(x = stats::runif(n, 0, cfg$screen[["width"]]),
                     y = stats::runif(n, 0, cfg$screen[["height"]]))
    k <- source == 2L
    if (any(k)) xy[k, ] <- density_delta_sample_template(
      design$templates$group[[t]], sum(k), fit$h_group, cfg$screen)
    k <- source == 3L
    if (any(k)) xy[k, ] <- density_delta_sample_template(
      design$templates$background[[t]], sum(k), fit$h_background, cfg$screen)
    xy$duration <- old$duration
    xy
  })
}

density_delta_generating_fits <- function(design, folds) {
  fits <- list()
  fold_of <- integer(nrow(design$trials))
  for (fold in folds) {
    in_p <- as.character(design$trials$participant) %in% as.character(fold$eval_participants)
    in_i <- design$trials$item %in% fold$eval_items
    fits[[as.character(fold$id)]] <- density_delta_fit(design, which(!in_p & !in_i))
    fold_of[in_p & in_i] <- fold$id
  }
  list(fits = fits, fold_of = fold_of)
}

run_density_delta_null1 <- function(
    config = density_delta_config(), cores = 6L,
    output_dir = density_delta_result_dir,
    replicates = config$null1_replicates, verify = TRUE) {
  if (verify) density_delta_verify_freeze(config = config)
  context <- density_delta_context()
  design <- density_delta_real_design(context, "old_lure", "primary",
                                      "combined", config)
  folds <- density_delta_real_folds(context)
  gen <- density_delta_generating_fits(design, folds)
  design$fold_of <- gen$fold_of
  started <- proc.time()[["elapsed"]]
  reps <- parallel::mclapply(seq_len(replicates), function(r) {
    seed <- config$seed + 500000L + r
    y <- density_delta_null1_draw(design, gen$fits, seed)
    null_design <- density_delta_redesign(design, y)
    cf <- density_delta_crossfit(null_design, folds)
    test <- density_delta_crossed_test(
      cf$scores, cf$scores$delta, "null1", config$simulation_bootstrap_draws,
      seed, config$alpha_one_sided
    )
    data.frame(replicate = r, estimate = test$estimate,
               lower_95 = test$lower_95, upper_95 = test$upper_95,
               reject = test$reject_one_sided,
               mean_own_weight = mean(cf$fits$w1_own))
  }, mc.cores = cores, mc.preschedule = FALSE)
  if (any(vapply(reps, inherits, logical(1), "try-error"))) {
    stop("Null (1) replicates failed.")
  }
  reps <- do.call(rbind, reps)
  ci <- density_delta_rate_interval(sum(reps$reject), nrow(reps))
  summary <- data.frame(
    null = "null1_simulated_group_background", replicates = nrow(reps),
    false_positive_rate = mean(reps$reject),
    rate_lower_95 = ci[["lower"]], rate_upper_95 = ci[["upper"]],
    mean_delta = mean(reps$estimate), sd_delta = stats::sd(reps$estimate),
    mean_own_weight = mean(reps$mean_own_weight),
    fpr_tolerance = config$null1_fpr_tolerance,
    # Amendment 1: judged on the decision rule's FPR only. E[Delta] under
    # this null is minus a KL divergence, not zero, so mean Delta is reported
    # descriptively and is not a gate.
    pass = mean(reps$reject) <= config$null1_fpr_tolerance,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    stringsAsFactors = FALSE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, file.path(output_dir, "null1-summary.csv"),
                   row.names = FALSE)
  utils::write.csv(reps, file.path(output_dir, "null1-replicates.csv"),
                   row.names = FALSE)
  invisible(summary)
}

# ---------------------------------------------------------------------------
# Real court: primary, nulls (2) and (3), secondary, sensitivity
# ---------------------------------------------------------------------------

density_delta_run_analysis <- function(context, label, family, variant, window,
                                       config, own_shift = 0L) {
  design <- density_delta_real_design(context, family, variant, window, config)
  cf <- density_delta_crossfit(design, density_delta_real_folds(context),
                               own_shift)
  cf$scores$probe_type <- design$trials$probe_type[cf$scores$trial]
  cf$scores$analysis <- label
  cf$fits$analysis <- label
  cf$support <- data.frame(
    analysis = label, trials = nrow(design$trials),
    min_group_donors = min(design$trials$group_donors),
    median_group_donors = stats::median(design$trials$group_donors),
    max_group_donors = max(design$trials$group_donors),
    template_item_differs = sum(design$trials$template_item != design$trials$item),
    stringsAsFactors = FALSE
  )
  cf
}

density_delta_align <- function(first, second) {
  key1 <- paste(first$participant, first$item)
  key2 <- paste(second$participant, second$item)
  common <- intersect(key1, key2)
  list(first = first[match(common, key1), ], second = second[match(common, key2), ])
}

density_delta_verdict <- function(summary, null1, config) {
  row <- function(label) summary[summary$contrast == label, , drop = FALSE]
  primary <- row("primary_old_lure")
  gates <- data.frame(
    gate = c(
      "primary_rejects", "null1_pass",
      "null2_wrong_item_within_tolerance", "null2_primary_exceeds",
      "null3_pseudo_newtest_within_tolerance", "null3_primary_exceeds",
      "attribution_own_minus_pseudo"
    ),
    pass = c(
      isTRUE(primary$reject_one_sided),
      isTRUE(null1$pass),
      isTRUE(row("null2_wrong_item")$upper_95 <= config$null_delta_tolerance) &&
        !isTRUE(row("null2_wrong_item")$reject_one_sided),
      isTRUE(row("primary_minus_null2_paired")$reject_one_sided),
      isTRUE(row("null3_pseudo_own_newtest")$upper_95 <= config$null_delta_tolerance) &&
        !isTRUE(row("null3_pseudo_own_newtest")$reject_one_sided),
      isTRUE(row("primary_minus_null3_unpaired")$reject_one_sided),
      isTRUE(row("attribution_own_minus_pseudo_paired")$reject_one_sided)
    ),
    stringsAsFactors = FALSE
  )
  gates$verdict <- if (all(gates$pass)) {
    "supports_participant_specific_correspondence"
  } else if (isTRUE(primary$reject_one_sided)) {
    "primary_positive_but_not_attributable"
  } else {
    "not_detected_at_protocol_sensitivity"
  }
  gates
}

run_density_delta_court <- function(config = density_delta_config(),
                                    output_dir = density_delta_result_dir,
                                    verify = TRUE, analyses = NULL) {
  if (verify) density_delta_verify_freeze(config = config)
  context <- density_delta_context()
  plan <- list(
    primary = list("old_lure", "primary", "combined", 0L),
    null2_wrong_item = list("old_lure", "wrong_item", "combined", 0L),
    null3_pseudo_own_newtest = list("newtest", "pseudo_own", "combined", 0L),
    attribution_pseudo_own_old_lure = list("old_lure", "pseudo_own", "combined", 0L),
    sensitivity_delay = list("old_lure", "primary", "delay", 0L),
    sensitivity_donor_support = list("old_lure", "donor_support", "combined", 0L),
    sensitivity_own_bandwidth_half = list(
      "old_lure", "primary", "combined", config$bandwidth_shifts[["half"]]),
    sensitivity_own_bandwidth_double = list(
      "old_lure", "primary", "combined", config$bandwidth_shifts[["double"]]),
    sensitivity_fixation_sum = list(
      "old_lure", "primary", "combined", 0L, "fixation_sum")
  )
  if (!is.null(analyses)) plan <- plan[analyses]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list()
  for (label in names(plan)) {
    p <- plan[[label]]
    path <- file.path(output_dir, paste0("analysis-", label, ".rds"))
    if (file.exists(path)) {
      results[[label]] <- readRDS(path)
      next
    }
    message("Density delta: ", label)
    started <- proc.time()[["elapsed"]]
    analysis_config <- config
    if (length(p) >= 5L) analysis_config$aggregation <- p[[5L]]
    results[[label]] <- density_delta_run_analysis(
      context, label, p[[1L]], p[[2L]], p[[3L]], analysis_config, p[[4L]]
    )
    results[[label]]$elapsed_seconds <- proc.time()[["elapsed"]] - started
    saveRDS(results[[label]], path, version = 3)
  }
  draws <- config$bootstrap_draws
  seed <- config$seed
  alpha <- config$alpha_one_sided
  test <- function(scores, label) {
    density_delta_crossed_test(scores, scores$delta, label, draws, seed, alpha)
  }
  rows <- list()
  primary <- results$primary$scores
  rows$primary <- test(primary, "primary_old_lure")
  for (probe in c("old", "lure")) {
    rows[[probe]] <- test(primary[primary$probe_type == probe, ],
                          paste0("primary_", probe))
  }
  if (!is.null(results$null2_wrong_item)) {
    null2 <- results$null2_wrong_item$scores
    rows$null2 <- test(null2, "null2_wrong_item")
    al <- density_delta_align(primary, null2)
    rows$null2_paired <- density_delta_crossed_test(
      al$first, al$first$delta - al$second$delta,
      "primary_minus_null2_paired", draws, seed, alpha)
  }
  if (!is.null(results$null3_pseudo_own_newtest)) {
    null3 <- results$null3_pseudo_own_newtest$scores
    rows$null3 <- test(null3, "null3_pseudo_own_newtest")
    rows$null3_contrast <- density_delta_contrast_test(
      primary, null3, "primary_minus_null3_unpaired", draws, seed, alpha)
  }
  if (!is.null(results$attribution_pseudo_own_old_lure)) {
    pseudo <- results$attribution_pseudo_own_old_lure$scores
    rows$pseudo <- test(pseudo, "attribution_pseudo_own_old_lure")
    al <- density_delta_align(primary, pseudo)
    rows$pseudo_paired <- density_delta_crossed_test(
      al$first, al$first$delta - al$second$delta,
      "attribution_own_minus_pseudo_paired", draws, seed, alpha)
  }
  for (label in grep("^sensitivity_", names(results), value = TRUE)) {
    rows[[label]] <- test(results[[label]]$scores, label)
  }
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  fits <- do.call(rbind, lapply(results, `[[`, "fits"))
  support <- do.call(rbind, lapply(results, `[[`, "support"))
  rownames(fits) <- rownames(support) <- NULL
  null1_path <- file.path(output_dir, "null1-summary.csv")
  null1 <- if (file.exists(null1_path)) utils::read.csv(null1_path) else NULL
  verdict <- if (is.null(analyses)) density_delta_verdict(summary, null1, config) else NULL
  utils::write.csv(summary, file.path(output_dir, "court-summary.csv"), row.names = FALSE)
  utils::write.csv(fits, file.path(output_dir, "fold-fits.csv"), row.names = FALSE)
  utils::write.csv(support, file.path(output_dir, "support.csv"), row.names = FALSE)
  if (!is.null(verdict)) {
    utils::write.csv(verdict, file.path(output_dir, "verdict.csv"), row.names = FALSE)
  }
  invisible(list(summary = summary, fits = fits, support = support,
                 verdict = verdict))
}
