# Frozen synthetic and invariance court for GazeWeave v2.
#
# Protocol: inst/validation/GAZEWEAVE-V2-COURT.md
# Run from a source checkout after devtools::load_all():
#   source("inst/validation/gaze-weave-v2-court.R")
#   run_gaze_weave_v2_court("inst/validation/gaze-weave-v2-results")

v2_court_path <- function(coords, duration = NULL, onset = NULL) {
  coords <- as.matrix(coords)
  if (is.null(duration)) duration <- rep(1, nrow(coords))
  if (is.null(onset)) onset <- c(0, head(cumsum(duration), -1L))
  fixation_group(
    x = coords[, 1], y = coords[, 2],
    duration = duration, onset = onset
  )
}

v2_court_permutations <- function(n_items, n_fixations) {
  base <- seq_len(n_fixations)
  lapply(seq_len(n_items), function(item) {
    shift <- (item - 1L) %% n_fixations
    shifted <- c(base[(shift + 1L):n_fixations],
                 if (shift > 0L) base[seq_len(shift)])
    if (item %% 2L == 0L && n_fixations >= 4L) {
      shifted[c(2L, 3L)] <- shifted[c(3L, 2L)]
    }
    shifted
  })
}

v2_court_templates <- function(n_items = 6L, n_fixations = 6L,
                               family, replicate = 1L, seed = 20260816L) {
  stopifnot(n_items >= 2L, n_fixations >= 4L)
  set.seed(seed + replicate * 1009L + match(
    family, c("transport_favouring", "replay_favouring",
              "neutral_misspecified")
  ) * 7919L)
  common <- rbind(
    c(5.0, 4.0), c(8.0, 13.5), c(11.0, 6.0),
    c(15.0, 14.0), c(19.0, 5.0), c(17.0, 10.0)
  )[seq_len(n_fixations), , drop = FALSE]
  permutations <- v2_court_permutations(n_items, n_fixations)

  lapply(seq_len(n_items), function(item) {
    if (identical(family, "neutral_misspecified")) {
      coords <- common[permutations[[item]], , drop = FALSE]
      duration <- rep(1, n_fixations)
    } else {
      theta <- seq(0, 2 * pi, length.out = n_fixations + 1L)[-1L]
      phase <- item * pi / (n_items + 1)
      radial <- 2.7 + ((seq_len(n_fixations) + item) %% 3) * 0.8
      center <- c(
        12 + 2.4 * cos(2 * pi * item / n_items),
        9 + 1.8 * sin(2 * pi * item / n_items)
      )
      coords <- cbind(
        center[[1]] + radial * cos(theta + phase),
        center[[2]] + 0.8 * radial * sin(theta + phase)
      )
      coords <- coords[permutations[[item]], , drop = FALSE]
      duration <- 0.65 + ((seq_len(n_fixations) + 2L * item) %% 5) / 4
    }
    v2_court_path(coords, duration)
  })
}

v2_inverse_warp <- function(coords, scale, translation) {
  sweep(as.matrix(coords), 2L, translation, FUN = "-") / scale
}

v2_court_signal <- function(reference, family, scale, translation, seed) {
  set.seed(seed)
  coords <- cbind(reference$x, reference$y)
  duration <- reference$duration
  n <- nrow(coords)

  if (identical(family, "transport_favouring")) {
    keep <- seq_len(n) != ((seed %% n) + 1L)
    coords <- coords[keep, , drop = FALSE]
    duration <- duration[keep]
    if (nrow(coords) >= 4L) coords[c(3L, 4L), ] <- coords[c(4L, 3L), ]
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 0.28), ncol = 2)
  } else if (identical(family, "replay_favouring")) {
    first <- seq_len(max(2L, floor(n / 2)))
    second <- seq.int(max(first) + 1L, n)
    order <- c(second, first)
    coords <- coords[order, , drop = FALSE]
    duration <- duration[order]
    background <- matrix(c(12, 9) + stats::rnorm(2L, sd = 2), nrow = 1L)
    coords <- rbind(coords, background)
    duration <- c(duration, stats::median(duration))
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 0.22), ncol = 2)
  } else {
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 0.32), ncol = 2)
    duration <- rep(1, nrow(coords))
  }

  raw <- v2_inverse_warp(coords, scale, translation)
  v2_court_path(raw, duration)
}

v2_court_null <- function(n_fixations, duration, scale, translation, seed) {
  set.seed(seed)
  angle <- stats::runif(n_fixations, 0, 2 * pi)
  radius <- abs(stats::rnorm(n_fixations, mean = 1.8, sd = 0.8))
  registered <- cbind(12 + radius * cos(angle), 9 + radius * sin(angle))
  raw <- v2_inverse_warp(registered, scale, translation)
  v2_court_path(raw, rep(duration, length.out = n_fixations))
}

v2_court_data <- function(n_participants = 1L, n_items = 6L,
                          n_replications = 2L, n_fixations = 6L,
                          seed = 20260816L) {
  families <- c(
    "transport_favouring", "replay_favouring", "neutral_misspecified"
  )
  reference_rows <- list()
  source_rows <- list()
  reference_index <- 1L
  source_index <- 1L
  truth <- list(scale = numeric(n_participants), translation = vector("list", n_participants))

  for (participant in seq_len(n_participants)) {
    scale <- 0.76 + 0.02 * participant
    translation <- c(1.1 + 0.15 * participant, -0.8 + 0.1 * participant)
    truth$scale[[participant]] <- scale
    truth$translation[[participant]] <- translation
    for (family in families) {
      for (replicate in seq_len(n_replications)) {
        templates <- v2_court_templates(
          n_items, n_fixations, family, replicate,
          seed + participant * 100003L
        )
        for (item in seq_len(n_items)) {
          id <- list(
            participant = paste0("p", participant),
            family = family,
            replicate = replicate,
            item = item
          )
          reference_rows[[reference_index]] <- c(
            id, list(fixgroup = templates[[item]])
          )
          reference_index <- reference_index + 1L
          signal <- v2_court_signal(
            templates[[item]], family, scale, translation,
            seed + participant * 1000003L + replicate * 10007L + item * 101L +
              match(family, families) * 1009L
          )
          null <- v2_court_null(
            nrow(signal), signal$duration, scale, translation,
            seed + participant * 2000003L + replicate * 20011L + item * 211L +
              match(family, families) * 2017L
          )
          source_rows[[source_index]] <- c(
            id, list(condition = "signal", fixgroup = signal)
          )
          source_rows[[source_index + 1L]] <- c(
            id, list(condition = "null", fixgroup = null)
          )
          source_index <- source_index + 2L
        }
      }
    }
  }

  to_tibble <- function(rows) {
    tibble::tibble(
      participant = vapply(rows, `[[`, character(1), "participant"),
      family = vapply(rows, `[[`, character(1), "family"),
      replicate = vapply(rows, `[[`, integer(1), "replicate"),
      item = vapply(rows, `[[`, integer(1), "item"),
      fixgroup = lapply(rows, `[[`, "fixgroup")
    )
  }
  reference <- to_tibble(reference_rows)
  source <- tibble::tibble(
    participant = vapply(source_rows, `[[`, character(1), "participant"),
    family = vapply(source_rows, `[[`, character(1), "family"),
    replicate = vapply(source_rows, `[[`, integer(1), "replicate"),
    item = vapply(source_rows, `[[`, integer(1), "item"),
    condition = vapply(source_rows, `[[`, character(1), "condition"),
    fixgroup = lapply(source_rows, `[[`, "fixgroup")
  )
  list(reference = reference, source = source, truth = truth)
}

v2_court_specs <- function(smoke = FALSE) {
  screen <- gaze_screen(24, 18, unit = "deg")
  warp <- gaze_warp_contraction(
    center = "screen", translation = TRUE, fit_by = "participant"
  )
  transport <- gaze_transport_v2_spec(
    spatial = gaze_gaussian_mixture(c(0.75, 1.5, 3), unit = "deg"),
    chronology = gaze_order_neighbours(2, unit = "deg"),
    coverage_grid = if (smoke) c(0.5, 1) else c(0.5, 0.75, 1),
    coverage_penalty_grid = c(0.25, 0.75, 1.5),
    temporal_weight = 1,
    entropy_schedule = if (smoke) 0.03 else c(0.05, 0.015),
    warp = warp, screen = screen,
    maxit = if (smoke) 60L else 80L,
    projection_maxit = if (smoke) 200L else 350L,
    tolerance = 5e-4,
    projection_tolerance = 1e-6,
    multistart = 1L
  )
  replay <- gaze_replay_spec(
    grid_size = if (smoke) 24L else 48L,
    max_skip = 2L,
    student_df = 4,
    scale_floor = 0.15,
    transition_grid = list(
      background = c(0.03, 0.10),
      restart = c(0.02, 0.10),
      advance = c(0.25, 0.55),
      background_stay = c(0.85, 0.95)
    ),
    warp = warp, screen = screen
  )
  baseline <- gaze_baseline_spec(
    screen = screen,
    density_sigmas = c(0.75, 1.5, 3),
    density_grid = if (smoke) 12L else 24L,
    warp = warp,
    lambda_grid = c(0.01, 0.1, 1, 10),
    inner_folds = 2L,
    elastic_radii = c(consensus = 2, rigidity = 10, matching = 1),
    elastic_maxit = if (smoke) 10L else 25L,
    elastic_tolerance = 1e-4
  )
  list(screen = screen, warp = warp, transport = transport,
       replay = replay, baseline = baseline)
}

v2_profile_fit <- function(label, expression) {
  invisible(gc(reset = TRUE))
  elapsed <- system.time(value <- force(expression))[["elapsed"]]
  memory <- gc()
  list(
    value = value,
    resources = data.frame(
      method = label,
      elapsed_seconds = unname(elapsed),
      gc_used_mb = sum(memory[, 2L]),
      gc_high_water_mb = sum(memory[, 7L]),
      fit_size_mb = as.numeric(utils::object.size(value)) / 1024^2,
      stringsAsFactors = FALSE
    )
  )
}

v2_candidate_ece <- function(candidate_tables, bins = 5L) {
  candidates <- do.call(rbind, lapply(candidate_tables, function(tab) {
    data.frame(
      posterior = tab$posterior,
      is_true = as.numeric(tab$is_true)
    )
  }))
  breaks <- seq(0, 1, length.out = bins + 1L)
  group <- cut(candidates$posterior, breaks, include.lowest = TRUE,
               labels = FALSE)
  total <- nrow(candidates)
  sum(vapply(split(seq_len(total), group), function(index) {
    length(index) / total * abs(
      mean(candidates$posterior[index]) - mean(candidates$is_true[index])
    )
  }, numeric(1)))
}

v2_engine_long <- function(fit, method) {
  results <- fit$results
  data.frame(
    participant = results$participant,
    family = results$family,
    replicate = results$replicate,
    item = results$item,
    condition = results$condition,
    method = method,
    calibrated = TRUE,
    gaze_info_bits = results$gaze_info_bits,
    log_loss = results$log_loss,
    brier_score = results$brier_score,
    template_rank = results$template_rank,
    top1_credit = results$top1_credit,
    posterior_true = results$posterior_true,
    candidates = I(results$candidates),
    stringsAsFactors = FALSE
  )
}

v2_baseline_long <- function(fit) {
  results <- fit$results
  keep <- results$status == "scored"
  results <- results[keep, , drop = FALSE]
  data.frame(
    participant = results$participant,
    family = results$family,
    replicate = results$replicate,
    item = results$item,
    condition = results$condition,
    method = results$method,
    calibrated = results$calibrated,
    gaze_info_bits = results$gaze_info_bits,
    log_loss = results$log_loss,
    brier_score = results$brier_score,
    template_rank = results$template_rank,
    top1_credit = results$top1_credit,
    posterior_true = results$posterior_true,
    compatibility_bits = results$compatibility_bits,
    compatibility_log_loss = results$compatibility_log_loss,
    candidates = I(results$candidates),
    stringsAsFactors = FALSE
  )
}

v2_method_summary <- function(scored) {
  supervised <- scored[scored$calibrated, , drop = FALSE]
  groups <- split(
    seq_len(nrow(supervised)),
    interaction(supervised$method, supervised$family,
                supervised$condition, drop = TRUE)
  )
  do.call(rbind, lapply(groups, function(index) {
    part <- supervised[index, , drop = FALSE]
    data.frame(
      method = part$method[[1]],
      family = part$family[[1]],
      condition = part$condition[[1]],
      n = nrow(part),
      mean_info_bits = mean(part$gaze_info_bits),
      mean_log_loss = mean(part$log_loss),
      mean_brier = mean(part$brier_score),
      mean_rank = mean(part$template_rank),
      top1_credit = mean(part$top1_credit),
      ece_5 = v2_candidate_ece(part$candidates),
      stringsAsFactors = FALSE
    )
  }))
}

v2_frozen_summary <- function(scored) {
  frozen <- scored[!scored$calibrated, , drop = FALSE]
  if (nrow(frozen) == 0L) return(data.frame())
  groups <- split(
    seq_len(nrow(frozen)),
    interaction(frozen$method, frozen$family, frozen$condition, drop = TRUE)
  )
  do.call(rbind, lapply(groups, function(index) {
    part <- frozen[index, , drop = FALSE]
    data.frame(
      method = part$method[[1]], family = part$family[[1]],
      condition = part$condition[[1]], n = nrow(part),
      mean_compatibility_bits = mean(part$compatibility_bits),
      mean_compatibility_log_loss = mean(part$compatibility_log_loss),
      mean_rank = mean(part$template_rank),
      top1_credit = mean(part$top1_credit),
      stringsAsFactors = FALSE
    )
  }))
}

v2_operating_characteristics <- function(scored) {
  supervised <- scored[scored$calibrated, , drop = FALSE]
  groups <- split(seq_len(nrow(supervised)), supervised$method)
  do.call(rbind, lapply(groups, function(index) {
    part <- supervised[index, , drop = FALSE]
    null <- sort(part$gaze_info_bits[part$condition == "null"])
    signal <- part$gaze_info_bits[part$condition == "signal"]
    threshold <- null[[ceiling(0.95 * length(null))]]
    credit <- function(value) mean(value > threshold) + 0.5 * mean(value == threshold)
    data.frame(
      method = part$method[[1]], threshold_bits = threshold,
      null_error = credit(null), power = credit(signal),
      stringsAsFactors = FALSE
    )
  }))
}

v2_extract_warp_scales <- function(fit, method, true_scale) {
  values <- unlist(lapply(fit$folds, function(fold) {
    if (is.null(fold)) return(numeric())
    vapply(fold$warp$groups, `[[`, numeric(1), "scale")
  }))
  data.frame(
    method = method, fold = seq_along(values), estimated_scale = values,
    true_scale = true_scale, absolute_error = abs(values - true_scale),
    stringsAsFactors = FALSE
  )
}

v2_alignment_diagnostics <- function(transport, replay) {
  transport_rows <- transport$results$condition == "signal"
  replay_rows <- replay$results$condition == "signal"
  transport_alignment <- transport$results$alignment[transport_rows]
  replay_alignment <- replay$results$alignment[replay_rows]
  data.frame(
    method = c("transport_v2", "replay"),
    mean_replay_coverage = c(
      mean(vapply(transport_alignment, function(x) x$diagnostics$replay_coverage,
                  numeric(1))),
      mean(vapply(replay_alignment, function(x) x$diagnostics$replay_coverage,
                  numeric(1)))
    ),
    mean_background_coverage = c(
      NA_real_,
      mean(vapply(replay_alignment, function(x) x$diagnostics$background_coverage,
                  numeric(1)))
    ),
    mean_spatial_rmse = c(
      mean(vapply(transport_alignment, function(x) x$diagnostics$spatial_rmse,
                  numeric(1))),
      mean(vapply(replay_alignment, function(x) x$diagnostics$spatial_rmse,
                  numeric(1)), na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

v2_score_engine_candidates <- function(results, true_key) {
  eyesim:::score_gaze_engine_results(
    results, true_key = true_key, candidate_pool_id = "v2-validation"
  )
}

v2_invariance_court <- function(seed = 20260816L) {
  set.seed(seed)
  coords <- rbind(c(4, 4), c(8, 13), c(13, 6), c(19, 14))
  reference <- v2_court_path(coords, c(1, 2, 1, 2))
  distractor <- v2_court_path(coords[c(1, 3, 2, 4), ] + c(0.7, -0.4), c(1, 2, 1, 2))
  source <- reference
  split_source <- v2_court_path(
    rbind(coords[1, ], coords[2, ], coords[2, ], coords[3:4, ]),
    c(1, 0.8, 1.2, 1, 2)
  )
  dilated <- v2_court_path(coords, c(3, 6, 3, 6), c(0, 3, 9, 12))
  detector <- v2_court_path(
    rbind(coords[1, ], coords[2, ] + c(-0.01, 0.01),
          coords[2, ] + c(0.01, -0.01), coords[3:4, ]),
    c(1, 1, 1, 1, 2)
  )

  transport_spec <- gaze_transport_v2_spec(
    gaze_gaussian_mixture(c(0.75, 1.5, 3), unit = "deg"),
    chronology = gaze_order_neighbours(2, coalesce_distance = 0.05, unit = "deg"),
    coverage_grid = c(0.5, 1), coverage_penalty_grid = 0.75,
    entropy_schedule = 0.03, maxit = 50, projection_maxit = 150,
    tolerance = 2e-4, projection_tolerance = 1e-6, multistart = 1
  )
  training_ref <- tibble::tibble(item = c("true", "wrong"),
                                 fixgroup = list(reference, distractor))
  training_source <- tibble::tibble(item = c("true", "wrong"),
                                    fixgroup = list(reference, distractor))
  replay_spec <- gaze_replay_spec(
    grid_size = 48, max_skip = 2,
    transition_grid = list(
      background = 0.03, restart = 0.03, advance = 0.4,
      background_stay = 0.9
    )
  )
  replay_model <- fit_gaze_replay_model(
    training_ref, training_source, "item", spec = replay_spec
  )
  score_transport <- function(path, spec = transport_spec,
                              references = list(reference, distractor)) {
    result <- Map(function(candidate, key) {
      gaze_transport_v2_align(candidate, path, spec, candidate_key = key,
                              coverage_penalty = 0.75)
    }, references, c("true", "wrong"))
    v2_score_engine_candidates(result, "true")$gaze_info_bits
  }
  score_replay <- function(path, model = replay_model,
                           references = list(reference, distractor)) {
    result <- Map(function(candidate, key) {
      gaze_replay_align(candidate, path, model, candidate_key = key)
    }, references, c("true", "wrong"))
    v2_score_engine_candidates(result, "true")$gaze_info_bits
  }
  base_transport <- score_transport(source)
  base_replay <- score_replay(source)

  px_factor <- 50
  px_path <- function(path) v2_court_path(
    cbind(path$x, path$y) * px_factor, path$duration, path$onset
  )
  px_transport <- gaze_transport_v2_spec(
    gaze_gaussian_mixture(c(0.75, 1.5, 3) * px_factor, unit = "px"),
    chronology = gaze_order_neighbours(
      2, coalesce_distance = 0.05 * px_factor, unit = "px"
    ),
    coverage_grid = c(0.5, 1), coverage_penalty_grid = 0.75,
    entropy_schedule = 0.03, maxit = 50, projection_maxit = 150,
    tolerance = 2e-4, projection_tolerance = 1e-6, multistart = 1
  )
  px_reference <- lapply(list(reference, distractor), px_path)
  px_transport_score <- score_transport(
    px_path(source), px_transport, px_reference
  )
  px_replay_spec <- replay_spec
  px_replay_spec$scale_floor <- replay_spec$scale_floor * px_factor
  px_model <- fit_gaze_replay_model(
    tibble::tibble(item = c("true", "wrong"), fixgroup = px_reference),
    tibble::tibble(item = c("true", "wrong"), fixgroup = px_reference),
    "item", spec = px_replay_spec
  )
  px_replay_score <- score_replay(px_path(source), px_model, px_reference)

  data.frame(
    check = c(
      "transport_split_merge", "replay_split_merge",
      "transport_time_dilation", "replay_time_dilation",
      "transport_detector_variant", "replay_detector_variant",
      "transport_unit_conversion", "replay_unit_conversion"
    ),
    absolute_difference = c(
      abs(score_transport(split_source) - base_transport),
      abs(score_replay(split_source) - base_replay),
      abs(score_transport(dilated) - base_transport),
      abs(score_replay(dilated) - base_replay),
      abs(score_transport(detector) - base_transport),
      abs(score_replay(detector) - base_replay),
      abs(px_transport_score - base_transport),
      abs(px_replay_score - base_replay)
    ),
    tolerance = c(1e-8, 1e-8, 1e-8, 1e-8, 0.15, 0.15, 1e-6, 1e-6),
    stringsAsFactors = FALSE
  )
}

v2_fractional_design <- function() {
  base <- as.matrix(expand.grid(A = 0:1, B = 0:1, C = 0:1, D = 0:1))
  xor <- function(...) Reduce(function(x, y) (x + y) %% 2L, list(...))
  design <- cbind(
    contraction = base[, 1], translation = base[, 2], noise = base[, 3],
    time_dilation = base[, 4],
    local_swaps = xor(base[, 1], base[, 2]),
    block_reorder = xor(base[, 1], base[, 3]),
    reversal = xor(base[, 1], base[, 4]),
    deletion = xor(base[, 2], base[, 3]),
    insertion = xor(base[, 2], base[, 4]),
    central_bias = xor(base[, 3], base[, 4]),
    fixation_count = xor(base[, 1], base[, 2], base[, 3]),
    duration_concentration = xor(base[, 1], base[, 2], base[, 4]),
    candidate_difficulty = xor(base[, 1], base[, 3], base[, 4])
  )
  as.data.frame(design)
}

v2_apply_perturbations <- function(path, row, seed) {
  set.seed(seed)
  coords <- cbind(path$x, path$y)
  duration <- path$duration
  if (row[["contraction"]] == 1) coords <- 0.72 * coords + 0.28 * c(12, 9)
  if (row[["translation"]] == 1) coords <- sweep(coords, 2, c(1.5, -1), "+")
  if (row[["noise"]] == 1) {
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 0.7), ncol = 2)
  }
  if (row[["local_swaps"]] == 1 && nrow(coords) >= 4) {
    coords[c(3, 4), ] <- coords[c(4, 3), ]
    duration[c(3, 4)] <- duration[c(4, 3)]
  }
  if (row[["block_reorder"]] == 1) {
    half <- floor(nrow(coords) / 2)
    order <- c(seq.int(half + 1L, nrow(coords)), seq_len(half))
    coords <- coords[order, , drop = FALSE]
    duration <- duration[order]
  }
  if (row[["reversal"]] == 1) {
    coords <- coords[nrow(coords):1L, , drop = FALSE]
    duration <- rev(duration)
  }
  if (row[["deletion"]] == 1 && nrow(coords) > 3) {
    keep <- seq_len(nrow(coords)) != ceiling(nrow(coords) / 2)
    coords <- coords[keep, , drop = FALSE]
    duration <- duration[keep]
  }
  if (row[["insertion"]] == 1) {
    coords <- rbind(coords, c(12, 9) + stats::rnorm(2, sd = 2.5))
    duration <- c(duration, stats::median(duration))
  }
  if (row[["central_bias"]] == 1) coords <- 0.65 * coords + 0.35 * c(12, 9)
  if (row[["fixation_count"]] == 1) {
    index <- which.max(duration)
    coords <- rbind(coords[seq_len(index), , drop = FALSE],
                    coords[index, , drop = FALSE],
                    if (index < nrow(coords)) coords[(index + 1L):nrow(coords), , drop = FALSE])
    duration <- c(duration[seq_len(index)], duration[index] / 2,
                  if (index < length(duration)) duration[(index + 1L):length(duration)])
    duration[index] <- duration[index] / 2
  }
  if (row[["duration_concentration"]] == 1) {
    duration[] <- 0.4
    duration[[1]] <- 4
  }
  dilation <- if (row[["time_dilation"]] == 1) 3 else 1
  v2_court_path(coords, duration * dilation)
}

v2_perturbation_court <- function(seed = 20260816L, smoke = FALSE) {
  templates <- v2_court_templates(4, 6, "replay_favouring", 1, seed)
  training <- tibble::tibble(item = paste0("i", 1:4), fixgroup = templates)
  replay_spec <- gaze_replay_spec(
    grid_size = if (smoke) 24 else 48, max_skip = 2,
    transition_grid = list(
      background = c(0.03, 0.1), restart = c(0.02, 0.1),
      advance = c(0.25, 0.55), background_stay = 0.9
    )
  )
  replay_model <- fit_gaze_replay_model(
    training, training, "item", spec = replay_spec
  )
  transport_spec <- gaze_transport_v2_spec(
    gaze_gaussian_mixture(c(0.75, 1.5, 3), unit = "deg"),
    chronology = gaze_order_neighbours(2, coalesce_distance = 0.05, unit = "deg"),
    coverage_grid = if (smoke) c(0.5, 1) else c(0.5, 0.75, 1),
    coverage_penalty_grid = 0.75,
    entropy_schedule = if (smoke) 0.03 else c(0.05, 0.015),
    maxit = if (smoke) 25 else 50, projection_maxit = 150,
    tolerance = if (smoke) 5e-4 else 2e-4,
    projection_tolerance = 1e-6, multistart = 1
  )
  density_spec <- gaze_baseline_spec(
    gaze_screen(24, 18, "deg"), c(0.75, 1.5, 3), density_grid = 20,
    methods = "density", inner_folds = 2
  )
  reference <- templates[[1]]
  easy <- templates[[2]]
  hard_coords <- 0.65 * cbind(easy$x, easy$y) +
    0.35 * cbind(reference$x, reference$y)
  hard <- v2_court_path(hard_coords, easy$duration)
  design <- v2_fractional_design()
  anchors <- as.data.frame(matrix(0L, nrow = ncol(design) + 1L,
                                  ncol = ncol(design)))
  names(anchors) <- names(design)
  for (i in seq_len(ncol(design))) anchors[i + 1L, i] <- 1L
  design <- rbind(anchors, design)
  design$run <- seq_len(nrow(design))
  design$type <- c("baseline", rep("one_factor", ncol(anchors)),
                   rep("fractional_factorial", nrow(design) - nrow(anchors)))

  rows <- lapply(seq_len(nrow(design)), function(i) {
    row <- design[i, , drop = FALSE]
    source <- v2_apply_perturbations(reference, row, seed + i * 101L)
    distractor <- if (row$candidate_difficulty[[1]] == 1) hard else easy
    transport_results <- list(
      gaze_transport_v2_align(reference, source, transport_spec,
                              candidate_key = "true", coverage_penalty = 0.75),
      gaze_transport_v2_align(distractor, source, transport_spec,
                              candidate_key = "wrong", coverage_penalty = 0.75)
    )
    replay_results <- list(
      gaze_replay_align(reference, source, replay_model, candidate_key = "true"),
      gaze_replay_align(distractor, source, replay_model, candidate_key = "wrong")
    )
    transport_evidence <- v2_score_engine_candidates(transport_results, "true")
    replay_evidence <- v2_score_engine_candidates(replay_results, "true")
    density <- rbind(
      eyesim:::baseline_density_features(reference, source, density_spec),
      eyesim:::baseline_density_features(distractor, source, density_spec)
    )
    density_evidence <- eyesim:::score_gaze_candidates(
      rowMeans(density), 1L, c("true", "wrong"),
      candidate_pool_id = "v2-perturbation"
    )
    data.frame(
      row[, setdiff(names(row), "type"), drop = FALSE],
      type = row$type,
      transport_info_bits = transport_evidence$gaze_info_bits,
      replay_info_bits = replay_evidence$gaze_info_bits,
      density_compatibility_bits = density_evidence$gaze_info_bits,
      transport_converged = all(vapply(transport_results, function(x) {
        x$convergence$converged
      }, logical(1))),
      replay_converged = all(vapply(replay_results, function(x) {
        x$convergence$converged
      }, logical(1))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

v2_stability_court <- function(seed = 20260816L) {
  templates <- v2_court_templates(3, 6, "replay_favouring", 1, seed)
  reference <- templates[[1]]
  source <- v2_court_signal(reference, "replay_favouring", 1, c(0, 0), seed + 1L)
  transport_spec <- gaze_transport_v2_spec(
    gaze_gaussian_mixture(c(0.75, 1.5, 3), unit = "deg"),
    chronology = gaze_order_neighbours(2, unit = "deg"),
    coverage_grid = c(0.5, 0.75, 1), coverage_penalty_grid = 0.75,
    entropy_schedule = c(0.05, 0.015), maxit = 60,
    projection_maxit = 250, tolerance = 2e-4,
    projection_tolerance = 1e-6, multistart = 2
  )
  transport <- lapply(seq_along(templates), function(i) {
    gaze_transport_v2_align(
      templates[[i]], source, transport_spec,
      candidate_key = paste0("i", i), coverage_penalty = 0.75
    )
  })
  forward <- v2_score_engine_candidates(transport, "i1")
  reverse <- v2_score_engine_candidates(rev(transport), "i1")
  forward_score <- forward$candidates$log_score[
    order(forward$candidates$candidate_key)
  ]
  reverse_score <- reverse$candidates$log_score[
    order(reverse$candidates$candidate_key)
  ]

  training <- tibble::tibble(item = paste0("i", 1:3), fixgroup = templates)
  normalized <- numeric(3)
  replay_scores <- list()
  resolution <- c(32L, 48L, 64L)
  for (g in seq_along(resolution)) {
    spec <- gaze_replay_spec(
      grid_size = resolution[[g]], max_skip = 2,
      scale_floor = 0.15,
      transition_grid = list(
        background = 0.03, restart = 0.03, advance = 0.4,
        background_stay = 0.9
      )
    )
    model <- fit_gaze_replay_model(training, training, "item", spec = spec)
    candidate <- lapply(seq_along(templates), function(i) {
      gaze_replay_align(
        templates[[i]], source, model, candidate_key = paste0("i", i)
      )
    })
    replay_scores[[g]] <- candidate
    candidate_score <- vapply(candidate, `[[`, numeric(1), "log_score")
    normalized[[g]] <- (
      candidate_score[[1]] - eyesim:::log_mean_exp(candidate_score[-1])
    ) / resolution[[g]]
  }
  replay_forward <- v2_score_engine_candidates(replay_scores[[2]], "i1")
  replay_reverse <- v2_score_engine_candidates(rev(replay_scores[[2]]), "i1")
  replay_forward_score <- replay_forward$candidates$log_score[
    order(replay_forward$candidates$candidate_key)
  ]
  replay_reverse_score <- replay_reverse$candidates$log_score[
    order(replay_reverse$candidates$candidate_key)
  ]
  true_ranks <- vapply(replay_scores, function(candidate) {
    v2_score_engine_candidates(candidate, "i1")$template_rank
  }, numeric(1))
  posterior_error <- max(vapply(replay_scores[[2]], function(result) {
    max(abs(rowSums(result$alignment$posterior) - 1))
  }, numeric(1)))
  true_transport <- transport[[1]]$diagnostics

  data.frame(
    check = c(
      "transport_start_scientific_spread",
      "transport_barycentric_start_rmse",
      "transport_candidate_order",
      "transport_stationarity_finite",
      "replay_candidate_order",
      "replay_posterior_normalization",
      "replay_resolution_per_bin",
      "replay_resolution_rank"
    ),
    value = c(
      true_transport$start_scientific_spread,
      true_transport$alignment_stability$value,
      max(abs(forward_score - reverse_score)),
      as.numeric(is.finite(true_transport$stationarity)),
      max(abs(replay_forward_score - replay_reverse_score)),
      posterior_error,
      max(normalized) - min(normalized),
      as.numeric(length(unique(true_ranks)) == 1L)
    ),
    comparator = c("le", "le", "le", "eq", "le", "le", "le", "eq"),
    tolerance = c(0.02, 0.25, 1e-8, 1, 1e-10, 1e-10, 0.15, 1),
    stringsAsFactors = FALSE
  )
}

v2_gate_verdict <- function(invariance, stability, summary, operating,
                            recovery, scored, alignment) {
  invariant_pass <- all(invariance$absolute_difference <= invariance$tolerance)
  stability_pass <- all(ifelse(
    stability$comparator == "le",
    stability$value <= stability$tolerance,
    stability$value == stability$tolerance
  ))
  engine_summary <- summary[
    summary$method %in% c("transport_v2", "replay") &
      summary$condition == "signal", , drop = FALSE
  ]
  positive <- aggregate(
    mean_info_bits ~ method, engine_summary,
    function(x) sum(x > 0)
  )
  names(positive)[[2]] <- "positive_families"
  top1 <- aggregate(top1_credit ~ method, engine_summary, mean)
  predictive_pass <- all(positive$positive_families >= 2L) &&
    all(top1$top1_credit > 1 / 6)
  engine_operating <- operating[operating$method %in% c("transport_v2", "replay"), ]
  null_pass <- all(engine_operating$null_error <= 0.075)
  recovery_pass <- all(recovery$absolute_error <= 0.08)
  convergence <- c(
    mean(scored$transport_v2$results$all_converged),
    mean(scored$replay$results$all_converged)
  )
  convergence_pass <- all(convergence >= 0.99)

  order <- summary[
    summary$family == "neutral_misspecified" &
      summary$condition == "signal" &
      summary$method %in% c(
        "transport_v2", "replay", "density_ridge_registered"
      ), , drop = FALSE
  ]
  density_loss <- order$mean_log_loss[
    order$method == "density_ridge_registered"
  ]
  engine_loss <- order$mean_log_loss[
    order$method %in% c("transport_v2", "replay")
  ]
  order_pass <- length(density_loss) == 1L && any(engine_loss < density_loss)

  baseline <- scored$baseline$results
  null_baseline <- baseline[baseline$condition == "null", , drop = FALSE]
  raw <- mean(null_baseline$top1_credit[
    null_baseline$method == "density_raw"
  ])
  registered <- mean(null_baseline$top1_credit[
    null_baseline$method == "density_registered"
  ])
  null_registration_pass <- is.finite(raw) && is.finite(registered) &&
    registered - raw <= 0.05

  directional <- summary[
    summary$family %in% c("replay_favouring", "neutral_misspecified") &
      summary$condition == "signal" &
      summary$method %in% c("transport_v2", "replay"), , drop = FALSE
  ]
  pooled_loss <- aggregate(mean_log_loss ~ method, directional, mean)
  replay_loss <- pooled_loss$mean_log_loss[pooled_loss$method == "replay"]
  transport_loss <- pooled_loss$mean_log_loss[pooled_loss$method == "transport_v2"]
  perturbation <- scored$perturbation
  local_row <- perturbation[
    perturbation$type == "one_factor" & perturbation$local_swaps == 1, , drop = FALSE
  ]
  local_default <- nrow(local_row) == 1L &&
    local_row$replay_info_bits > local_row$transport_info_bits
  default_replay <- isTRUE(replay_loss <= transport_loss + 0.05) && local_default

  gates <- data.frame(
    gate = c(
      "representation_and_stability", "null_error", "predictive_signal",
      "fixed_density_order", "contraction_recovery",
      "null_registration", "engine_convergence"
    ),
    passed = c(
      invariant_pass && stability_pass, null_pass, predictive_pass,
      order_pass, recovery_pass, null_registration_pass, convergence_pass
    ),
    stringsAsFactors = FALSE
  )
  list(
    gates = gates,
    advance = all(gates$passed),
    provisional_default = if (all(gates$passed) && default_replay) {
      "replay"
    } else if (all(gates$passed)) {
      "none_directional"
    } else {
      "gate_failed"
    },
    pooled_directional_log_loss = pooled_loss,
    convergence_rate = data.frame(
      method = c("transport_v2", "replay"), rate = convergence
    ),
    alignment = alignment
  )
}

v2_write_results <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  csv <- list(
    "supervised-summary.csv" = result$supervised_summary,
    "frozen-summary.csv" = result$frozen_summary,
    "operating-characteristics.csv" = result$operating,
    "parameter-recovery.csv" = result$parameter_recovery,
    "alignment-diagnostics.csv" = result$alignment_diagnostics,
    "invariance.csv" = result$invariance,
    "perturbations.csv" = result$perturbation,
    "stability.csv" = result$stability,
    "resources.csv" = result$resources,
    "gate-verdict.csv" = result$verdict$gates,
    "convergence.csv" = result$verdict$convergence_rate,
    "directional-log-loss.csv" = result$verdict$pooled_directional_log_loss
  )
  for (name in names(csv)) {
    utils::write.csv(csv[[name]], file.path(output_dir, name), row.names = FALSE)
  }
  config <- data.frame(
    protocol_version = 1L, final_seed = result$config$seed,
    smoke = result$config$smoke,
    n_participants = result$config$n_participants,
    n_items = result$config$n_items,
    n_replications = result$config$n_replications,
    advance = result$verdict$advance,
    provisional_default = result$verdict$provisional_default
  )
  utils::write.csv(config, file.path(output_dir, "configuration.csv"), row.names = FALSE)
  saveRDS(
    result[c(
      "supervised_summary", "frozen_summary", "operating",
      "parameter_recovery", "alignment_diagnostics", "invariance",
      "perturbation", "stability", "verdict", "config"
    )],
    file.path(output_dir, "scientific-results.rds"), version = 3
  )
  utils::capture.output(
    utils::sessionInfo(), file = file.path(output_dir, "session-info.txt")
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[basename(files) != "manifest-md5.csv"]
  manifest <- data.frame(
    file = basename(files), md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    manifest, file.path(output_dir, "manifest-md5.csv"), row.names = FALSE
  )
  invisible(result)
}

run_gaze_weave_v2_court <- function(
    output_dir = NULL, seed = 20260816L, smoke = FALSE,
    n_participants = if (smoke) 1L else 1L,
    n_items = if (smoke) 4L else 6L,
    n_replications = 2L) {
  if (!smoke && seed != 20260816L) {
    stop("The final court is frozen to seed 20260816; use smoke = TRUE for debugging.")
  }
  data <- v2_court_data(
    n_participants, n_items, n_replications,
    n_fixations = if (smoke) 5L else 6L, seed = seed
  )
  specs <- v2_court_specs(smoke)
  match_on <- c("participant", "family", "replicate", "item")
  contrast_on <- c("participant", "family")
  outer_split_on <- c("participant", "family", "replicate")
  fit_filter <- function(tab) tab$condition == "signal"
  detected_cores <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 2L
  transport_workers <- if (smoke) 2L else min(4L, detected_cores)

  transport_profile <- v2_profile_fit("transport_v2", {
    gaze_transport_v2_cv(
      data$reference, data$source, match_on, contrast_on,
      spec = specs$transport, split_on = outer_split_on, n_folds = 2L,
      seed = seed, fit_source_filter = fit_filter,
      workers = transport_workers
    )
  })
  replay_profile <- v2_profile_fit("replay", {
    gaze_replay_cv(
      data$reference, data$source, match_on, contrast_on,
      spec = specs$replay, split_on = outer_split_on, n_folds = 2L,
      seed = seed, fit_source_filter = fit_filter
    )
  })
  baseline_profile <- v2_profile_fit("baselines", {
    gaze_baseline_cv(
      data$reference, data$source, match_on, contrast_on,
      spec = specs$baseline, split_on = outer_split_on, n_folds = 2L,
      seed = seed, fit_source_filter = fit_filter,
      inner_split_on = match_on
    )
  })
  transport <- transport_profile$value
  replay <- replay_profile$value
  baseline <- baseline_profile$value
  scored <- dplyr::bind_rows(
    v2_engine_long(transport, "transport_v2"),
    v2_engine_long(replay, "replay"),
    v2_baseline_long(baseline)
  )
  supervised_summary <- v2_method_summary(scored)
  frozen_summary <- v2_frozen_summary(scored)
  operating <- v2_operating_characteristics(scored)
  true_scale <- mean(data$truth$scale)
  parameter_recovery <- rbind(
    v2_extract_warp_scales(transport, "transport_v2", true_scale),
    v2_extract_warp_scales(replay, "replay", true_scale),
    v2_extract_warp_scales(baseline, "baselines", true_scale)
  )
  alignment_diagnostics <- v2_alignment_diagnostics(transport, replay)
  invariance <- v2_invariance_court(seed)
  perturbation <- v2_perturbation_court(seed, smoke)
  stability <- v2_stability_court(seed)
  scored_objects <- list(
    transport_v2 = transport, replay = replay, baseline = baseline,
    perturbation = perturbation
  )
  verdict <- v2_gate_verdict(
    invariance, stability, supervised_summary, operating,
    parameter_recovery, scored_objects, alignment_diagnostics
  )
  result <- list(
    supervised_summary = supervised_summary,
    frozen_summary = frozen_summary,
    operating = operating,
    parameter_recovery = parameter_recovery,
    alignment_diagnostics = alignment_diagnostics,
    invariance = invariance,
    perturbation = perturbation,
    stability = stability,
    resources = rbind(
      transport_profile$resources, replay_profile$resources,
      baseline_profile$resources
    ),
    verdict = verdict,
    fits = scored_objects,
    config = list(
      protocol_version = 1L, seed = seed, smoke = smoke,
      n_participants = n_participants, n_items = n_items,
      n_replications = n_replications
    )
  )
  if (!is.null(output_dir)) v2_write_results(result, output_dir)
  result
}
