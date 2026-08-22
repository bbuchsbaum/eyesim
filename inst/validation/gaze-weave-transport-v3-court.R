# Frozen public simulation court for GazeWeave Transport v3.
#
# Run from a source checkout after installing eyesim:
#   library(eyesim)
#   source("inst/validation/gaze-weave-transport-v3-court.R")
#   run_gaze_weave_transport_v3_court(
#     "inst/validation/gaze-weave-transport-v3-court-results"
#   )

v3_court_path <- function(coords, duration = NULL) {
  if (is.null(duration)) duration <- rep(1, nrow(coords))
  fixation_group(
    x = coords[, 1], y = coords[, 2], duration = duration,
    onset = c(0, head(cumsum(duration), -1L))
  )
}

v3_court_rotate <- function(coords, angle, center) {
  rotation <- matrix(
    c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2,
    byrow = TRUE
  )
  sweep(sweep(coords, 2, center, FUN = "-") %*% t(rotation),
        2, center, FUN = "+")
}

v3_court_data <- function(n_items = 8L, seed = 20260822L) {
  set.seed(seed)
  center <- c(600, 450)
  motif <- rbind(
    c(-180, -70), c(-95, 155), c(35, -135),
    c(175, 105), c(135, -185), c(-145, 55)
  )
  duration <- c(1, 2, 0.8, 1.6, 1.1, 0.7)
  reference <- lapply(seq_len(n_items), function(item) {
    coords <- v3_court_rotate(
      sweep(motif, 2, center, FUN = "+"),
      2 * pi * (item - 1L) / n_items,
      center
    )
    v3_court_path(coords, duration)
  })
  scale <- 0.8
  translation <- c(52, -38)
  signal <- lapply(seq_len(n_items), function(item) {
    coords <- cbind(reference[[item]]$x, reference[[item]]$y)
    raw <- sweep(coords, 2, translation, FUN = "-") / scale
    v3_court_path(raw, duration * c(1, 1.2, 0.9, 1.1, 1, 0.8))
  })
  null_registered <- matrix(center, nrow = nrow(motif), ncol = 2,
                            byrow = TRUE)
  null_raw <- sweep(null_registered, 2, translation, FUN = "-") / scale
  null <- lapply(seq_len(n_items), function(item) {
    v3_court_path(null_raw, duration)
  })
  list(
    reference = tibble::tibble(
      participant = "public", item = seq_len(n_items),
      prior_weight = 1, fixgroup = reference
    ),
    source = tibble::tibble(
      participant = "public",
      item = rep(seq_len(n_items), each = 2L),
      condition = rep(c("signal", "null"), n_items),
      fixgroup = unlist(
        Map(function(signal_path, null_path) {
          list(signal_path, null_path)
        }, signal, null),
        recursive = FALSE
      )
    ),
    truth = list(scale = scale, translation = translation, center = center)
  )
}

v3_court_specs <- function(smoke = FALSE) {
  screen <- gaze_screen(1200, 900, unit = "px")
  warp <- gaze_warp_contraction(
    center = "screen", translation = TRUE, fit_by = "participant"
  )
  transport_v3 <- gaze_transport_spec(
    spatial = gaze_gaussian_mixture(c(40, 80, 160), unit = "px"),
    chronology = gaze_order_neighbours(
      2, coalesce_distance = 2, unit = "px"
    ),
    coverage_nodes = if (smoke) 2L else 12L,
    temporal_weight = 2,
    entropy_schedule = if (smoke) 0.03 else c(0.05, 0.015),
    warp = warp, screen = screen,
    maxit = 1000L,
    tolerance = 5e-5,
    projection_maxit = if (smoke) 300L else 1000L,
    projection_tolerance = 1e-7,
    backend = "optimized", reliability = "effective_fixations",
    calibration_folds = 2L, calibration_seed = 20260822L
  )
  transport_v2 <- gaze_transport_v2_spec(
    spatial = gaze_gaussian_mixture(c(40, 80, 160), unit = "px"),
    chronology = gaze_order_neighbours(2, coalesce_distance = 2, unit = "px"),
    coverage_grid = c(0.5, 0.75, 1),
    coverage_penalty_grid = c(0.25, 0.75, 1.5),
    temporal_weight = 1,
    entropy_schedule = if (smoke) 0.03 else c(0.05, 0.015),
    warp = warp, screen = screen,
    maxit = if (smoke) 60L else 100L,
    projection_maxit = if (smoke) 250L else 400L,
    tolerance = 5e-4, projection_tolerance = 1e-6,
    multistart = 1L
  )
  replay <- gaze_replay_spec(
    grid_size = if (smoke) 24L else 48L,
    max_skip = 2L, student_df = 4, scale_floor = 8,
    transition_grid = list(
      background = c(0.03, 0.10), restart = c(0.02, 0.10),
      advance = c(0.25, 0.55), background_stay = c(0.85, 0.95)
    ),
    warp = warp, screen = screen,
    reliability = "effective_fixations"
  )
  baseline <- gaze_baseline_spec(
    screen = screen,
    density_sigmas = c(80, 160),
    density_grid = if (smoke) 8L else 24L,
    methods = c("multimatch", "density"),
    warp = warp,
    lambda_grid = c(0.01, 0.1, 1, 10), inner_folds = 2L
  )
  list(
    screen = screen, warp = warp, transport_v3 = transport_v3,
    transport_v2 = transport_v2, replay = replay, baseline = baseline
  )
}

v3_court_profile <- function(label, expression) {
  invisible(gc(reset = TRUE))
  elapsed <- system.time(value <- force(expression))[["elapsed"]]
  memory <- gc()
  list(
    value = value,
    resources = data.frame(
      method = label, elapsed_seconds = unname(elapsed),
      gc_used_mb = sum(memory[, 2L]),
      gc_high_water_mb = sum(memory[, 7L]),
      fit_size_mb = as.numeric(utils::object.size(value)) / 1024^2,
      stringsAsFactors = FALSE
    )
  )
}

v3_court_engine_row <- function(evidence, item, condition, method,
                                converged, candidates) {
  data.frame(
    item = item, condition = condition, method = method,
    gaze_info_bits = evidence$gaze_info_bits,
    log_loss = evidence$log_loss,
    brier_score = evidence$brier_score,
    template_rank = evidence$template_rank,
    top1_credit = evidence$top1_credit,
    posterior_true = evidence$posterior_true,
    prior_true = evidence$prior_true,
    candidate_count = evidence$candidate_count,
    converged = converged,
    candidates = I(list(candidates)),
    stringsAsFactors = FALSE
  )
}

v3_court_score_common_fold <- function(
    data, specs, fold, fold_id, seed) {
  source <- data$source
  reference <- data$reference
  eval_items <- which(fold_id == fold)
  train_items <- which(fold_id != fold)
  train_source <- source[
    source$item %in% train_items & source$condition == "signal", ,
    drop = FALSE
  ]
  ref_train <- reference[reference$item %in% train_items, , drop = FALSE]
  eval_source <- source[source$item %in% eval_items, , drop = FALSE]

  v3_warp <- eyesim:::fit_transport_v3_warp(
    ref_train, train_source, "item", "fixgroup", "fixgroup", NULL,
    specs$transport_v3
  )
  train_source$..gaze_row_id <- match(train_source$item, source$item)
  v3_inner <- eyesim:::fit_transport_v3_inner_calibration(
    reference, train_source, "item", NULL, "fixgroup", "fixgroup",
    NULL, "prior_weight", specs$transport_v3
  )
  v2_model <- eyesim:::fit_gaze_transport_v2_model(
    ref_train, train_source, "item", NULL, "fixgroup", "fixgroup",
    specs$transport_v2, workers = 1L
  )
  replay_model <- fit_gaze_replay_model(
    ref_train, train_source, "item", contrast_on = NULL,
    spec = specs$replay
  )
  baseline_model <- eyesim:::fit_gaze_baseline_model(
    ref_train, train_source, "item", NULL, "item",
    "fixgroup", "fixgroup", specs$baseline, seed + fold
  )
  baseline_sets <- eyesim:::baseline_build_feature_sets(
    reference, eval_source, "item", NULL, "fixgroup", "fixgroup",
    specs$baseline, baseline_model$warp, baseline_model$availability
  )
  baseline_scored <- eyesim:::score_gaze_baseline_sets(
    baseline_sets, baseline_model
  )

  rows <- list()
  row_index <- 0L
  keep_baselines <- c(
    "density_ridge_registered",
    "density_sigma_80_raw_calibrated",
    "density_sigma_160_raw_calibrated",
    "multimatch_ridge_registered",
    "multimatch_mm_vector_raw_calibrated",
    "multimatch_mm_direction_raw_calibrated",
    "multimatch_mm_length_raw_calibrated",
    "multimatch_mm_position_raw_calibrated",
    "multimatch_mm_duration_raw_calibrated",
    "multimatch_mm_position_emd_raw_calibrated"
  )
  for (index in seq_len(nrow(eval_source))) {
    source_row <- eval_source[index, , drop = FALSE]
    item <- source_row$item[[1L]]
    condition <- source_row$condition[[1L]]
    v3 <- eyesim:::score_transport_v3_cv_row(
      source_row, reference, "item", NULL, "fixgroup", "fixgroup",
      NULL, "prior_weight", specs$transport_v3, v3_warp,
      temperature = v3_inner$calibration$temperature,
      kappa = v3_inner$calibration$kappa
    )
    row_index <- row_index + 1L
    rows[[row_index]] <- v3_court_engine_row(
      v3$evidence, item, condition, "transport_v3", v3$all_converged,
      v3$evidence$candidates
    )
    v2 <- eyesim:::score_gaze_transport_v2_row(
      source_row, reference, "item", NULL, "fixgroup", "fixgroup",
      v2_model, workers = 1L
    )
    row_index <- row_index + 1L
    rows[[row_index]] <- v3_court_engine_row(
      v2$evidence, item, condition, "transport_v2", v2$all_converged,
      v2$evidence$candidates
    )
    replay <- eyesim:::score_gaze_replay_row(
      source_row, reference, "item", NULL, "fixgroup", "fixgroup",
      replay_model
    )
    row_index <- row_index + 1L
    rows[[row_index]] <- v3_court_engine_row(
      replay$evidence, item, condition, "replay", replay$all_converged,
      replay$evidence$candidates
    )
    for (method in intersect(keep_baselines, names(baseline_scored[[index]]))) {
      result <- baseline_scored[[index]][[method]]
      if (!identical(result$status, "scored") || !result$calibrated) next
      evidence <- list(
        gaze_info_bits = result$gaze_info_bits,
        log_loss = result$log_loss,
        brier_score = result$brier_score,
        template_rank = result$template_rank,
        top1_credit = result$top1_credit,
        posterior_true = result$posterior_true,
        prior_true = 1 / result$candidate_count,
        candidate_count = result$candidate_count
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- v3_court_engine_row(
        evidence, item, condition, method, TRUE, result$candidates
      )
    }
  }
  list(
    scored = dplyr::bind_rows(rows),
    receipt = list(
      fold = fold, train_items = train_items, eval_items = eval_items,
      overlap_item_n = length(intersect(train_items, eval_items)),
      candidate_keys = as.character(reference$item),
      candidate_prior = rep(1 / nrow(reference), nrow(reference)),
      v3_inner = v3_inner$receipts,
      warp = list(
        transport_v3 = v3_warp$info,
        transport_v2 = v2_model$warp$info,
        replay = replay_model$warp$info,
        baselines = baseline_model$warp$info
      ),
      baseline_inner = baseline_model$inner_splits
    ),
    models = list(
      transport_v3 = list(warp = v3_warp, calibration = v3_inner$calibration),
      transport_v2 = v2_model, replay = replay_model,
      baselines = baseline_model
    )
  )
}

v3_court_predictive <- function(data, specs, seed = 20260822L) {
  set.seed(seed)
  item_order <- sample(seq_len(nrow(data$reference)))
  fold_id <- integer(length(item_order))
  fold_id[item_order] <- rep(1:2, length.out = length(item_order))
  profiles <- lapply(1:2, function(fold) {
    v3_court_profile(
      paste0("shared_fold_", fold),
      v3_court_score_common_fold(data, specs, fold, fold_id, seed)
    )
  })
  list(
    scored = dplyr::bind_rows(lapply(profiles, function(x) x$value$scored)),
    receipts = lapply(profiles, function(x) x$value$receipt),
    models = lapply(profiles, function(x) x$value$models),
    resources = do.call(rbind, lapply(profiles, `[[`, "resources")),
    fold_id = fold_id
  )
}

v3_court_ece <- function(candidate_tables, bins = 10L) {
  evidence <- do.call(rbind, lapply(candidate_tables, function(table) {
    data.frame(posterior = table$posterior, is_true = table$is_true)
  }))
  group <- cut(
    evidence$posterior, seq(0, 1, length.out = bins + 1L),
    include.lowest = TRUE, labels = FALSE
  )
  sum(vapply(split(seq_len(nrow(evidence)), group), function(index) {
    length(index) / nrow(evidence) * abs(
      mean(evidence$posterior[index]) - mean(evidence$is_true[index])
    )
  }, numeric(1)))
}

v3_court_summaries <- function(scored) {
  groups <- split(
    seq_len(nrow(scored)),
    interaction(scored$method, scored$condition, drop = TRUE)
  )
  summary <- do.call(rbind, lapply(groups, function(index) {
    part <- scored[index, , drop = FALSE]
    data.frame(
      method = part$method[[1L]], condition = part$condition[[1L]],
      n = nrow(part), mean_info_bits = mean(part$gaze_info_bits),
      mean_log_loss = mean(part$log_loss),
      mean_brier = mean(part$brier_score),
      mean_rank = mean(part$template_rank),
      top1_credit = mean(part$top1_credit),
      candidate_count = unique(part$candidate_count),
      ece_10 = v3_court_ece(part$candidates),
      convergence_rate = mean(part$converged),
      stringsAsFactors = FALSE
    )
  }))
  operating <- do.call(rbind, lapply(
    split(seq_len(nrow(scored)), scored$method), function(index) {
      part <- scored[index, , drop = FALSE]
      null <- sort(part$gaze_info_bits[part$condition == "null"])
      signal <- part$gaze_info_bits[part$condition == "signal"]
      threshold <- null[[ceiling(0.95 * length(null))]]
      credit <- function(value) {
        mean(value > threshold) + 0.5 * mean(value == threshold)
      }
      data.frame(
        method = part$method[[1L]], threshold_bits = threshold,
        false_positive_rate = credit(null), power = credit(signal),
        stringsAsFactors = FALSE
      )
    }
  ))
  list(summary = summary, operating = operating)
}

v3_court_bootstrap <- function(scored, seed = 20260822L,
                               replicates = 500L) {
  set.seed(seed)
  items <- sort(unique(scored$item))
  draws <- replicate(
    replicates, sample(items, length(items), replace = TRUE),
    simplify = FALSE
  )
  groups <- unique(scored[c("method", "condition")])
  do.call(rbind, lapply(seq_len(nrow(groups)), function(index) {
    method <- groups$method[[index]]
    condition <- groups$condition[[index]]
    part <- scored[
      scored$method == method & scored$condition == condition, ,
      drop = FALSE
    ]
    estimate <- vapply(draws, function(draw) {
      mean(vapply(draw, function(item) {
        mean(part$gaze_info_bits[part$item == item])
      }, numeric(1)))
    }, numeric(1))
    data.frame(
      method = method, condition = condition,
      mean_info_bits = mean(part$gaze_info_bits),
      lower_95 = unname(stats::quantile(estimate, 0.025, type = 8)),
      upper_95 = unname(stats::quantile(estimate, 0.975, type = 8)),
      seed = seed, replicates = replicates,
      stringsAsFactors = FALSE
    )
  }))
}

v3_court_score_candidates <- function(references, source, spec,
                                      batch_size = length(references)) {
  results <- gaze_transport_align_batch(
    references, source, spec, batch_size = batch_size
  )
  eyesim:::score_gaze_engine_results(
    results, true_key = names(references)[[1L]],
    candidate_pool_id = "transport-v3-invariance"
  )
}

v3_court_invariance <- function(seed = 20260822L) {
  set.seed(seed)
  coords <- rbind(
    c(100, 120), c(280, 470), c(520, 180), c(760, 590), c(980, 240)
  )
  references <- list(
    true = v3_court_path(coords, c(1, 2, 1, 2, 1)),
    wrong_a = v3_court_path(coords[c(1, 3, 2, 5, 4), ] + c(25, -15)),
    wrong_b = v3_court_path(coords[c(5, 4, 3, 2, 1), ] + c(-30, 20))
  )
  source <- references$true
  split <- v3_court_path(
    rbind(coords[1:2, ], coords[2, ], coords[3:5, ]),
    c(1, 0.8, 1.2, 1, 2, 1)
  )
  dilated <- v3_court_path(
    coords, c(3, 6, 3, 6, 3)
  )
  make_spec <- function(nodes = 12L, scale = 1) {
    gaze_transport_spec(
      spatial = gaze_gaussian_mixture(c(40, 80, 160) * scale),
      chronology = gaze_order_neighbours(
        2, coalesce_distance = 2 * scale
      ),
      coverage_nodes = nodes, entropy_schedule = 0.03,
      maxit = 1000, tolerance = 5e-5,
      projection_maxit = 300, projection_tolerance = 1e-7,
      backend = "optimized", reliability = "none"
    )
  }
  spec <- make_spec()
  base <- v3_court_score_candidates(references, source, spec)
  reverse_order <- v3_court_score_candidates(
    rev(references), source, spec, batch_size = 1L
  )
  refined <- v3_court_score_candidates(references, source, make_spec(24L))
  factor <- 0.01
  scale_path <- function(path) {
    v3_court_path(cbind(path$x, path$y) * factor, path$duration)
  }
  unit <- v3_court_score_candidates(
    lapply(references, scale_path), scale_path(source),
    make_spec(12L, factor)
  )
  checks <- data.frame(
    check = c(
      "split_merge", "time_dilation", "unit_conversion",
      "candidate_order", "batch_size", "coverage_refinement_energy",
      "coverage_refinement_information"
    ),
    value = c(
      abs(v3_court_score_candidates(references, split, spec)$gaze_info_bits -
            base$gaze_info_bits),
      abs(v3_court_score_candidates(references, dilated, spec)$gaze_info_bits -
            base$gaze_info_bits),
      abs(unit$gaze_info_bits - base$gaze_info_bits),
      max(abs(
        base$candidates$log_score[order(base$candidates$candidate_key)] -
          reverse_order$candidates$log_score[
            order(reverse_order$candidates$candidate_key)
          ]
      )),
      max(abs(
        v3_court_score_candidates(references, source, spec, 1L)$candidates$log_score -
          v3_court_score_candidates(references, source, spec, 3L)$candidates$log_score
      )),
      max(abs(base$candidates$log_score - refined$candidates$log_score)),
      abs(base$gaze_info_bits - refined$gaze_info_bits)
    ),
    tolerance = c(1e-8, 1e-8, 1e-8, 1e-10, 1e-10, 1e-3, 0.01),
    stringsAsFactors = FALSE
  )
  checks$passed <- checks$value <= checks$tolerance
  checks
}

v3_court_perturb_path <- function(path, factor, center, seed) {
  set.seed(seed)
  coords <- cbind(path$x, path$y)
  duration <- path$duration
  if (factor == "contraction") coords <- 0.75 * coords + 0.25 * center
  if (factor == "translation") coords <- sweep(coords, 2, c(70, -50), "+")
  if (factor == "spatial_noise") {
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 35), ncol = 2)
  }
  if (factor == "duration_heterogeneity") {
    duration <- c(4, rep(0.4, length(duration) - 1L))
  }
  if (factor == "global_time_dilation") duration <- duration * 4
  if (factor == "local_swap") {
    coords[c(3, 4), ] <- coords[c(4, 3), ]
    duration[c(3, 4)] <- duration[c(4, 3)]
  }
  if (factor == "block_reorder") {
    order <- c(4:6, 1:3)
    coords <- coords[order, , drop = FALSE]
    duration <- duration[order]
  }
  if (factor == "reversal") {
    coords <- coords[nrow(coords):1L, , drop = FALSE]
    duration <- rev(duration)
  }
  if (factor == "deletion") {
    coords <- coords[-3L, , drop = FALSE]
    duration <- duration[-3L]
  }
  if (factor == "insertion") {
    coords <- rbind(coords, center + c(120, -90))
    duration <- c(duration, stats::median(duration))
  }
  if (factor == "central_bias") coords <- 0.6 * coords + 0.4 * center
  if (factor == "fixation_count") {
    coords <- rbind(coords[1:2, ], coords[2, ], coords[3:6, ])
    duration <- c(duration[1], duration[2] / 2, duration[2] / 2,
                  duration[3:6])
  }
  if (factor == "coherent_short_subsequence") {
    coords <- coords[1:3, , drop = FALSE]
    duration <- duration[1:3]
  }
  v3_court_path(coords, duration)
}

v3_court_perturbations <- function(data, seed = 20260822L) {
  references <- stats::setNames(data$reference$fixgroup, data$reference$item)
  reference <- references[[1L]]
  center <- data$truth$center
  factors <- c(
    "intact", "contraction", "translation", "spatial_noise",
    "duration_heterogeneity", "global_time_dilation", "local_swap",
    "block_reorder", "reversal", "deletion", "insertion",
    "central_bias", "fixation_count", "candidate_difficulty",
    "coherent_short_subsequence"
  )
  spec <- gaze_transport_spec(
    spatial = gaze_gaussian_mixture(c(40, 80, 160)),
    chronology = gaze_order_neighbours(2, coalesce_distance = 2),
    coverage_nodes = 12, entropy_schedule = 0.03,
    maxit = 1000, tolerance = 5e-5,
    projection_maxit = 300, projection_tolerance = 1e-7,
    backend = "optimized", reliability = "none"
  )
  baseline_spec <- gaze_baseline_spec(
    gaze_screen(1200, 900, "px"), c(80, 160), density_grid = 12,
    methods = c("multimatch", "density"), inner_folds = 2
  )
  availability <- gaze_baseline_availability(baseline_spec)
  rows <- lapply(seq_along(factors), function(index) {
    factor <- factors[[index]]
    source <- if (factor == "intact" || factor == "candidate_difficulty") {
      reference
    } else {
      v3_court_perturb_path(reference, factor, center, seed + index * 101L)
    }
    candidates <- references
    if (factor == "candidate_difficulty") {
      hard <- candidates[[2L]]
      hard$x <- 0.75 * reference$x + 0.25 * hard$x
      hard$y <- 0.75 * reference$y + 0.25 * hard$y
      candidates[[2L]] <- hard
    }
    evidence <- v3_court_score_candidates(candidates, source, spec)
    pair <- lapply(candidates, function(candidate) {
      eyesim:::baseline_pair_features(
        candidate, source, source, baseline_spec, availability
      )$features
    })
    feature <- do.call(rbind, pair)
    density <- eyesim:::score_gaze_candidates(
      rowMeans(feature[, c(
        "raw_density_sigma_80", "raw_density_sigma_160"
      ), drop = FALSE]),
      1L, names(candidates), candidate_pool_id = "v3-perturbation"
    )
    mm <- eyesim:::score_gaze_candidates(
      feature[, "raw_mm_vector"], 1L, names(candidates),
      candidate_pool_id = "v3-perturbation"
    )
    data.frame(
      factor = factor,
      transport_v3_info_bits = evidence$gaze_info_bits,
      density_info_bits = density$gaze_info_bits,
      multimatch_vector_info_bits = mm$gaze_info_bits,
      transport_v3_rank = evidence$template_rank,
      transport_v3_converged = all(is.finite(evidence$candidates$log_score)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

v3_court_registration <- function(data, specs) {
  signal <- data$source[data$source$condition == "signal", , drop = FALSE]
  warp <- eyesim:::fit_transport_v3_warp(
    data$reference, signal, "item", "fixgroup", "fixgroup", NULL,
    specs$transport_v3
  )
  identity <- eyesim:::identity_gaze_warp_model()
  improvements <- lapply(seq_len(nrow(signal)), function(index) {
    source <- signal[index, , drop = FALSE]
    raw <- eyesim:::score_transport_v3_cv_row(
      source, data$reference, "item", NULL, "fixgroup", "fixgroup",
      NULL, "prior_weight", specs$transport_v3, identity
    )
    registered <- eyesim:::score_transport_v3_cv_row(
      source, data$reference, "item", NULL, "fixgroup", "fixgroup",
      NULL, "prior_weight", specs$transport_v3, warp
    )
    raw_score <- raw$evidence$candidates$log_score
    registered_score <- registered$evidence$candidates$log_score
    true <- raw$evidence$candidates$is_true
    data.frame(
      item = source$item[[1L]],
      true_improvement = registered_score[true] - raw_score[true],
      wrong_improvement = mean(
        registered_score[!true] - raw_score[!true]
      ),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, improvements)
  result$passed <- result$wrong_improvement <= result$true_improvement + 1e-10
  result
}

v3_court_ledger <- function(summary, operating, invariance, perturbation,
                            registration) {
  v3_signal <- summary[
    summary$method == "transport_v3" & summary$condition == "signal", ]
  v3_null <- summary[
    summary$method == "transport_v3" & summary$condition == "null", ]
  v3_operating <- operating[operating$method == "transport_v3", ]
  value <- function(factor) {
    perturbation$transport_v3_info_bits[perturbation$factor == factor]
  }
  ledger <- data.frame(
    gate = c(
      "representation_invariance", "null_mean_information",
      "null_top1_prior", "heldout_ece", "matched_false_positive_rate",
      "finite_converged_rate", "harmless_time_dilation",
      "fixed_density_chronology_destruction", "partial_coherent_intermediate",
      "registration_wrong_candidate_guard"
    ),
    value = c(
      as.numeric(all(invariance$passed)), abs(v3_null$mean_info_bits),
      abs(v3_null$top1_credit - 1 / v3_null$candidate_count), v3_null$ece_10,
      v3_operating$false_positive_rate, v3_signal$convergence_rate,
      abs(value("global_time_dilation") - value("intact")),
      value("intact") - value("reversal"),
      as.numeric(
        value("coherent_short_subsequence") < value("intact") &&
          value("coherent_short_subsequence") > value("reversal")
      ),
      as.numeric(all(registration$passed))
    ),
    comparator = c("eq", "le", "le", "le", "le", "ge", "le",
                   "gt", "eq", "eq"),
    threshold = c(1, 0.02, 0.02, 0.10, 0.075, 0.99, 1e-8, 0, 1, 1),
    stringsAsFactors = FALSE
  )
  ledger$passed <- mapply(function(value, comparator, threshold) {
    switch(comparator, le = value <= threshold, ge = value >= threshold,
           gt = value > threshold, eq = value == threshold)
  }, ledger$value, ledger$comparator, ledger$threshold)
  ledger$failure_action <- ifelse(
    ledger$passed, "none",
    "revise_algorithm_before_study_or_retrieval_scoring"
  )
  ledger
}

v3_court_write <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(
    "predictive-summary.csv" = result$summary,
    "operating-characteristics.csv" = result$operating,
    "bootstrap-information.csv" = result$bootstrap,
    "invariance.csv" = result$invariance,
    "perturbations.csv" = result$perturbation,
    "registration.csv" = result$registration,
    "resources.csv" = result$resources,
    "gate-ledger.csv" = result$ledger,
    "fold-receipts.csv" = result$fold_receipts,
    "court-manifest.csv" = result$manifest
  )
  for (name in names(tables)) {
    utils::write.csv(tables[[name]], file.path(output_dir, name),
                     row.names = FALSE)
  }
  saveRDS(
    result[c(
      "summary", "operating", "bootstrap", "invariance", "perturbation",
      "registration", "resources", "ledger", "fold_receipts", "manifest"
    )],
    file.path(output_dir, "scientific-results.rds"), version = 3
  )
  utils::capture.output(
    utils::sessionInfo(), file = file.path(output_dir, "session-info.txt")
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[basename(files) != "manifest-md5.csv"]
  hashes <- data.frame(
    file = basename(files), md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    hashes, file.path(output_dir, "manifest-md5.csv"), row.names = FALSE
  )
  invisible(result)
}

run_gaze_weave_transport_v3_court <- function(
    output_dir = NULL, seed = 20260822L, smoke = FALSE,
    n_items = if (smoke) 4L else 8L) {
  if (!smoke && seed != 20260822L) {
    stop("The frozen final court uses seed 20260822; use smoke for debugging.")
  }
  data <- v3_court_data(n_items, seed)
  specs <- v3_court_specs(smoke)
  predictive <- v3_court_predictive(data, specs, seed)
  summaries <- v3_court_summaries(predictive$scored)
  invariance <- v3_court_invariance(seed)
  perturbation <- v3_court_perturbations(data, seed)
  registration <- v3_court_registration(data, specs)
  ledger <- v3_court_ledger(
    summaries$summary, summaries$operating, invariance,
    perturbation, registration
  )
  fold_receipts <- do.call(rbind, lapply(predictive$receipts, function(receipt) {
    data.frame(
      fold = receipt$fold,
      train_items = paste(receipt$train_items, collapse = ";"),
      eval_items = paste(receipt$eval_items, collapse = ";"),
      overlap_item_n = receipt$overlap_item_n,
      candidate_count = length(receipt$candidate_keys),
      candidate_prior_sum = sum(receipt$candidate_prior),
      v3_inner_overlap_max = max(vapply(
        receipt$v3_inner, `[[`, integer(1), "overlap_match_n"
      )),
      baseline_inner_overlap_max = if (length(receipt$baseline_inner)) {
        max(vapply(
          receipt$baseline_inner, `[[`, integer(1), "overlap_match_n"
        ))
      } else {
        NA_integer_
      },
      stringsAsFactors = FALSE
    )
  }))
  manifest <- data.frame(
    protocol_version = "3.0.2",
    seed = seed,
    smoke = smoke,
    n_items = n_items,
    candidate_policy = "shared_exhaustive",
    prior_policy = "shared_actual_uniform_design",
    outer_folds = 2L,
    bootstrap_seed = 20260822L,
    bootstrap_replicates = if (smoke) 100L else 500L,
    retrieval_fields_read = FALSE,
    private_data_read = FALSE,
    passed = all(ledger$passed),
    stringsAsFactors = FALSE
  )
  result <- list(
    summary = summaries$summary,
    operating = summaries$operating,
    bootstrap = v3_court_bootstrap(
      predictive$scored, 20260822L, if (smoke) 100L else 500L
    ),
    invariance = invariance,
    perturbation = perturbation,
    registration = registration,
    resources = predictive$resources,
    ledger = ledger,
    fold_receipts = fold_receipts,
    manifest = manifest,
    scored = predictive$scored,
    fits = predictive$models
  )
  if (!is.null(output_dir)) v3_court_write(result, output_dir)
  result
}
