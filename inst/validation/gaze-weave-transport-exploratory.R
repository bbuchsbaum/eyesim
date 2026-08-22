# Post-court exploratory analysis of canonical Transport recognition effects.
#
# This analysis consumes the frozen GW-18.10 out-of-fold score checkpoints. It
# does not alter the Transport estimator or the frozen retrieval conclusion.

transport_exploratory_seed <- 20260826L
transport_exploratory_draws <- 2000L
transport_exploratory_kappa <- c(0, 0.25, 0.5, 1, 2, 4, 8, 16)
transport_exploratory_result_dir <- file.path(
  "inst", "validation", "gaze-weave-transport-v3-retrieval-results",
  "exploratory"
)

transport_exploratory_alignment_grid <- function() {
  data.frame(
    spec_id = c(
      "default", "temporal_0", "temporal_1", "temporal_4",
      "neighbours_1", "neighbours_3", "spatial_0.67", "spatial_1.5",
      "coverage_uniform", "coverage_concentrated"
    ),
    temporal_weight = c(2, 0, 1, 4, 2, 2, 2, 2, 2, 2),
    neighbours = c(2, 2, 2, 2, 1, 3, 2, 2, 2, 2),
    spatial_scale = c(1, 1, 1, 1, 1, 1, 0.67, 1.5, 1, 1),
    coverage_a = c(2, 2, 2, 2, 2, 2, 2, 2, 1, 4),
    coverage_b = c(2, 2, 2, 2, 2, 2, 2, 2, 1, 4),
    stringsAsFactors = FALSE
  )
}

transport_exploratory_weighted_mean <- function(values, weights) {
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  if (!any(keep)) return(NA_real_)
  sum(values[keep] * weights[keep]) / sum(weights[keep])
}

transport_exploratory_interval <- function(values, probability = 0.95) {
  values <- values[is.finite(values)]
  if (!length(values)) return(c(lower = NA_real_, upper = NA_real_))
  alpha <- 1 - probability
  stats::setNames(
    unname(stats::quantile(
      values, c(alpha / 2, 1 - alpha / 2), type = 8
    )),
    c("lower", "upper")
  )
}

transport_exploratory_bootstrap_plan <- function(
    tab, draws = transport_exploratory_draws,
    seed = transport_exploratory_seed) {
  tab <- tab[order(tab$participant, tab$item), , drop = FALSE]
  participants <- sort(unique(as.character(tab$participant)))
  items <- sort(unique(as.integer(tab$item)))
  set.seed(seed)
  weights <- vapply(seq_len(draws), function(draw) {
    participant_frequency <- table(sample(
      participants, length(participants), replace = TRUE
    ))
    item_frequency <- table(sample(items, length(items), replace = TRUE))
    participant_weight <- as.numeric(
      participant_frequency[as.character(tab$participant)]
    )
    item_weight <- as.numeric(item_frequency[as.character(tab$item)])
    participant_weight[is.na(participant_weight)] <- 0
    item_weight[is.na(item_weight)] <- 0
    participant_weight * item_weight
  }, numeric(nrow(tab)))
  if (draws == 1L) weights <- matrix(weights, ncol = 1L)
  list(
    weights = weights,
    trial_keys = paste(tab$participant, tab$item, sep = ":"),
    participants = participants,
    items = items,
    seed = seed,
    draws = draws
  )
}

transport_exploratory_align_weights <- function(tab, plan) {
  keys <- paste(tab$participant, tab$item, sep = ":")
  position <- match(keys, plan$trial_keys)
  if (anyNA(position) || anyDuplicated(position) || length(position) != nrow(tab)) {
    stop("The crossed bootstrap plan does not match the trial support.")
  }
  plan$weights[position, , drop = FALSE]
}

transport_exploratory_prepare <- function(scores, metadata) {
  required_scores <- c(
    "participant", "item", "outer_fold", "gaze_info_bits",
    "effective_fixations", "total_duration", "candidates"
  )
  required_metadata <- c(
    "participant", "item", "probe_type", "degradation", "Accuracy"
  )
  if (length(setdiff(required_scores, names(scores))) ||
      length(setdiff(required_metadata, names(metadata)))) {
    stop("Transport scores or recognition metadata have an unexpected schema.")
  }
  metadata <- unique(metadata[required_metadata])
  if (anyDuplicated(metadata[c("participant", "item")])) {
    stop("Recognition metadata are not unique by participant and item.")
  }
  tab <- merge(
    scores, metadata, by = c("participant", "item"),
    all.x = TRUE, sort = FALSE
  )
  if (anyNA(tab$Accuracy) || any(!tab$Accuracy %in% 0:1) ||
      any(!tab$probe_type %in% c("old", "lure")) ||
      any(!tab$degradation %in% c(20, 40, 60, 80, 100))) {
    stop("Recognition outcomes violate the full-cohort value contract.")
  }
  tab$accuracy <- as.integer(tab$Accuracy)
  tab$saliency_z <- (as.numeric(tab$degradation) - 60) / 20
  tab$correct_ec <- tab$accuracy - 0.5
  tab$probe_type_ec <- ifelse(tab$probe_type == "old", 0.5, -0.5)
  tab$z_quality <- as.numeric(scale(log(tab$effective_fixations)))
  tab$z_duration <- as.numeric(scale(tab$total_duration))
  tab$trial_key <- paste(tab$participant, tab$item, sep = ":")
  tab[order(tab$participant, tab$item), , drop = FALSE]
}

transport_exploratory_design <- function(tab) {
  stats::model.matrix(
    ~ saliency_z * correct_ec * probe_type_ec + z_quality + z_duration,
    data = tab
  )
}

transport_exploratory_fit_association <- function(tab, weights = NULL) {
  design <- transport_exploratory_design(tab)
  if (is.null(weights)) weights <- rep(1, nrow(tab))
  keep <- is.finite(weights) & weights > 0
  if (sum(keep) < ncol(design)) {
    return(list(status = "insufficient", coefficients = NULL, rank = 0L))
  }
  fit <- stats::lm.wfit(
    design[keep, , drop = FALSE], tab$gaze_info_bits[keep], weights[keep]
  )
  coefficients <- stats::setNames(
    as.numeric(fit$coefficients), colnames(design)
  )
  valid <- fit$rank == ncol(design) && all(is.finite(coefficients))
  list(
    status = if (valid) "scored" else "rank_deficient",
    coefficients = if (valid) coefficients else NULL,
    rank = fit$rank
  )
}

transport_exploratory_contrasts <- function(coefficients) {
  terms <- c(
    saliency_20_to_100 = "saliency_z",
    correct_at_60 = "correct_ec",
    old_minus_lure_at_60 = "probe_type_ec",
    saliency_by_correct_per_20 = "saliency_z:correct_ec",
    saliency_by_probe_per_20 = "saliency_z:probe_type_ec",
    correct_by_probe_at_60 = "correct_ec:probe_type_ec",
    three_way_per_20 = "saliency_z:correct_ec:probe_type_ec"
  )
  if (is.null(coefficients) || !all(terms %in% names(coefficients))) {
    return(stats::setNames(rep(NA_real_, length(terms)), names(terms)))
  }
  result <- coefficients[terms]
  names(result) <- names(terms)
  result[["saliency_20_to_100"]] <- 4 * result[["saliency_20_to_100"]]
  result
}

transport_exploratory_association <- function(tab, plan) {
  fixed <- transport_exploratory_fit_association(tab)
  if (!identical(fixed$status, "scored")) {
    stop("The full Transport interaction design is rank deficient.")
  }
  estimate <- transport_exploratory_contrasts(fixed$coefficients)
  weights <- transport_exploratory_align_weights(tab, plan)
  bootstrap <- vapply(seq_len(ncol(weights)), function(draw) {
    transport_exploratory_contrasts(
      transport_exploratory_fit_association(tab, weights[, draw])$coefficients
    )
  }, numeric(length(estimate)))
  rownames(bootstrap) <- names(estimate)
  participants <- sort(unique(tab$participant))
  lopo <- vapply(participants, function(participant) {
    transport_exploratory_contrasts(
      transport_exploratory_fit_association(
        tab[tab$participant != participant, , drop = FALSE]
      )$coefficients
    )
  }, numeric(length(estimate)))
  rownames(lopo) <- names(estimate)
  family_probability <- 1 - 0.05 / length(estimate)
  effects <- do.call(rbind, lapply(names(estimate), function(name) {
    ordinary <- transport_exploratory_interval(bootstrap[name, ], 0.95)
    family <- transport_exploratory_interval(
      bootstrap[name, ], family_probability
    )
    data.frame(
      contrast = name,
      estimate = estimate[[name]],
      lower_95 = ordinary[["lower"]],
      upper_95 = ordinary[["upper"]],
      lower_family = family[["lower"]],
      upper_family = family[["upper"]],
      valid_draws = sum(is.finite(bootstrap[name, ])),
      lopo_sign_n = sum(
        sign(lopo[name, ]) == sign(estimate[[name]]), na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
  list(
    effects = effects,
    coefficients = fixed$coefficients,
    bootstrap = bootstrap,
    lopo = lopo,
    design_kappa = kappa(transport_exploratory_design(tab), exact = TRUE)
  )
}

transport_exploratory_log_loss <- function(response, probability) {
  probability <- pmin(pmax(probability, 1e-8), 1 - 1e-8)
  -(response * log(probability) + (1 - response) * log1p(-probability))
}

transport_exploratory_standardize <- function(train, evaluate, source, target) {
  center <- mean(train[[source]])
  scale_value <- stats::sd(train[[source]])
  if (!is.finite(center) || !is.finite(scale_value) || scale_value <= 0) {
    stop("Invalid training scale for ", source, ".")
  }
  train[[target]] <- (train[[source]] - center) / scale_value
  evaluate[[target]] <- (evaluate[[source]] - center) / scale_value
  list(train = train, evaluate = evaluate)
}

transport_exploratory_fit_behavior_fold <- function(train, evaluate) {
  if (length(intersect(train$participant, evaluate$participant)) ||
      length(intersect(train$item, evaluate$item))) {
    stop("Behavior training and evaluation support overlaps.")
  }
  train$log_effective_fixations <- log(train$effective_fixations)
  evaluate$log_effective_fixations <- log(evaluate$effective_fixations)
  for (mapping in list(
    c("gaze_info_bits", "z_metric"),
    c("log_effective_fixations", "z_quality_fold"),
    c("total_duration", "z_duration_fold")
  )) {
    value <- transport_exploratory_standardize(
      train, evaluate, mapping[[1L]], mapping[[2L]]
    )
    train <- value$train
    evaluate <- value$evaluate
  }
  base_formula <- accuracy ~ saliency_z * probe_type_ec +
    z_quality_fold + z_duration_fold
  additive_formula <- stats::update.formula(base_formula, ". ~ . + z_metric")
  heterogeneous_formula <- stats::update.formula(
    base_formula, ". ~ . + z_metric * saliency_z * probe_type_ec"
  )
  fits <- lapply(
    list(base_formula, additive_formula, heterogeneous_formula),
    function(formula) stats::glm(
      formula, family = stats::binomial(), data = train
    )
  )
  probability <- lapply(
    fits, stats::predict, newdata = evaluate, type = "response"
  )
  loss <- lapply(probability, function(value) {
    transport_exploratory_log_loss(evaluate$accuracy, value)
  })
  data.frame(
    participant = evaluate$participant,
    item = evaluate$item,
    trial_key = evaluate$trial_key,
    outer_fold = evaluate$outer_fold,
    accuracy = evaluate$accuracy,
    probe_type = evaluate$probe_type,
    degradation = evaluate$degradation,
    gaze_info_bits = evaluate$gaze_info_bits,
    loss_base = loss[[1L]],
    loss_additive = loss[[2L]],
    loss_heterogeneous = loss[[3L]],
    gain_additive_bits = (loss[[1L]] - loss[[2L]]) / log(2),
    gain_heterogeneous_bits = (loss[[1L]] - loss[[3L]]) / log(2),
    gain_heterogeneity_bits = (loss[[2L]] - loss[[3L]]) / log(2),
    stringsAsFactors = FALSE
  )
}

transport_exploratory_behavior_crossfit <- function(tab) {
  rows <- lapply(sort(unique(tab$outer_fold)), function(fold_id) {
    evaluate <- tab[tab$outer_fold == fold_id, , drop = FALSE]
    train <- tab[
      !tab$participant %in% evaluate$participant &
        !tab$item %in% evaluate$item, , drop = FALSE
    ]
    transport_exploratory_fit_behavior_fold(train, evaluate)
  })
  result <- do.call(rbind, rows)
  if (nrow(result) != nrow(tab) || anyDuplicated(result$trial_key)) {
    stop("Behavior cross-fitting did not return one row per trial.")
  }
  result[order(result$participant, result$item), , drop = FALSE]
}

transport_exploratory_behavior_summary <- function(predictions, plan) {
  weights <- transport_exploratory_align_weights(predictions, plan)
  metrics <- c(
    additive_vs_base = "gain_additive_bits",
    heterogeneous_vs_base = "gain_heterogeneous_bits",
    heterogeneity_vs_additive = "gain_heterogeneity_bits"
  )
  family_probability <- 1 - 0.05 / length(metrics)
  do.call(rbind, lapply(names(metrics), function(name) {
    values <- predictions[[metrics[[name]]]]
    bootstrap <- vapply(seq_len(ncol(weights)), function(draw) {
      transport_exploratory_weighted_mean(values, weights[, draw])
    }, numeric(1))
    ordinary <- transport_exploratory_interval(bootstrap, 0.95)
    family <- transport_exploratory_interval(bootstrap, family_probability)
    data.frame(
      contrast = name,
      mean_information_bits = mean(values),
      lower_95 = ordinary[["lower"]],
      upper_95 = ordinary[["upper"]],
      lower_family = family[["lower"]],
      upper_family = family[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
}

transport_exploratory_kappa_scores <- function(tab, kappa) {
  if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
      kappa < 0) {
    stop("kappa must be one finite non-negative value.")
  }
  base_true <- vapply(tab$candidates, function(candidate) {
    candidate$base_posterior[candidate$is_true][[1L]]
  }, numeric(1))
  prior_true <- vapply(tab$candidates, function(candidate) {
    candidate$prior[candidate$is_true][[1L]]
  }, numeric(1))
  reliability <- tab$effective_fixations / (tab$effective_fixations + kappa)
  posterior_true <- reliability * base_true + (1 - reliability) * prior_true
  result <- tab
  result$gaze_info_bits <- log2(posterior_true / prior_true)
  result$kappa <- kappa
  result
}

transport_exploratory_verify_scores <- function(result_dir) {
  manifest <- utils::read.csv(
    file.path(result_dir, "checkpoint-manifest.csv"),
    stringsAsFactors = FALSE
  )
  paths <- file.path(result_dir, manifest$file)
  if (!all(file.exists(paths)) ||
      !identical(unname(tools::md5sum(paths)), manifest$md5)) {
    stop("Frozen Transport checkpoint hashes do not verify.")
  }
  checkpoint_scores <- do.call(rbind, lapply(paths, function(path) {
    readRDS(path)$scored
  }))
  measurement <- readRDS(file.path(result_dir, "measurement-scores.rds"))
  core <- c(
    "participant", "item", "outer_fold", "gaze_info_bits", "log_loss"
  )
  order_core <- function(tab) {
    tab[order(tab$participant, tab$item), core, drop = FALSE]
  }
  if (!identical(order_core(checkpoint_scores), order_core(measurement))) {
    stop("Measurement scores differ from the hashed fold checkpoints.")
  }
  measurement
}

transport_exploratory_read_data <- function(
    result_dir = dirname(transport_exploratory_result_dir)) {
  scores <- transport_exploratory_verify_scores(result_dir)
  full_file <- system.file(
    "validation", "gaze-weave-recognition-full-cohort.R", package = "eyesim"
  )
  if (!nzchar(full_file)) {
    full_file <- file.path(
      "inst", "validation", "gaze-weave-recognition-full-cohort.R"
    )
  }
  if (!file.exists(full_file)) {
    stop("The full-cohort recognition input reader is unavailable.")
  }
  full_file <- normalizePath(full_file, mustWork = TRUE)
  source_root <- normalizePath(
    file.path(dirname(full_file), "..", ".."), mustWork = TRUE
  )
  previous_directory <- setwd(source_root)
  on.exit(setwd(previous_directory), add = TRUE)
  environment <- new.env(parent = globalenv())
  sys.source(full_file, envir = environment)
  raw <- environment$full_recognition_read_inputs(verify = TRUE)
  metadata <- unique(raw$retrieval[c(
    "participant", "item", "probe_type", "degradation", "Accuracy"
  )])
  tab <- transport_exploratory_prepare(scores, metadata)
  if (nrow(tab) != 1295L || length(unique(tab$participant)) != 36L ||
      length(unique(tab$item)) != 120L) {
    stop("The exploratory cohort differs from frozen GW-18.10 support.")
  }
  tab
}

transport_exploratory_support <- function(tab) {
  counts <- as.data.frame(table(
    saliency = tab$degradation,
    accuracy = tab$accuracy,
    probe_type = tab$probe_type
  ), stringsAsFactors = FALSE)
  counts[counts$Freq > 0L, , drop = FALSE]
}

transport_exploratory_write <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  write(result$support, "support-counts.csv")
  write(result$association, "association-effects.csv")
  write(result$behavior, "behavior-prediction.csv")
  write(result$kappa, "kappa-sensitivity.csv")
  write(transport_exploratory_alignment_grid(), "alignment-grid.csv")
  write(result$configuration, "configuration.csv")
  saveRDS(result$predictions, file.path(output_dir, "predictions.rds"), version = 3)
  saveRDS(
    result[c("support", "association", "behavior", "kappa", "configuration")],
    file.path(output_dir, "scientific-results.rds"), version = 3
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[basename(files) != "manifest-md5.csv"]
  write(data.frame(
    file = basename(files),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  ), "manifest-md5.csv")
  invisible(result)
}

run_gaze_weave_transport_exploratory <- function(
    output_dir = transport_exploratory_result_dir,
    draws = transport_exploratory_draws,
    seed = transport_exploratory_seed,
    kappa_grid = transport_exploratory_kappa) {
  tab <- transport_exploratory_read_data()
  plan <- transport_exploratory_bootstrap_plan(tab, draws, seed)
  association <- transport_exploratory_association(tab, plan)
  predictions <- transport_exploratory_behavior_crossfit(tab)
  behavior <- transport_exploratory_behavior_summary(predictions, plan)

  kappa_rows <- lapply(kappa_grid, function(kappa) {
    message("Exploratory reliability kappa: ", kappa)
    candidate <- transport_exploratory_kappa_scores(tab, kappa)
    candidate_association <- transport_exploratory_association(candidate, plan)
    candidate_predictions <- transport_exploratory_behavior_crossfit(candidate)
    candidate_behavior <- transport_exploratory_behavior_summary(
      candidate_predictions, plan
    )
    three_way <- candidate_association$effects[
      candidate_association$effects$contrast == "three_way_per_20",
    ]
    heterogeneous <- candidate_behavior[
      candidate_behavior$contrast == "heterogeneous_vs_base",
    ]
    data.frame(
      kappa = kappa,
      mean_gaze_info_bits = mean(candidate$gaze_info_bits),
      three_way_estimate = three_way$estimate,
      three_way_lower_95 = three_way$lower_95,
      three_way_upper_95 = three_way$upper_95,
      behavior_gain_bits = heterogeneous$mean_information_bits,
      behavior_lower_95 = heterogeneous$lower_95,
      behavior_upper_95 = heterogeneous$upper_95,
      stringsAsFactors = FALSE
    )
  })
  result <- list(
    support = transport_exploratory_support(tab),
    association = association$effects,
    behavior = behavior,
    kappa = do.call(rbind, kappa_rows),
    predictions = predictions,
    configuration = data.frame(
      status = "post_court_exploratory",
      seed = seed,
      bootstrap_draws = draws,
      trials = nrow(tab),
      participants = length(unique(tab$participant)),
      items = length(unique(tab$item)),
      design_kappa = association$design_kappa,
      frozen_primary_unchanged = TRUE,
      local_only = TRUE,
      stringsAsFactors = FALSE
    )
  )
  transport_exploratory_write(result, output_dir)
}
