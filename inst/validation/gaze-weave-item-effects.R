# Known-item/new-participant item reinstatement court.

item_effect_seed <- 20260828L
item_effect_draws <- 2000L
item_effect_result_dir <- file.path(
  "inst", "validation", "gaze-weave-transport-v3-retrieval-results",
  "item-effects"
)
item_effect_study_dir <- file.path(
  "inst", "validation", "gaze-weave-transport-v3-repeated-viewing-results"
)

item_effect_dependency <- system.file(
  "validation", "gaze-weave-transport-density-diagnostic.R",
  package = "eyesim"
)
if (!nzchar(item_effect_dependency)) {
  item_effect_dependency <- file.path(
    "inst", "validation", "gaze-weave-transport-density-diagnostic.R"
  )
}
if (!file.exists(item_effect_dependency)) {
  stop("The Transport-density diagnostic helpers are unavailable.")
}
source(item_effect_dependency, local = TRUE)

item_effect_shrink <- function(residual, item) {
  residual <- as.numeric(residual)
  item <- as.character(item)
  if (length(residual) != length(item) || !length(residual) ||
      any(!is.finite(residual)) || anyNA(item)) {
    stop("Residuals and item labels must be finite, aligned, and non-empty.")
  }
  groups <- split(residual, item)
  item_n <- vapply(groups, length, integer(1))
  raw_mean <- vapply(groups, mean, numeric(1))
  within_ss <- sum(vapply(groups, function(value) {
    sum((value - mean(value))^2)
  }, numeric(1)))
  within_df <- length(residual) - length(groups)
  sigma2 <- if (within_df > 0L) within_ss / within_df else 0
  sampling_variance <- sigma2 / item_n
  tau2 <- max(0, stats::var(raw_mean) - mean(sampling_variance))
  weight <- if (tau2 <= 0) {
    rep(0, length(item_n))
  } else {
    tau2 / (tau2 + sampling_variance)
  }
  data.frame(
    item = names(groups),
    source_n = as.integer(item_n),
    raw_mean = as.numeric(raw_mean),
    shrinkage_weight = as.numeric(weight),
    item_propensity = as.numeric(weight * raw_mean),
    tau2 = tau2,
    sigma2 = sigma2,
    stringsAsFactors = FALSE
  )
}

item_effect_nuisance_fit <- function(train, outcome, density = NULL) {
  train$log_effective_fixations <- log(train$effective_fixations)
  candidates <- c(
    "saliency_z", "probe_type", "log_effective_fixations",
    "total_duration", density
  )
  informative <- vapply(candidates, function(variable) {
    value <- train[[variable]]
    length(unique(value[!is.na(value)])) > 1L
  }, logical(1))
  formula <- stats::reformulate(
    candidates[informative], response = outcome
  )
  fit <- stats::lm(formula, data = train)
  if (fit$rank != length(stats::coef(fit)) ||
      any(!is.finite(stats::coef(fit)))) {
    stop("The item-effect nuisance model is rank deficient.")
  }
  fit
}

item_effect_nuisance_predict <- function(fit, tab) {
  tab$log_effective_fixations <- log(tab$effective_fixations)
  as.numeric(stats::predict(fit, newdata = tab))
}

item_effect_assign <- function(tab, propensity) {
  position <- match(as.character(tab$item), propensity$item)
  if (anyNA(position)) stop("A held-out item has no training-participant source.")
  data.frame(
    item_propensity = propensity$item_propensity[position],
    source_n = propensity$source_n[position],
    tau2 = propensity$tau2[position],
    sigma2 = propensity$sigma2[position],
    stringsAsFactors = FALSE
  )
}

item_effect_crossparticipant <- function(tab, outcome, density = NULL) {
  participants <- sort(unique(as.character(tab$participant)))
  rows <- lapply(participants, function(participant) {
    train <- tab[as.character(tab$participant) != participant, , drop = FALSE]
    evaluate <- tab[as.character(tab$participant) == participant, , drop = FALSE]
    fit <- item_effect_nuisance_fit(train, outcome, density)
    train_fixed <- item_effect_nuisance_predict(fit, train)
    evaluate_fixed <- item_effect_nuisance_predict(fit, evaluate)
    train_residual <- train[[outcome]] - train_fixed
    evaluate_residual <- evaluate[[outcome]] - evaluate_fixed
    propensity <- item_effect_shrink(train_residual, train$item)
    assigned <- item_effect_assign(evaluate, propensity)
    if (participant %in% unique(as.character(train$participant))) {
      stop("Held-out participant leaked into the item propensity source.")
    }
    data.frame(
      participant = evaluate$participant,
      item = evaluate$item,
      trial_key = evaluate$trial_key,
      probe_type = evaluate$probe_type,
      degradation = evaluate$degradation,
      observed = evaluate[[outcome]],
      fixed_prediction = evaluate_fixed,
      residual = evaluate_residual,
      item_propensity = assigned$item_propensity,
      item_prediction = evaluate_fixed + assigned$item_propensity,
      squared_gain = evaluate_residual^2 -
        (evaluate_residual - assigned$item_propensity)^2,
      source_n = assigned$source_n,
      tau2 = assigned$tau2,
      sigma2 = assigned$sigma2,
      source_participants = length(unique(train$participant)),
      self_participant_overlap = 0L,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  if (nrow(result) != nrow(tab) || anyDuplicated(result$trial_key) ||
      any(result$self_participant_overlap != 0L) ||
      any(result$source_n < 1L)) {
    stop("Cross-participant item predictions violate their support contract.")
  }
  result[order(result$participant, result$item), , drop = FALSE]
}

item_effect_weighted_r2 <- function(prediction, weights) {
  numerator <- transport_exploratory_weighted_mean(
    prediction$squared_gain, weights
  )
  denominator <- transport_exploratory_weighted_mean(
    prediction$residual^2, weights
  )
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  numerator / denominator
}

item_effect_score_summary <- function(predictions, draws, seed) {
  reference <- predictions[["transport"]]
  plan <- transport_exploratory_bootstrap_plan(reference, draws, seed)
  weights <- transport_exploratory_align_weights(reference, plan)
  bootstrap <- list()
  rows <- lapply(names(predictions), function(method) {
    part <- predictions[[method]]
    part <- part[match(plan$trial_keys, part$trial_key), ]
    values <- vapply(seq_len(ncol(weights)), function(draw) {
      item_effect_weighted_r2(part, weights[, draw])
    }, numeric(1))
    bootstrap[[method]] <<- values
    interval <- transport_exploratory_interval(values, 0.95)
    data.frame(
      method = method,
      trials = nrow(part),
      participants = length(unique(part$participant)),
      items = length(unique(part$item)),
      min_item_sources = min(part$source_n),
      mean_item_sources = mean(part$source_n),
      mean_tau2 = mean(part$tau2),
      mean_squared_gain = mean(part$squared_gain),
      cv_r2 = sum(part$squared_gain) / sum(part$residual^2),
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      pearson = stats::cor(part$item_propensity, part$residual),
      spearman = stats::cor(
        part$item_propensity, part$residual, method = "spearman"
      ),
      stringsAsFactors = FALSE
    )
  })
  densities <- setdiff(names(predictions), c("transport", "transport_specific"))
  contrast_rows <- lapply(densities, function(method) {
    values <- bootstrap[["transport"]] - bootstrap[[method]]
    interval <- transport_exploratory_interval(values, 0.95)
    data.frame(
      contrast = paste("transport_minus", method, sep = "_"),
      cv_r2_difference = mean(bootstrap[["transport"]], na.rm = TRUE) -
        mean(bootstrap[[method]], na.rm = TRUE),
      lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
      stringsAsFactors = FALSE
    )
  })
  list(
    methods = do.call(rbind, rows),
    contrasts = do.call(rbind, contrast_rows),
    bootstrap = bootstrap
  )
}

item_effect_read_study <- function(
    study_dir = item_effect_study_dir,
    retrieval_reference = NULL) {
  tasks <- c(
    "study_p1_p2", "study_p1_p3", "study_p1_p4",
    "study_p2_p3", "study_p3_p4"
  )
  paths <- unlist(lapply(tasks, function(task) {
    file.path(study_dir, sprintf(
      "checkpoint-%s-fold-%02d.rds", task, 1:4
    ))
  }), use.names = FALSE)
  if (length(paths) != 20L || !all(file.exists(paths))) {
    stop("The complete 20-checkpoint repeated-viewing court is unavailable.")
  }
  rows <- lapply(paths, function(path) {
    value <- readRDS(path)$scored
    value[value$control == "intact", c(
      "task", "participant", "item", "gaze_info_bits",
      "effective_fixations", "outer_fold"
    )]
  })
  scored <- do.call(rbind, rows)
  if (nrow(scored) != 5L * 1295L ||
      any(table(scored$task) != 1295L) ||
      anyDuplicated(scored[c("task", "participant", "item")])) {
    stop("Repeated-viewing intact support differs from the frozen court.")
  }
  aggregate_mean <- function(value) mean(value)
  study <- stats::aggregate(
    cbind(
      study_info_bits = scored$gaze_info_bits,
      study_effective_fixations = scored$effective_fixations
    ) ~ participant + item,
    data = scored, FUN = aggregate_mean
  )
  names(study)[names(study) == "study_info_bits"] <- "score_study"
  if (is.null(retrieval_reference)) return(study)
  metadata <- unique(retrieval_reference[c(
    "participant", "item", "probe_type", "degradation", "saliency_z",
    "trial_key"
  )])
  study <- merge(study, metadata, by = c("participant", "item"),
                 all.x = TRUE, sort = FALSE)
  if (nrow(study) != 1295L || anyNA(study$saliency_z)) {
    stop("Study and retrieval support do not align.")
  }
  study$effective_fixations <- study$study_effective_fixations
  study$total_duration <- 1
  study[order(study$participant, study$item), , drop = FALSE]
}

item_effect_study_propensity <- function(study_source) {
  study_source$log_effective_fixations <- log(
    study_source$study_effective_fixations
  )
  fit <- stats::lm(
    score_study ~ saliency_z + probe_type + log_effective_fixations,
    data = study_source
  )
  residual <- study_source$score_study -
    stats::predict(fit, newdata = study_source)
  item_effect_shrink(residual, study_source$item)
}

item_effect_study_transfer <- function(retrieval, study) {
  participants <- sort(unique(as.character(retrieval$participant)))
  rows <- lapply(participants, function(participant) {
    retrieval_train <- retrieval[
      as.character(retrieval$participant) != participant, , drop = FALSE
    ]
    retrieval_evaluate <- retrieval[
      as.character(retrieval$participant) == participant, , drop = FALSE
    ]
    study_train <- study[
      as.character(study$participant) != participant, , drop = FALSE
    ]
    retrieval_fit <- item_effect_nuisance_fit(
      retrieval_train, "score_transport"
    )
    retrieval_fixed <- item_effect_nuisance_predict(
      retrieval_fit, retrieval_evaluate
    )
    retrieval_train_fixed <- item_effect_nuisance_predict(
      retrieval_fit, retrieval_train
    )
    retrieval_residual <- retrieval_evaluate$score_transport - retrieval_fixed
    retrieval_train_residual <- retrieval_train$score_transport -
      retrieval_train_fixed

    propensity <- item_effect_study_propensity(study_train)
    assigned <- item_effect_assign(retrieval_evaluate, propensity)
    train_propensity <- numeric(nrow(retrieval_train))
    for (source_participant in unique(
      as.character(retrieval_train$participant)
    )) {
      retrieval_index <-
        as.character(retrieval_train$participant) == source_participant
      source <- study_train[
        as.character(study_train$participant) != source_participant,
        , drop = FALSE
      ]
      source_propensity <- item_effect_study_propensity(source)
      train_propensity[retrieval_index] <- item_effect_assign(
        retrieval_train[retrieval_index, , drop = FALSE], source_propensity
      )$item_propensity
    }
    calibration <- stats::lm(retrieval_train_residual ~ train_propensity)
    calibrated <- as.numeric(stats::predict(
      calibration,
      newdata = data.frame(train_propensity = assigned$item_propensity)
    ))
    data.frame(
      participant = retrieval_evaluate$participant,
      item = retrieval_evaluate$item,
      trial_key = retrieval_evaluate$trial_key,
      probe_type = retrieval_evaluate$probe_type,
      degradation = retrieval_evaluate$degradation,
      observed = retrieval_evaluate$score_transport,
      fixed_prediction = retrieval_fixed,
      residual = retrieval_residual,
      source_propensity = assigned$item_propensity,
      item_propensity = calibrated,
      item_prediction = retrieval_fixed + calibrated,
      squared_gain = retrieval_residual^2 -
        (retrieval_residual - calibrated)^2,
      source_n = assigned$source_n,
      tau2 = assigned$tau2,
      sigma2 = assigned$sigma2,
      source_participants = length(unique(study_train$participant)),
      self_participant_overlap = 0L,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  if (nrow(result) != nrow(retrieval) || anyDuplicated(result$trial_key)) {
    stop("Study-to-retrieval transfer support is incomplete.")
  }
  result[order(result$participant, result$item), , drop = FALSE]
}

item_effect_full_item_table <- function(tab, outcome, density = NULL,
                                        label = outcome) {
  fit <- item_effect_nuisance_fit(tab, outcome, density)
  residual <- tab[[outcome]] - item_effect_nuisance_predict(fit, tab)
  result <- item_effect_shrink(residual, tab$item)
  names(result)[names(result) == "raw_mean"] <- label
  result
}

item_effect_transfer_summary <- function(prediction, draws, seed) {
  plan <- transport_exploratory_bootstrap_plan(prediction, draws, seed)
  weights <- transport_exploratory_align_weights(prediction, plan)
  values <- vapply(seq_len(ncol(weights)), function(draw) {
    item_effect_weighted_r2(prediction, weights[, draw])
  }, numeric(1))
  interval <- transport_exploratory_interval(values, 0.95)
  data.frame(
    source = "study_repeated_viewing",
    target = "retrieval_transport",
    trials = nrow(prediction),
    cv_r2 = sum(prediction$squared_gain) / sum(prediction$residual^2),
    lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
    pearson = stats::cor(prediction$item_propensity, prediction$residual),
    spearman = stats::cor(
      prediction$item_propensity, prediction$residual, method = "spearman"
    ),
    stringsAsFactors = FALSE
  )
}

item_effect_propensity_for_confidence <- function(
    train, evaluate, outcome, density = NULL) {
  fit <- item_effect_nuisance_fit(train, outcome, density)
  train_residual <- train[[outcome]] - item_effect_nuisance_predict(fit, train)
  evaluate_propensity <- item_effect_assign(
    evaluate, item_effect_shrink(train_residual, train$item)
  )$item_propensity
  train_propensity <- numeric(nrow(train))
  for (participant in unique(as.character(train$participant))) {
    index <- as.character(train$participant) == participant
    source <- train[!index, , drop = FALSE]
    source_fit <- item_effect_nuisance_fit(source, outcome, density)
    source_residual <- source[[outcome]] -
      item_effect_nuisance_predict(source_fit, source)
    propensity <- item_effect_shrink(
      source_residual, source$item
    )
    train_propensity[index] <- item_effect_assign(
      train[index, , drop = FALSE], propensity
    )$item_propensity
  }
  list(train = train_propensity, evaluate = evaluate_propensity)
}

item_effect_study_propensity_for_confidence <- function(
    study_train, retrieval_train, retrieval_evaluate) {
  evaluate <- item_effect_assign(
    retrieval_evaluate, item_effect_study_propensity(study_train)
  )$item_propensity
  train <- numeric(nrow(retrieval_train))
  for (participant in unique(as.character(retrieval_train$participant))) {
    retrieval_index <- as.character(retrieval_train$participant) == participant
    study_index <- as.character(study_train$participant) == participant
    source <- study_train[!study_index, , drop = FALSE]
    propensity <- item_effect_study_propensity(source)
    train[retrieval_index] <- item_effect_assign(
      retrieval_train[retrieval_index, , drop = FALSE], propensity
    )$item_propensity
  }
  list(train = train, evaluate = evaluate)
}

item_effect_safe_standardize <- function(train, evaluate, source, target) {
  center <- mean(train[[source]])
  scale_value <- stats::sd(train[[source]])
  if (!is.finite(center) || !is.finite(scale_value)) {
    stop("A confidence predictor has an invalid training distribution.")
  }
  if (scale_value <= 1e-12) {
    train[[target]] <- 0
    evaluate[[target]] <- 0
  } else {
    train[[target]] <- (train[[source]] - center) / scale_value
    evaluate[[target]] <- (evaluate[[source]] - center) / scale_value
  }
  list(train = train, evaluate = evaluate)
}

item_effect_confidence_fold <- function(
    retrieval, study, participant, probe_type) {
  retrieval_train <- retrieval[
    as.character(retrieval$participant) != participant, , drop = FALSE
  ]
  retrieval_evaluate <- retrieval[
    as.character(retrieval$participant) == participant, , drop = FALSE
  ]
  study_train <- study[as.character(study$participant) != participant, ]

  transport <- item_effect_propensity_for_confidence(
    retrieval_train, retrieval_evaluate, "score_transport"
  )
  density <- item_effect_propensity_for_confidence(
    retrieval_train, retrieval_evaluate, "score_density_sigma_80"
  )
  specific <- item_effect_propensity_for_confidence(
    retrieval_train, retrieval_evaluate, "score_transport",
    density = "score_density_sigma_80"
  )
  study_propensity <- item_effect_study_propensity_for_confidence(
    study_train, retrieval_train, retrieval_evaluate
  )
  retrieval_train$prop_transport <- transport$train
  retrieval_evaluate$prop_transport <- transport$evaluate
  retrieval_train$prop_density <- density$train
  retrieval_evaluate$prop_density <- density$evaluate
  retrieval_train$prop_transport_specific <- specific$train
  retrieval_evaluate$prop_transport_specific <- specific$evaluate
  retrieval_train$prop_study <- study_propensity$train
  retrieval_evaluate$prop_study <- study_propensity$evaluate

  train <- retrieval_train[
    retrieval_train$probe_type == probe_type & !is.na(retrieval_train$oldness),
  ]
  evaluate <- retrieval_evaluate[
    retrieval_evaluate$probe_type == probe_type &
      !is.na(retrieval_evaluate$oldness),
  ]
  if (!nrow(evaluate)) return(NULL)
  train$log_effective_fixations <- log(train$effective_fixations)
  evaluate$log_effective_fixations <- log(evaluate$effective_fixations)
  mappings <- list(
    c("log_effective_fixations", "z_quality"),
    c("total_duration", "z_duration"),
    c("score_transport", "z_individual_transport"),
    c("prop_transport", "z_prop_transport"),
    c("prop_density", "z_prop_density"),
    c("prop_transport_specific", "z_prop_transport_specific"),
    c("prop_study", "z_prop_study")
  )
  for (mapping in mappings) {
    value <- item_effect_safe_standardize(
      train, evaluate, mapping[[1L]], mapping[[2L]]
    )
    train <- value$train
    evaluate <- value$evaluate
  }
  base <- oldness ~ saliency_z + z_quality + z_duration
  formulas <- list(
    base = base,
    individual_transport = stats::update.formula(
      base, ". ~ . + z_individual_transport"
    ),
    retrieval_item = stats::update.formula(base, ". ~ . + z_prop_transport"),
    density_item = stats::update.formula(base, ". ~ . + z_prop_density"),
    transport_specific_item = stats::update.formula(
      base, ". ~ . + z_prop_transport_specific"
    ),
    study_item = stats::update.formula(base, ". ~ . + z_prop_study"),
    retrieval_and_study_item = stats::update.formula(
      base, ". ~ . + z_prop_transport + z_prop_study"
    )
  )
  fits <- lapply(formulas, transport_density_ordinal_fit, data = train)
  losses <- lapply(fits, transport_density_ordinal_loss, evaluate = evaluate)
  result <- data.frame(
    participant = evaluate$participant,
    item = evaluate$item,
    trial_key = evaluate$trial_key,
    probe_type = evaluate$probe_type,
    response = evaluate$Response,
    oldness = as.integer(evaluate$oldness),
    source_participant_overlap = 0L,
    stringsAsFactors = FALSE
  )
  for (name in names(losses)) result[[paste0("loss_", name)]] <- losses[[name]]
  result
}

item_effect_confidence_crossfit <- function(retrieval, study, probe_type) {
  participants <- sort(unique(as.character(retrieval$participant)))
  rows <- lapply(participants, function(participant) {
    message("Known-item confidence: ", probe_type, ", participant ", participant)
    item_effect_confidence_fold(retrieval, study, participant, probe_type)
  })
  result <- do.call(rbind, rows)
  expected <- retrieval[
    retrieval$probe_type == probe_type & !is.na(retrieval$oldness),
  ]
  if (nrow(result) != nrow(expected) || anyDuplicated(result$trial_key) ||
      any(result$source_participant_overlap != 0L)) {
    stop("Known-item confidence predictions violate the leakage contract.")
  }
  result[order(result$participant, result$item), , drop = FALSE]
}

item_effect_confidence_summary <- function(predictions, draws, seed) {
  model_names <- sub("^loss_", "", grep(
    "^loss_", names(predictions), value = TRUE
  ))
  model_names <- setdiff(model_names, "base")
  plan <- transport_exploratory_bootstrap_plan(predictions, draws, seed)
  weights <- transport_exploratory_align_weights(predictions, plan)
  family_probability <- 1 - 0.05 / length(model_names)
  do.call(rbind, lapply(model_names, function(model) {
    gain <- (
      predictions$loss_base - predictions[[paste0("loss_", model)]]
    ) / log(2)
    bootstrap <- vapply(seq_len(ncol(weights)), function(draw) {
      transport_exploratory_weighted_mean(gain, weights[, draw])
    }, numeric(1))
    interval <- transport_exploratory_interval(bootstrap, 0.95)
    family <- transport_exploratory_interval(bootstrap, family_probability)
    data.frame(
      probe_type = predictions$probe_type[[1L]],
      model = model,
      trials = nrow(predictions),
      information_gain_bits = mean(gain),
      lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
      lower_family = family[["lower"]],
      upper_family = family[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
}

item_effect_write <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  write(result$score_summary$methods, "score-prediction.csv")
  write(result$score_summary$contrasts, "score-contrasts.csv")
  write(result$transfer_summary, "study-transfer.csv")
  write(result$item_correlations, "item-correlations.csv")
  write(result$confidence_summary, "newness-prediction.csv")
  write(result$leakage_audit, "leakage-audit.csv")
  write(result$configuration, "configuration.csv")
  saveRDS(
    result[c(
      "score_summary", "transfer_summary", "item_correlations",
      "confidence_summary", "leakage_audit", "configuration"
    )], file.path(output_dir, "scientific-results.rds"), version = 3
  )
  saveRDS(
    result$local_predictions,
    file.path(output_dir, "local-predictions.rds"), version = 3
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[basename(files) != "manifest-md5.csv"]
  write(data.frame(
    file = basename(files), md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  ), "manifest-md5.csv")
  invisible(result)
}

run_gaze_weave_item_effects <- function(
    output_dir = item_effect_result_dir,
    draws = item_effect_draws,
    seed = item_effect_seed) {
  panel <- transport_density_read_panel()
  retrieval <- transport_density_wide(panel)
  retrieval$probe_type_ec <- ifelse(retrieval$probe_type == "old", 0.5, -0.5)
  methods <- names(transport_density_methods)
  predictions <- lapply(methods, function(method) {
    message("Cross-participant item score: ", method)
    item_effect_crossparticipant(retrieval, paste0("score_", method))
  })
  names(predictions) <- methods
  predictions$transport_specific <- item_effect_crossparticipant(
    retrieval, "score_transport", density = "score_density_sigma_80"
  )
  score_summary <- item_effect_score_summary(predictions, draws, seed)

  study <- item_effect_read_study(retrieval_reference = retrieval)
  study_prediction <- item_effect_crossparticipant(study, "score_study")
  study_transfer <- item_effect_study_transfer(retrieval, study)
  transfer_summary <- rbind(
    transform(
      item_effect_score_summary(
        list(transport = study_prediction), draws, seed
      )$methods,
      source = "study_repeated_viewing",
      target = "study_repeated_viewing"
    )[c("source", "target", "trials", "cv_r2", "lower_95", "upper_95",
        "pearson", "spearman")],
    item_effect_transfer_summary(study_transfer, draws, seed)
  )

  retrieval_items <- item_effect_full_item_table(
    retrieval, "score_transport", label = "retrieval_transport"
  )
  density_items <- item_effect_full_item_table(
    retrieval, "score_density_sigma_80", label = "retrieval_density"
  )
  specific_items <- item_effect_full_item_table(
    retrieval, "score_transport", density = "score_density_sigma_80",
    label = "retrieval_transport_specific"
  )
  study_items <- item_effect_full_item_table(
    study, "score_study", label = "study_transport"
  )
  item_table <- Reduce(function(x, y) merge(x, y, by = "item"), list(
    retrieval_items[c("item", "retrieval_transport")],
    density_items[c("item", "retrieval_density")],
    specific_items[c("item", "retrieval_transport_specific")],
    study_items[c("item", "study_transport")]
  ))
  item_correlations <- as.data.frame(as.table(stats::cor(
    item_table[-1L], method = "spearman"
  )), stringsAsFactors = FALSE)
  names(item_correlations) <- c("first", "second", "spearman")

  confidence_predictions <- lapply(c("old", "lure"), function(probe) {
    item_effect_confidence_crossfit(retrieval, study, probe)
  })
  names(confidence_predictions) <- c("old", "lure")
  confidence_summary <- do.call(rbind, lapply(confidence_predictions, function(x) {
    item_effect_confidence_summary(x, draws, seed)
  }))
  leakage_audit <- data.frame(
    artifact = c(
      names(predictions), "study_within", "study_to_retrieval",
      "newness_old", "newness_lure"
    ),
    rows = c(
      vapply(predictions, nrow, integer(1)), nrow(study_prediction),
      nrow(study_transfer), vapply(confidence_predictions, nrow, integer(1))
    ),
    self_participant_overlap = 0L,
    stringsAsFactors = FALSE
  )
  result <- list(
    score_summary = score_summary,
    transfer_summary = transfer_summary,
    item_correlations = item_correlations,
    confidence_summary = confidence_summary,
    leakage_audit = leakage_audit,
    configuration = data.frame(
      status = "known_item_new_participant_exploratory",
      seed = seed, bootstrap_draws = draws,
      retrieval_trials = nrow(retrieval),
      study_trials = nrow(study),
      participants = length(unique(retrieval$participant)),
      items = length(unique(retrieval$item)),
      response_scale = "newness_1_to_4",
      frozen_item_disjoint_primary_unchanged = TRUE,
      local_only = TRUE,
      stringsAsFactors = FALSE
    ),
    local_predictions = list(
      retrieval = predictions,
      study = study_prediction,
      study_transfer = study_transfer,
      confidence = confidence_predictions,
      item_table = item_table
    )
  )
  item_effect_write(result, output_dir)
}
