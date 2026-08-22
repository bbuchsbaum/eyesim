# Frozen behavioral court for the multi-presentation Replay score.
#
# Protocol: inst/validation/GAZEWEAVE-PCMRI-MULTI-REPLAY-BEHAVIOR.md
# Run from a source checkout after devtools::load_all(). Outputs remain local
# under the Git-ignored multi-replay results directory.

multi_replay_behavior_protocol_version <- 1L
multi_replay_behavior_seed <- 20260829L
multi_replay_behavior_draws <- 2000L
multi_replay_behavior_methods <- c(
  "replay_all4_shrunk_exhaustive",
  "replay_all4_unshrunk_exhaustive",
  "replay_p4_exhaustive",
  "density_sigma_80_all4_exhaustive"
)
multi_replay_behavior_primary <- "replay_all4_shrunk_exhaustive"
multi_replay_behavior_conditional <- "replay_beyond_density80"
multi_replay_behavior_expected <- c(
  trials = 989L,
  said_new = 97L,
  said_old = 892L
)

multi_replay_behavior_measurement_file <- system.file(
  "validation", "gaze-weave-pcmri-multi-replay.R", package = "eyesim"
)
if (!nzchar(multi_replay_behavior_measurement_file)) {
  multi_replay_behavior_measurement_file <- file.path(
    "inst", "validation", "gaze-weave-pcmri-multi-replay.R"
  )
}
if (!file.exists(multi_replay_behavior_measurement_file)) {
  stop("Run the behavioral court from the eyesim source root.")
}
source(multi_replay_behavior_measurement_file, local = TRUE)

multi_replay_behavior_response <- function(response) {
  response <- suppressWarnings(as.integer(response))
  ifelse(
    response == 0L | is.na(response), NA_integer_,
    as.integer(response %in% c(1L, 2L))
  )
}

multi_replay_behavior_verify_measurement <- function(output_dir) {
  manifest_path <- file.path(output_dir, "manifest-md5.csv")
  verdict_path <- file.path(output_dir, "measurement-verdict.csv")
  if (!file.exists(manifest_path) || !file.exists(verdict_path)) {
    stop("The frozen measurement court has not been finalized.")
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  paths <- file.path(output_dir, manifest$file)
  if (!all(file.exists(paths)) ||
      !identical(unname(tools::md5sum(paths)), manifest$md5)) {
    stop("The frozen measurement manifest failed verification.")
  }
  verdict <- utils::read.csv(verdict_path, stringsAsFactors = FALSE)
  if (nrow(verdict) != 1L || !isTRUE(verdict$measurement_gates_pass) ||
      !isTRUE(verdict$freeze_for_behavior) ||
      isTRUE(verdict$behavior_tested)) {
    stop("The measurement verdict does not authorize behavioral testing.")
  }
  invisible(verdict)
}

multi_replay_behavior_measurement_scores <- function(output_dir) {
  checkpoints <- unlist(lapply(1:4, function(fold_id) {
    lapply(c("p4", "all4", "density"), function(method_group) {
      path <- multi_replay_checkpoint(
        output_dir, method_group, fold_id, smoke = FALSE
      )
      if (!file.exists(path)) stop("Missing measurement checkpoint: ", path)
      readRDS(path)
    })
  }), recursive = FALSE)
  scored <- dplyr::bind_rows(lapply(checkpoints, `[[`, "scored"))
  scored <- scored[
    scored$method %in% multi_replay_behavior_methods, , drop = FALSE
  ]
  support <- table(scored$method)
  if (!all(multi_replay_behavior_methods %in% names(support)) ||
      any(support[multi_replay_behavior_methods] != 2055L)) {
    stop("Behavioral methods do not share frozen measurement support.")
  }
  scored[c(
    "participant", "item", "outer_fold", "method", "gaze_info_bits",
    "candidate_count", "effective_fixations", "total_duration"
  )]
}

multi_replay_behavior_metadata <- function(verify = TRUE) {
  raw <- full_recognition_read_inputs(verify = verify)
  cohort <- multi_replay_select_cohort(raw)
  response <- raw$retrieval[c(
    "participant", "item", "probe_type", "degradation", "Response"
  )]
  response <- unique(response)
  key <- paste(response$participant, response$item, sep = ":")
  if (anyDuplicated(key)) stop("Behavior metadata is not unique by trial.")
  metadata <- merge(
    cohort$pairs[c("participant", "item", "probe_type", "degradation")],
    response[c("participant", "item", "Response")],
    by = c("participant", "item"), all.x = TRUE, sort = FALSE
  )
  metadata$said_old <- multi_replay_behavior_response(metadata$Response)
  metadata$sal10 <- (as.numeric(metadata$degradation) - 60) / 10
  metadata
}

multi_replay_behavior_scores <- function(measurement_dir, verify = TRUE) {
  if (verify) multi_replay_behavior_verify_measurement(measurement_dir)
  scored <- multi_replay_behavior_measurement_scores(measurement_dir)
  metadata <- multi_replay_behavior_metadata(verify = verify)
  scored <- merge(
    scored, metadata,
    by = c("participant", "item"), all.x = TRUE, sort = FALSE
  )
  scored <- scored[
    scored$probe_type == "old" & !is.na(scored$said_old), , drop = FALSE
  ]
  support <- table(scored$method)
  if (any(support[multi_replay_behavior_methods] !=
          multi_replay_behavior_expected[["trials"]])) {
    stop("Usable old-item behavioral support differs across methods.")
  }
  primary <- scored[scored$method == multi_replay_behavior_primary, ]
  observed <- c(
    trials = nrow(primary),
    said_new = sum(primary$said_old == 0L),
    said_old = sum(primary$said_old == 1L)
  )
  if (!identical(as.integer(observed),
                 as.integer(multi_replay_behavior_expected))) {
    stop("Old-item response support differs from the frozen protocol.")
  }
  scored[order(scored$method, scored$participant, scored$item), ]
}

multi_replay_behavior_log_loss <- function(response, probability) {
  probability <- pmin(pmax(as.numeric(probability), 1e-8), 1 - 1e-8)
  -(response * log(probability) + (1 - response) * log1p(-probability))
}

multi_replay_behavior_standardize <- function(train, evaluate,
                                              source, target) {
  center <- mean(train[[source]])
  scale <- stats::sd(train[[source]])
  if (!is.finite(center) || !is.finite(scale) || scale <= 0) {
    stop("Behavioral training predictor has zero or invalid scale: ", source)
  }
  train[[target]] <- (train[[source]] - center) / scale
  evaluate[[target]] <- (evaluate[[source]] - center) / scale
  list(train = train, evaluate = evaluate, center = center, scale = scale)
}

multi_replay_behavior_prepare_fold <- function(train, evaluate,
                                               metric = "gaze_info_bits") {
  train$log_effective_fixations <- log(train$effective_fixations)
  evaluate$log_effective_fixations <- log(evaluate$effective_fixations)
  train$log_candidate_count <- log(train$candidate_count)
  evaluate$log_candidate_count <- log(evaluate$candidate_count)
  variables <- c(
    metric = metric,
    log_effective_fixations = "log_effective_fixations",
    total_duration = "total_duration",
    log_candidate_count = "log_candidate_count"
  )
  targets <- c(
    metric = "z_metric",
    log_effective_fixations = "z_log_effective_fixations",
    total_duration = "z_total_duration",
    log_candidate_count = "z_log_candidate_count"
  )
  scaling <- list()
  for (name in names(variables)) {
    transformed <- multi_replay_behavior_standardize(
      train, evaluate, variables[[name]], targets[[name]]
    )
    train <- transformed$train
    evaluate <- transformed$evaluate
    scaling[[name]] <- transformed[c("center", "scale")]
  }
  list(train = train, evaluate = evaluate, scaling = scaling)
}

multi_replay_behavior_base_formula <- stats::as.formula(
  "said_old ~ sal10 + z_log_effective_fixations + z_total_duration + z_log_candidate_count"
)
multi_replay_behavior_full_formula <- stats::update.formula(
  multi_replay_behavior_base_formula, ". ~ . + z_metric"
)

multi_replay_behavior_fit_fold <- function(train, evaluate, method, fold_id) {
  if (length(intersect(train$participant, evaluate$participant)) > 0L ||
      length(intersect(train$item, evaluate$item)) > 0L) {
    stop("Behavioral training and evaluation support overlaps.")
  }
  prepared <- multi_replay_behavior_prepare_fold(train, evaluate)
  train <- prepared$train
  evaluate <- prepared$evaluate
  base <- stats::glm(
    multi_replay_behavior_base_formula,
    family = stats::binomial(), data = train
  )
  full <- stats::glm(
    multi_replay_behavior_full_formula,
    family = stats::binomial(), data = train
  )
  probability_base <- stats::predict(
    base, newdata = evaluate, type = "response"
  )
  probability_full <- stats::predict(
    full, newdata = evaluate, type = "response"
  )
  loss_base <- multi_replay_behavior_log_loss(
    evaluate$said_old, probability_base
  )
  loss_full <- multi_replay_behavior_log_loss(
    evaluate$said_old, probability_full
  )
  data.frame(
    participant = as.character(evaluate$participant),
    item = as.integer(evaluate$item),
    trial_key = paste(evaluate$participant, evaluate$item, sep = ":"),
    method = method,
    outer_fold = as.integer(fold_id),
    said_old = as.integer(evaluate$said_old),
    sal10 = evaluate$sal10,
    gaze_info_bits = evaluate$gaze_info_bits,
    effective_fixations = evaluate$effective_fixations,
    total_duration = evaluate$total_duration,
    candidate_count = evaluate$candidate_count,
    probability_base = probability_base,
    probability_full = probability_full,
    loss_base = loss_base,
    loss_full = loss_full,
    behavior_info_bits = (loss_base - loss_full) / log(2),
    coefficient = unname(stats::coef(full)[["z_metric"]]),
    train_n = nrow(train),
    eval_n = nrow(evaluate),
    train_participants = length(unique(train$participant)),
    eval_participants = length(unique(evaluate$participant)),
    train_items = length(unique(train$item)),
    eval_items = length(unique(evaluate$item)),
    stringsAsFactors = FALSE
  )
}

multi_replay_behavior_crossfit <- function(scores) {
  rows <- list()
  index <- 1L
  for (method in multi_replay_behavior_methods) {
    method_scores <- scores[scores$method == method, , drop = FALSE]
    for (fold_id in sort(unique(method_scores$outer_fold))) {
      evaluate <- method_scores[
        method_scores$outer_fold == fold_id, , drop = FALSE
      ]
      train <- method_scores[
        !method_scores$participant %in% unique(evaluate$participant) &
          !method_scores$item %in% unique(evaluate$item),
        , drop = FALSE
      ]
      if (!nrow(train) || !nrow(evaluate)) {
        stop("A behavioral cross-fit fold has empty support.")
      }
      rows[[index]] <- multi_replay_behavior_fit_fold(
        train, evaluate, method, fold_id
      )
      index <- index + 1L
    }
  }
  result <- dplyr::bind_rows(rows)
  expected <- NULL
  for (method in multi_replay_behavior_methods) {
    keys <- sort(result$trial_key[result$method == method])
    if (is.null(expected)) expected <- keys
    if (!identical(keys, expected) || anyDuplicated(keys)) {
      stop("Behavioral methods do not share unique evaluation support.")
    }
  }
  result
}

multi_replay_behavior_conditional_data <- function(scores) {
  replay <- scores[
    scores$method == multi_replay_behavior_primary, , drop = FALSE
  ]
  density <- scores[
    scores$method == "density_sigma_80_all4_exhaustive",
    c("participant", "item", "gaze_info_bits"), drop = FALSE
  ]
  names(density)[names(density) == "gaze_info_bits"] <- "density_info_bits"
  result <- merge(
    replay, density, by = c("participant", "item"),
    all = FALSE, sort = FALSE
  )
  if (nrow(result) != multi_replay_behavior_expected[["trials"]]) {
    stop("Replay and density do not share conditional behavioral support.")
  }
  result
}

multi_replay_behavior_fit_conditional_fold <- function(train, evaluate,
                                                       fold_id) {
  if (length(intersect(train$participant, evaluate$participant)) > 0L ||
      length(intersect(train$item, evaluate$item)) > 0L) {
    stop("Conditional behavioral training and evaluation support overlaps.")
  }
  prepared <- multi_replay_behavior_prepare_fold(train, evaluate)
  train <- prepared$train
  evaluate <- prepared$evaluate
  density <- multi_replay_behavior_standardize(
    train, evaluate, "density_info_bits", "z_density"
  )
  train <- density$train
  evaluate <- density$evaluate
  base_formula <- stats::update.formula(
    multi_replay_behavior_base_formula, ". ~ . + z_density"
  )
  full_formula <- stats::update.formula(
    base_formula, ". ~ . + z_metric"
  )
  base <- stats::glm(base_formula, family = stats::binomial(), data = train)
  full <- stats::glm(full_formula, family = stats::binomial(), data = train)
  probability_base <- stats::predict(
    base, newdata = evaluate, type = "response"
  )
  probability_full <- stats::predict(
    full, newdata = evaluate, type = "response"
  )
  loss_base <- multi_replay_behavior_log_loss(
    evaluate$said_old, probability_base
  )
  loss_full <- multi_replay_behavior_log_loss(
    evaluate$said_old, probability_full
  )
  data.frame(
    participant = as.character(evaluate$participant),
    item = as.integer(evaluate$item),
    trial_key = paste(evaluate$participant, evaluate$item, sep = ":"),
    method = multi_replay_behavior_conditional,
    outer_fold = as.integer(fold_id),
    said_old = as.integer(evaluate$said_old),
    sal10 = evaluate$sal10,
    gaze_info_bits = evaluate$gaze_info_bits,
    effective_fixations = evaluate$effective_fixations,
    total_duration = evaluate$total_duration,
    candidate_count = evaluate$candidate_count,
    probability_base = probability_base,
    probability_full = probability_full,
    loss_base = loss_base,
    loss_full = loss_full,
    behavior_info_bits = (loss_base - loss_full) / log(2),
    coefficient = unname(stats::coef(full)[["z_metric"]]),
    train_n = nrow(train),
    eval_n = nrow(evaluate),
    train_participants = length(unique(train$participant)),
    eval_participants = length(unique(evaluate$participant)),
    train_items = length(unique(train$item)),
    eval_items = length(unique(evaluate$item)),
    stringsAsFactors = FALSE
  )
}

multi_replay_behavior_conditional_crossfit <- function(scores) {
  data <- multi_replay_behavior_conditional_data(scores)
  rows <- lapply(sort(unique(data$outer_fold)), function(fold_id) {
    evaluate <- data[data$outer_fold == fold_id, , drop = FALSE]
    train <- data[
      !data$participant %in% unique(evaluate$participant) &
        !data$item %in% unique(evaluate$item),
      , drop = FALSE
    ]
    multi_replay_behavior_fit_conditional_fold(train, evaluate, fold_id)
  })
  dplyr::bind_rows(rows)
}

multi_replay_behavior_weighted_mean <- function(value, weights) {
  numerator <- colSums(weights * value)
  denominator <- colSums(weights)
  result <- numerator / denominator
  result[!is.finite(result)] <- NA_real_
  result
}

multi_replay_behavior_interval <- function(value, probability = 0.95) {
  value <- value[is.finite(value)]
  if (!length(value)) return(c(lower = NA_real_, upper = NA_real_))
  alpha <- 1 - probability
  stats::setNames(
    unname(stats::quantile(
      value, c(alpha / 2, 1 - alpha / 2), type = 8
    )),
    c("lower", "upper")
  )
}

multi_replay_behavior_bootstrap_p <- function(value) {
  value <- value[is.finite(value)]
  if (!length(value)) return(NA_real_)
  min(1, 2 * min(mean(value <= 0), mean(value >= 0)))
}

multi_replay_behavior_prediction_summary <- function(
    predictions, draws = multi_replay_behavior_draws,
    seed = multi_replay_behavior_seed) {
  methods <- c(
    multi_replay_behavior_methods, multi_replay_behavior_conditional
  )
  reference <- predictions[
    predictions$method == multi_replay_behavior_primary, , drop = FALSE
  ]
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws = draws, seed = seed)
  method_draws <- list()
  method_values <- list()
  rows <- lapply(methods, function(method) {
    part <- predictions[predictions$method == method, , drop = FALSE]
    part <- part[match(plan$trial_keys, part$trial_key), , drop = FALSE]
    if (anyNA(part$trial_key)) stop("Could not align behavioral predictions.")
    values <- multi_replay_behavior_weighted_mean(
      part$behavior_info_bits, plan$weights
    )
    method_draws[[method]] <<- values
    method_values[[method]] <<- part$behavior_info_bits
    interval <- multi_replay_behavior_interval(values)
    data.frame(
      method = method,
      role = if (method == multi_replay_behavior_primary) {
        "primary"
      } else if (method == multi_replay_behavior_conditional) {
        "secondary_conditional"
      } else {
        "comparator"
      },
      trials = nrow(part),
      mean_behavior_info_bits = mean(part$behavior_info_bits),
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      bootstrap_p = multi_replay_behavior_bootstrap_p(values),
      mean_full_log_loss = mean(part$loss_full),
      mean_base_log_loss = mean(part$loss_base),
      mean_fold_coefficient = mean(part$coefficient),
      stringsAsFactors = FALSE
    )
  })
  comparison_methods <- setdiff(
    multi_replay_behavior_methods, multi_replay_behavior_primary
  )
  comparisons <- lapply(comparison_methods, function(comparator) {
    difference <- method_draws[[multi_replay_behavior_primary]] -
      method_draws[[comparator]]
    interval <- multi_replay_behavior_interval(difference)
    data.frame(
      contrast = paste0("primary_minus_", comparator),
      estimate = mean(
        method_values[[multi_replay_behavior_primary]] -
          method_values[[comparator]]
      ),
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      bootstrap_p = multi_replay_behavior_bootstrap_p(difference),
      stringsAsFactors = FALSE
    )
  })
  list(
    methods = dplyr::bind_rows(rows),
    comparisons = dplyr::bind_rows(comparisons),
    draws = method_draws,
    plan = plan
  )
}

multi_replay_behavior_association_data <- function(part) {
  part$log_effective_fixations <- log(part$effective_fixations)
  part$log_candidate_count <- log(part$candidate_count)
  variables <- c(
    gaze_info_bits = "z_metric",
    log_effective_fixations = "z_log_effective_fixations",
    total_duration = "z_total_duration",
    log_candidate_count = "z_log_candidate_count"
  )
  for (source in names(variables)) {
    part[[variables[[source]]]] <- as.numeric(scale(part[[source]]))
  }
  part
}

multi_replay_behavior_association_bootstrap <- function(
    scores, draws = multi_replay_behavior_draws,
    seed = multi_replay_behavior_seed + 10000L) {
  reference <- scores[
    scores$method == multi_replay_behavior_primary, , drop = FALSE
  ]
  reference$trial_key <- paste(reference$participant, reference$item, sep = ":")
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws = draws, seed = seed)
  rows <- list()
  index <- 1L
  for (method in multi_replay_behavior_methods) {
    part <- scores[scores$method == method, , drop = FALSE]
    part$trial_key <- paste(part$participant, part$item, sep = ":")
    part <- part[match(plan$trial_keys, part$trial_key), , drop = FALSE]
    part <- multi_replay_behavior_association_data(part)
    design <- stats::model.matrix(
      multi_replay_behavior_full_formula, data = part
    )
    response <- part$said_old
    point <- stats::glm.fit(
      design, response, family = stats::binomial()
    )
    coefficient_draws <- vapply(seq_len(ncol(plan$weights)), function(draw) {
      fit <- suppressWarnings(tryCatch(
        stats::glm.fit(
          design, response, weights = plan$weights[, draw],
          family = stats::binomial()
        ),
        error = function(error) NULL
      ))
      if (is.null(fit) || !fit$converged || fit$rank != ncol(design) ||
          any(!is.finite(fit$coefficients))) {
        return(rep(NA_real_, ncol(design)))
      }
      fit$coefficients
    }, numeric(ncol(design)))
    rownames(coefficient_draws) <- colnames(design)
    terms <- if (method == multi_replay_behavior_primary) {
      c(
        "z_metric", "sal10", "z_log_effective_fixations",
        "z_total_duration", "z_log_candidate_count"
      )
    } else {
      "z_metric"
    }
    for (term in terms) {
      values <- coefficient_draws[term, ]
      interval <- multi_replay_behavior_interval(values)
      rows[[index]] <- data.frame(
        method = method,
        term = term,
        estimate = unname(point$coefficients[[term]]),
        lower_95 = interval[["lower"]],
        upper_95 = interval[["upper"]],
        bootstrap_p = multi_replay_behavior_bootstrap_p(values),
        valid_draws = sum(is.finite(values)),
        valid_fraction = mean(is.finite(values)),
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  dplyr::bind_rows(rows)
}

multi_replay_behavior_glmer_audit <- function(scores) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(data.frame(
      status = "missing_lme4", singular = NA, converged = FALSE,
      term = "z_metric", estimate = NA_real_, standard_error = NA_real_,
      z_value = NA_real_, p_value = NA_real_, warnings = NA_character_,
      messages = NA_character_, stringsAsFactors = FALSE
    ))
  }
  part <- scores[
    scores$method == multi_replay_behavior_primary, , drop = FALSE
  ]
  part <- multi_replay_behavior_association_data(part)
  part$participant <- factor(part$participant)
  part$item <- factor(part$item)
  warnings <- character()
  messages <- character()
  fit <- tryCatch(
    withCallingHandlers(
      lme4::glmer(
        said_old ~ sal10 + z_log_effective_fixations + z_total_duration +
          z_log_candidate_count + z_metric +
          (1 | participant) + (1 | item),
        family = stats::binomial(), data = part,
        control = lme4::glmerControl(
          optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)
        )
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      },
      message = function(message) {
        messages <<- c(messages, conditionMessage(message))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(error) error
  )
  if (inherits(fit, "error")) {
    return(data.frame(
      status = "error", singular = NA, converged = FALSE,
      term = "z_metric", estimate = NA_real_, standard_error = NA_real_,
      z_value = NA_real_, p_value = NA_real_,
      warnings = paste(warnings, collapse = " | "),
      messages = conditionMessage(fit), stringsAsFactors = FALSE
    ))
  }
  coefficients <- summary(fit)$coefficients
  convergence_messages <- fit@optinfo$conv$lme4$messages
  converged <- is.null(convergence_messages) &&
    identical(fit@optinfo$conv$opt, 0L)
  data.frame(
    status = "scored",
    singular = lme4::isSingular(fit),
    converged = converged,
    term = "z_metric",
    estimate = coefficients["z_metric", "Estimate"],
    standard_error = coefficients["z_metric", "Std. Error"],
    z_value = coefficients["z_metric", "z value"],
    p_value = coefficients["z_metric", "Pr(>|z|)"],
    warnings = paste(warnings, collapse = " | "),
    messages = paste(c(messages, convergence_messages), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

multi_replay_behavior_write_results <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  write(result$predictions, "behavior-predictions.csv")
  write(result$predictive$methods, "predictive-summary.csv")
  write(result$predictive$comparisons, "predictive-comparisons.csv")
  write(result$association, "association-bootstrap.csv")
  write(result$glmer, "glmer-audit.csv")
  write(result$verdict, "behavior-verdict.csv")
  configuration <- data.frame(
    protocol_version = multi_replay_behavior_protocol_version,
    seed = multi_replay_behavior_seed,
    bootstrap_draws = multi_replay_behavior_draws,
    usable_old_trials = multi_replay_behavior_expected[["trials"]],
    primary_method = multi_replay_behavior_primary,
    base_model = paste(deparse(multi_replay_behavior_base_formula), collapse = ""),
    recognition_window_ms = "[0,3000)",
    measurement_frozen_before_behavior = TRUE,
    local_only = TRUE,
    stringsAsFactors = FALSE
  )
  write(configuration, "configuration.csv")
  saveRDS(
    list(
      predictive = result$predictive[c("methods", "comparisons")],
      association = result$association,
      glmer = result$glmer,
      verdict = result$verdict,
      configuration = configuration
    ),
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
  write(manifest, "manifest-md5.csv")
  invisible(result)
}

run_gaze_weave_pcmri_multi_replay_behavior <- function(
    measurement_dir = file.path(
      "inst", "validation",
      "gaze-weave-pcmri-catalog-replication-results", "multi-replay"
    ),
    output_dir = file.path(measurement_dir, "behavior"),
    draws = multi_replay_behavior_draws,
    seed = multi_replay_behavior_seed,
    verify = TRUE) {
  if (draws != multi_replay_behavior_draws ||
      seed != multi_replay_behavior_seed) {
    stop("The full behavioral court uses the frozen draws and seed.")
  }
  scores <- multi_replay_behavior_scores(
    measurement_dir, verify = verify
  )
  predictions <- multi_replay_behavior_crossfit(scores)
  conditional <- multi_replay_behavior_conditional_crossfit(scores)
  predictions <- dplyr::bind_rows(predictions, conditional)
  predictive <- multi_replay_behavior_prediction_summary(
    predictions, draws = draws, seed = seed
  )
  association <- multi_replay_behavior_association_bootstrap(
    scores, draws = draws, seed = seed + 10000L
  )
  glmer <- multi_replay_behavior_glmer_audit(scores)
  primary <- predictive$methods[
    predictive$methods$method == multi_replay_behavior_primary,
    , drop = FALSE
  ]
  conditional_row <- predictive$methods[
    predictive$methods$method == multi_replay_behavior_conditional,
    , drop = FALSE
  ]
  verdict <- data.frame(
    primary_predictive_gain = primary$mean_behavior_info_bits,
    primary_lower_95 = primary$lower_95,
    primary_upper_95 = primary$upper_95,
    primary_supported = primary$lower_95 > 0,
    replay_beyond_density_gain = conditional_row$mean_behavior_info_bits,
    replay_beyond_density_lower_95 = conditional_row$lower_95,
    replay_beyond_density_upper_95 = conditional_row$upper_95,
    replay_beyond_density_supported = conditional_row$lower_95 > 0,
    measurement_was_frozen = TRUE,
    exploratory_search_performed = FALSE,
    stringsAsFactors = FALSE
  )
  result <- list(
    scores = scores,
    predictions = predictions,
    predictive = predictive,
    association = association,
    glmer = glmer,
    verdict = verdict,
    config = list(draws = draws, seed = seed)
  )
  multi_replay_behavior_write_results(result, output_dir)
  result
}
