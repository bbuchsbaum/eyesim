# Cross-fitted old-response prediction for the pcmri GazeWeave comparison.
#
# This is the confirmatory follow-up to the exploratory mixed-model coefficient
# reported by gaze-weave-pcmri-catalog-replication.R. It uses already held-out
# gaze scores, then cross-fits the behavioral association across disjoint
# participant and item folds. The primary value is improvement over a
# saliency-only model in held-out log loss, measured in bits per trial.

pcmri_behavior_protocol_version <- 1L
pcmri_behavior_seed <- 20260821L
pcmri_behavior_draws <- 2000L
pcmri_behavior_methods <- c(
  "replay", "transport_v2", "density_sigma_80_raw_calibrated",
  "density_ridge_registered", "multimatch_ridge_registered",
  "multimatch_mm_position_raw_calibrated"
)
pcmri_behavior_output_dir <- file.path(
  "inst", "validation", "gaze-weave-pcmri-catalog-replication-results"
)

pcmri_behavior_catalog_file <- system.file(
  "validation", "gaze-weave-pcmri-catalog-replication.R", package = "eyesim"
)
pcmri_behavior_this_file <- tryCatch(
  normalizePath(sys.frame(1L)$ofile, mustWork = TRUE),
  error = function(error) ""
)
if (!nzchar(pcmri_behavior_catalog_file) &&
    nzchar(pcmri_behavior_this_file)) {
  pcmri_behavior_catalog_file <- file.path(
    dirname(pcmri_behavior_this_file),
    "gaze-weave-pcmri-catalog-replication.R"
  )
}
if (!nzchar(pcmri_behavior_catalog_file)) {
  pcmri_behavior_catalog_file <- file.path(
    "inst", "validation", "gaze-weave-pcmri-catalog-replication.R"
  )
}
if (!file.exists(pcmri_behavior_catalog_file)) {
  stop("Run the behavioral prediction court from the eyesim source root.")
}
source(pcmri_behavior_catalog_file, local = TRUE)

pcmri_behavior_log_loss <- function(response, probability) {
  probability <- pmin(pmax(as.numeric(probability), 1e-8), 1 - 1e-8)
  -(response * log(probability) + (1 - response) * log1p(-probability))
}

pcmri_behavior_fit_fold <- function(train, evaluate, method, fold_id) {
  if (length(intersect(train$participant, evaluate$participant)) > 0L ||
      length(intersect(train$item, evaluate$item)) > 0L) {
    stop("Behavioral training and evaluation support overlaps.")
  }
  center <- mean(train$gaze_info_bits)
  scale <- stats::sd(train$gaze_info_bits)
  if (!is.finite(scale) || scale <= 0) {
    stop("Training gaze information has zero or invalid scale.")
  }
  train$z_metric <- (train$gaze_info_bits - center) / scale
  evaluate$z_metric <- (evaluate$gaze_info_bits - center) / scale
  null <- stats::glm(
    said_old ~ sal10, family = stats::binomial(), data = train
  )
  full <- stats::glm(
    said_old ~ sal10 + z_metric,
    family = stats::binomial(), data = train
  )
  probability_null <- stats::predict(
    null, newdata = evaluate, type = "response"
  )
  probability_full <- stats::predict(
    full, newdata = evaluate, type = "response"
  )
  loss_null <- pcmri_behavior_log_loss(
    evaluate$said_old, probability_null
  )
  loss_full <- pcmri_behavior_log_loss(
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
    probability_null = probability_null,
    probability_full = probability_full,
    loss_null = loss_null,
    loss_full = loss_full,
    behavior_info_bits = (loss_null - loss_full) / log(2),
    train_n = nrow(train),
    eval_n = nrow(evaluate),
    train_participants = length(unique(train$participant)),
    eval_participants = length(unique(evaluate$participant)),
    train_items = length(unique(train$item)),
    eval_items = length(unique(evaluate$item)),
    coefficient = unname(stats::coef(full)[["z_metric"]]),
    stringsAsFactors = FALSE
  )
}

pcmri_behavior_crossfit <- function(scores, methods = pcmri_behavior_methods) {
  scores <- scores[
    scores$method %in% methods & scores$probe_type == "old" &
      !is.na(scores$said_old) & is.finite(scores$gaze_info_bits),
    , drop = FALSE
  ]
  rows <- list()
  index <- 1L
  for (method in methods) {
    method_scores <- scores[scores$method == method, , drop = FALSE]
    for (fold_id in sort(unique(method_scores$outer_fold))) {
      evaluate <- method_scores[
        method_scores$outer_fold == fold_id, , drop = FALSE
      ]
      eval_participants <- unique(evaluate$participant)
      eval_items <- unique(evaluate$item)
      train <- method_scores[
        !method_scores$participant %in% eval_participants &
          !method_scores$item %in% eval_items,
        , drop = FALSE
      ]
      if (!nrow(train) || !nrow(evaluate)) {
        stop("A behavioral cross-fit fold has empty support.")
      }
      rows[[index]] <- pcmri_behavior_fit_fold(
        train, evaluate, method, fold_id
      )
      index <- index + 1L
    }
  }
  result <- dplyr::bind_rows(rows)
  expected <- NULL
  for (method in methods) {
    keys <- sort(result$trial_key[result$method == method])
    if (is.null(expected)) expected <- keys
    if (!identical(keys, expected) || anyDuplicated(keys)) {
      stop("Behavioral methods do not share unique evaluation support.")
    }
  }
  result
}

pcmri_behavior_weighted_mean <- function(value, weights) {
  numerator <- colSums(weights * value)
  denominator <- colSums(weights)
  result <- numerator / denominator
  result[!is.finite(result)] <- NA_real_
  result
}

pcmri_behavior_interval <- function(value, probability = 0.95) {
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

pcmri_behavior_bootstrap_p <- function(value) {
  value <- value[is.finite(value)]
  if (!length(value)) return(NA_real_)
  min(1, 2 * min(mean(value <= 0), mean(value >= 0)))
}

pcmri_behavior_prediction_summary <- function(
    predictions, draws = pcmri_behavior_draws,
    seed = pcmri_behavior_seed) {
  reference <- predictions[
    predictions$method == pcmri_behavior_methods[[1L]], , drop = FALSE
  ]
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws = draws, seed = seed)
  method_draws <- list()
  method_values <- list()
  rows <- lapply(pcmri_behavior_methods, function(method) {
    part <- predictions[predictions$method == method, , drop = FALSE]
    part <- part[match(plan$trial_keys, part$trial_key), , drop = FALSE]
    if (anyNA(part$trial_key)) stop("Could not align behavioral predictions.")
    values <- pcmri_behavior_weighted_mean(
      part$behavior_info_bits, plan$weights
    )
    method_draws[[method]] <<- values
    method_values[[method]] <<- part$behavior_info_bits
    interval <- pcmri_behavior_interval(values)
    data.frame(
      method = method,
      n = nrow(part),
      mean_behavior_info_bits = mean(part$behavior_info_bits),
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      bootstrap_p = pcmri_behavior_bootstrap_p(values),
      full_log_loss = mean(part$loss_full),
      null_log_loss = mean(part$loss_null),
      mean_coefficient = mean(part$coefficient),
      stringsAsFactors = FALSE
    )
  })
  comparison_n <- length(pcmri_behavior_methods) - 1L
  family_probability <- 1 - 0.05 / comparison_n
  comparisons <- lapply(setdiff(pcmri_behavior_methods, "replay"), function(
      comparator) {
    difference <- method_draws[["replay"]] - method_draws[[comparator]]
    ordinary <- pcmri_behavior_interval(difference)
    family <- pcmri_behavior_interval(difference, family_probability)
    observed <- mean(
      method_values[["replay"]] - method_values[[comparator]]
    )
    data.frame(
      contrast = paste("replay_minus", comparator, sep = "_"),
      estimate = observed,
      lower_95 = ordinary[["lower"]],
      upper_95 = ordinary[["upper"]],
      lower_family = family[["lower"]],
      upper_family = family[["upper"]],
      bootstrap_p = pcmri_behavior_bootstrap_p(difference),
      valid_draws = sum(is.finite(difference)),
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

pcmri_behavior_association_bootstrap <- function(
    scores, draws = pcmri_behavior_draws,
    seed = pcmri_behavior_seed + 10000L) {
  reference <- scores[
    scores$method == pcmri_behavior_methods[[1L]] &
      scores$probe_type == "old" & !is.na(scores$said_old) &
      is.finite(scores$gaze_info_bits), , drop = FALSE
  ]
  reference$trial_key <- paste(reference$participant, reference$item, sep = ":")
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws = draws, seed = seed)
  rows <- lapply(pcmri_behavior_methods, function(method) {
    part <- scores[
      scores$method == method & scores$probe_type == "old" &
        !is.na(scores$said_old) & is.finite(scores$gaze_info_bits),
      , drop = FALSE
    ]
    part$trial_key <- paste(part$participant, part$item, sep = ":")
    part <- part[match(plan$trial_keys, part$trial_key), , drop = FALSE]
    part$z_metric <- as.numeric(scale(part$gaze_info_bits))
    design <- stats::model.matrix(~ z_metric + sal10, data = part)
    response <- part$said_old
    estimate_fit <- stats::glm.fit(
      design, response, family = stats::binomial()
    )
    values <- vapply(seq_len(ncol(plan$weights)), function(draw) {
      fit <- suppressWarnings(tryCatch(
        stats::glm.fit(
          design, response, weights = plan$weights[, draw],
          family = stats::binomial()
        ),
        error = function(error) NULL
      ))
      if (is.null(fit) || !fit$converged) return(NA_real_)
      unname(fit$coefficients[["z_metric"]])
    }, numeric(1))
    interval <- pcmri_behavior_interval(values)
    data.frame(
      method = method,
      estimate = unname(estimate_fit$coefficients[["z_metric"]]),
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      bootstrap_p = pcmri_behavior_bootstrap_p(values),
      valid_draws = sum(is.finite(values)),
      valid_fraction = mean(is.finite(values)),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

pcmri_behavior_run <- function(
    output_dir = pcmri_behavior_output_dir,
    draws = pcmri_behavior_draws, seed = pcmri_behavior_seed) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  scores <- pcmri_catalog_common_support()
  predictions <- pcmri_behavior_crossfit(scores)
  predictive <- pcmri_behavior_prediction_summary(
    predictions, draws = draws, seed = seed
  )
  association <- pcmri_behavior_association_bootstrap(
    scores, draws = draws, seed = seed + 10000L
  )
  utils::write.csv(
    predictions,
    file.path(output_dir, "behavior-predictions.csv"), row.names = FALSE
  )
  utils::write.csv(
    predictive$methods,
    file.path(output_dir, "behavior-predictive-summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    predictive$comparisons,
    file.path(output_dir, "behavior-predictive-comparisons.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    association,
    file.path(output_dir, "behavior-association-bootstrap.csv"),
    row.names = FALSE
  )
  invisible(list(
    predictions = predictions,
    predictive = predictive,
    association = association
  ))
}
