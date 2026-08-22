# Metric-natural-support recovery gate for the pcmri recognition data.
#
# Protocol: inst/validation/GAZEWEAVE-PCMRI-SIGNAL-RECOVERY.md
# Run from the eyesim source root after devtools::load_all(). Raw scores,
# checkpoints, predictions, and bootstrap draws remain in a Git-ignored folder.

expanded_support_protocol_version <- 1L
expanded_support_seed <- 20260827L
expanded_support_draws <- 2000L
expanded_support_methods <- c(
  "replay", "density_sigma_80_raw_calibrated"
)
expanded_support_expected <- c(
  raw_participants = 46L,
  eligible_participants = 46L,
  eligible_pairs = 2225L,
  retained_participants = 45L,
  retained_pairs = 2222L,
  items = 120L
)
expanded_support_output_dir <- file.path(
  "inst", "validation", "gaze-weave-pcmri-catalog-replication-results",
  "expanded-support"
)

expanded_support_behavior_file <- system.file(
  "validation", "gaze-weave-pcmri-behavior-prediction.R", package = "eyesim"
)
if (!nzchar(expanded_support_behavior_file)) {
  expanded_support_behavior_file <- file.path(
    "inst", "validation", "gaze-weave-pcmri-behavior-prediction.R"
  )
}
if (!file.exists(expanded_support_behavior_file)) {
  stop("Run the expanded-support gate from the eyesim source root.")
}
source(expanded_support_behavior_file, local = TRUE)

expanded_support_eligible_pairs <- function(raw) {
  study <- raw$study[
    raw$study$presentation == 4L & raw$study$nfix >= 1L,
    c("participant", "item"), drop = FALSE
  ]
  study <- unique(study)
  retrieval <- raw$retrieval[
    raw$retrieval$probe_type %in% c("old", "lure") &
      raw$retrieval$nfix >= 1L,
    c(
      "participant", "item", "probe_type", "degradation", "Accuracy"
    ), drop = FALSE
  ]
  names(retrieval)[names(retrieval) == "Accuracy"] <- "accuracy"
  retrieval_key <- paste(retrieval$participant, retrieval$item, sep = ":")
  if (anyDuplicated(retrieval_key)) {
    stop("Combined retrieval paths are not unique by participant and item.")
  }
  pairs <- merge(
    study, retrieval, by = c("participant", "item"),
    all = FALSE, sort = FALSE
  )
  pairs[order(pairs$participant, pairs$item), , drop = FALSE]
}

expanded_support_select_cohort <- function(
    raw, item_seed = full_recognition_item_seed,
    candidate_n = full_recognition_candidate_n) {
  eligible <- expanded_support_eligible_pairs(raw)
  set.seed(item_seed)
  item_order <- sample(sort(unique(eligible$item)))
  item_map <- data.frame(
    item = item_order,
    item_fold = rep(1:2, length.out = length(item_order)),
    stringsAsFactors = FALSE
  )
  eligible <- merge(eligible, item_map, by = "item", sort = FALSE)
  counts <- table(
    factor(
      eligible$participant,
      levels = sort(unique(eligible$participant))
    ),
    factor(eligible$item_fold, levels = 1:2)
  )
  keep_participant <- rownames(counts)[
    counts[, 1L] >= candidate_n & counts[, 2L] >= candidate_n
  ]
  pairs <- eligible[
    eligible$participant %in% keep_participant, , drop = FALSE
  ]
  pairs <- pairs[order(pairs$participant, pairs$item_fold, pairs$item), ]
  rownames(pairs) <- NULL
  list(
    participants = sort(unique(pairs$participant)),
    pairs = pairs,
    item_map = item_map[order(item_map$item), ],
    eligible_pairs = eligible,
    raw_participant_n = length(unique(raw$trials$study$Subject)),
    eligible_participant_n = length(unique(eligible$participant)),
    eligible_pair_n = nrow(eligible),
    candidate_n = candidate_n,
    balance = list(
      probe_type = table(pairs$probe_type),
      degradation = table(pairs$degradation),
      accuracy = table(pairs$accuracy)
    )
  )
}

expanded_support_validate_design <- function(cohort, candidates, folds,
                                             strict = TRUE) {
  design <- full_recognition_validate_design(
    cohort, candidates, folds, strict = FALSE
  )
  if (strict && !identical(
    as.integer(design$observed), as.integer(expanded_support_expected)
  )) {
    stop("Expanded support differs from the frozen protocol.")
  }
  if (strict && (!identical(
    as.integer(cohort$balance$probe_type), c(1120L, 1102L)
  ) || !identical(
    as.integer(cohort$balance$degradation),
    c(439L, 451L, 446L, 440L, 446L)
  ))) {
    stop("Expanded condition balance differs from the frozen protocol.")
  }
  design
}

expanded_support_specs <- function(smoke = FALSE) {
  specs <- real_specs(smoke)
  specs$baseline <- gaze_baseline_spec(
    screen = specs$screen,
    density_sigmas = 80,
    density_grid = if (smoke) 12L else 24L,
    methods = "density",
    warp = specs$warp,
    lambda_grid = c(0.01, 0.1, 1, 10),
    inner_folds = 2L
  )
  specs
}

expanded_support_checkpoint <- function(output_dir, method, fold_id) {
  file.path(
    output_dir,
    sprintf("checkpoint-%s-fold-%02d.rds", method, fold_id)
  )
}

expanded_support_attach_behavior <- function(scored, raw) {
  response <- unique(raw$trials$retrieval[c(
    "Subject", "ImageNumber", "Response"
  )])
  names(response) <- c("participant", "item", "response")
  response$participant <- as.character(response$participant)
  response$item <- as.integer(response$item)
  result <- merge(
    scored, response, by = c("participant", "item"),
    all.x = TRUE, sort = FALSE
  )
  result$said_old <- pcmri_catalog_response(result$response)
  result$sal10 <- (as.numeric(result$degradation) - 60) / 10
  result[order(result$method, result$participant, result$item), ]
}

expanded_support_metric_summary <- function(
    scores, draws = expanded_support_draws,
    seed = expanded_support_seed) {
  methods <- expanded_support_methods
  reference <- scores[
    scores$method == methods[[1L]] & is.finite(scores$gaze_info_bits),
    , drop = FALSE
  ]
  reference$trial_key <- paste(
    reference$participant, reference$item, sep = ":"
  )
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws = draws, seed = seed)
  method_draws <- list()
  rows <- lapply(methods, function(method) {
    part <- scores[
      scores$method == method & is.finite(scores$gaze_info_bits),
      , drop = FALSE
    ]
    part$trial_key <- paste(part$participant, part$item, sep = ":")
    part <- part[match(plan$trial_keys, part$trial_key), , drop = FALSE]
    if (anyNA(part$trial_key) || anyDuplicated(part$trial_key)) {
      stop("Expanded methods do not share unique trial support.")
    }
    values <- pcmri_behavior_weighted_mean(
      part$gaze_info_bits, plan$weights
    )
    method_draws[[method]] <<- values
    interval <- pcmri_behavior_interval(values)
    data.frame(
      method = method,
      n = nrow(part),
      participants = length(unique(part$participant)),
      mean_gaze_info_bits = mean(part$gaze_info_bits),
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      bootstrap_p = pcmri_behavior_bootstrap_p(values),
      mean_rank = mean(part$template_rank),
      top1_credit = mean(part$top1_credit),
      mean_log_loss = mean(part$log_loss),
      stringsAsFactors = FALSE
    )
  })
  difference <- method_draws[[methods[[1L]]]] -
    method_draws[[methods[[2L]]]]
  interval <- pcmri_behavior_interval(difference)
  comparison <- data.frame(
    contrast = "replay_minus_density_sigma_80",
    estimate = mean(
      scores$gaze_info_bits[scores$method == methods[[1L]]] -
        scores$gaze_info_bits[scores$method == methods[[2L]]]
    ),
    lower_95 = interval[["lower"]],
    upper_95 = interval[["upper"]],
    bootstrap_p = pcmri_behavior_bootstrap_p(difference),
    stringsAsFactors = FALSE
  )
  list(
    methods = dplyr::bind_rows(rows),
    comparison = comparison,
    draws = method_draws,
    plan = plan
  )
}

expanded_support_behavior_summary <- function(
    predictions, draws = expanded_support_draws,
    seed = expanded_support_seed + 1000L) {
  methods <- expanded_support_methods
  reference <- predictions[
    predictions$method == methods[[1L]], , drop = FALSE
  ]
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws = draws, seed = seed)
  method_draws <- list()
  method_values <- list()
  rows <- lapply(methods, function(method) {
    part <- predictions[predictions$method == method, , drop = FALSE]
    part <- part[match(plan$trial_keys, part$trial_key), , drop = FALSE]
    if (anyNA(part$trial_key)) {
      stop("Could not align expanded behavioral predictions.")
    }
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
  difference <- method_draws[[methods[[1L]]]] -
    method_draws[[methods[[2L]]]]
  interval <- pcmri_behavior_interval(difference)
  comparison <- data.frame(
    contrast = "replay_minus_density_sigma_80",
    estimate = mean(
      method_values[[methods[[1L]]]] - method_values[[methods[[2L]]]]
    ),
    lower_95 = interval[["lower"]],
    upper_95 = interval[["upper"]],
    bootstrap_p = pcmri_behavior_bootstrap_p(difference),
    stringsAsFactors = FALSE
  )
  list(
    methods = dplyr::bind_rows(rows),
    comparison = comparison,
    draws = method_draws,
    plan = plan
  )
}

expanded_support_association_bootstrap <- function(
    scores, draws = expanded_support_draws,
    seed = expanded_support_seed + 2000L) {
  methods <- expanded_support_methods
  reference <- scores[
    scores$method == methods[[1L]] & scores$probe_type == "old" &
      !is.na(scores$said_old) & is.finite(scores$gaze_info_bits),
    , drop = FALSE
  ]
  reference$trial_key <- paste(
    reference$participant, reference$item, sep = ":"
  )
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws = draws, seed = seed)
  rows <- lapply(methods, function(method) {
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

expanded_support_candidate_scores <- function(scored) {
  full_recognition_candidate_scores(scored)
}

expanded_support_configuration <- function(cohort, folds, draws, seed) {
  data.frame(
    protocol_version = expanded_support_protocol_version,
    seed = seed,
    item_seed = full_recognition_item_seed,
    bootstrap_draws = draws,
    raw_participants = cohort$raw_participant_n,
    eligible_participants = cohort$eligible_participant_n,
    eligible_pairs = cohort$eligible_pair_n,
    retained_participants = length(cohort$participants),
    retained_pairs = nrow(cohort$pairs),
    old_pairs = unname(cohort$balance$probe_type[["old"]]),
    lure_pairs = unname(cohort$balance$probe_type[["lure"]]),
    items = nrow(cohort$item_map),
    candidate_n = full_recognition_candidate_n,
    outer_folds = length(folds$folds),
    study_presentation = 4L,
    window = "[0,3000)",
    minimum_path_fixations = 1L,
    methods = paste(expanded_support_methods, collapse = ";"),
    local_only = TRUE,
    stringsAsFactors = FALSE
  )
}

expanded_support_write <- function(result, output_dir) {
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  score_columns <- setdiff(names(result$scores), "candidates")
  write(result$scores[score_columns], "trial-scores.csv")
  write(
    expanded_support_candidate_scores(result$scores),
    "candidate-scores.csv"
  )
  write(result$metric$methods, "metric-summary.csv")
  write(result$metric$comparison, "metric-comparison.csv")
  write(result$behavior$methods, "behavior-predictive-summary.csv")
  write(result$behavior$comparison, "behavior-predictive-comparison.csv")
  write(result$association, "behavior-association-bootstrap.csv")
  write(result$predictions, "behavior-predictions.csv")
  write(result$config$design$fold_audit, "design-fold-audit.csv")
  write(result$config$folds$participant_map, "participant-folds.csv")
  write(result$config$folds$item_map, "item-folds.csv")
  write(result$config$candidates, "candidate-plan.csv")
  write(result$resources, "resources.csv")
  write(result$warps, "warp-audit.csv")
  write(result$fold_audit, "scoring-fold-audit.csv")
  write(result$config$configuration, "configuration.csv")
  saveRDS(
    result[c("metric", "behavior", "association", "config")],
    file.path(output_dir, "scientific-results.rds"), version = 3
  )
  utils::capture.output(
    utils::sessionInfo(), file = file.path(output_dir, "session-info.txt")
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[basename(files) != "manifest-md5.csv"]
  manifest <- data.frame(
    file = basename(files),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  write(manifest, "manifest-md5.csv")
  invisible(result)
}

run_gaze_weave_pcmri_expanded_support <- function(
    output_dir = expanded_support_output_dir,
    seed = expanded_support_seed,
    draws = expanded_support_draws,
    smoke = FALSE,
    method_groups = c("replay", "baselines"),
    fold_ids = NULL,
    workers = NULL,
    resume = TRUE,
    finalize = !smoke) {
  if (!smoke && (seed != expanded_support_seed ||
                 draws != expanded_support_draws)) {
    stop("Non-smoke expanded-support runs use the frozen seed and draws.")
  }
  method_groups <- match.arg(
    method_groups, c("replay", "baselines"), several.ok = TRUE
  )
  raw <- full_recognition_read_inputs(verify = !smoke)
  cohort <- expanded_support_select_cohort(raw)
  if (smoke) cohort <- full_recognition_smoke_cohort(cohort)
  candidates <- full_recognition_candidate_plan(cohort, seed)
  folds <- full_recognition_fold_plan(cohort, seed)
  design <- expanded_support_validate_design(
    cohort, candidates, folds, strict = !smoke
  )
  tables <- full_recognition_task_tables(raw, cohort)
  specs <- expanded_support_specs(smoke)
  if (is.null(workers)) {
    detected <- parallel::detectCores(logical = FALSE)
    if (!is.finite(detected) || detected < 1L) detected <- 1L
    workers <- if (smoke) min(2L, detected) else min(4L, detected)
  }
  workers <- as.integer(workers)
  selected_folds <- folds$folds
  if (!is.null(fold_ids)) {
    selected_folds <- selected_folds[vapply(
      selected_folds, function(fold) fold$id %in% fold_ids, logical(1)
    )]
  }
  if (!length(selected_folds)) stop("fold_ids selected no outer folds.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (fold in selected_folds) {
    for (method in method_groups) {
      path <- expanded_support_checkpoint(output_dir, method, fold$id)
      if (resume && file.exists(path)) {
        checkpoint <- readRDS(path)
        if (!identical(
          checkpoint$expanded_support_protocol_version,
          expanded_support_protocol_version
        ) || !identical(checkpoint$smoke, smoke) ||
            !identical(checkpoint$expanded_support_seed, seed)) {
          stop("A checkpoint uses another expanded protocol: ", path)
        }
        message("Using checkpoint: ", basename(path))
        next
      }
      message(
        "Expanded-support scoring: fold ", fold$id,
        ", method ", method
      )
      checkpoint <- full_recognition_score_fold_method(
        tables, candidates, fold, specs, method,
        workers = workers, seed = seed, smoke = smoke
      )
      checkpoint$expanded_support_protocol_version <-
        expanded_support_protocol_version
      checkpoint$expanded_support_seed <- seed
      saveRDS(checkpoint, path, version = 3)
    }
  }
  if (!finalize) {
    return(invisible(list(
      cohort = cohort, candidates = candidates, folds = folds,
      design = design
    )))
  }
  required <- expand.grid(
    method = c("replay", "baselines"), fold = 1:4,
    stringsAsFactors = FALSE
  )
  paths <- mapply(
    expanded_support_checkpoint,
    MoreArgs = list(output_dir = output_dir),
    method = required$method, fold_id = required$fold,
    USE.NAMES = FALSE
  )
  if (!all(file.exists(paths))) {
    stop("Expanded-support finalization requires all eight checkpoints.")
  }
  checkpoints <- lapply(paths, readRDS)
  scored <- dplyr::bind_rows(lapply(checkpoints, `[[`, "scored"))
  scores <- scored[scored$method %in% expanded_support_methods, ]
  support <- table(scores$method)
  if (!all(expanded_support_methods %in% names(support)) ||
      any(support[expanded_support_methods] != nrow(cohort$pairs)) ||
      any(scores$status != "scored") || any(!scores$calibrated) ||
      any(!is.finite(scores$gaze_info_bits)) ||
      any(scores$candidate_count != full_recognition_candidate_n)) {
    stop("Expanded methods failed the frozen score/support contract.")
  }
  scores <- expanded_support_attach_behavior(scores, raw)
  metric <- expanded_support_metric_summary(scores, draws, seed)
  predictions <- pcmri_behavior_crossfit(
    scores, methods = expanded_support_methods
  )
  behavior <- expanded_support_behavior_summary(
    predictions, draws, seed + 1000L
  )
  association <- expanded_support_association_bootstrap(
    scores, draws, seed + 2000L
  )
  configuration <- expanded_support_configuration(cohort, folds, draws, seed)
  result <- list(
    scores = scores,
    metric = metric,
    predictions = predictions,
    behavior = behavior,
    association = association,
    resources = dplyr::bind_rows(lapply(checkpoints, `[[`, "resources")),
    warps = dplyr::bind_rows(lapply(checkpoints, `[[`, "warp")),
    fold_audit = dplyr::bind_rows(lapply(checkpoints, `[[`, "audit")),
    config = list(
      seed = seed, draws = draws, smoke = smoke, workers = workers,
      cohort = cohort, candidates = candidates, folds = folds,
      design = design, configuration = configuration
    )
  )
  expanded_support_write(result, output_dir)
  result
}
