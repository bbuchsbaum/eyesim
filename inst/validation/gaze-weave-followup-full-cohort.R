# Full-cohort study-repeat and phase-conditional follow-up for GazeWeave.
#
# Protocol: inst/validation/GAZEWEAVE-FOLLOWUP-FULL-COHORT.md
# Run from the eyesim source root after devtools::load_all().

followup_protocol_version <- 1L
followup_seed <- 20260826L
followup_bootstrap_draws <- 2000L
followup_study_contrasts <- data.frame(
  task = c("study_p1_p2", "study_p1_p3", "study_p1_p4"),
  reference_presentation = 1L,
  source_presentation = 2:4,
  stringsAsFactors = FALSE
)
followup_phases <- c("probe", "delay", "early_delay", "late_delay")
followup_mandatory_methods <- c(
  "transport_v2", "replay",
  "density_sigma_80_raw_calibrated",
  "density_sigma_160_raw_calibrated",
  "density_ridge_registered"
)

followup_full_file <- system.file(
  "validation", "gaze-weave-recognition-full-cohort.R", package = "eyesim"
)
if (!nzchar(followup_full_file)) {
  followup_full_file <- file.path(
    "inst", "validation", "gaze-weave-recognition-full-cohort.R"
  )
}
if (!file.exists(followup_full_file)) {
  stop("Run the follow-up court from the eyesim source root.")
}
source(followup_full_file, local = TRUE)

followup_response_class <- function(probe_type, accuracy) {
  probe_type <- as.character(probe_type)
  accuracy <- as.integer(accuracy)
  if (length(probe_type) != length(accuracy) ||
      any(!probe_type %in% c("old", "lure")) ||
      any(!accuracy %in% 0:1)) {
    stop("Response class requires old/lure probe type and binary accuracy.")
  }
  value <- ifelse(
    probe_type == "old",
    ifelse(accuracy == 1L, "hit", "miss"),
    ifelse(accuracy == 1L, "correct_rejection", "false_alarm")
  )
  factor(
    value,
    levels = c("hit", "miss", "correct_rejection", "false_alarm")
  )
}

followup_metadata <- function(raw, cohort) {
  pairs <- cohort$pairs[c(
    "participant", "item", "item_fold", "probe_type", "degradation",
    "accuracy"
  )]
  combined <- raw$retrieval[c(
    "participant", "age", "item", "cue_duration",
    "retrieval_image_version"
  )]
  combined <- combined[!duplicated(combined[c("participant", "item")]), ]
  merge(
    pairs, combined, by = c("participant", "item"),
    all.x = TRUE, sort = FALSE
  )
}

followup_study_tables <- function(raw, cohort, reference_presentation,
                                  source_presentation) {
  reference_presentation <- as.integer(reference_presentation)
  source_presentation <- as.integer(source_presentation)
  if (!reference_presentation %in% 1:4 ||
      !source_presentation %in% 1:4 ||
      reference_presentation >= source_presentation) {
    stop("Study contrasts require an earlier reference presentation.")
  }
  metadata <- followup_metadata(raw, cohort)
  study <- merge(
    raw$study,
    metadata[c("participant", "item", "item_fold")],
    by = c("participant", "item"), all = FALSE, sort = FALSE
  )
  make_part <- function(presentation) {
    part <- study[
      study$presentation == presentation,
      c(
        "participant", "age", "item", "item_fold", "study_image_version",
        "nfix", "fixgroup"
      )
    ]
    part <- merge(
      part,
      metadata[setdiff(names(metadata), c("age", "item_fold"))],
      by = c("participant", "item"), all.x = TRUE, sort = FALSE
    )
    part
  }
  reference <- make_part(reference_presentation)
  signal <- make_part(source_presentation)
  task <- sprintf(
    "study_p%d_p%d", reference_presentation, source_presentation
  )
  reference$task <- signal$task <- task
  reference$condition <- signal$condition <- "signal"
  reference$cue_duration <- signal$cue_duration <- 2500
  reference$presentation <- reference_presentation
  signal$presentation <- source_presentation
  reference <- reference[order(reference$participant, reference$item), ]
  signal <- signal[order(signal$participant, signal$item), ]
  rownames(reference) <- rownames(signal) <- NULL
  if (nrow(reference) != nrow(cohort$pairs) ||
      nrow(signal) != nrow(reference) ||
      any(reference$nfix < 3L) || any(signal$nfix < 3L)) {
    stop("Study-repeat tables violate the frozen support contract.")
  }
  list(
    reference = tibble::as_tibble(reference),
    signal = tibble::as_tibble(signal)
  )
}

followup_read_retrieval_phase <- function(phase) {
  phase <- match.arg(phase, followup_phases)
  files <- probe_delay_data_files()
  retrieval <- probe_delay_read_csv(files[["retrieval"]], "retrieval")
  probe_delay_retrieval_paths(retrieval, phase)
}

followup_phase_tables <- function(raw, cohort, phase,
                                  phase_paths = NULL) {
  phase <- match.arg(phase, followup_phases)
  if (is.null(phase_paths)) {
    phase_paths <- followup_read_retrieval_phase(phase)
  }
  metadata <- cohort$pairs[c("participant", "item", "item_fold")]
  reference <- merge(
    raw$study[raw$study$presentation == 4L, ], metadata,
    by = c("participant", "item"), all = FALSE, sort = FALSE
  )
  signal <- merge(
    phase_paths[phase_paths$probe_type %in% c("old", "lure"), ], metadata,
    by = c("participant", "item"), all = FALSE, sort = FALSE
  )
  signal$accuracy <- as.integer(signal$Accuracy)
  signal$response_class <- followup_response_class(
    signal$probe_type, signal$accuracy
  )
  signal$saliency <- as.integer(signal$degradation)
  reference$task <- signal$task <- phase
  reference$condition <- signal$condition <- "signal"
  reference <- reference[order(reference$participant, reference$item), ]
  signal <- signal[order(signal$participant, signal$item), ]
  rownames(reference) <- rownames(signal) <- NULL
  if (anyDuplicated(signal[c("participant", "item")]) ||
      any(signal$nfix < 1L) || nrow(signal) > nrow(cohort$pairs)) {
    stop("Phase tables violate the score-blind eligibility contract.")
  }
  list(
    reference = tibble::as_tibble(reference),
    signal = tibble::as_tibble(signal)
  )
}

followup_available_candidate_plan <- function(candidate_plan, signal) {
  keys <- paste(signal$participant, signal$item, sep = ":")
  plan <- candidate_plan[candidate_plan$candidate_set_id %in% keys, ]
  if (nrow(plan) != 5L * nrow(signal) ||
      any(table(plan$candidate_set_id) != 5L) ||
      any(tapply(plan$is_true, plan$candidate_set_id, sum) != 1L)) {
    stop("Available phase candidates violate the five-candidate contract.")
  }
  plan
}

followup_specs <- function(task, smoke = FALSE) {
  specs <- full_recognition_specs(smoke)
  if (task %in% followup_phases) {
    # Density is the mandatory phase comparator. MultiMatch is undefined for
    # most probe paths and elastic matching is not part of the frozen primary
    # phase family, so avoid paying their training cost on phase-only runs.
    specs$baseline$methods <- "density"
  } else if (task %in% followup_study_contrasts$task) {
    # Study sensitivity needs the requested density and MultiMatch landscape;
    # elastic matching is outside this follow-up question and is already
    # represented in the completed GW-13 court.
    specs$baseline$methods <- c("density", "multimatch")
  }
  specs
}

followup_fold_subsets <- function(tables, candidate_plan, fold) {
  eval_target <- with(
    tables$signal,
    participant %in% fold$eval_participants & item %in% fold$eval_items
  )
  train_target <- with(
    tables$signal,
    !participant %in% fold$eval_participants & !item %in% fold$eval_items
  )
  source_train <- tables$signal[train_target, , drop = FALSE]
  source_eval <- tables$signal[eval_target, , drop = FALSE]

  train_key <- paste(
    source_train$participant, source_train$item, sep = ":"
  )
  reference_key <- paste(
    tables$reference$participant, tables$reference$item, sep = ":"
  )
  ref_train <- tables$reference[reference_key %in% train_key, , drop = FALSE]
  if (nrow(ref_train) != nrow(source_train)) {
    stop("Training phase support does not have one study template per source.")
  }

  source_eval$candidate_set_id <- paste(
    source_eval$participant, source_eval$item, sep = ":"
  )
  plan <- candidate_plan[
    candidate_plan$candidate_set_id %in% source_eval$candidate_set_id,
    , drop = FALSE
  ]
  reference_lookup <- tables$reference
  names(reference_lookup)[names(reference_lookup) == "item"] <-
    "candidate_item"
  ref_eval <- merge(
    plan, reference_lookup,
    by = c("participant", "candidate_item"), all.x = TRUE, sort = FALSE
  )
  ref_eval$item <- ref_eval$candidate_item
  ref_eval <- ref_eval[order(
    ref_eval$candidate_set_id, ref_eval$candidate_position
  ), , drop = FALSE]
  source_eval <- source_eval[order(source_eval$candidate_set_id), ]
  if (nrow(ref_eval) != nrow(source_eval) * 5L ||
      anyNA(ref_eval$fixgroup)) {
    stop("Could not expand support-aware follow-up candidate sets.")
  }
  list(
    ref_train = ref_train,
    source_train = source_train,
    ref_eval = tibble::as_tibble(ref_eval),
    source_eval = tibble::as_tibble(source_eval)
  )
}

followup_score_fold_method <- function(
    tables, candidate_plan, fold, specs, method,
    workers = 1L, seed = followup_seed, smoke = FALSE) {
  value <- full_recognition_score_fold_method(
    tables, candidate_plan, fold, specs, method,
    workers = workers, seed = seed, smoke = smoke,
    subset_function = followup_fold_subsets
  )
  task <- unique(as.character(tables$signal$task))
  if (length(task) != 1L) stop("A follow-up scoring call must contain one task.")
  value$resources$task <- task
  value$warp$task <- task
  value$protocol_version <- followup_protocol_version
  value
}

followup_checkpoint <- function(output_dir, task, method, fold_id) {
  file.path(
    output_dir,
    sprintf("checkpoint-%s-%s-fold-%02d.rds", task, method, fold_id)
  )
}

followup_task_tables <- function(raw, cohort, task) {
  study <- followup_study_contrasts[
    followup_study_contrasts$task == task, , drop = FALSE
  ]
  if (nrow(study) == 1L) {
    return(followup_study_tables(
      raw, cohort, study$reference_presentation,
      study$source_presentation
    ))
  }
  followup_phase_tables(raw, cohort, task)
}

followup_score_tasks <- function(
    tasks = c(followup_study_contrasts$task, followup_phases),
    output_dir = file.path(
      "inst", "validation", "gaze-weave-followup-full-cohort-results"
    ),
    method_groups = c("replay", "baselines", "transport_v2"),
    fold_ids = NULL, workers = NULL, smoke = FALSE, resume = TRUE) {
  tasks <- match.arg(
    tasks,
    c(followup_study_contrasts$task, followup_phases),
    several.ok = TRUE
  )
  method_groups <- match.arg(
    method_groups, c("replay", "baselines", "transport_v2"),
    several.ok = TRUE
  )
  raw <- full_recognition_read_inputs(verify = !smoke)
  cohort <- full_recognition_select_cohort(raw)
  if (smoke) cohort <- full_recognition_smoke_cohort(cohort)
  candidate_plan <- full_recognition_candidate_plan(cohort, followup_seed)
  fold_plan <- full_recognition_fold_plan(cohort, followup_seed)
  folds <- fold_plan$folds
  if (!is.null(fold_ids)) {
    folds <- folds[vapply(
      folds, function(fold) fold$id %in% fold_ids, logical(1)
    )]
  }
  if (length(folds) == 0L) stop("fold_ids selected no outer folds.")
  if (is.null(workers)) {
    detected <- parallel::detectCores(logical = FALSE)
    if (!is.finite(detected) || detected < 1L) detected <- 1L
    workers <- if (smoke) min(2L, detected) else min(4L, detected)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (task in tasks) {
    tables <- followup_task_tables(raw, cohort, task)
    plan <- followup_available_candidate_plan(candidate_plan, tables$signal)
    specs <- followup_specs(task, smoke)
    for (fold in folds) {
      for (method in method_groups) {
        path <- followup_checkpoint(output_dir, task, method, fold$id)
        if (resume && file.exists(path)) {
          checkpoint <- readRDS(path)
          if (!identical(
            checkpoint$protocol_version, followup_protocol_version
          ) || !identical(checkpoint$smoke, smoke) ||
              !identical(checkpoint$seed, followup_seed)) {
            stop("A checkpoint uses another follow-up protocol: ", path)
          }
          message("Using checkpoint: ", basename(path))
          next
        }
        message(
          "Follow-up scoring: ", task, ", fold ", fold$id,
          ", method ", method
        )
        checkpoint <- followup_score_fold_method(
          tables, plan, fold, specs, method,
          workers = workers, seed = followup_seed, smoke = smoke
        )
        saveRDS(checkpoint, path, version = 3)
      }
    }
  }
  invisible(list(
    cohort = cohort, candidate_plan = candidate_plan,
    fold_plan = fold_plan, tasks = tasks
  ))
}

followup_read_checkpoints <- function(
    output_dir = file.path(
      "inst", "validation", "gaze-weave-followup-full-cohort-results"
    ),
    tasks = c(followup_study_contrasts$task, followup_phases),
    method_groups = c("replay", "baselines", "transport_v2"),
    fold_ids = 1:4, strict = TRUE) {
  paths <- unlist(lapply(tasks, function(task) {
    unlist(lapply(method_groups, function(method) {
      vapply(
        fold_ids, function(fold) {
          followup_checkpoint(output_dir, task, method, fold)
        }, character(1)
      )
    }), use.names = FALSE)
  }), use.names = FALSE)
  if (strict && !all(file.exists(paths))) {
    stop("One or more requested follow-up checkpoints are missing.")
  }
  paths <- paths[file.exists(paths)]
  checkpoints <- lapply(paths, readRDS)
  if (any(vapply(
    checkpoints,
    function(value) !identical(
      value$protocol_version, followup_protocol_version
    ), logical(1)
  ))) {
    stop("A follow-up checkpoint has another protocol version.")
  }
  list(
    scored = dplyr::bind_rows(lapply(checkpoints, `[[`, "scored")),
    resources = dplyr::bind_rows(lapply(checkpoints, `[[`, "resources")),
    warps = dplyr::bind_rows(lapply(checkpoints, `[[`, "warp")),
    audit = dplyr::bind_rows(lapply(checkpoints, `[[`, "audit")),
    checkpoints = checkpoints,
    paths = paths
  )
}

followup_weighted_mean <- function(value, weight) {
  keep <- is.finite(value) & is.finite(weight) & weight > 0
  if (!any(keep) || sum(weight[keep]) <= 0) return(NA_real_)
  stats::weighted.mean(value[keep], weight[keep])
}

followup_bootstrap_plan <- function(cohort, draws = followup_bootstrap_draws,
                                    seed = followup_seed) {
  tab <- cohort$pairs[c("participant", "item")]
  tab$trial_key <- paste(tab$participant, tab$item, sep = ":")
  recognition_bootstrap_plan(tab, draws, seed)
}

followup_score_weights <- function(tab, plan) {
  tab$trial_key <- paste(tab$participant, tab$item, sep = ":")
  position <- match(tab$trial_key, plan$trial_keys)
  if (anyNA(position)) {
    stop("Follow-up scores do not map to the bootstrap support.")
  }
  plan$weights[position, , drop = FALSE]
}

followup_study_summary <- function(scored, cohort,
                                   draws = followup_bootstrap_draws,
                                   seed = followup_seed) {
  keep <- scored$task %in% followup_study_contrasts$task &
    scored$calibrated & scored$status == "scored"
  tab <- scored[keep, , drop = FALSE]
  plan <- followup_bootstrap_plan(cohort, draws, seed)
  groups <- split(
    seq_len(nrow(tab)), interaction(tab$task, tab$method, drop = TRUE)
  )
  dplyr::bind_rows(lapply(groups, function(index) {
    part <- tab[index, , drop = FALSE]
    weights <- followup_score_weights(part, plan)
    boot <- vapply(seq_len(ncol(weights)), function(draw) {
      followup_weighted_mean(part$gaze_info_bits, weights[, draw])
    }, numeric(1))
    interval <- recognition_interval(boot, 0.95)
    data.frame(
      task = part$task[[1L]], method = part$method[[1L]],
      trials = nrow(part),
      mean_info_bits = mean(part$gaze_info_bits),
      lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
      mean_log_loss = mean(part$log_loss),
      mean_rank = mean(part$template_rank),
      top1_credit = mean(part$top1_credit),
      convergence = mean(part$converged),
      stringsAsFactors = FALSE
    )
  }))
}

followup_prepare_phase_scores <- function(scored) {
  tab <- scored[scored$task %in% followup_phases, , drop = FALSE]
  tab$saliency <- as.integer(tab$degradation)
  tab$response_class <- followup_response_class(
    tab$probe_type, tab$accuracy
  )
  tab
}

followup_phase_contrasts <- function(tab, weight = NULL) {
  if (is.null(weight)) weight <- rep(1, nrow(tab))
  cell <- function(phase, saliency, response_class) {
    keep <- tab$task == phase & tab$saliency == saliency &
      as.character(tab$response_class) == response_class
    followup_weighted_mean(tab$gaze_info_bits[keep], weight[keep])
  }
  old_difference <- function(phase, saliency) {
    cell(phase, saliency, "hit") - cell(phase, saliency, "miss")
  }
  c(
    delay_hit_minus_miss_20 = old_difference("delay", 20L),
    late_delay_hit_minus_miss_20 = old_difference("late_delay", 20L),
    delay_minus_probe_hit_minus_miss_20 =
      old_difference("delay", 20L) - old_difference("probe", 20L),
    delay_false_alarm_minus_correct_rejection_20 =
      cell("delay", 20L, "false_alarm") -
      cell("delay", 20L, "correct_rejection"),
    delay_hit_minus_miss_20_minus_100 =
      old_difference("delay", 20L) - old_difference("delay", 100L)
  )
}

followup_phase_effects <- function(scored, cohort,
                                   draws = followup_bootstrap_draws,
                                   seed = followup_seed) {
  tab <- followup_prepare_phase_scores(scored)
  tab <- tab[tab$calibrated & tab$status == "scored", , drop = FALSE]
  plan <- followup_bootstrap_plan(cohort, draws, seed)
  methods <- intersect(
    c("transport_v2", "replay", full_recognition_calibrated_comparators),
    unique(tab$method)
  )
  dplyr::bind_rows(lapply(methods, function(method) {
    part <- tab[tab$method == method, , drop = FALSE]
    estimate <- followup_phase_contrasts(part)
    weights <- followup_score_weights(part, plan)
    bootstrap <- vapply(seq_len(ncol(weights)), function(draw) {
      followup_phase_contrasts(part, weights[, draw])
    }, numeric(length(estimate)))
    if (ncol(weights) == 1L) bootstrap <- matrix(bootstrap, ncol = 1L)
    rownames(bootstrap) <- names(estimate)
    family_probability <- if (method %in% full_recognition_engines) {
      1 - 0.05 / 10
    } else {
      0.95
    }
    dplyr::bind_rows(lapply(names(estimate), function(name) {
      ordinary <- recognition_interval(bootstrap[name, ], 0.95)
      family <- recognition_interval(
        bootstrap[name, ], family_probability
      )
      data.frame(
        method = method, contrast = name, estimate = estimate[[name]],
        lower_95 = ordinary[["lower"]], upper_95 = ordinary[["upper"]],
        lower_family = family[["lower"]],
        upper_family = family[["upper"]],
        valid_draws = sum(is.finite(bootstrap[name, ])),
        valid_fraction = mean(is.finite(bootstrap[name, ])),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

followup_phase_cells <- function(scored) {
  tab <- followup_prepare_phase_scores(scored)
  keep <- tab$calibrated & tab$status == "scored"
  tab <- tab[keep, , drop = FALSE]
  groups <- split(
    seq_len(nrow(tab)),
    interaction(
      tab$method, tab$task, tab$saliency, tab$response_class,
      drop = TRUE
    )
  )
  dplyr::bind_rows(lapply(groups, function(index) {
    part <- tab[index, , drop = FALSE]
    data.frame(
      method = part$method[[1L]], phase = part$task[[1L]],
      saliency = part$saliency[[1L]],
      response_class = as.character(part$response_class[[1L]]),
      trials = nrow(part), mean_info_bits = mean(part$gaze_info_bits),
      mean_rank = mean(part$template_rank),
      top1_credit = mean(part$top1_credit),
      stringsAsFactors = FALSE
    )
  }))
}
