# Frozen Transport-v3 0-3000 ms retrieval court.
#
# Measurement is scored and checksummed without behavioral columns. The
# behavioral functions refuse to merge outcomes until the measurement manifest
# verifies byte for byte. Participant-linked artifacts remain Git-ignored.

v3_retrieval_protocol_version <- 1L
v3_retrieval_seed <- 20260825L
v3_retrieval_bootstrap_draws <- 2000L
v3_retrieval_output_dir <- file.path(
  "inst", "validation", "gaze-weave-transport-v3-retrieval-results"
)
v3_retrieval_freeze_dir <- file.path(
  "inst", "validation", "gaze-weave-transport-v3-retrieval-freeze"
)
v3_retrieval_forbidden_score_columns <- c(
  "Response", "response", "Accuracy", "accuracy", "probe_type",
  "degradation", "saliency", "said_old", "confidence"
)
v3_retrieval_methods <- c(
  "transport_v3", "transport_v2", "replay",
  "density_sigma_80_raw_calibrated",
  "density_sigma_160_raw_calibrated",
  "density_ridge_registered",
  "multimatch_mm_vector_raw_calibrated",
  "multimatch_mm_direction_raw_calibrated",
  "multimatch_mm_length_raw_calibrated",
  "multimatch_mm_position_raw_calibrated",
  "multimatch_mm_duration_raw_calibrated",
  "multimatch_mm_position_emd_raw_calibrated",
  "multimatch_ridge_registered"
)

v3_retrieval_repeat_file <- system.file(
  "validation", "gaze-weave-transport-v3-repeated-viewing.R",
  package = "eyesim"
)
if (!nzchar(v3_retrieval_repeat_file)) {
  v3_retrieval_repeat_file <- file.path(
    "inst", "validation", "gaze-weave-transport-v3-repeated-viewing.R"
  )
}
if (!file.exists(v3_retrieval_repeat_file)) {
  stop("Run the retrieval court from the eyesim source root.")
}
source(v3_retrieval_repeat_file, local = TRUE)

v3_retrieval_spec <- function(smoke = FALSE) v3_repeated_spec(smoke)

v3_retrieval_score_blind_tables <- function(raw, cohort) {
  pair_keys <- cohort$pairs[c("participant", "item", "item_fold")]
  reference <- merge(
    raw$study, pair_keys, by = c("participant", "item"),
    all = FALSE, sort = FALSE
  )
  reference <- reference[c(
    "participant", "item", "item_fold", "presentation", "nfix", "fixgroup"
  )]
  names(reference)[names(reference) == "presentation"] <- "episode_id"
  source <- merge(
    raw$retrieval, pair_keys, by = c("participant", "item"),
    all = FALSE, sort = FALSE
  )
  source <- source[c("participant", "item", "item_fold", "nfix", "fixgroup")]
  source$task <- "retrieval_0_3000"
  source$condition <- "score_blind"
  source <- source[order(source$participant, source$item), ]
  reference <- reference[order(
    reference$participant, reference$item, reference$episode_id
  ), ]
  if (nrow(source) != nrow(cohort$pairs) ||
      nrow(reference) != 4L * nrow(source) ||
      any(table(paste(reference$participant, reference$item, sep = ":")) != 4L) ||
      any(v3_retrieval_forbidden_score_columns %in% names(source)) ||
      any(v3_retrieval_forbidden_score_columns %in% names(reference))) {
    stop("Score-blind retrieval tables violate the frozen four-episode support.")
  }
  list(
    reference = tibble::as_tibble(reference),
    source = tibble::as_tibble(source)
  )
}

v3_retrieval_expand_pool <- function(reference, source, candidate_plan) {
  source$candidate_set_id <- paste(source$participant, source$item, sep = ":")
  plan <- candidate_plan[
    candidate_plan$candidate_set_id %in% source$candidate_set_id,
    , drop = FALSE
  ]
  lookup <- reference
  names(lookup)[names(lookup) == "item"] <- "candidate_item"
  pool <- merge(
    plan, lookup, by = c("participant", "candidate_item"),
    all.x = TRUE, sort = FALSE
  )
  pool$item <- pool$candidate_item
  pool$prior_weight <- 1
  pool <- pool[order(
    pool$candidate_set_id, pool$candidate_position, pool$episode_id
  ), ]
  source <- source[order(source$candidate_set_id), ]
  counts <- table(pool$candidate_set_id)
  episode_counts <- table(
    paste(pool$candidate_set_id, pool$candidate_item, sep = "|")
  )
  if (nrow(pool) != 20L * nrow(source) || any(counts != 20L) ||
      any(episode_counts != 4L) || anyNA(pool$fixgroup) ||
      any(tapply(pool$is_true, pool$candidate_set_id, sum) != 4L)) {
    stop("Retrieval candidate expansion violates K=5 by four episodes.")
  }
  list(reference = tibble::as_tibble(pool), source = tibble::as_tibble(source))
}

v3_retrieval_fold_data <- function(tables, candidate_plan, fold) {
  eval_mask <- with(
    tables$source,
    participant %in% fold$eval_participants & item %in% fold$eval_items
  )
  train_mask <- with(
    tables$source,
    !participant %in% fold$eval_participants & !item %in% fold$eval_items
  )
  source_train <- tables$source[train_mask, , drop = FALSE]
  source_eval <- tables$source[eval_mask, , drop = FALSE]
  train <- v3_retrieval_expand_pool(
    tables$reference, source_train, candidate_plan
  )
  eval <- v3_retrieval_expand_pool(
    tables$reference, source_eval, candidate_plan
  )
  train_key <- paste(source_train$participant, source_train$item, sep = ":")
  reference_key <- paste(
    tables$reference$participant, tables$reference$item, sep = ":"
  )
  ref_true_train <- tables$reference[
    reference_key %in% train_key, , drop = FALSE
  ]
  if (nrow(ref_true_train) != 4L * nrow(source_train)) {
    stop("Outer training support lacks four true study episodes.")
  }
  list(
    ref_true_train = ref_true_train,
    source_train = train$source,
    ref_train_pool = train$reference,
    source_eval = eval$source,
    ref_eval_pool = eval$reference
  )
}

v3_retrieval_score_row <- function(source_row, reference_pool, spec, warp,
                                   calibration) {
  scored <- eyesim:::score_transport_v3_cv_row(
    source_row, reference_pool,
    c("participant", "item"), "candidate_set_id",
    "fixgroup", "fixgroup", "episode_id", "prior_weight",
    spec, warp, calibration$temperature, calibration$kappa
  )
  evidence <- scored$evidence
  candidate_converged <- vapply(scored$candidates, function(candidate) {
    isTRUE(candidate$convergence$converged)
  }, logical(1))
  path <- source_row$fixgroup[[1L]]
  data.frame(
    participant = source_row$participant[[1L]],
    item = source_row$item[[1L]],
    method = "transport_v3", gaze_info_bits = evidence$gaze_info_bits,
    log_loss = evidence$log_loss, brier_score = evidence$brier_score,
    posterior_true = evidence$posterior_true,
    prior_true = evidence$prior_true,
    template_rank = evidence$template_rank,
    top1_credit = evidence$top1_credit,
    candidate_count = evidence$candidate_count,
    effective_fixations = scored$quality$effective_fixations,
    total_duration = sum(path$duration), fixation_count = nrow(path),
    reliability = evidence$reliability,
    row_converged = scored$all_converged,
    converged_candidates = sum(candidate_converged),
    candidate_alignments = length(candidate_converged),
    common_episode_count = length(scored$common_episode_ids),
    candidates = I(list(evidence$candidates)),
    stringsAsFactors = FALSE
  )
}

v3_retrieval_score_fold <- function(tables, candidate_plan, fold, spec) {
  started <- proc.time()[["elapsed"]]
  data <- v3_retrieval_fold_data(tables, candidate_plan, fold)
  warp <- eyesim:::fit_transport_v3_warp(
    data$ref_true_train, data$source_train,
    c("participant", "item"), "fixgroup", "fixgroup", "episode_id", spec
  )
  data$source_train$..gaze_row_id <- seq_len(nrow(data$source_train))
  inner <- eyesim:::fit_transport_v3_inner_calibration(
    data$ref_train_pool, data$source_train,
    c("participant", "item"), "candidate_set_id",
    "fixgroup", "fixgroup", "episode_id", "prior_weight", spec,
    fold_contrast_on = NULL
  )
  rows <- lapply(seq_len(nrow(data$source_eval)), function(index) {
    v3_retrieval_score_row(
      data$source_eval[index, , drop = FALSE], data$ref_eval_pool,
      spec, warp, inner$calibration
    )
  })
  scored <- dplyr::bind_rows(rows)
  scored$outer_fold <- fold$id
  elapsed <- proc.time()[["elapsed"]] - started
  list(
    scored = scored,
    resources = data.frame(
      outer_fold = fold$id, train_trials = nrow(data$source_train),
      eval_trials = nrow(data$source_eval), elapsed_seconds = elapsed,
      seconds_per_trial = elapsed / nrow(scored), stringsAsFactors = FALSE
    ),
    audit = data.frame(
      outer_fold = fold$id, train_trials = nrow(data$source_train),
      eval_trials = nrow(data$source_eval),
      participant_overlap = length(intersect(
        unique(data$source_train$participant), unique(data$source_eval$participant)
      )),
      item_overlap = length(intersect(
        unique(data$source_train$item), unique(data$source_eval$item)
      )),
      episode_count_min = min(scored$common_episode_count),
      episode_count_max = max(scored$common_episode_count),
      row_convergence = mean(scored$row_converged),
      alignment_convergence = sum(scored$converged_candidates) /
        sum(scored$candidate_alignments),
      stringsAsFactors = FALSE
    ),
    calibration = inner$calibration, warp = warp$info,
    protocol_version = v3_retrieval_protocol_version,
    seed = v3_retrieval_seed
  )
}

v3_retrieval_checkpoint <- function(output_dir, fold_id) {
  file.path(output_dir, sprintf("checkpoint-transport-v3-fold-%02d.rds", fold_id))
}

v3_retrieval_source_files <- function() {
  sort(unique(c(
    Sys.glob(file.path("R", "gaze_weave_transport_v3*.R")),
    "src/transport_v3.cpp",
    "inst/validation/gaze-weave-transport-v3-retrieval.R",
    "inst/validation/gaze-weave-transport-v3-manifest.json",
    "inst/validation/GAZEWEAVE-TRANSPORT-V3-PROTOCOL.md",
    "inst/validation/GAZEWEAVE-TRANSPORT-V3-DECISIONS.md"
  )))
}

v3_retrieval_freeze <- function(
    output_dir = v3_retrieval_output_dir,
    freeze_dir = v3_retrieval_freeze_dir) {
  raw <- full_recognition_read_inputs(verify = TRUE)
  cohort <- full_recognition_select_cohort(raw)
  candidate_plan <- full_recognition_candidate_plan(cohort, v3_retrieval_seed)
  fold_plan <- full_recognition_fold_plan(cohort, v3_retrieval_seed)
  design <- full_recognition_validate_design(
    cohort, candidate_plan, fold_plan, strict = TRUE
  )
  tables <- v3_retrieval_score_blind_tables(raw, cohort)
  spec <- v3_retrieval_spec(FALSE)
  private <- file.path(output_dir, "freeze-private")
  dir.create(private, recursive = TRUE, showWarnings = FALSE)
  dir.create(freeze_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(candidate_plan, file.path(private, "candidate-plan.csv"),
                   row.names = FALSE)
  utils::write.csv(fold_plan$participant_map,
                   file.path(private, "participant-folds.csv"), row.names = FALSE)
  utils::write.csv(fold_plan$item_map,
                   file.path(private, "item-folds.csv"), row.names = FALSE)
  utils::write.csv(design$fold_audit,
                   file.path(private, "fold-audit.csv"), row.names = FALSE)
  saveRDS(spec, file.path(private, "spec.rds"), version = 3)
  saveRDS(tables, file.path(private, "score-blind-tables.rds"), version = 3)
  private_files <- list.files(private, full.names = TRUE)
  design_manifest <- data.frame(
    artifact = basename(private_files),
    md5 = unname(tools::md5sum(private_files)), stringsAsFactors = FALSE
  )
  utils::write.csv(
    design_manifest, file.path(freeze_dir, "design-manifest.csv"),
    row.names = FALSE
  )
  source_files <- v3_retrieval_source_files()
  source_manifest <- data.frame(
    file = source_files, md5 = unname(tools::md5sum(source_files)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    source_manifest, file.path(freeze_dir, "source-manifest.csv"),
    row.names = FALSE
  )
  thresholds <- data.frame(
    field = c(
      "protocol", "window_ms", "candidate_count", "episode_count",
      "coverage_nodes", "alignment_convergence_min", "outer_seed",
      "calibration_seed", "bootstrap_draws"
    ),
    value = c(
      "gazeweave-transport-v3/3.0.2", "[0,3000)", "5", "4", "12",
      "0.99", as.character(v3_retrieval_seed), "20260822",
      as.character(v3_retrieval_bootstrap_draws)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    thresholds, file.path(freeze_dir, "thresholds.csv"), row.names = FALSE
  )
  writeLines(
    capture.output(utils::sessionInfo()),
    file.path(freeze_dir, "session-info.txt")
  )
  receipt <- data.frame(
    frozen_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    protocol_version = v3_retrieval_protocol_version,
    score_rows = nrow(tables$source), study_episode_rows = nrow(tables$reference),
    candidates_per_trial = 5L, episodes_per_candidate = 4L,
    behavior_opened = FALSE, stringsAsFactors = FALSE
  )
  utils::write.csv(receipt, file.path(freeze_dir, "freeze-receipt.csv"),
                   row.names = FALSE)
  invisible(receipt)
}

v3_retrieval_verify_freeze <- function(
    output_dir = v3_retrieval_output_dir,
    freeze_dir = v3_retrieval_freeze_dir) {
  source_manifest <- utils::read.csv(
    file.path(freeze_dir, "source-manifest.csv"), stringsAsFactors = FALSE
  )
  if (!all(file.exists(source_manifest$file)) ||
      !identical(unname(tools::md5sum(source_manifest$file)),
                 source_manifest$md5)) {
    stop("The frozen Transport-v3 source manifest has drifted.")
  }
  design_manifest <- utils::read.csv(
    file.path(freeze_dir, "design-manifest.csv"), stringsAsFactors = FALSE
  )
  private <- file.path(output_dir, "freeze-private")
  design_paths <- file.path(private, design_manifest$artifact)
  if (!all(file.exists(design_paths)) ||
      !identical(unname(tools::md5sum(design_paths)), design_manifest$md5)) {
    stop("The frozen private design manifest has drifted.")
  }
  invisible(TRUE)
}

run_gaze_weave_transport_v3_retrieval_measurement <- function(
    output_dir = v3_retrieval_output_dir, fold_ids = NULL,
    smoke = FALSE, resume = TRUE) {
  if (!smoke) v3_retrieval_verify_freeze(output_dir)
  raw <- full_recognition_read_inputs(verify = !smoke)
  cohort <- full_recognition_select_cohort(raw)
  if (smoke) cohort <- full_recognition_smoke_cohort(cohort)
  candidate_plan <- full_recognition_candidate_plan(cohort, v3_retrieval_seed)
  fold_plan <- full_recognition_fold_plan(cohort, v3_retrieval_seed)
  tables <- v3_retrieval_score_blind_tables(raw, cohort)
  spec <- v3_retrieval_spec(smoke)
  folds <- fold_plan$folds
  if (!is.null(fold_ids)) {
    folds <- folds[vapply(folds, function(x) x$id %in% fold_ids, logical(1))]
  }
  if (!length(folds)) stop("fold_ids selected no retrieval folds.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (fold in folds) {
    path <- v3_retrieval_checkpoint(output_dir, fold$id)
    if (resume && file.exists(path)) {
      value <- readRDS(path)
      if (!identical(value$protocol_version, v3_retrieval_protocol_version) ||
          !identical(value$seed, v3_retrieval_seed)) {
        stop("A retrieval checkpoint belongs to another protocol.")
      }
      message("Using checkpoint: ", basename(path))
      next
    }
    message("Transport v3 frozen retrieval measurement: fold ", fold$id)
    saveRDS(
      v3_retrieval_score_fold(tables, candidate_plan, fold, spec),
      path, version = 3
    )
  }
  invisible(list(cohort = cohort, fold_plan = fold_plan,
                 candidate_plan = candidate_plan))
}

v3_retrieval_read_measurement <- function(
    output_dir = v3_retrieval_output_dir, strict = TRUE) {
  paths <- vapply(1:4, function(fold) {
    v3_retrieval_checkpoint(output_dir, fold)
  }, character(1))
  if (strict && !all(file.exists(paths))) stop("Retrieval checkpoints incomplete.")
  paths <- paths[file.exists(paths)]
  values <- lapply(paths, readRDS)
  list(
    scored = dplyr::bind_rows(lapply(values, `[[`, "scored")),
    resources = dplyr::bind_rows(lapply(values, `[[`, "resources")),
    audit = dplyr::bind_rows(lapply(values, `[[`, "audit")),
    checkpoints = values, paths = paths
  )
}

v3_retrieval_finalize_measurement <- function(
    output_dir = v3_retrieval_output_dir) {
  v3_retrieval_verify_freeze(output_dir)
  value <- v3_retrieval_read_measurement(output_dir, strict = TRUE)
  if (nrow(value$scored) != 1295L || any(value$audit$participant_overlap != 0L) ||
      any(value$audit$item_overlap != 0L) ||
      any(value$audit$episode_count_min != 4L) ||
      any(value$audit$episode_count_max != 4L)) {
    stop("The retrieval measurement violates frozen support or episode gates.")
  }
  convergence <- sum(value$scored$converged_candidates) /
    sum(value$scored$candidate_alignments)
  if (convergence < 0.99) stop("Retrieval alignment convergence is below 0.99.")
  checkpoint_manifest <- data.frame(
    file = basename(value$paths), md5 = unname(tools::md5sum(value$paths)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    checkpoint_manifest, file.path(output_dir, "checkpoint-manifest.csv"),
    row.names = FALSE
  )
  saveRDS(value$scored, file.path(output_dir, "measurement-scores.rds"),
          version = 3)
  verdict <- data.frame(
    measurement_frozen = TRUE, behavior_opened = FALSE,
    trials = nrow(value$scored), candidate_count = 5L, episode_count = 4L,
    mean_info_bits = mean(value$scored$gaze_info_bits),
    mean_log_loss = mean(value$scored$log_loss),
    top1_credit = mean(value$scored$top1_credit),
    alignment_convergence = convergence,
    stringsAsFactors = FALSE
  )
  utils::write.csv(verdict, file.path(output_dir, "measurement-verdict.csv"),
                   row.names = FALSE)
  invisible(verdict)
}

v3_retrieval_verify_measurement <- function(output_dir = v3_retrieval_output_dir) {
  v3_retrieval_verify_freeze(output_dir)
  manifest <- utils::read.csv(
    file.path(output_dir, "checkpoint-manifest.csv"), stringsAsFactors = FALSE
  )
  paths <- file.path(output_dir, manifest$file)
  verdict <- utils::read.csv(
    file.path(output_dir, "measurement-verdict.csv"), stringsAsFactors = FALSE
  )
  if (!all(file.exists(paths)) ||
      !identical(unname(tools::md5sum(paths)), manifest$md5) ||
      nrow(verdict) != 1L || !isTRUE(verdict$measurement_frozen) ||
      isTRUE(verdict$behavior_opened)) {
    stop("The retrieval measurement is not immutable and behavior-sealed.")
  }
  invisible(TRUE)
}

v3_retrieval_comparator_scores <- function(
    comparator_dir = file.path(
      "inst", "validation", "gaze-weave-recognition-full-cohort-results"
    )) {
  manifest <- utils::read.csv(
    file.path(comparator_dir, "manifest-md5.csv"), stringsAsFactors = FALSE
  )
  paths <- file.path(comparator_dir, manifest$file)
  if (!all(file.exists(paths)) ||
      !identical(unname(tools::md5sum(paths)), manifest$md5)) {
    stop("The frozen full-cohort comparator manifest failed verification.")
  }
  checkpoints <- unlist(lapply(1:4, function(fold) {
    lapply(c("replay", "baselines", "transport_v2"), function(method) {
      readRDS(full_recognition_checkpoint(comparator_dir, method, fold))
    })
  }), recursive = FALSE)
  scored <- dplyr::bind_rows(lapply(checkpoints, `[[`, "scored"))
  scored <- scored[scored$method %in% setdiff(v3_retrieval_methods,
                                               "transport_v3"), ]
  scored[c(
    "participant", "item", "outer_fold", "method", "gaze_info_bits",
    "log_loss", "template_rank", "top1_credit", "candidate_count"
  )]
}

v3_retrieval_measurement_panel <- function(
    output_dir = v3_retrieval_output_dir) {
  v3_retrieval_verify_measurement(output_dir)
  v3 <- readRDS(file.path(output_dir, "measurement-scores.rds"))
  quality <- v3[c(
    "participant", "item", "effective_fixations", "total_duration",
    "fixation_count"
  )]
  comparators <- merge(
    v3_retrieval_comparator_scores(), quality,
    by = c("participant", "item"), all.x = TRUE, sort = FALSE
  )
  columns <- c(
    "participant", "item", "outer_fold", "method", "gaze_info_bits",
    "log_loss", "template_rank", "top1_credit", "candidate_count",
    "effective_fixations", "total_duration", "fixation_count"
  )
  v3 <- v3[columns]
  result <- dplyr::bind_rows(v3, comparators[columns])
  support <- table(result$method)
  if (!all(v3_retrieval_methods %in% names(support)) ||
      any(support[v3_retrieval_methods] != 1295L)) {
    stop("Retrieval methods do not share the 1,295-row measurement support.")
  }
  result
}

v3_retrieval_response <- function(response) {
  response <- suppressWarnings(as.integer(response))
  ifelse(response == 0L | is.na(response), NA_integer_,
         as.integer(response %in% c(1L, 2L)))
}

v3_retrieval_behavior_data <- function(output_dir = v3_retrieval_output_dir) {
  panel <- v3_retrieval_measurement_panel(output_dir)
  raw <- full_recognition_read_inputs(verify = TRUE)
  metadata <- unique(raw$retrieval[c(
    "participant", "item", "probe_type", "degradation", "Response", "Accuracy"
  )])
  if (anyDuplicated(metadata[c("participant", "item")])) {
    stop("Behavior metadata is not unique by trial.")
  }
  metadata$said_old <- v3_retrieval_response(metadata$Response)
  metadata$sal10 <- (as.numeric(metadata$degradation) - 60) / 10
  result <- merge(panel, metadata, by = c("participant", "item"),
                  all.x = TRUE, sort = FALSE)
  result
}

v3_retrieval_log_loss <- function(response, probability) {
  probability <- pmin(pmax(probability, 1e-8), 1 - 1e-8)
  -(response * log(probability) + (1 - response) * log1p(-probability))
}

v3_retrieval_standardize <- function(train, evaluate, source, target) {
  center <- mean(train[[source]])
  scale <- stats::sd(train[[source]])
  if (!is.finite(center) || !is.finite(scale) || scale <= 0) {
    stop("Behavior predictor has invalid training scale: ", source)
  }
  train[[target]] <- (train[[source]] - center) / scale
  evaluate[[target]] <- (evaluate[[source]] - center) / scale
  list(train = train, evaluate = evaluate)
}

v3_retrieval_fit_behavior_fold <- function(train, evaluate, method, fold_id,
                                           family_label) {
  if (length(intersect(train$participant, evaluate$participant)) ||
      length(intersect(train$item, evaluate$item))) {
    stop("Behavior training and evaluation support overlaps.")
  }
  train$log_effective_fixations <- log(train$effective_fixations)
  evaluate$log_effective_fixations <- log(evaluate$effective_fixations)
  for (mapping in list(
    c("gaze_info_bits", "z_metric"),
    c("log_effective_fixations", "z_log_effective_fixations"),
    c("total_duration", "z_total_duration")
  )) {
    value <- v3_retrieval_standardize(
      train, evaluate, mapping[[1L]], mapping[[2L]]
    )
    train <- value$train
    evaluate <- value$evaluate
  }
  base_formula <- said_old ~ sal10 + z_log_effective_fixations + z_total_duration
  full_formula <- stats::update.formula(base_formula, ". ~ . + z_metric")
  base <- stats::glm(base_formula, family = stats::binomial(), data = train)
  full <- stats::glm(full_formula, family = stats::binomial(), data = train)
  probability_base <- stats::predict(base, evaluate, type = "response")
  probability_full <- stats::predict(full, evaluate, type = "response")
  loss_base <- v3_retrieval_log_loss(evaluate$said_old, probability_base)
  loss_full <- v3_retrieval_log_loss(evaluate$said_old, probability_full)
  data.frame(
    participant = evaluate$participant, item = evaluate$item,
    trial_key = paste(evaluate$participant, evaluate$item, sep = ":"),
    method = method, family = family_label, outer_fold = fold_id,
    said_old = evaluate$said_old, response = evaluate$Response,
    gaze_info_bits = evaluate$gaze_info_bits,
    probability_base = probability_base, probability_full = probability_full,
    loss_base = loss_base, loss_full = loss_full,
    behavior_info_bits = (loss_base - loss_full) / log(2),
    coefficient = unname(stats::coef(full)[["z_metric"]]),
    stringsAsFactors = FALSE
  )
}

v3_retrieval_behavior_crossfit <- function(scores, probe_type = "old") {
  scores <- scores[
    scores$probe_type == probe_type & !is.na(scores$said_old), , drop = FALSE
  ]
  rows <- list()
  index <- 1L
  for (method in v3_retrieval_methods) {
    part <- scores[scores$method == method, , drop = FALSE]
    for (fold_id in sort(unique(part$outer_fold))) {
      evaluate <- part[part$outer_fold == fold_id, , drop = FALSE]
      train <- part[
        !part$participant %in% evaluate$participant &
          !part$item %in% evaluate$item, , drop = FALSE
      ]
      rows[[index]] <- v3_retrieval_fit_behavior_fold(
        train, evaluate, method, fold_id, probe_type
      )
      index <- index + 1L
    }
  }
  result <- dplyr::bind_rows(rows)
  support <- table(result$method)
  if (length(unique(support)) != 1L) stop("Behavior support differs by method.")
  result
}

v3_retrieval_behavior_summary <- function(
    predictions, draws = v3_retrieval_bootstrap_draws,
    seed = v3_retrieval_seed) {
  reference <- predictions[predictions$method == "transport_v3", ]
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws, seed)
  dplyr::bind_rows(lapply(v3_retrieval_methods, function(method) {
    part <- predictions[predictions$method == method, ]
    part <- part[match(plan$trial_keys, part$trial_key), ]
    values <- vapply(seq_len(ncol(plan$weights)), function(draw) {
      followup_weighted_mean(part$behavior_info_bits, plan$weights[, draw])
    }, numeric(1))
    interval <- recognition_interval(values, 0.95)
    data.frame(
      method = method, family = part$family[[1L]], trials = nrow(part),
      mean_behavior_info_bits = mean(part$behavior_info_bits),
      lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
      mean_full_log_loss = mean(part$loss_full),
      mean_base_log_loss = mean(part$loss_base),
      mean_fold_coefficient = mean(part$coefficient),
      stringsAsFactors = FALSE
    )
  }))
}

v3_retrieval_glmm <- function(scores) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(data.frame(status = "missing_lme4", estimate = NA_real_,
                      standard_error = NA_real_, p_value = NA_real_))
  }
  part <- scores[
    scores$method == "transport_v3" & scores$probe_type == "old" &
      !is.na(scores$said_old), , drop = FALSE
  ]
  part$z_metric <- as.numeric(scale(part$gaze_info_bits))
  part$z_quality <- as.numeric(scale(log(part$effective_fixations)))
  part$z_duration <- as.numeric(scale(part$total_duration))
  part$participant <- factor(part$participant)
  part$item <- factor(part$item)
  fit <- tryCatch(
    lme4::glmer(
      said_old ~ sal10 + z_quality + z_duration + z_metric +
        (1 | participant) + (1 | item),
      family = stats::binomial(), data = part,
      control = lme4::glmerControl(
        optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)
      )
    ), error = function(error) error
  )
  if (inherits(fit, "error")) {
    return(data.frame(status = "error", estimate = NA_real_,
                      standard_error = NA_real_, p_value = NA_real_))
  }
  coefficient <- summary(fit)$coefficients["z_metric", ]
  data.frame(
    status = "scored", estimate = coefficient[["Estimate"]],
    standard_error = coefficient[["Std. Error"]],
    p_value = coefficient[["Pr(>|z|)"]],
    singular = lme4::isSingular(fit), stringsAsFactors = FALSE
  )
}

run_gaze_weave_transport_v3_retrieval_behavior <- function(
    output_dir = v3_retrieval_output_dir) {
  v3_retrieval_verify_measurement(output_dir)
  scores <- v3_retrieval_behavior_data(output_dir)
  old <- v3_retrieval_behavior_crossfit(scores, "old")
  lure <- v3_retrieval_behavior_crossfit(scores, "lure")
  old_summary <- v3_retrieval_behavior_summary(old)
  lure_summary <- v3_retrieval_behavior_summary(lure)
  glmm <- v3_retrieval_glmm(scores)
  confidence <- scores[
    scores$method == "transport_v3" & !is.na(scores$said_old), , drop = FALSE
  ]
  confidence_summary <- data.frame(
    role = "secondary_descriptive",
    trials = nrow(confidence),
    spearman_response_information = suppressWarnings(stats::cor(
      as.numeric(confidence$Response), confidence$gaze_info_bits,
      method = "spearman", use = "complete.obs"
    )),
    stringsAsFactors = FALSE
  )
  behavior_dir <- file.path(output_dir, "behavior")
  dir.create(behavior_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(old, file.path(behavior_dir, "old-predictions.rds"), version = 3)
  saveRDS(lure, file.path(behavior_dir, "lure-predictions.rds"), version = 3)
  utils::write.csv(old_summary, file.path(behavior_dir, "old-summary.csv"),
                   row.names = FALSE)
  utils::write.csv(lure_summary, file.path(behavior_dir, "lure-summary.csv"),
                   row.names = FALSE)
  utils::write.csv(glmm, file.path(behavior_dir, "glmm-diagnostic.csv"),
                   row.names = FALSE)
  utils::write.csv(confidence_summary,
                   file.path(behavior_dir, "confidence-secondary.csv"),
                   row.names = FALSE)
  verdict_path <- file.path(output_dir, "measurement-verdict.csv")
  verdict <- utils::read.csv(verdict_path, stringsAsFactors = FALSE)
  verdict$behavior_opened <- TRUE
  utils::write.csv(verdict, verdict_path, row.names = FALSE)
  list(old = old_summary, lure = lure_summary, glmm = glmm,
       confidence = confidence_summary)
}
