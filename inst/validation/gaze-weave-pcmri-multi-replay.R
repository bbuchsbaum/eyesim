# Frozen multi-presentation Replay measurement court for the private
# probe-delay recognition data.
#
# Protocol: inst/validation/GAZEWEAVE-PCMRI-MULTI-REPLAY.md
# Run from a source checkout after devtools::load_all(). Raw paths, trial-level
# scores, and checkpoints remain under the Git-ignored private results folder.

multi_replay_protocol_version <- 1L
multi_replay_seed <- 20260828L
multi_replay_item_seed <- 20260820L
multi_replay_bootstrap_draws <- 2000L
multi_replay_expected <- c(
  raw_participants = 46L,
  eligible_participants = 46L,
  eligible_pairs = 2056L,
  retained_participants = 45L,
  retained_pairs = 2055L,
  items = 120L,
  old = 1011L,
  lure = 1044L
)
multi_replay_method_groups <- c("p4", "all4", "density", "study_control")
multi_replay_methods <- c(
  "replay_p4_exhaustive",
  "replay_all4_unshrunk_exhaustive",
  "replay_all4_shrunk_exhaustive",
  "density_sigma_80_all4_exhaustive",
  "density_sigma_160_all4_exhaustive",
  "study_replay_intact",
  "study_replay_reversed"
)

multi_replay_parent_file <- system.file(
  "validation", "gaze-weave-recognition-full-cohort.R", package = "eyesim"
)
if (!nzchar(multi_replay_parent_file)) {
  multi_replay_parent_file <- file.path(
    "inst", "validation", "gaze-weave-recognition-full-cohort.R"
  )
}
if (!file.exists(multi_replay_parent_file)) {
  stop("Run the multi-presentation court from the eyesim source root.")
}
source(multi_replay_parent_file, local = TRUE)

multi_replay_eligible_pairs <- function(raw) {
  study_groups <- split(
    seq_len(nrow(raw$study)),
    interaction(raw$study$participant, raw$study$item, drop = TRUE)
  )
  study_pairs <- dplyr::bind_rows(lapply(study_groups, function(index) {
    part <- raw$study[index, , drop = FALSE]
    data.frame(
      participant = part$participant[[1L]],
      item = part$item[[1L]],
      four_presentations = identical(
        sort(unique(part$presentation)), 1:4
      ),
      study_min_nfix = min(part$nfix),
      stringsAsFactors = FALSE
    )
  }))
  study_pairs <- study_pairs[
    study_pairs$four_presentations & study_pairs$study_min_nfix >= 1L,
    c("participant", "item"), drop = FALSE
  ]
  retrieval <- raw$retrieval[
    raw$retrieval$probe_type %in% c("old", "lure") &
      raw$retrieval$nfix >= 1L, , drop = FALSE
  ]
  retrieval_key <- paste(retrieval$participant, retrieval$item, sep = ":")
  if (anyDuplicated(retrieval_key)) {
    stop("Combined retrieval mappings are not unique by participant and item.")
  }
  retrieval_pairs <- retrieval[c(
    "participant", "item", "probe_type", "degradation", "Accuracy"
  )]
  names(retrieval_pairs)[names(retrieval_pairs) == "Accuracy"] <- "accuracy"
  pairs <- merge(
    study_pairs, retrieval_pairs,
    by = c("participant", "item"), all = FALSE, sort = FALSE
  )
  pairs[order(pairs$participant, pairs$item), , drop = FALSE]
}

multi_replay_select_cohort <- function(
    raw, item_seed = multi_replay_item_seed) {
  eligible <- multi_replay_eligible_pairs(raw)
  set.seed(item_seed)
  item_order <- sample(sort(unique(eligible$item)))
  item_map <- data.frame(
    item = item_order,
    item_fold = rep(1:2, length.out = length(item_order)),
    stringsAsFactors = FALSE
  )
  eligible <- merge(eligible, item_map, by = "item", sort = FALSE)
  counts <- table(
    factor(eligible$participant, levels = sort(unique(eligible$participant))),
    factor(eligible$item_fold, levels = 1:2)
  )
  keep_participant <- rownames(counts)[counts[, 1L] >= 2L & counts[, 2L] >= 2L]
  pairs <- eligible[eligible$participant %in% keep_participant, , drop = FALSE]
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
    balance = list(
      probe_type = table(pairs$probe_type),
      degradation = table(pairs$degradation),
      accuracy = table(pairs$accuracy)
    )
  )
}

multi_replay_smoke_cohort <- function(cohort, participant_n = 4L) {
  participant_n <- as.integer(participant_n)
  if (participant_n < 4L || participant_n > length(cohort$participants)) {
    stop("Smoke participant_n must be between four and the cohort size.")
  }
  participants <- cohort$participants[seq_len(participant_n)]
  smoke <- cohort
  smoke$participants <- participants
  smoke$pairs <- cohort$pairs[
    cohort$pairs$participant %in% participants, , drop = FALSE
  ]
  smoke
}

multi_replay_candidate_plan <- function(cohort) {
  groups <- split(
    seq_len(nrow(cohort$pairs)),
    interaction(
      cohort$pairs$participant, cohort$pairs$item_fold,
      drop = TRUE, lex.order = TRUE
    )
  )
  plan <- dplyr::bind_rows(lapply(groups, function(index) {
    part <- cohort$pairs[index, , drop = FALSE]
    participant <- part$participant[[1L]]
    item_fold <- part$item_fold[[1L]]
    candidates <- sort(unique(part$item))
    dplyr::bind_rows(lapply(candidates, function(target) {
      data.frame(
        participant = participant,
        target_item = target,
        item_fold = item_fold,
        candidate_set_id = paste(participant, target, sep = ":"),
        candidate_position = seq_along(candidates),
        candidate_item = candidates,
        is_true = candidates == target,
        stringsAsFactors = FALSE
      )
    }))
  }))
  plan <- plan[order(
    plan$participant, plan$target_item, plan$candidate_position
  ), ]
  rownames(plan) <- NULL
  plan
}

multi_replay_validate_design <- function(cohort, candidate_plan, fold_plan,
                                         strict = TRUE) {
  observed <- c(
    raw_participants = cohort$raw_participant_n,
    eligible_participants = cohort$eligible_participant_n,
    eligible_pairs = cohort$eligible_pair_n,
    retained_participants = length(cohort$participants),
    retained_pairs = nrow(cohort$pairs),
    items = nrow(cohort$item_map),
    old = sum(cohort$pairs$probe_type == "old"),
    lure = sum(cohort$pairs$probe_type == "lure")
  )
  if (strict && !identical(
    as.integer(observed), as.integer(multi_replay_expected)
  )) {
    stop("The multi-presentation support differs from the frozen protocol.")
  }
  set_sizes <- table(candidate_plan$candidate_set_id)
  pool_sizes <- vapply(split(
    candidate_plan$candidate_item,
    interaction(
      candidate_plan$participant, candidate_plan$item_fold,
      drop = TRUE, lex.order = TRUE
    )
  ), function(value) length(unique(value)), integer(1))
  truth_sizes <- tapply(
    candidate_plan$is_true, candidate_plan$candidate_set_id, sum
  )
  if (length(set_sizes) != nrow(cohort$pairs) || any(truth_sizes != 1L) ||
      min(set_sizes) < 2L) {
    stop("Exhaustive candidate sets violate the truth contract.")
  }
  if (strict && (
      min(pool_sizes) != 5L || max(pool_sizes) != 33L ||
      stats::median(pool_sizes) != 25 ||
      abs(mean(pool_sizes) - 22.83) > 0.01)) {
    stop("Exhaustive candidate counts differ from the frozen protocol.")
  }

  evaluated <- character()
  audits <- lapply(fold_plan$folds, function(fold) {
    eval <- cohort$pairs$participant %in% fold$eval_participants &
      cohort$pairs$item %in% fold$eval_items
    train <- !cohort$pairs$participant %in% fold$eval_participants &
      !cohort$pairs$item %in% fold$eval_items
    evaluated <<- c(
      evaluated,
      paste(cohort$pairs$participant[eval], cohort$pairs$item[eval], sep = ":")
    )
    data.frame(
      outer_fold = fold$id,
      eval_trials = sum(eval),
      train_trials = sum(train),
      participant_overlap = length(intersect(
        cohort$pairs$participant[eval], cohort$pairs$participant[train]
      )),
      item_overlap = length(intersect(
        cohort$pairs$item[eval], cohort$pairs$item[train]
      )),
      stringsAsFactors = FALSE
    )
  })
  expected_keys <- paste(cohort$pairs$participant, cohort$pairs$item, sep = ":")
  if (!identical(sort(evaluated), sort(expected_keys)) ||
      anyDuplicated(evaluated)) {
    stop("Outer folds do not evaluate every retained trial exactly once.")
  }
  list(
    observed = observed,
    candidate_counts = data.frame(
      candidate_set_id = names(set_sizes),
      candidate_count = as.integer(set_sizes),
      stringsAsFactors = FALSE
    ),
    fold_audit = dplyr::bind_rows(audits)
  )
}

multi_replay_task_tables <- function(raw, cohort) {
  pair_keys <- cohort$pairs[c("participant", "item", "item_fold")]
  study <- merge(
    raw$study, pair_keys, by = c("participant", "item"),
    all = FALSE, sort = FALSE
  )
  retrieval <- merge(
    raw$retrieval, pair_keys, by = c("participant", "item"),
    all = FALSE, sort = FALSE
  )
  retrieval <- retrieval[retrieval$probe_type %in% c("old", "lure"), ]
  metadata <- retrieval[c(
    "participant", "age", "item", "probe_type", "degradation",
    "cue_duration", "Accuracy", "retrieval_image_version", "phase"
  )]
  names(metadata)[names(metadata) == "Accuracy"] <- "accuracy"
  reference <- study[c(
    "participant", "age", "item", "item_fold", "presentation",
    "study_image_version", "nfix", "fixgroup"
  )]
  reference <- merge(
    reference, metadata, by = c("participant", "age", "item"),
    all.x = TRUE, sort = FALSE
  )
  signal <- retrieval
  names(signal)[names(signal) == "Accuracy"] <- "accuracy"
  reference$task <- signal$task <- "combined"
  reference$condition <- signal$condition <- "signal"
  reference <- reference[order(
    reference$participant, reference$item, reference$presentation
  ), ]
  signal <- signal[order(signal$participant, signal$item), ]
  rownames(reference) <- rownames(signal) <- NULL
  reference_count <- table(paste(
    reference$participant, reference$item, sep = ":"
  ))
  if (nrow(signal) != nrow(cohort$pairs) ||
      nrow(reference) != 4L * nrow(signal) ||
      any(reference_count != 4L)) {
    stop("Episode tables do not contain four study paths per retrieval trial.")
  }
  list(
    reference = tibble::as_tibble(reference),
    signal = tibble::as_tibble(signal)
  )
}

multi_replay_fold_subsets <- function(tables, candidate_plan, fold) {
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
  train_key <- paste(source_train$participant, source_train$item, sep = ":")
  reference_key <- paste(
    tables$reference$participant, tables$reference$item, sep = ":"
  )
  ref_train <- tables$reference[reference_key %in% train_key, , drop = FALSE]

  eval_key <- paste(source_eval$participant, source_eval$item, sep = ":")
  ref_eval <- tables$reference[reference_key %in% eval_key, , drop = FALSE]
  source_eval$candidate_pool_id <- paste(
    source_eval$participant, source_eval$item_fold, sep = ":fold:"
  )
  source_eval <- source_eval[order(
    source_eval$participant, source_eval$item_fold, source_eval$item
  ), ]
  ref_eval <- ref_eval[order(
    ref_eval$participant, ref_eval$item_fold, ref_eval$item,
    ref_eval$presentation
  ), ]
  if (nrow(ref_train) != 4L * nrow(source_train) ||
      nrow(ref_eval) != 4L * nrow(source_eval) || anyNA(ref_eval$fixgroup)) {
    stop("Could not expand four-path training or exhaustive evaluation sets.")
  }
  list(
    ref_train = tibble::as_tibble(ref_train),
    source_train = tibble::as_tibble(source_train),
    ref_eval = tibble::as_tibble(ref_eval),
    source_eval = tibble::as_tibble(source_eval)
  )
}

multi_replay_specs <- function(smoke = FALSE) {
  screen <- gaze_screen(800, 600, unit = "px")
  warp <- gaze_warp_contraction(
    center = "screen", translation = TRUE, fit_by = NULL
  )
  replay_args <- list(
    grid_size = if (smoke) 16L else 48L,
    max_skip = 2L,
    student_df = 4,
    scale_floor = 6,
    transition_grid = list(
      background = c(0.03, 0.10),
      restart = c(0.02, 0.10),
      advance = c(0.25, 0.55),
      background_stay = c(0.85, 0.95)
    ),
    warp = warp,
    screen = screen
  )
  list(
    screen = screen,
    replay_p4 = do.call(
      gaze_replay_spec, c(replay_args, list(reliability = "none"))
    ),
    replay_all4 = do.call(
      gaze_replay_spec,
      c(replay_args, list(reliability = "effective_fixations"))
    ),
    density = gaze_baseline_spec(
      screen = screen,
      density_sigmas = c(80, 160),
      density_grid = if (smoke) 12L else 32L,
      methods = "density",
      warp = gaze_warp_none(),
      inner_folds = 2L
    )
  )
}

multi_replay_candidate_references <- function(ref_eval, source_row) {
  keep <- ref_eval$participant == source_row$participant[[1L]] &
    ref_eval$item_fold == source_row$item_fold[[1L]]
  result <- ref_eval[keep, , drop = FALSE]
  if (nrow(result) == 0L ||
      !source_row$item[[1L]] %in% result$item) {
    stop("A retrieval trial has no exhaustive candidate reference set.")
  }
  result
}

multi_replay_slim_candidates <- function(candidates) {
  columns <- intersect(c(
    "participant", "item", "candidate_key", "log_score", "prior",
    "base_log_posterior", "base_posterior", "log_posterior",
    "posterior", "is_true", "template_count"
  ), names(candidates))
  candidates[columns]
}

multi_replay_measurement_row <- function(source_row, scored, evidence,
                                         method, outer_fold) {
  diagnostics <- scored$alignment$diagnostics
  quality <- scored$quality
  tibble::tibble(
    participant = source_row$participant,
    item = source_row$item,
    item_fold = source_row$item_fold,
    outer_fold = outer_fold,
    candidate_pool_id = source_row$candidate_pool_id,
    method = method,
    gaze_info_bits = evidence$gaze_info_bits,
    base_gaze_info_bits = evidence$base_gaze_info_bits,
    log_loss = evidence$log_loss,
    base_log_loss = evidence$base_log_loss,
    posterior_true = evidence$posterior_true,
    prior_true = evidence$prior_true,
    template_rank = evidence$template_rank,
    top1_credit = evidence$top1_credit,
    candidate_count = evidence$candidate_count,
    temperature = evidence$temperature,
    reliability = evidence$reliability,
    raw_fixation_count = quality$raw_fixation_count,
    coalesced_fixation_count = quality$coalesced_fixation_count,
    effective_fixations = quality$effective_fixations,
    total_duration = quality$total_duration,
    duration_concentration = quality$duration_concentration,
    spatial_dispersion = quality$spatial_dispersion,
    replay_coverage = real_value(diagnostics, "replay_coverage"),
    background_coverage = real_value(diagnostics, "background_coverage"),
    expected_restarts = real_value(diagnostics, "expected_restarts"),
    restart_rate = real_value(diagnostics, "restart_rate"),
    spatial_rmse_px = real_value(diagnostics, "spatial_rmse"),
    template_count = real_value(diagnostics, "template_count"),
    template_effective_count = real_value(
      diagnostics, "template_effective_count"
    ),
    template_weight_max = real_value(diagnostics, "template_weight_max"),
    modal_template = as.character(
      if (is.null(diagnostics$modal_template)) NA else
        diagnostics$modal_template
    ),
    converged = scored$all_converged,
    candidates = list(multi_replay_slim_candidates(evidence$candidates))
  )
}

multi_replay_component_rows <- function(source_row, scored, outer_fold,
                                        method) {
  component <- scored$candidate_components
  true_key <- real_ns("gaze_key")(
    source_row, c("participant", "item"), "match_on"
  )
  component$is_true_candidate <-
    as.character(component$candidate_key) == true_key
  component$participant <- source_row$participant[[1L]]
  component$target_item <- source_row$item[[1L]]
  component$item_fold <- source_row$item_fold[[1L]]
  component$outer_fold <- outer_fold
  component$method <- method
  component
}

multi_replay_score_replay_fold <- function(subset, fold, spec,
                                           method_group) {
  if (!method_group %in% c("p4", "all4")) {
    stop("Replay method_group must be p4 or all4.")
  }
  is_all4 <- identical(method_group, "all4")
  ref_train <- if (is_all4) subset$ref_train else
    subset$ref_train[subset$ref_train$presentation == 4L, , drop = FALSE]
  ref_eval <- if (is_all4) subset$ref_eval else
    subset$ref_eval[subset$ref_eval$presentation == 4L, , drop = FALSE]
  template_on <- if (is_all4) "presentation" else NULL
  profile <- real_profile(method_group, "combined", fold$id, {
    model <- fit_gaze_replay_model(
      ref_train,
      subset$source_train,
      match_on = c("participant", "item"),
      contrast_on = "participant",
      spec = spec,
      template_on = template_on
    )
    scored <- lapply(seq_len(nrow(subset$source_eval)), function(i) {
      source_row <- subset$source_eval[i, , drop = FALSE]
      candidate_ref <- multi_replay_candidate_references(ref_eval, source_row)
      real_ns("score_gaze_replay_row")(
        source_row,
        candidate_ref,
        match_on = c("participant", "item"),
        contrast_on = NULL,
        refvar = "fixgroup",
        sourcevar = "fixgroup",
        model = model
      )
    })
    list(model = model, scored = scored)
  })
  model <- profile$value$model
  scored <- profile$value$scored
  rows <- list()
  components <- list()
  row_index <- 1L
  component_index <- 1L
  for (i in seq_along(scored)) {
    source_row <- subset$source_eval[i, , drop = FALSE]
    value <- scored[[i]]
    if (is_all4) {
      temperature_only <- model$calibration$temperature_only_temperature
      unshrunk <- real_ns("score_gaze_candidates")(
        log_score = value$candidates$log_score,
        true_index = which(value$candidates$is_true),
        candidate_key = value$candidates$item,
        temperature = temperature_only,
        reliability = 1,
        candidate_pool_id = source_row$candidate_pool_id[[1L]]
      )
      rows[[row_index]] <- multi_replay_measurement_row(
        source_row, value, unshrunk,
        "replay_all4_unshrunk_exhaustive", fold$id
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- multi_replay_measurement_row(
        source_row, value, value$evidence,
        "replay_all4_shrunk_exhaustive", fold$id
      )
      components[[component_index]] <- multi_replay_component_rows(
        source_row, value, fold$id,
        "replay_all4_shrunk_exhaustive"
      )
      component_index <- component_index + 1L
    } else {
      rows[[row_index]] <- multi_replay_measurement_row(
        source_row, value, value$evidence,
        "replay_p4_exhaustive", fold$id
      )
    }
    row_index <- row_index + 1L
  }
  calibration <- model$calibration
  calibration_gate <- data.frame(
    outer_fold = fold$id,
    method_group = method_group,
    temperature = model$temperature,
    reliability_kappa = model$reliability_kappa,
    calibration_log_loss = if (is.null(calibration)) NA_real_ else
      calibration$log_loss,
    temperature_only_log_loss = if (
      is.null(calibration$temperature_only_log_loss)
    ) NA_real_ else calibration$temperature_only_log_loss,
    shrinkage_not_worse = if (
      is.null(calibration$temperature_only_log_loss) ||
        !is.finite(calibration$temperature_only_log_loss)
    ) TRUE else calibration$log_loss <=
      calibration$temperature_only_log_loss + 1e-8,
    stringsAsFactors = FALSE
  )
  list(
    scored = dplyr::bind_rows(rows),
    components = dplyr::bind_rows(components),
    resources = profile$resource,
    calibration = calibration_gate,
    model_audit = list(
      parameters = model$parameters,
      calibration = model$calibration,
      template_on = model$template_on,
      template_count = model$training$template_count,
      warp = model$warp$info
    )
  )
}

multi_replay_reverse_path <- function(path) {
  n <- nrow(path)
  if (n < 2L) return(path)
  index <- rev(seq_len(n))
  duration <- path$duration[index]
  onset <- c(0, head(cumsum(duration), -1L))
  fixation_group(
    x = path$x[index],
    y = path$y[index],
    duration = duration,
    onset = onset
  )
}

multi_replay_study_control_subset <- function(subset) {
  build <- function(reference, role) {
    source <- reference[reference$presentation == 4L, , drop = FALSE]
    episode <- reference[reference$presentation %in% 1:3, , drop = FALSE]
    if (nrow(episode) != 3L * nrow(source)) {
      stop("Study control requires P4 sources and three earlier templates.")
    }
    if (identical(role, "eval")) {
      source$candidate_pool_id <- paste(
        source$participant, source$item_fold, sep = ":fold:"
      )
    }
    list(source = source, reference = episode)
  }
  train <- build(subset$ref_train, "train")
  eval <- build(subset$ref_eval, "eval")
  list(
    ref_train = train$reference,
    source_train = train$source,
    ref_eval = eval$reference,
    source_eval = eval$source
  )
}

multi_replay_score_study_control_fold <- function(subset, fold, spec) {
  control <- multi_replay_study_control_subset(subset)
  profile <- real_profile("study_control", "study", fold$id, {
    model <- fit_gaze_replay_model(
      control$ref_train,
      control$source_train,
      match_on = c("participant", "item"),
      contrast_on = "participant",
      spec = spec,
      template_on = "presentation"
    )
    intact <- vector("list", nrow(control$source_eval))
    shuffled <- vector("list", nrow(control$source_eval))
    for (i in seq_len(nrow(control$source_eval))) {
      source_row <- control$source_eval[i, , drop = FALSE]
      candidate_ref <- multi_replay_candidate_references(
        control$ref_eval, source_row
      )
      intact[[i]] <- real_ns("score_gaze_replay_row")(
        source_row, candidate_ref,
        match_on = c("participant", "item"), contrast_on = NULL,
        refvar = "fixgroup", sourcevar = "fixgroup", model = model
      )
      shuffled_row <- source_row
      shuffled_row$fixgroup <- list(multi_replay_reverse_path(
        source_row$fixgroup[[1L]]
      ))
      shuffled[[i]] <- real_ns("score_gaze_replay_row")(
        shuffled_row, candidate_ref,
        match_on = c("participant", "item"), contrast_on = NULL,
        refvar = "fixgroup", sourcevar = "fixgroup", model = model
      )
    }
    list(model = model, intact = intact, shuffled = shuffled)
  })
  rows <- list()
  components <- list()
  index <- 1L
  for (i in seq_len(nrow(control$source_eval))) {
    source_row <- control$source_eval[i, , drop = FALSE]
    rows[[index]] <- multi_replay_measurement_row(
      source_row, profile$value$intact[[i]],
      profile$value$intact[[i]]$evidence,
      "study_replay_intact", fold$id
    )
    index <- index + 1L
    shuffled_row <- source_row
    shuffled_row$fixgroup <- list(multi_replay_reverse_path(
      source_row$fixgroup[[1L]]
    ))
    rows[[index]] <- multi_replay_measurement_row(
      shuffled_row, profile$value$shuffled[[i]],
      profile$value$shuffled[[i]]$evidence,
      "study_replay_reversed", fold$id
    )
    index <- index + 1L
    components[[i]] <- multi_replay_component_rows(
      source_row, profile$value$intact[[i]], fold$id,
      "study_replay_intact"
    )
  }
  model <- profile$value$model
  calibration <- model$calibration
  calibration_gate <- data.frame(
    outer_fold = fold$id,
    method_group = "study_control",
    temperature = model$temperature,
    reliability_kappa = model$reliability_kappa,
    calibration_log_loss = calibration$log_loss,
    temperature_only_log_loss = calibration$temperature_only_log_loss,
    shrinkage_not_worse = calibration$log_loss <=
      calibration$temperature_only_log_loss + 1e-8,
    stringsAsFactors = FALSE
  )
  list(
    scored = dplyr::bind_rows(rows),
    components = dplyr::bind_rows(components),
    resources = profile$resource,
    calibration = calibration_gate,
    model_audit = list(
      parameters = model$parameters,
      calibration = model$calibration,
      template_on = model$template_on,
      template_count = model$training$template_count,
      warp = model$warp$info
    )
  )
}

multi_replay_density_episode_table <- function(reference, spec) {
  signatures <- lapply(
    reference$fixgroup,
    real_ns("baseline_density_signature"),
    spec = spec
  )
  key <- interaction(
    reference$participant, reference$item,
    drop = TRUE, lex.order = TRUE
  )
  rows <- lapply(split(seq_len(nrow(reference)), key), function(index) {
    if (length(index) != 4L) {
      stop("Density episode references require four study paths.")
    }
    averaged <- lapply(seq_along(spec$density_sigmas), function(sigma_index) {
      Reduce(
        `+`, lapply(signatures[index], `[[`, sigma_index)
      ) / length(index)
    })
    data.frame(
      participant = reference$participant[[index[[1L]]]],
      item = reference$item[[index[[1L]]]],
      item_fold = reference$item_fold[[index[[1L]]]],
      signature = I(list(averaged)),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

multi_replay_density_source_signatures <- function(source, spec) {
  lapply(
    source$fixgroup,
    real_ns("baseline_density_signature"),
    spec = spec
  )
}

multi_replay_density_score_set <- function(source_row, source_signature,
                                           episode, sigma_index) {
  candidate <- episode[
    episode$participant == source_row$participant[[1L]], , drop = FALSE
  ]
  candidate <- candidate[order(candidate$item), , drop = FALSE]
  true_index <- match(source_row$item[[1L]], candidate$item)
  if (is.na(true_index) || nrow(candidate) < 2L) {
    stop("A density row lacks a true candidate or nonmatch.")
  }
  score <- vapply(candidate$signature, function(value) {
    real_ns("baseline_cosine")(
      value[[sigma_index]], source_signature[[sigma_index]]
    )
  }, numeric(1))
  list(
    score = score,
    true_index = true_index,
    candidate_item = candidate$item
  )
}

multi_replay_density_row <- function(source_row, score_set, evidence,
                                     method, outer_fold) {
  quality <- real_ns("gaze_replay_path_quality")(
    source_row$fixgroup[[1L]]
  )
  candidates <- evidence$candidates
  candidates$item <- score_set$candidate_item
  tibble::tibble(
    participant = source_row$participant,
    item = source_row$item,
    item_fold = source_row$item_fold,
    outer_fold = outer_fold,
    candidate_pool_id = source_row$candidate_pool_id,
    method = method,
    gaze_info_bits = evidence$gaze_info_bits,
    base_gaze_info_bits = evidence$base_gaze_info_bits,
    log_loss = evidence$log_loss,
    base_log_loss = evidence$base_log_loss,
    posterior_true = evidence$posterior_true,
    prior_true = evidence$prior_true,
    template_rank = evidence$template_rank,
    top1_credit = evidence$top1_credit,
    candidate_count = evidence$candidate_count,
    temperature = evidence$temperature,
    reliability = 1,
    raw_fixation_count = quality$raw_fixation_count,
    coalesced_fixation_count = quality$coalesced_fixation_count,
    effective_fixations = quality$effective_fixations,
    total_duration = quality$total_duration,
    duration_concentration = quality$duration_concentration,
    spatial_dispersion = quality$spatial_dispersion,
    replay_coverage = NA_real_,
    background_coverage = NA_real_,
    expected_restarts = NA_real_,
    restart_rate = NA_real_,
    spatial_rmse_px = NA_real_,
    template_count = 4,
    template_effective_count = NA_real_,
    template_weight_max = NA_real_,
    modal_template = NA_character_,
    converged = TRUE,
    candidates = list(multi_replay_slim_candidates(candidates))
  )
}

multi_replay_score_density_fold <- function(subset, fold, spec) {
  profile <- real_profile("density", "combined", fold$id, {
    train_episode <- multi_replay_density_episode_table(
      subset$ref_train, spec
    )
    eval_episode <- multi_replay_density_episode_table(
      subset$ref_eval, spec
    )
    train_signature <- multi_replay_density_source_signatures(
      subset$source_train, spec
    )
    eval_signature <- multi_replay_density_source_signatures(
      subset$source_eval, spec
    )
    fits <- vector("list", length(spec$density_sigmas))
    evaluated <- vector("list", length(spec$density_sigmas))
    for (sigma_index in seq_along(spec$density_sigmas)) {
      train_sets <- lapply(seq_len(nrow(subset$source_train)), function(i) {
        multi_replay_density_score_set(
          subset$source_train[i, , drop = FALSE],
          train_signature[[i]], train_episode, sigma_index
        )
      })
      fits[[sigma_index]] <- real_ns("fit_gaze_temperature")(
        lapply(train_sets, `[[`, "score"),
        vapply(train_sets, `[[`, integer(1), "true_index"),
        bounds = c(0.05, 100)
      )
      evaluated[[sigma_index]] <- lapply(
        seq_len(nrow(subset$source_eval)), function(i) {
          source_row <- subset$source_eval[i, , drop = FALSE]
          candidates <- eval_episode[
            eval_episode$participant == source_row$participant[[1L]] &
              eval_episode$item_fold == source_row$item_fold[[1L]],
            , drop = FALSE
          ]
          score_set <- multi_replay_density_score_set(
            source_row, eval_signature[[i]], candidates, sigma_index
          )
          evidence <- real_ns("score_gaze_candidates")(
            score_set$score,
            true_index = score_set$true_index,
            candidate_key = score_set$candidate_item,
            temperature = fits[[sigma_index]]$temperature,
            candidate_pool_id = source_row$candidate_pool_id[[1L]]
          )
          list(score_set = score_set, evidence = evidence)
        }
      )
    }
    list(fits = fits, evaluated = evaluated)
  })
  rows <- list()
  index <- 1L
  for (sigma_index in seq_along(spec$density_sigmas)) {
    method <- paste0(
      "density_sigma_", format(spec$density_sigmas[[sigma_index]], trim = TRUE),
      "_all4_exhaustive"
    )
    for (i in seq_along(profile$value$evaluated[[sigma_index]])) {
      value <- profile$value$evaluated[[sigma_index]][[i]]
      rows[[index]] <- multi_replay_density_row(
        subset$source_eval[i, , drop = FALSE],
        value$score_set, value$evidence, method, fold$id
      )
      index <- index + 1L
    }
  }
  calibration <- dplyr::bind_rows(lapply(
    seq_along(spec$density_sigmas), function(i) {
      fit <- profile$value$fits[[i]]
      data.frame(
        outer_fold = fold$id,
        method_group = paste0("density_sigma_", spec$density_sigmas[[i]]),
        temperature = fit$temperature,
        reliability_kappa = 0,
        calibration_log_loss = fit$log_loss,
        temperature_only_log_loss = fit$log_loss,
        shrinkage_not_worse = TRUE,
        stringsAsFactors = FALSE
      )
    }
  ))
  list(
    scored = dplyr::bind_rows(rows),
    components = tibble::tibble(),
    resources = profile$resource,
    calibration = calibration,
    model_audit = profile$value$fits
  )
}

multi_replay_score_fold_method <- function(tables, candidate_plan, fold,
                                           specs, method_group) {
  subset <- multi_replay_fold_subsets(tables, candidate_plan, fold)
  scored <- if (identical(method_group, "p4")) {
    multi_replay_score_replay_fold(
      subset, fold, specs$replay_p4, "p4"
    )
  } else if (identical(method_group, "all4")) {
    multi_replay_score_replay_fold(
      subset, fold, specs$replay_all4, "all4"
    )
  } else if (identical(method_group, "density")) {
    multi_replay_score_density_fold(subset, fold, specs$density)
  } else if (identical(method_group, "study_control")) {
    multi_replay_score_study_control_fold(
      subset, fold, specs$replay_all4
    )
  } else {
    stop("Unknown multi-presentation method group: ", method_group)
  }
  scored$audit <- data.frame(
    outer_fold = fold$id,
    method_group = method_group,
    train_trials = nrow(subset$source_train),
    eval_trials = nrow(subset$source_eval),
    train_reference_paths = nrow(subset$ref_train),
    eval_reference_paths = nrow(subset$ref_eval),
    train_participants = length(unique(subset$source_train$participant)),
    eval_participants = length(unique(subset$source_eval$participant)),
    train_items = length(unique(subset$source_train$item)),
    eval_items = length(unique(subset$source_eval$item)),
    participant_overlap = length(intersect(
      unique(subset$source_train$participant),
      unique(subset$source_eval$participant)
    )),
    item_overlap = length(intersect(
      unique(subset$source_train$item), unique(subset$source_eval$item)
    )),
    stringsAsFactors = FALSE
  )
  scored$protocol_version <- multi_replay_protocol_version
  scored$seed <- multi_replay_seed
  scored$method_group <- method_group
  scored$outer_fold <- fold$id
  scored
}

multi_replay_checkpoint <- function(output_dir, method_group, fold_id,
                                    smoke = FALSE) {
  suffix <- if (smoke) "-smoke" else ""
  file.path(
    output_dir,
    sprintf(
      "checkpoint-%s-fold-%02d%s.rds", method_group, fold_id, suffix
    )
  )
}

multi_replay_collect_checkpoints <- function(output_dir, method_groups,
                                             fold_ids, smoke = FALSE) {
  paths <- unlist(lapply(fold_ids, function(fold_id) {
    vapply(method_groups, function(method_group) {
      multi_replay_checkpoint(
        output_dir, method_group, fold_id, smoke = smoke
      )
    }, character(1))
  }), use.names = FALSE)
  if (!all(file.exists(paths))) {
    stop("Finalization requires every selected method-fold checkpoint.")
  }
  checkpoints <- lapply(paths, readRDS)
  valid <- vapply(checkpoints, function(value) {
    identical(value$protocol_version, multi_replay_protocol_version) &&
      identical(value$seed, multi_replay_seed)
  }, logical(1))
  if (!all(valid)) stop("A checkpoint belongs to another protocol or seed.")
  checkpoints
}

multi_replay_crossed_bootstrap <- function(value, participant, item,
                                           draws, seed) {
  keep <- is.finite(value) & !is.na(participant) & !is.na(item)
  value <- value[keep]
  participant <- as.character(participant[keep])
  item <- as.character(item[keep])
  participant_levels <- sort(unique(participant))
  item_levels <- sort(unique(item))
  set.seed(seed)
  result <- numeric(draws)
  for (draw in seq_len(draws)) {
    participant_count <- table(factor(
      sample(participant_levels, replace = TRUE),
      levels = participant_levels
    ))
    item_count <- table(factor(
      sample(item_levels, replace = TRUE), levels = item_levels
    ))
    weight <- as.numeric(participant_count[participant]) *
      as.numeric(item_count[item])
    result[[draw]] <- if (sum(weight) > 0) {
      stats::weighted.mean(value, weight)
    } else {
      NA_real_
    }
  }
  result
}

multi_replay_interval <- function(value, participant, item, draws, seed) {
  bootstrap <- multi_replay_crossed_bootstrap(
    value, participant, item, draws, seed
  )
  data.frame(
    estimate = mean(value, na.rm = TRUE),
    bootstrap_se = stats::sd(bootstrap, na.rm = TRUE),
    lower_95 = stats::quantile(
      bootstrap, 0.025, na.rm = TRUE, names = FALSE
    ),
    upper_95 = stats::quantile(
      bootstrap, 0.975, na.rm = TRUE, names = FALSE
    ),
    bootstrap_p = full_recognition_bootstrap_p(bootstrap),
    stringsAsFactors = FALSE
  )
}

multi_replay_method_summary <- function(scored, draws, seed) {
  rows <- lapply(seq_along(multi_replay_methods), function(index) {
    method <- multi_replay_methods[[index]]
    tab <- scored[scored$method == method, , drop = FALSE]
    if (nrow(tab) == 0L) return(NULL)
    interval <- multi_replay_interval(
      tab$gaze_info_bits, tab$participant, tab$item,
      draws, seed + index * 1009L
    )
    data.frame(
      method = method,
      trials = nrow(tab),
      mean_gaze_info_bits = interval$estimate,
      lower_95_bits = interval$lower_95,
      upper_95_bits = interval$upper_95,
      bits_bootstrap_p = interval$bootstrap_p,
      mean_log_loss = mean(tab$log_loss),
      prior_log_loss = mean(-log(tab$prior_true)),
      mean_top1_credit = mean(tab$top1_credit),
      mean_template_rank = mean(tab$template_rank),
      mean_candidate_count = mean(tab$candidate_count),
      convergence_rate = mean(tab$converged),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

multi_replay_paired_comparison <- function(scored, first, second,
                                           label, draws, seed) {
  first_tab <- scored[scored$method == first, c(
    "participant", "item", "gaze_info_bits", "log_loss"
  )]
  second_tab <- scored[scored$method == second, c(
    "participant", "item", "gaze_info_bits", "log_loss"
  )]
  names(first_tab)[3:4] <- c("first_bits", "first_loss")
  names(second_tab)[3:4] <- c("second_bits", "second_loss")
  paired <- merge(
    first_tab, second_tab,
    by = c("participant", "item"), all = FALSE, sort = FALSE
  )
  if (nrow(paired) != nrow(first_tab) || nrow(paired) != nrow(second_tab)) {
    stop("Paired measurement methods do not share trial support: ", label)
  }
  bits <- multi_replay_interval(
    paired$second_bits - paired$first_bits,
    paired$participant, paired$item, draws, seed
  )
  loss <- multi_replay_interval(
    paired$first_loss - paired$second_loss,
    paired$participant, paired$item, draws, seed + 1L
  )
  data.frame(
    comparison = label,
    first = first,
    second = second,
    trials = nrow(paired),
    bits_gain = bits$estimate,
    bits_gain_lower_95 = bits$lower_95,
    bits_gain_upper_95 = bits$upper_95,
    bits_gain_p = bits$bootstrap_p,
    log_loss_reduction = loss$estimate,
    log_loss_reduction_lower_95 = loss$lower_95,
    log_loss_reduction_upper_95 = loss$upper_95,
    log_loss_reduction_p = loss$bootstrap_p,
    stringsAsFactors = FALSE
  )
}

multi_replay_fixation_bands <- function(scored) {
  tab <- scored[scored$method %in% c(
    "replay_all4_unshrunk_exhaustive",
    "replay_all4_shrunk_exhaustive"
  ), , drop = FALSE]
  tab$fixation_band <- cut(
    tab$effective_fixations,
    breaks = c(1, 2, 3, 5, Inf),
    right = FALSE,
    include.lowest = TRUE,
    labels = c("[1,2)", "[2,3)", "[3,5)", "5+")
  )
  dplyr::bind_rows(lapply(split(
    seq_len(nrow(tab)), interaction(tab$method, tab$fixation_band, drop = TRUE)
  ), function(index) {
    part <- tab[index, , drop = FALSE]
    data.frame(
      method = part$method[[1L]],
      fixation_band = as.character(part$fixation_band[[1L]]),
      trials = nrow(part),
      mean_effective_fixations = mean(part$effective_fixations),
      mean_reliability = mean(part$reliability),
      mean_gaze_info_bits = mean(part$gaze_info_bits),
      mean_log_loss = mean(part$log_loss),
      stringsAsFactors = FALSE
    )
  }))
}

multi_replay_weighted_slope <- function(x, y, weight) {
  keep <- is.finite(x) & is.finite(y) & is.finite(weight) & weight > 0
  x <- x[keep]
  y <- y[keep]
  weight <- weight[keep]
  weight <- weight / sum(weight)
  x_centered <- x - sum(weight * x)
  y_centered <- y - sum(weight * y)
  denominator <- sum(weight * x_centered^2)
  if (denominator <= 0) return(NA_real_)
  sum(weight * x_centered * y_centered) / denominator
}

multi_replay_fixation_association <- function(scored, draws, seed) {
  tab <- scored[
    scored$method == "replay_all4_unshrunk_exhaustive", , drop = FALSE
  ]
  x <- as.numeric(scale(log(tab$effective_fixations)))
  y <- as.numeric(scale(tab$gaze_info_bits))
  participant <- as.character(tab$participant)
  item <- as.character(tab$item)
  participant_levels <- sort(unique(participant))
  item_levels <- sort(unique(item))
  set.seed(seed)
  bootstrap <- numeric(draws)
  for (draw in seq_len(draws)) {
    participant_count <- table(factor(
      sample(participant_levels, replace = TRUE),
      levels = participant_levels
    ))
    item_count <- table(factor(
      sample(item_levels, replace = TRUE), levels = item_levels
    ))
    weight <- as.numeric(participant_count[participant]) *
      as.numeric(item_count[item])
    bootstrap[[draw]] <- multi_replay_weighted_slope(x, y, weight)
  }
  data.frame(
    predictor = "log_effective_fixations",
    outcome = "unshrunk_gaze_info_bits",
    standardized_slope = multi_replay_weighted_slope(
      x, y, rep(1, length(x))
    ),
    lower_95 = stats::quantile(
      bootstrap, 0.025, na.rm = TRUE, names = FALSE
    ),
    upper_95 = stats::quantile(
      bootstrap, 0.975, na.rm = TRUE, names = FALSE
    ),
    bootstrap_p = full_recognition_bootstrap_p(bootstrap),
    spearman_rho = stats::cor(
      tab$effective_fixations, tab$gaze_info_bits,
      method = "spearman", use = "complete.obs"
    ),
    stringsAsFactors = FALSE
  )
}

multi_replay_presentation_summary <- function(components) {
  tab <- components[
    components$method == "replay_all4_shrunk_exhaustive" &
      components$is_true_candidate, , drop = FALSE
  ]
  presentation <- suppressWarnings(as.integer(tab$template_key))
  if (anyNA(presentation)) presentation <- tab$template_key
  summary <- dplyr::bind_rows(lapply(split(
    seq_len(nrow(tab)), presentation
  ), function(index) {
    part <- tab[index, , drop = FALSE]
    data.frame(
      presentation = as.character(part$template_key[[1L]]),
      trials = nrow(part),
      mean_posterior_responsibility = mean(part$template_posterior),
      mean_effective_fixations = mean(part$template_effective_fixations),
      mean_raw_fixation_count = mean(part$template_raw_fixation_count),
      mean_total_duration = mean(part$template_total_duration),
      stringsAsFactors = FALSE
    )
  }))
  episode <- interaction(
    tab$participant, tab$target_item, drop = TRUE, lex.order = TRUE
  )
  centered_effective <- tab$template_effective_fixations - ave(
    tab$template_effective_fixations, episode
  )
  centered_raw <- tab$template_raw_fixation_count - ave(
    tab$template_raw_fixation_count, episode
  )
  centered_posterior <- tab$template_posterior - ave(
    tab$template_posterior, episode
  )
  association <- data.frame(
    quantity = "within-episode presentation responsibility",
    within_episode_spearman_effective_fixations = stats::cor(
      centered_effective,
      centered_posterior,
      method = "spearman", use = "complete.obs"
    ),
    within_episode_spearman_raw_fixation_count = stats::cor(
      centered_raw,
      centered_posterior,
      method = "spearman", use = "complete.obs"
    ),
    stringsAsFactors = FALSE
  )
  list(summary = summary, association = association)
}

multi_replay_panel_variance <- function(scored, panel_count = 8L) {
  tab <- scored[
    scored$method == "replay_all4_shrunk_exhaustive", , drop = FALSE
  ]
  rows <- lapply(seq_len(nrow(tab)), function(i) {
    candidates <- tab$candidates[[i]]
    true_index <- which(candidates$is_true)
    nonmatch <- setdiff(seq_len(nrow(candidates)), true_index)
    panel_n <- min(panel_count, length(nonmatch))
    if (length(nonmatch) < 4L || panel_n < 1L) return(NULL)
    nonmatch <- nonmatch[order(candidates$candidate_key[nonmatch])]
    panel_bits <- vapply(seq_len(panel_n), function(panel) {
      position <- ((panel - 1L + 0:3) %% length(nonmatch)) + 1L
      selected <- c(true_index, nonmatch[position])
      evidence <- real_ns("score_gaze_candidates")(
        candidates$log_score[selected],
        true_index = 1L,
        candidate_key = candidates$candidate_key[selected],
        temperature = tab$temperature[[i]],
        reliability = tab$reliability[[i]],
        candidate_pool_id = paste0(tab$candidate_pool_id[[i]], ":panel:", panel)
      )
      evidence$gaze_info_bits
    }, numeric(1))
    data.frame(
      participant = tab$participant[[i]],
      item = tab$item[[i]],
      candidate_count = tab$candidate_count[[i]],
      exhaustive_bits = tab$gaze_info_bits[[i]],
      panel_count = panel_n,
      panel_mean_bits = mean(panel_bits),
      panel_sd_bits = stats::sd(panel_bits),
      panel_min_bits = min(panel_bits),
      panel_max_bits = max(panel_bits),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

multi_replay_generic_quality_summary <- function(scored) {
  tab <- scored[
    scored$method == "replay_all4_shrunk_exhaustive", , drop = FALSE
  ]
  quantities <- c(
    "raw_fixation_count", "coalesced_fixation_count",
    "effective_fixations", "total_duration", "duration_concentration",
    "spatial_dispersion", "reliability", "replay_coverage",
    "background_coverage", "spatial_rmse_px", "template_effective_count"
  )
  dplyr::bind_rows(lapply(quantities, function(quantity) {
    value <- tab[[quantity]]
    data.frame(
      quantity = quantity,
      observations = sum(is.finite(value)),
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      q05 = stats::quantile(value, 0.05, na.rm = TRUE, names = FALSE),
      median = stats::median(value, na.rm = TRUE),
      q95 = stats::quantile(value, 0.95, na.rm = TRUE, names = FALSE),
      stringsAsFactors = FALSE
    )
  }))
}

multi_replay_abstention_gate <- function(scored, tolerance = 1e-12) {
  tab <- scored[
    scored$method == "replay_all4_shrunk_exhaustive", , drop = FALSE
  ]
  valid <- vapply(seq_len(nrow(tab)), function(i) {
    candidate <- tab$candidates[[i]]
    all(abs(candidate$posterior - candidate$prior) <=
          abs(candidate$base_posterior - candidate$prior) + tolerance)
  }, logical(1))
  all(valid)
}

multi_replay_analyse <- function(scored, components, calibration, audit,
                                 draws = multi_replay_bootstrap_draws,
                                 seed = multi_replay_seed,
                                 expected_trials = 2055L) {
  forbidden <- intersect(
    c("accuracy", "probe_type", "degradation", "response"), names(scored)
  )
  if (length(forbidden) > 0L) {
    stop("Measurement analysis contains behavior-bearing columns: ",
         paste(forbidden, collapse = ", "))
  }
  support <- table(scored$method)
  if (!all(multi_replay_methods %in% names(support)) ||
      any(support[multi_replay_methods] != expected_trials)) {
    stop("Full measurement methods do not share the frozen trial support.")
  }
  if (any(!is.finite(scored$gaze_info_bits)) ||
      any(!is.finite(scored$log_loss)) ||
      any(!scored$converged)) {
    stop("A measurement method violates the finite-score contract.")
  }
  summary <- multi_replay_method_summary(scored, draws, seed)
  comparisons <- dplyr::bind_rows(
    multi_replay_paired_comparison(
      scored,
      "replay_p4_exhaustive",
      "replay_all4_unshrunk_exhaustive",
      "all4_unshrunk_minus_p4", draws, seed + 10000L
    ),
    multi_replay_paired_comparison(
      scored,
      "replay_all4_unshrunk_exhaustive",
      "replay_all4_shrunk_exhaustive",
      "reliability_shrinkage_minus_unshrunk", draws, seed + 20000L
    ),
    multi_replay_paired_comparison(
      scored,
      "study_replay_reversed",
      "study_replay_intact",
      "study_intact_minus_reversed", draws, seed + 30000L
    )
  )
  fixation_bands <- multi_replay_fixation_bands(scored)
  fixation_association <- multi_replay_fixation_association(
    scored, draws, seed + 40000L
  )
  presentation <- multi_replay_presentation_summary(components)
  panels <- multi_replay_panel_variance(scored)
  generic_quality <- multi_replay_generic_quality_summary(scored)
  study_summary <- summary[summary$method %in% c(
    "study_replay_intact", "study_replay_reversed"
  ), ]
  study_comparison <- comparisons[
    comparisons$comparison == "study_intact_minus_reversed", , drop = FALSE
  ]
  integrity <- all(audit$participant_overlap == 0L) &&
    all(audit$item_overlap == 0L) &&
    all(scored$candidate_count >= 2L) &&
    all(scored$candidate_count <= 33L)
  calibration_ok <- all(calibration$shrinkage_not_worse)
  abstention_ok <- multi_replay_abstention_gate(scored)
  positive_control_ok <- study_summary$lower_95_bits[
    study_summary$method == "study_replay_intact"
  ] > 0
  reversal_ok <- study_comparison$bits_gain_lower_95 > 0
  gates <- data.frame(
    gate = c(
      "score_and_fold_integrity",
      "calibration_not_worse",
      "abstention_toward_prior",
      "repeated_viewing_positive",
      "complete_reversal_loses_information",
      "measurement_reported_response_blind"
    ),
    passed = c(
      integrity,
      calibration_ok,
      abstention_ok,
      positive_control_ok,
      reversal_ok,
      TRUE
    ),
    stringsAsFactors = FALSE
  )
  verdict <- data.frame(
    measurement_gates_pass = all(gates$passed),
    freeze_for_behavior = all(gates$passed),
    behavior_tested = FALSE,
    protocol_version = multi_replay_protocol_version,
    stringsAsFactors = FALSE
  )
  list(
    method_summary = summary,
    paired_comparisons = comparisons,
    fixation_bands = fixation_bands,
    fixation_association = fixation_association,
    presentation_summary = presentation$summary,
    presentation_association = presentation$association,
    panel_variance = panels,
    generic_quality = generic_quality,
    gates = gates,
    verdict = verdict
  )
}

multi_replay_write_results <- function(result, output_dir) {
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  analysis <- result$analysis
  write(analysis$method_summary, "measurement-summary.csv")
  write(analysis$paired_comparisons, "paired-comparisons.csv")
  write(analysis$fixation_bands, "fixation-bands.csv")
  write(analysis$fixation_association, "fixation-association.csv")
  write(analysis$presentation_summary, "presentation-summary.csv")
  write(analysis$presentation_association, "presentation-association.csv")
  write(analysis$panel_variance, "candidate-panel-variance.csv")
  write(analysis$generic_quality, "generic-gaze-quality.csv")
  write(analysis$gates, "measurement-gates.csv")
  write(analysis$verdict, "measurement-verdict.csv")
  write(result$calibration, "calibration-audit.csv")
  write(result$audit, "fold-audit.csv")
  write(result$resources, "resources.csv")
  configuration <- data.frame(
    protocol_version = multi_replay_protocol_version,
    seed = multi_replay_seed,
    item_seed = multi_replay_item_seed,
    bootstrap_draws = multi_replay_bootstrap_draws,
    retained_participants = multi_replay_expected[["retained_participants"]],
    retained_pairs = multi_replay_expected[["retained_pairs"]],
    candidate_policy = "exhaustive_within_participant_item_fold",
    study_template_policy = "equal_mixture_separate_presentations",
    recognition_window_ms = "[0,3000)",
    behavior_tested = FALSE,
    local_only = TRUE,
    stringsAsFactors = FALSE
  )
  write(configuration, "configuration.csv")
  saveRDS(
    list(
      analysis = analysis,
      calibration = result$calibration,
      audit = result$audit,
      resources = result$resources,
      configuration = configuration
    ),
    file.path(output_dir, "scientific-results.rds"), version = 3
  )
  utils::capture.output(
    utils::sessionInfo(), file = file.path(output_dir, "session-info.txt")
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[
    !grepl("^checkpoint-", basename(files)) &
      basename(files) != "manifest-md5.csv"
  ]
  manifest <- data.frame(
    file = basename(files),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  write(manifest, "manifest-md5.csv")
  invisible(result)
}

run_gaze_weave_pcmri_multi_replay <- function(
    output_dir = file.path(
      "inst", "validation",
      "gaze-weave-pcmri-catalog-replication-results", "multi-replay"
    ),
    smoke = FALSE,
    smoke_participants = 4L,
    method_groups = multi_replay_method_groups,
    fold_ids = NULL,
    resume = TRUE,
    finalize = !smoke) {
  method_groups <- match.arg(
    method_groups, multi_replay_method_groups, several.ok = TRUE
  )
  raw <- full_recognition_read_inputs(verify = !smoke)
  cohort <- multi_replay_select_cohort(raw)
  if (smoke) {
    cohort <- multi_replay_smoke_cohort(cohort, smoke_participants)
  }
  candidate_plan <- multi_replay_candidate_plan(cohort)
  fold_plan <- full_recognition_fold_plan(cohort, multi_replay_seed)
  design <- multi_replay_validate_design(
    cohort, candidate_plan, fold_plan, strict = !smoke
  )
  tables <- multi_replay_task_tables(raw, cohort)
  specs <- multi_replay_specs(smoke)
  folds <- fold_plan$folds
  if (!is.null(fold_ids)) {
    folds <- folds[vapply(
      folds, function(fold) fold$id %in% fold_ids, logical(1)
    )]
  }
  if (length(folds) == 0L) stop("fold_ids selected no outer folds.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (fold in folds) {
    for (method_group in method_groups) {
      path <- multi_replay_checkpoint(
        output_dir, method_group, fold$id, smoke = smoke
      )
      if (resume && file.exists(path)) {
        checkpoint <- readRDS(path)
        if (!identical(
          checkpoint$protocol_version, multi_replay_protocol_version
        ) || !identical(checkpoint$seed, multi_replay_seed)) {
          stop("A checkpoint uses another protocol version: ", path)
        }
        message("Using checkpoint: ", basename(path))
        next
      }
      message(
        "Multi-presentation scoring: fold ", fold$id,
        ", method ", method_group
      )
      checkpoint <- multi_replay_score_fold_method(
        tables, candidate_plan, fold, specs, method_group
      )
      saveRDS(checkpoint, path, version = 3)
    }
  }
  if (!finalize) {
    return(invisible(list(
      cohort = cohort,
      candidate_plan = candidate_plan,
      fold_plan = fold_plan,
      design = design,
      output_dir = output_dir
    )))
  }
  selected_folds <- vapply(folds, `[[`, integer(1), "id")
  if (!smoke && (!setequal(method_groups, multi_replay_method_groups) ||
                 !identical(sort(selected_folds), 1:4))) {
    stop("Full finalization requires all frozen method groups and outer folds.")
  }
  checkpoints <- multi_replay_collect_checkpoints(
    output_dir, method_groups, selected_folds, smoke = smoke
  )
  result <- list(
    scored = dplyr::bind_rows(lapply(checkpoints, `[[`, "scored")),
    components = dplyr::bind_rows(lapply(checkpoints, `[[`, "components")),
    resources = dplyr::bind_rows(lapply(checkpoints, `[[`, "resources")),
    calibration = dplyr::bind_rows(lapply(checkpoints, `[[`, "calibration")),
    audit = dplyr::bind_rows(lapply(checkpoints, `[[`, "audit")),
    model_audit = lapply(checkpoints, `[[`, "model_audit"),
    config = list(
      cohort = cohort,
      candidate_plan = candidate_plan,
      fold_plan = fold_plan,
      design = design,
      smoke = smoke,
      method_groups = method_groups,
      fold_ids = selected_folds,
      seed = multi_replay_seed
    )
  )
  if (!smoke) {
    result$analysis <- multi_replay_analyse(
      result$scored,
      result$components,
      result$calibration,
      result$audit,
      draws = multi_replay_bootstrap_draws,
      seed = multi_replay_seed
    )
    multi_replay_write_results(result, output_dir)
  }
  result
}
