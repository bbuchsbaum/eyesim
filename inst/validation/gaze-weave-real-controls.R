# Frozen real-data court for GazeWeave v2.
#
# Protocol: inst/validation/GAZEWEAVE-REAL-CONTROLS.md
# Run from a source checkout after devtools::load_all():
#   source("inst/validation/gaze-weave-real-controls.R")
#   run_gaze_weave_real_court("inst/validation/gaze-weave-real-results")

real_ns <- function(name) getFromNamespace(name, "eyesim")

real_path <- function(rows, order_column) {
  rows <- rows[order(rows[[order_column]], seq_len(nrow(rows))), , drop = FALSE]
  fixation_group(
    x = rows$FixX - 112,
    y = rows$FixY - 84,
    duration = rows$FixDuration,
    onset = seq_len(nrow(rows)) - 1L
  )
}

real_split_paths <- function(tab, keys, order_column) {
  key <- do.call(interaction, c(tab[keys], list(drop = TRUE, lex.order = TRUE)))
  groups <- split(seq_len(nrow(tab)), key)
  rows <- lapply(groups, function(index) {
    part <- tab[index, , drop = FALSE]
    metadata <- part[1L, keys, drop = FALSE]
    metadata$fixgroup <- list(real_path(part, order_column))
    metadata
  })
  dplyr::bind_rows(rows)
}

real_read_wynn <- function(blank_window = c(0, 3000)) {
  study_file <- file.path("test_data", "study_fixations_all.csv")
  recall_file <- file.path("test_data", "testdelay_fixations.csv")
  if (!file.exists(study_file) || !file.exists(recall_file)) {
    stop("Run the real court from the eyesim source root; Wynn CSV files are missing.")
  }
  study <- utils::read.csv(study_file, stringsAsFactors = FALSE)
  recall <- utils::read.csv(recall_file, stringsAsFactors = FALSE)
  shared <- intersect(unique(study$Subject), unique(recall$Subject))
  valid_spatial <- function(tab) {
    is.finite(tab$FixX) & is.finite(tab$FixY) &
      tab$FixX >= 112 & tab$FixX <= 912 &
      tab$FixY >= 84 & tab$FixY <= 684
  }
  study <- study[
    study$Subject %in% shared & study$Image != "." &
      valid_spatial(study) & is.finite(study$FixDuration) &
      study$FixDuration > 80,
    , drop = FALSE
  ]
  recall <- recall[
    recall$Subject %in% shared & recall$Image != "." &
      valid_spatial(recall) & is.finite(recall$FixDuration) &
      recall$FixDuration > 80 & is.finite(recall$FixOffset),
    , drop = FALSE
  ]

  trial_key <- interaction(
    study$Subject, study$ImageNumber, drop = TRUE, lex.order = TRUE
  )
  study$repetition <- ave(study$Trial, trial_key, FUN = function(value) {
    match(value, sort(unique(value)))
  })
  trial_count <- ave(study$Trial, trial_key, FUN = function(value) {
    length(unique(value))
  })
  study <- study[trial_count == 4L & study$repetition %in% c(1L, 4L), ]

  if (!is.null(blank_window)) {
    stopifnot(length(blank_window) == 2L, blank_window[[1]] < blank_window[[2]])
    recall <- recall[
      recall$FixOffset >= blank_window[[1]] &
        recall$FixOffset < blank_window[[2]],
      , drop = FALSE
    ]
    recall$FixDuration <- pmin(
      recall$FixDuration,
      blank_window[[2]] - recall$FixOffset
    )
    recall <- recall[recall$FixDuration > 80, , drop = FALSE]
  }

  study_paths <- real_split_paths(
    study, c("Subject", "Age", "ImageNumber", "repetition"), "FixStartTime"
  )
  recall_paths <- real_split_paths(
    recall,
    c(
      "Subject", "Age", "ImageNumber", "ImageRepetition", "Saliency",
      "Duration", "Accuracy"
    ),
    "FixOffset"
  )
  names(study_paths)[names(study_paths) == "Subject"] <- "participant"
  names(study_paths)[names(study_paths) == "Age"] <- "age"
  names(study_paths)[names(study_paths) == "ImageNumber"] <- "item"
  names(recall_paths)[names(recall_paths) == "Subject"] <- "participant"
  names(recall_paths)[names(recall_paths) == "Age"] <- "age"
  names(recall_paths)[names(recall_paths) == "ImageNumber"] <- "item"
  names(recall_paths)[names(recall_paths) == "ImageRepetition"] <- "probe_type"
  names(recall_paths)[names(recall_paths) == "Saliency"] <- "degradation"
  names(recall_paths)[names(recall_paths) == "Duration"] <- "cue_duration"
  names(recall_paths)[names(recall_paths) == "Accuracy"] <- "accuracy"
  recall_paths$probe_type[recall_paths$probe_type == "new"] <- "lure"
  study_paths$participant <- as.character(study_paths$participant)
  recall_paths$participant <- as.character(recall_paths$participant)
  study_paths$item <- as.integer(study_paths$item)
  recall_paths$item <- as.integer(recall_paths$item)
  list(study = study_paths, recall = recall_paths)
}

real_select_cohort <- function(raw, seed = 20260817L,
                               n_participants_per_age = 4L,
                               n_items = 8L) {
  study <- raw$study
  recall <- raw$recall
  study$nfix <- vapply(study$fixgroup, nrow, integer(1))
  recall$nfix <- vapply(recall$fixgroup, nrow, integer(1))
  first <- study[study$repetition == 1L & study$nfix >= 3L,
                 c("participant", "age", "item")]
  fourth <- study[study$repetition == 4L & study$nfix >= 3L,
                  c("participant", "age", "item")]
  retrieval <- recall[recall$nfix >= 3L, c("participant", "age", "item")]
  complete <- merge(merge(first, fourth, by = c("participant", "age", "item")),
                    retrieval, by = c("participant", "age", "item"))
  counts <- aggregate(item ~ participant + age, complete, function(x) {
    length(unique(x))
  })
  eligible <- counts[counts$item >= 100L, , drop = FALSE]
  if (length(unique(eligible$age)) < 2L) {
    stop("The frozen age-stratified cohort cannot be formed.")
  }
  set.seed(seed)
  chosen_participants <- unlist(lapply(sort(unique(eligible$age)), function(group) {
    ids <- sort(eligible$participant[eligible$age == group])
    if (length(ids) < n_participants_per_age) {
      stop("Too few eligible participants in age group ", group, ".")
    }
    sample(ids, n_participants_per_age)
  }), use.names = FALSE)
  complete_selected <- complete[complete$participant %in% chosen_participants, ]
  item_counts <- table(complete_selected$item)
  eligible_items <- as.integer(names(item_counts)[
    item_counts == length(chosen_participants)
  ])
  if (length(eligible_items) < n_items) {
    stop("Too few items are complete for the frozen participant sample.")
  }
  chosen_items <- sample(sort(eligible_items), n_items)
  list(
    participants = sort(as.character(chosen_participants)),
    items = sort(as.integer(chosen_items)),
    eligibility = list(
      participants_per_age = table(eligible$age),
      common_item_n = length(eligible_items)
    )
  )
}

real_generic_path <- function(paths) {
  coords <- do.call(rbind, lapply(paths, function(path) cbind(path$x, path$y)))
  weight <- unlist(lapply(paths, `[[`, "duration"), use.names = FALSE)
  weight <- weight / sum(weight)
  center <- colSums(coords * weight)
  centered <- sweep(coords, 2L, center, "-")
  covariance <- t(centered * weight) %*% centered
  eig <- eigen(covariance, symmetric = TRUE)
  transform <- eig$vectors %*% diag(sqrt(pmax(eig$values, 1)), 2L)
  angle <- seq(0, 2 * pi, length.out = 9L)[-9L]
  standard <- sqrt(2) * cbind(cos(angle), sin(angle))
  generic <- sweep(standard %*% t(transform), 2L, center, "+")
  generic[, 1] <- pmin(pmax(generic[, 1], 0), 800)
  generic[, 2] <- pmin(pmax(generic[, 2], 0), 600)
  fixation_group(
    generic[, 1], generic[, 2],
    duration = rep(stats::median(weight) * 1000, nrow(generic)),
    onset = seq_len(nrow(generic)) - 1L
  )
}

real_shuffle_path <- function(path, seed) {
  set.seed(seed)
  order <- sample(seq_len(nrow(path)))
  fixation_group(
    path$x[order], path$y[order], path$duration[order],
    onset = seq_len(nrow(path)) - 1L
  )
}

real_task_tables <- function(raw, cohort, task, seed = 20260817L) {
  study <- raw$study[
    raw$study$participant %in% cohort$participants &
      raw$study$item %in% cohort$items,
    , drop = FALSE
  ]
  recall <- raw$recall[
    raw$recall$participant %in% cohort$participants &
      raw$recall$item %in% cohort$items,
    , drop = FALSE
  ]
  metadata <- recall[
    !duplicated(recall[c("participant", "item")]),
    c(
      "participant", "age", "item", "probe_type", "degradation",
      "cue_duration", "accuracy"
    )
  ]
  reference_rep <- if (task == "repeated_viewing") 1L else 4L
  source_rep <- 4L
  reference <- study[study$repetition == reference_rep,
                     c("participant", "age", "item", "fixgroup")]
  reference <- merge(reference, metadata, by = c("participant", "age", "item"),
                     all.x = TRUE, sort = FALSE)
  if (task == "repeated_viewing") {
    signal <- study[study$repetition == source_rep,
                    c("participant", "age", "item", "fixgroup")]
    signal <- merge(signal, metadata, by = c("participant", "age", "item"),
                    all.x = TRUE, sort = FALSE)
  } else if (task == "blank_screen_imagery") {
    signal <- recall
  } else {
    stop("Unknown real-control task: ", task)
  }
  reference$task <- task
  signal$task <- task
  reference <- reference[order(reference$participant, reference$item), ]
  signal <- signal[order(signal$participant, signal$item), ]
  rownames(reference) <- rownames(signal) <- NULL

  signal_lookup <- split(
    seq_len(nrow(signal)), interaction(signal$participant, signal$item, drop = TRUE)
  )
  generic <- lapply(split(signal$fixgroup, signal$participant), real_generic_path)
  item_position <- stats::setNames(seq_along(cohort$items), cohort$items)
  controls <- lapply(seq_len(nrow(signal)), function(i) {
    row <- signal[i, , drop = FALSE]
    participant <- row$participant[[1]]
    item <- row$item[[1]]
    next_item <- cohort$items[(item_position[[as.character(item)]] %%
                                length(cohort$items)) + 1L]
    wrong_key <- interaction(participant, next_item, drop = TRUE)
    wrong_row <- signal_lookup[[as.character(wrong_key)]]
    if (is.null(wrong_row) || length(wrong_row) != 1L) {
      stop("Could not construct the frozen wrong-item control.")
    }
    paths <- list(
      signal = row$fixgroup[[1]],
      shuffled_order = real_shuffle_path(
        row$fixgroup[[1]], seed + as.integer(participant) * 1009L + item * 101L +
          if (task == "blank_screen_imagery") 1L else 0L
      ),
      wrong_item = signal$fixgroup[[wrong_row]],
      generic_gaze = generic[[participant]]
    )
    dplyr::bind_rows(lapply(names(paths), function(condition) {
      out <- row
      out$condition <- condition
      out$fixgroup <- list(paths[[condition]])
      out
    }))
  })
  list(reference = tibble::as_tibble(reference),
       signal = tibble::as_tibble(signal),
       source = dplyr::bind_rows(controls))
}

real_fold_plan <- function(cohort, participant_age, seed = 20260817L) {
  set.seed(seed)
  participant_map <- do.call(rbind, lapply(sort(unique(participant_age$age)), function(group) {
    ids <- sample(sort(participant_age$participant[participant_age$age == group]))
    data.frame(participant = ids, participant_fold = rep(1:2, length.out = length(ids)))
  }))
  items <- sample(cohort$items)
  item_map <- data.frame(item = items, item_fold = rep(1:2, length.out = length(items)))
  folds <- list()
  index <- 1L
  for (participant_fold in 1:2) {
    for (item_fold in 1:2) {
      folds[[index]] <- list(
        id = index,
        participant_fold = participant_fold,
        item_fold = item_fold,
        eval_participants = participant_map$participant[
          participant_map$participant_fold == participant_fold
        ],
        eval_items = item_map$item[item_map$item_fold == item_fold]
      )
      index <- index + 1L
    }
  }
  list(folds = folds, participant_map = participant_map, item_map = item_map)
}

real_specs <- function(smoke = FALSE) {
  screen <- gaze_screen(800, 600, unit = "px")
  warp <- gaze_warp_contraction(
    center = "screen", translation = TRUE, fit_by = NULL
  )
  transport <- gaze_transport_v2_spec(
    spatial = gaze_gaussian_mixture(c(30, 60, 120), unit = "px"),
    chronology = gaze_order_neighbours(2, coalesce_distance = 2, unit = "px"),
    coverage_grid = if (smoke) c(0.5, 1) else c(0.5, 0.75, 1),
    coverage_penalty_grid = c(0.25, 0.75, 1.5),
    temporal_weight = 1,
    entropy_schedule = if (smoke) 0.03 else c(0.05, 0.015),
    warp = warp, screen = screen,
    maxit = if (smoke) 30L else 80L,
    projection_maxit = if (smoke) 150L else 350L,
    tolerance = 5e-4, projection_tolerance = 1e-6,
    multistart = 1L
  )
  replay <- gaze_replay_spec(
    grid_size = if (smoke) 24L else 48L,
    max_skip = 2L, student_df = 4, scale_floor = 6,
    transition_grid = list(
      background = c(0.03, 0.10), restart = c(0.02, 0.10),
      advance = c(0.25, 0.55), background_stay = c(0.85, 0.95)
    ),
    warp = warp, screen = screen
  )
  baseline <- gaze_baseline_spec(
    screen = screen, density_sigmas = c(30, 60, 120),
    density_grid = if (smoke) 12L else 24L,
    warp = warp, lambda_grid = c(0.01, 0.1, 1, 10), inner_folds = 2L,
    elastic_radii = c(consensus = 80, rigidity = 400, matching = 40),
    elastic_maxit = if (smoke) 10L else 25L,
    elastic_tolerance = 1e-4
  )
  list(screen = screen, warp = warp, transport = transport,
       replay = replay, baseline = baseline)
}

real_value <- function(x, name, default = NA_real_) {
  value <- x[[name]]
  if (is.null(value) || length(value) == 0L) default else value[[1L]]
}

real_profile <- function(label, task, fold, expression) {
  invisible(gc(reset = TRUE))
  elapsed <- system.time(value <- force(expression))[["elapsed"]]
  memory <- gc()
  list(
    value = value,
    resource = data.frame(
      method = label, task = task, outer_fold = fold,
      elapsed_seconds = unname(elapsed),
      gc_used_mb = sum(memory[, 2L]),
      gc_high_water_mb = sum(memory[, 7L]),
      object_size_mb = as.numeric(utils::object.size(value)) / 1024^2,
      stringsAsFactors = FALSE
    )
  )
}

real_engine_row <- function(source_row, scored, method, outer_fold) {
  evidence <- scored$evidence
  diagnostics <- scored$alignment$diagnostics
  tibble::tibble(
    task = source_row$task,
    condition = source_row$condition,
    participant = source_row$participant,
    age = source_row$age,
    item = source_row$item,
    probe_type = source_row$probe_type,
    degradation = source_row$degradation,
    cue_duration = source_row$cue_duration,
    accuracy = source_row$accuracy,
    outer_fold = outer_fold,
    method = method,
    calibrated = TRUE,
    status = "scored",
    compatibility_bits = evidence$gaze_info_bits,
    gaze_info_bits = evidence$gaze_info_bits,
    log_loss = evidence$log_loss,
    brier_score = evidence$brier_score,
    posterior_true = evidence$posterior_true,
    template_rank = evidence$template_rank,
    top1_credit = evidence$top1_credit,
    candidate_count = evidence$candidate_count,
    converged = scored$all_converged,
    replay_coverage = real_value(diagnostics, "replay_coverage"),
    background_coverage = real_value(diagnostics, "background_coverage"),
    spatial_rmse_px = real_value(diagnostics, "spatial_rmse"),
    local_order_error = real_value(diagnostics, "local_order_error"),
    expected_restarts = real_value(diagnostics, "expected_restarts"),
    candidates = list(evidence$candidates)
  )
}

real_baseline_row <- function(source_row, evidence, outer_fold) {
  tibble::tibble(
    task = source_row$task,
    condition = source_row$condition,
    participant = source_row$participant,
    age = source_row$age,
    item = source_row$item,
    probe_type = source_row$probe_type,
    degradation = source_row$degradation,
    cue_duration = source_row$cue_duration,
    accuracy = source_row$accuracy,
    outer_fold = outer_fold,
    method = evidence$method,
    calibrated = evidence$calibrated,
    status = evidence$status,
    compatibility_bits = evidence$compatibility_bits,
    gaze_info_bits = evidence$gaze_info_bits,
    log_loss = evidence$log_loss,
    brier_score = evidence$brier_score,
    posterior_true = evidence$posterior_true,
    template_rank = evidence$template_rank,
    top1_credit = evidence$top1_credit,
    candidate_count = evidence$candidate_count,
    converged = identical(evidence$status, "scored"),
    replay_coverage = NA_real_, background_coverage = NA_real_,
    spatial_rmse_px = NA_real_, local_order_error = NA_real_,
    expected_restarts = NA_real_,
    candidates = list(evidence$candidates)
  )
}

real_fold_subsets <- function(tables, fold) {
  eval_reference <- with(
    tables$reference,
    participant %in% fold$eval_participants & item %in% fold$eval_items
  )
  train_reference <- with(
    tables$reference,
    !participant %in% fold$eval_participants & !item %in% fold$eval_items
  )
  eval_source <- with(
    tables$source,
    participant %in% fold$eval_participants & item %in% fold$eval_items
  )
  train_source <- with(
    tables$signal,
    !participant %in% fold$eval_participants & !item %in% fold$eval_items
  )
  list(
    ref_train = tables$reference[train_reference, , drop = FALSE],
    source_train = tables$signal[train_source, , drop = FALSE],
    ref_eval = tables$reference[eval_reference, , drop = FALSE],
    source_eval = tables$source[eval_source, , drop = FALSE]
  )
}

real_warp_audit <- function(info, task, fold, method) {
  groups <- info$groups
  if (!is.list(groups) || length(groups) == 0L) return(data.frame())
  do.call(rbind, lapply(groups, function(group) {
    data.frame(
      task = task, outer_fold = fold, method = method,
      group = group$group, matched_pairs = group$matched_pairs,
      scale = group$scale,
      translation_x = group$translation[[1]],
      translation_y = group$translation[[2]],
      reference_radius = group$reference_radius,
      source_radius = group$source_radius,
      stringsAsFactors = FALSE
    )
  }))
}

real_score_fold <- function(tables, fold, specs, workers = 1L,
                            engines = c("transport_v2", "replay", "baselines"),
                            conditions = c(
                              "signal", "shuffled_order", "wrong_item",
                              "generic_gaze"
                            ),
                            seed = 20260817L) {
  subset <- real_fold_subsets(tables, fold)
  subset$source_eval <- subset$source_eval[
    subset$source_eval$condition %in% conditions, , drop = FALSE
  ]
  match_on <- c("participant", "item")
  contrast_on <- "participant"
  rows <- list()
  resources <- list()
  warps <- list()
  model_audit <- list()
  result_index <- 1L
  resource_index <- 1L

  if ("transport_v2" %in% engines) {
    profile <- real_profile("transport_v2", tables$signal$task[[1]], fold$id, {
      model <- real_ns("fit_gaze_transport_v2_model")(
        subset$ref_train, subset$source_train, match_on, contrast_on,
        "fixgroup", "fixgroup", specs$transport, workers = workers
      )
      scored <- lapply(seq_len(nrow(subset$source_eval)), function(i) {
        real_ns("score_gaze_transport_v2_row")(
          subset$source_eval[i, , drop = FALSE], subset$ref_eval,
          match_on, contrast_on, "fixgroup", "fixgroup", model,
          workers = workers
        )
      })
      list(model = model, scored = scored)
    })
    resources[[resource_index]] <- profile$resource
    resource_index <- resource_index + 1L
    model <- profile$value$model
    for (i in seq_along(profile$value$scored)) {
      rows[[result_index]] <- real_engine_row(
        subset$source_eval[i, , drop = FALSE], profile$value$scored[[i]],
        "transport_v2", fold$id
      )
      result_index <- result_index + 1L
    }
    warps[["transport_v2"]] <- real_warp_audit(
      model$warp$info, tables$signal$task[[1]], fold$id, "transport_v2"
    )
    model_audit$transport_v2 <- list(
      coverage_penalty = model$coverage_policy$coverage_penalty,
      temperature = model$coverage_policy$temperature
    )
  }

  if ("replay" %in% engines) {
    profile <- real_profile("replay", tables$signal$task[[1]], fold$id, {
      model <- fit_gaze_replay_model(
        subset$ref_train, subset$source_train, match_on,
        contrast_on = contrast_on,
        spec = specs$replay
      )
      scored <- lapply(seq_len(nrow(subset$source_eval)), function(i) {
        real_ns("score_gaze_replay_row")(
          subset$source_eval[i, , drop = FALSE], subset$ref_eval,
          match_on, contrast_on, "fixgroup", "fixgroup", model
        )
      })
      list(model = model, scored = scored)
    })
    resources[[resource_index]] <- profile$resource
    resource_index <- resource_index + 1L
    model <- profile$value$model
    for (i in seq_along(profile$value$scored)) {
      rows[[result_index]] <- real_engine_row(
        subset$source_eval[i, , drop = FALSE], profile$value$scored[[i]],
        "replay", fold$id
      )
      result_index <- result_index + 1L
    }
    warps[["replay"]] <- real_warp_audit(
      model$warp$info, tables$signal$task[[1]], fold$id, "replay"
    )
    model_audit$replay <- list(
      parameters = model$parameters,
      calibration = model$calibration,
      emission_models = model$emission_models
    )
  }

  if ("baselines" %in% engines) {
    profile <- real_profile("baselines", tables$signal$task[[1]], fold$id, {
      model <- real_ns("fit_gaze_baseline_model")(
        subset$ref_train, subset$source_train, match_on, contrast_on,
        c("participant", "item"), "fixgroup", "fixgroup",
        specs$baseline, seed + fold$id
      )
      sets <- real_ns("baseline_build_feature_sets")(
        subset$ref_eval, subset$source_eval, match_on, contrast_on,
        "fixgroup", "fixgroup", specs$baseline, model$warp,
        model$availability
      )
      list(model = model, scored = real_ns("score_gaze_baseline_sets")(sets, model))
    })
    resources[[resource_index]] <- profile$resource
    model <- profile$value$model
    keep_methods <- c(
      "multimatch_ridge_registered", "density_ridge_registered",
      "elastic_ridge_registered", "density_raw", "density_registered"
    )
    for (i in seq_along(profile$value$scored)) {
      for (method in intersect(keep_methods, names(profile$value$scored[[i]]))) {
        rows[[result_index]] <- real_baseline_row(
          subset$source_eval[i, , drop = FALSE],
          profile$value$scored[[i]][[method]], fold$id
        )
        result_index <- result_index + 1L
      }
    }
    warps[["baselines"]] <- real_warp_audit(
      model$warp$info, tables$signal$task[[1]], fold$id, "baselines"
    )
    model_audit$baselines <- lapply(model$composites, function(composite) {
      if (!identical(composite$status, "scored")) return(composite)
      list(
        status = composite$status,
        lambda = composite$selection$lambda,
        coefficients = composite$model$coefficients,
        inner_log_loss = composite$selection$mean_log_loss,
        inner_overlap = vapply(model$inner_splits, `[[`, integer(1), "overlap_match_n")
      )
    })
  }

  train_participants <- sort(unique(subset$source_train$participant))
  train_items <- sort(unique(subset$source_train$item))
  audit <- list(
    task = tables$signal$task[[1]], fold = fold$id,
    train_participants = train_participants,
    eval_participants = sort(unique(subset$source_eval$participant)),
    train_items = train_items,
    eval_items = sort(unique(subset$source_eval$item)),
    participant_overlap = length(intersect(
      train_participants, unique(subset$source_eval$participant)
    )),
    item_overlap = length(intersect(
      train_items, unique(subset$source_eval$item)
    )),
    training_conditions = "signal",
    model = model_audit
  )
  list(
    scored = dplyr::bind_rows(rows),
    resources = dplyr::bind_rows(resources),
    warps = dplyr::bind_rows(warps),
    audit = audit
  )
}

real_candidate_ece <- function(candidate_tables, bins = 5L) {
  candidates <- do.call(rbind, lapply(candidate_tables, function(tab) {
    data.frame(
      posterior = tab$posterior,
      is_true = as.numeric(tab$is_true)
    )
  }))
  breaks <- seq(0, 1, length.out = bins + 1L)
  bin <- cut(candidates$posterior, breaks, include.lowest = TRUE, labels = FALSE)
  total <- nrow(candidates)
  sum(vapply(split(seq_len(total), bin), function(index) {
    length(index) / total * abs(
      mean(candidates$posterior[index]) - mean(candidates$is_true[index])
    )
  }, numeric(1)))
}

real_method_summary <- function(scored) {
  tab <- scored[scored$calibrated & scored$status == "scored", , drop = FALSE]
  groups <- split(
    seq_len(nrow(tab)),
    interaction(tab$task, tab$condition, tab$method, drop = TRUE)
  )
  dplyr::bind_rows(lapply(groups, function(index) {
    part <- tab[index, , drop = FALSE]
    data.frame(
      task = part$task[[1]], condition = part$condition[[1]],
      method = part$method[[1]], n = nrow(part),
      mean_info_bits = mean(part$gaze_info_bits),
      mean_log_loss = mean(part$log_loss),
      mean_brier = mean(part$brier_score),
      mean_rank = mean(part$template_rank),
      top1_credit = mean(part$top1_credit),
      ece_5 = real_candidate_ece(part$candidates),
      stringsAsFactors = FALSE
    )
  }))
}

real_frozen_summary <- function(scored) {
  tab <- scored[!scored$calibrated & scored$status == "scored", , drop = FALSE]
  groups <- split(
    seq_len(nrow(tab)),
    interaction(tab$task, tab$condition, tab$method, drop = TRUE)
  )
  dplyr::bind_rows(lapply(groups, function(index) {
    part <- tab[index, , drop = FALSE]
    data.frame(
      task = part$task[[1]], condition = part$condition[[1]],
      method = part$method[[1]], n = nrow(part),
      mean_compatibility_bits = mean(part$compatibility_bits),
      mean_rank = mean(part$template_rank),
      top1_credit = mean(part$top1_credit),
      stringsAsFactors = FALSE
    )
  }))
}

real_crossed_draw <- function(tab, value, draws, seed) {
  participants <- sort(unique(tab$participant))
  items <- sort(unique(tab$item))
  set.seed(seed)
  vapply(seq_len(draws), function(draw) {
    participant_frequency <- table(sample(
      participants, length(participants), replace = TRUE
    ))
    item_frequency <- table(sample(items, length(items), replace = TRUE))
    participant_weight <- as.numeric(participant_frequency[tab$participant])
    item_weight <- as.numeric(item_frequency[as.character(tab$item)])
    participant_weight[is.na(participant_weight)] <- 0
    item_weight[is.na(item_weight)] <- 0
    weight <- participant_weight * item_weight
    stats::weighted.mean(tab[[value]], weight)
  }, numeric(1))
}

real_bootstrap_summary <- function(scored, draws = 1000L, seed = 20260817L) {
  tab <- scored[scored$calibrated & scored$status == "scored", , drop = FALSE]
  groups <- split(
    seq_len(nrow(tab)),
    interaction(tab$task, tab$condition, tab$method, drop = TRUE)
  )
  metrics <- c("gaze_info_bits", "log_loss", "top1_credit")
  row_index <- 1L
  rows <- list()
  for (index in groups) {
    part <- tab[index, , drop = FALSE]
    for (metric in metrics) {
      values <- real_crossed_draw(
        part, metric, draws,
        seed + row_index * 1009L + match(metric, metrics) * 101L
      )
      rows[[row_index]] <- data.frame(
        task = part$task[[1]], condition = part$condition[[1]],
        method = part$method[[1]], metric = metric,
        estimate = mean(part[[metric]]),
        lower_95 = unname(stats::quantile(values, 0.025, type = 8)),
        upper_95 = unname(stats::quantile(values, 0.975, type = 8)),
        draws = draws, stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  dplyr::bind_rows(rows)
}

real_operating_characteristics <- function(scored) {
  tab <- scored[scored$calibrated & scored$status == "scored", , drop = FALSE]
  groups <- split(seq_len(nrow(tab)), interaction(tab$task, tab$method, drop = TRUE))
  dplyr::bind_rows(lapply(groups, function(index) {
    part <- tab[index, , drop = FALSE]
    null <- sort(part$gaze_info_bits[
      part$condition %in% c("wrong_item", "generic_gaze")
    ])
    signal <- part$gaze_info_bits[part$condition == "signal"]
    threshold <- null[[ceiling(0.95 * length(null))]]
    credit <- function(value) mean(value > threshold) + 0.5 * mean(value == threshold)
    data.frame(
      task = part$task[[1]], method = part$method[[1]],
      threshold_bits = threshold, null_error = credit(null),
      power = credit(signal), null_n = length(null), signal_n = length(signal),
      stringsAsFactors = FALSE
    )
  }))
}

real_order_checks <- function(scored) {
  methods <- unique(scored$method)
  dplyr::bind_rows(lapply(unique(scored$task), function(task) {
    dplyr::bind_rows(lapply(methods, function(method) {
      part <- scored[scored$task == task & scored$method == method &
                       scored$condition %in% c("signal", "shuffled_order"), ]
      if (nrow(part) == 0L ||
          !all(c("signal", "shuffled_order") %in% part$condition)) return(NULL)
      wide <- merge(
        part[part$condition == "signal", c(
          "participant", "item", "gaze_info_bits", "compatibility_bits"
        )],
        part[part$condition == "shuffled_order", c(
          "participant", "item", "gaze_info_bits", "compatibility_bits"
        )],
        by = c("participant", "item"), suffixes = c("_signal", "_shuffled")
      )
      if (nrow(wide) == 0L) return(NULL)
      calibrated <- part$calibrated[[1]]
      signal <- if (calibrated) wide$gaze_info_bits_signal else
        wide$compatibility_bits_signal
      shuffled <- if (calibrated) wide$gaze_info_bits_shuffled else
        wide$compatibility_bits_shuffled
      data.frame(
        task = task, method = method, calibrated = calibrated,
        mean_signal_minus_shuffled = mean(signal - shuffled),
        max_absolute_pair_difference = max(abs(signal - shuffled)),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

real_diagnostic_summary <- function(scored) {
  tab <- scored[
    scored$condition == "signal" & scored$method %in% c("transport_v2", "replay"),
    , drop = FALSE
  ]
  groups <- split(seq_len(nrow(tab)), interaction(tab$task, tab$method, drop = TRUE))
  dplyr::bind_rows(lapply(groups, function(index) {
    part <- tab[index, , drop = FALSE]
    data.frame(
      task = part$task[[1]], method = part$method[[1]], n = nrow(part),
      mean_replay_coverage = mean(part$replay_coverage, na.rm = TRUE),
      mean_background_coverage = if (all(is.na(part$background_coverage))) {
        NA_real_
      } else mean(part$background_coverage, na.rm = TRUE),
      mean_spatial_rmse_px = mean(part$spatial_rmse_px, na.rm = TRUE),
      mean_local_order_error = if (all(is.na(part$local_order_error))) {
        NA_real_
      } else mean(part$local_order_error, na.rm = TRUE),
      mean_expected_restarts = if (all(is.na(part$expected_restarts))) {
        NA_real_
      } else mean(part$expected_restarts, na.rm = TRUE),
      convergence_rate = mean(part$converged),
      stringsAsFactors = FALSE
    )
  }))
}

real_covariate_summary <- function(scored) {
  tab <- scored[
    scored$calibrated & scored$status == "scored" & scored$condition == "signal",
    , drop = FALSE
  ]
  variables <- c("age", "probe_type", "degradation", "cue_duration")
  rows <- list()
  index <- 1L
  for (variable in variables) {
    groups <- split(
      seq_len(nrow(tab)),
      interaction(tab$task, tab$method, tab[[variable]], drop = TRUE)
    )
    for (group in groups) {
      part <- tab[group, , drop = FALSE]
      rows[[index]] <- data.frame(
        task = part$task[[1]], method = part$method[[1]],
        variable = variable, level = as.character(part[[variable]][[1]]),
        n = nrow(part), mean_info_bits = mean(part$gaze_info_bits),
        mean_log_loss = mean(part$log_loss),
        top1_credit = mean(part$top1_credit), stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  dplyr::bind_rows(rows)
}

real_paired_comparison <- function(scored, first, second,
                                   draws = 1000L, seed = 20260817L) {
  columns <- c("task", "participant", "item", "gaze_info_bits", "log_loss",
               "top1_credit")
  one <- scored[
    scored$calibrated & scored$status == "scored" &
      scored$condition == "signal" & scored$method == first,
    columns, drop = FALSE
  ]
  two <- scored[
    scored$calibrated & scored$status == "scored" &
      scored$condition == "signal" & scored$method == second,
    columns, drop = FALSE
  ]
  paired <- merge(one, two, by = c("task", "participant", "item"),
                  suffixes = c("_first", "_second"))
  paired$delta_info <- paired$gaze_info_bits_first - paired$gaze_info_bits_second
  paired$delta_log_loss <- paired$log_loss_second - paired$log_loss_first
  paired$delta_top1 <- paired$top1_credit_first - paired$top1_credit_second
  metrics <- c("delta_info", "delta_log_loss", "delta_top1")
  dplyr::bind_rows(lapply(seq_along(metrics), function(i) {
    metric <- metrics[[i]]
    values <- real_crossed_draw(paired, metric, draws, seed + i * 10007L)
    data.frame(
      first = first, second = second, metric = metric,
      estimate = mean(paired[[metric]]),
      lower_95 = unname(stats::quantile(values, 0.025, type = 8)),
      upper_95 = unname(stats::quantile(values, 0.975, type = 8)),
      draws = draws, stringsAsFactors = FALSE
    )
  }))
}

real_fold_audit_table <- function(audits) {
  dplyr::bind_rows(lapply(audits, function(audit) {
    data.frame(
      task = audit$task, outer_fold = audit$fold,
      train_participants = paste(audit$train_participants, collapse = ";"),
      eval_participants = paste(audit$eval_participants, collapse = ";"),
      train_items = paste(audit$train_items, collapse = ";"),
      eval_items = paste(audit$eval_items, collapse = ";"),
      participant_overlap = audit$participant_overlap,
      item_overlap = audit$item_overlap,
      training_conditions = audit$training_conditions,
      stringsAsFactors = FALSE
    )
  }))
}

real_gate_verdict <- function(scored, summary, operating, order_checks,
                              fold_audit, comparison) {
  engine_summary <- summary[
    summary$condition == "signal" &
      summary$method %in% c("transport_v2", "replay"),
    , drop = FALSE
  ]
  repeated <- engine_summary[engine_summary$task == "repeated_viewing", ]
  imagery <- engine_summary[engine_summary$task == "blank_screen_imagery", ]
  repeated_pass <- nrow(repeated) == 2L &&
    all(repeated$mean_info_bits > 0) && all(repeated$top1_credit > 0.25)
  imagery_pass <- any(imagery$mean_info_bits > 0 & imagery$top1_credit > 0.25)
  engine_order <- order_checks[
    order_checks$task == "repeated_viewing" &
      order_checks$method %in% c("transport_v2", "replay"),
  ]
  density_order <- order_checks[
    order_checks$task == "repeated_viewing" &
      order_checks$method == "density_registered",
  ]
  order_pass <- any(engine_order$mean_signal_minus_shuffled > 0) &&
    nrow(density_order) == 1L && density_order$max_absolute_pair_difference <= 1e-10
  engine_operating <- operating[
    operating$method %in% c("transport_v2", "replay"), , drop = FALSE
  ]
  null_pass <- all(engine_operating$null_error <= 0.075)
  density_null <- scored[
    !scored$calibrated & scored$condition != "signal" &
      scored$method %in% c("density_raw", "density_registered"),
  ]
  raw_top1 <- mean(density_null$top1_credit[density_null$method == "density_raw"])
  registered_top1 <- mean(
    density_null$top1_credit[density_null$method == "density_registered"]
  )
  registration_pass <- registered_top1 - raw_top1 <= 0.05
  convergence_pass <- all(vapply(c("transport_v2", "replay"), function(method) {
    mean(scored$converged[scored$method == method]) >= 0.99
  }, logical(1)))
  leakage_pass <- all(fold_audit$participant_overlap == 0L) &&
    all(fold_audit$item_overlap == 0L)
  gates <- data.frame(
    gate = c(
      "repeated_viewing", "blank_screen_imagery", "local_order",
      "null_error", "null_registration", "engine_convergence",
      "two_way_crossfit"
    ),
    passed = c(
      repeated_pass, imagery_pass, order_pass, null_pass, registration_pass,
      convergence_pass, leakage_pass
    ),
    stringsAsFactors = FALSE
  )

  signal <- summary[summary$condition == "signal", ]
  pooled <- aggregate(
    cbind(mean_log_loss, top1_credit) ~ method, signal, mean
  )
  replay <- pooled[pooled$method == "replay", ]
  transport <- pooled[pooled$method == "transport_v2", ]
  baselines <- pooled[pooled$method %in% c(
    "multimatch_ridge_registered", "density_ridge_registered",
    "elastic_ridge_registered"
  ), ]
  best <- baselines[which.min(baselines$mean_log_loss), ]
  competitive <- nrow(replay) == 1L && nrow(transport) == 1L && nrow(best) == 1L &&
    replay$mean_log_loss <= transport$mean_log_loss + 0.05 &&
    replay$top1_credit >= transport$top1_credit - 0.05 &&
    replay$mean_log_loss <= best$mean_log_loss + 0.10 &&
    replay$top1_credit >= best$top1_credit - 0.10
  comparison_info <- comparison[comparison$metric == "delta_info", ]
  comparison_loss <- comparison[comparison$metric == "delta_log_loss", ]
  superiority <- nrow(comparison_info) == 1L && nrow(comparison_loss) == 1L &&
    comparison_info$lower_95 > 0 && comparison_loss$lower_95 > 0
  list(
    gates = gates,
    advance = all(gates$passed),
    provisional_default = if (all(gates$passed) && competitive) {
      "replay"
    } else {
      "no_real_data_default"
    },
    strongest_baseline = best$method,
    superiority_supported = superiority,
    pooled_signal = pooled,
    density_null_top1 = data.frame(
      method = c("density_raw", "density_registered"),
      top1_credit = c(raw_top1, registered_top1)
    )
  )
}

real_write_results <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  scored <- result$scored
  scored$candidates <- NULL
  csv <- list(
    "method-summary.csv" = result$method_summary,
    "frozen-summary.csv" = result$frozen_summary,
    "bootstrap-intervals.csv" = result$bootstrap,
    "operating-characteristics.csv" = result$operating,
    "order-controls.csv" = result$order_checks,
    "alignment-diagnostics.csv" = result$diagnostics,
    "warp-estimates.csv" = result$warps,
    "fold-audit.csv" = result$fold_audit,
    "covariate-sensitivity.csv" = result$covariates,
    "paired-comparison.csv" = result$comparison,
    "gate-verdict.csv" = result$verdict$gates,
    "pooled-signal.csv" = result$verdict$pooled_signal,
    "density-null-registration.csv" = result$verdict$density_null_top1,
    "resources.csv" = result$resources,
    "sensitivity-full-window.csv" = result$sensitivity,
    "trial-scores.csv" = scored
  )
  for (name in names(csv)) {
    utils::write.csv(csv[[name]], file.path(output_dir, name), row.names = FALSE)
  }
  configuration <- data.frame(
    protocol_version = 1L, seed = result$config$seed,
    smoke = result$config$smoke,
    participants = paste(result$config$cohort$participants, collapse = ";"),
    items = paste(result$config$cohort$items, collapse = ";"),
    advance = result$verdict$advance,
    provisional_default = result$verdict$provisional_default,
    strongest_baseline = result$verdict$strongest_baseline,
    superiority_supported = result$verdict$superiority_supported,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    configuration, file.path(output_dir, "configuration.csv"), row.names = FALSE
  )
  saveRDS(
    result[c(
      "method_summary", "frozen_summary", "bootstrap", "operating",
      "order_checks", "diagnostics", "warps", "fold_audit", "covariates",
      "comparison", "verdict", "resources", "sensitivity", "config"
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

run_gaze_weave_real_court <- function(
    output_dir = NULL, seed = 20260817L, smoke = FALSE,
    bootstrap_draws = if (smoke) 100L else 1000L) {
  if (!smoke && seed != 20260817L) {
    stop("The final real court is frozen to seed 20260817; use smoke = TRUE for debugging.")
  }
  raw <- real_read_wynn(c(0, 3000))
  cohort <- real_select_cohort(raw, seed)
  expected_participants <- c("104", "121", "124", "128", "315", "317", "319", "327")
  expected_items <- c(5L, 18L, 69L, 73L, 81L, 82L, 88L, 105L)
  if (!smoke && (!identical(cohort$participants, expected_participants) ||
                 !identical(cohort$items, expected_items))) {
    stop("The data-only cohort no longer matches the frozen protocol.")
  }
  tasks <- lapply(
    c("repeated_viewing", "blank_screen_imagery"),
    function(task) real_task_tables(raw, cohort, task, seed)
  )
  names(tasks) <- c("repeated_viewing", "blank_screen_imagery")
  participant_age <- unique(tasks[[1]]$signal[c("participant", "age")])
  plan <- real_fold_plan(cohort, participant_age, seed)
  specs <- real_specs(smoke)
  detected <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected) || detected < 1L) detected <- 2L
  workers <- if (smoke) min(2L, detected) else min(4L, detected)

  fold_results <- list()
  index <- 1L
  for (task in names(tasks)) {
    for (fold in plan$folds) {
      message("Real court: ", task, ", outer fold ", fold$id, "/4")
      fold_results[[index]] <- real_score_fold(
        tasks[[task]], fold, specs, workers = workers
      )
      index <- index + 1L
    }
  }
  scored <- dplyr::bind_rows(lapply(fold_results, `[[`, "scored"))
  resources <- dplyr::bind_rows(lapply(fold_results, `[[`, "resources"))
  warps <- dplyr::bind_rows(lapply(fold_results, `[[`, "warps"))
  audits <- lapply(fold_results, `[[`, "audit")
  fold_audit <- real_fold_audit_table(audits)
  method_summary <- real_method_summary(scored)
  frozen_summary <- real_frozen_summary(scored)
  bootstrap <- real_bootstrap_summary(scored, bootstrap_draws, seed)
  operating <- real_operating_characteristics(scored)
  order_checks <- real_order_checks(scored)
  diagnostics <- real_diagnostic_summary(scored)
  covariates <- real_covariate_summary(scored)
  baseline_signal <- method_summary[
    method_summary$condition == "signal" &
      method_summary$method %in% c(
        "multimatch_ridge_registered", "density_ridge_registered",
        "elastic_ridge_registered"
      ),
  ]
  baseline_pooled <- aggregate(mean_log_loss ~ method, baseline_signal, mean)
  strongest_baseline <- baseline_pooled$method[
    which.min(baseline_pooled$mean_log_loss)
  ]
  comparison <- real_paired_comparison(
    scored, "replay", strongest_baseline, bootstrap_draws, seed
  )
  verdict <- real_gate_verdict(
    scored, method_summary, operating, order_checks, fold_audit, comparison
  )

  raw_full <- real_read_wynn(NULL)
  full_task <- real_task_tables(raw_full, cohort, "blank_screen_imagery", seed)
  sensitivity_rows <- list()
  for (fold in plan$folds) {
    sensitivity_rows[[fold$id]] <- real_score_fold(
      full_task, fold, specs, workers = 1L, engines = "replay",
      conditions = "signal"
    )$scored
  }
  sensitivity_scored <- dplyr::bind_rows(sensitivity_rows)
  sensitivity <- real_method_summary(sensitivity_scored)
  sensitivity$window <- "all_deposited_posttest_fixations"

  result <- list(
    scored = scored, method_summary = method_summary,
    frozen_summary = frozen_summary, bootstrap = bootstrap,
    operating = operating, order_checks = order_checks,
    diagnostics = diagnostics, warps = warps, fold_audit = fold_audit,
    covariates = covariates, comparison = comparison, verdict = verdict,
    resources = resources, sensitivity = sensitivity,
    config = list(
      protocol_version = 1L, seed = seed, smoke = smoke,
      bootstrap_draws = bootstrap_draws, cohort = cohort,
      fold_plan = plan, workers = workers
    ),
    audits = audits
  )
  if (!is.null(output_dir)) real_write_results(result, output_dir)
  result
}
