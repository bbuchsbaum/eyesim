# Maximal-cohort recognition sensitivity court for GazeWeave.
#
# Protocol: inst/validation/GAZEWEAVE-RECOGNITION-FULL-COHORT.md
# Run from the eyesim source root after devtools::load_all().

full_recognition_protocol_version <- 1L
full_recognition_seed <- 20260825L
full_recognition_item_seed <- 20260820L
full_recognition_bootstrap_draws <- 2000L
full_recognition_candidate_n <- 5L
full_recognition_expected <- c(
  raw_participants = 46L,
  eligible_participants = 45L,
  eligible_pairs = 1348L,
  retained_participants = 36L,
  retained_pairs = 1295L,
  items = 120L
)
full_recognition_mm_metrics <- c(
  "mm_vector", "mm_direction", "mm_length", "mm_position",
  "mm_duration", "mm_position_emd"
)
full_recognition_engines <- c("transport_v2", "replay")
full_recognition_native_methods <- c(
  "density_sigma_80_raw", "density_sigma_160_raw",
  paste0("multimatch_", full_recognition_mm_metrics, "_raw")
)
full_recognition_calibrated_comparators <- c(
  "density_sigma_80_raw_calibrated",
  "density_sigma_160_raw_calibrated",
  paste0(
    "multimatch_", full_recognition_mm_metrics, "_raw_calibrated"
  ),
  "density_ridge_registered",
  "multimatch_ridge_registered",
  "elastic_ridge_registered"
)
full_recognition_all_methods <- c(
  full_recognition_engines,
  full_recognition_native_methods,
  full_recognition_calibrated_comparators
)

full_recognition_probe_file <- system.file(
  "validation", "gaze-weave-probe-delay.R", package = "eyesim"
)
full_recognition_model_file <- system.file(
  "validation", "gaze-weave-recognition-sensitivity.R", package = "eyesim"
)
if (!nzchar(full_recognition_probe_file)) {
  full_recognition_probe_file <- file.path(
    "inst", "validation", "gaze-weave-probe-delay.R"
  )
}
if (!nzchar(full_recognition_model_file)) {
  full_recognition_model_file <- file.path(
    "inst", "validation", "gaze-weave-recognition-sensitivity.R"
  )
}
if (!file.exists(full_recognition_probe_file) ||
    !file.exists(full_recognition_model_file)) {
  stop("Run the full recognition court from the eyesim source root.")
}
source(full_recognition_probe_file, local = TRUE)
source(full_recognition_model_file, local = TRUE)

full_recognition_read_inputs <- function(verify = TRUE) {
  data_dir <- file.path("test_data", "wynn_probe_delay")
  if (verify) probe_delay_verify_inputs(data_dir)
  files <- probe_delay_data_files(data_dir)
  study_raw <- probe_delay_read_csv(files[["study"]], "study")
  retrieval_raw <- probe_delay_read_csv(files[["retrieval"]], "retrieval")
  list(
    study = probe_delay_study_paths(study_raw),
    retrieval = probe_delay_retrieval_paths(retrieval_raw, "combined"),
    trials = probe_delay_trial_tables(study_raw, retrieval_raw)
  )
}

full_recognition_eligible_pairs <- function(raw) {
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
    study_pairs$four_presentations & study_pairs$study_min_nfix >= 3L,
    c("participant", "item"), drop = FALSE
  ]
  retrieval <- raw$retrieval[
    raw$retrieval$probe_type %in% c("old", "lure") &
      raw$retrieval$nfix >= 3L, , drop = FALSE
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

full_recognition_select_cohort <- function(
    raw, item_seed = full_recognition_item_seed,
    candidate_n = full_recognition_candidate_n) {
  eligible <- full_recognition_eligible_pairs(raw)
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
  keep_participant <- rownames(counts)[
    counts[, 1L] >= candidate_n & counts[, 2L] >= candidate_n
  ]
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
    candidate_n = candidate_n,
    balance = list(
      probe_type = table(pairs$probe_type),
      degradation = table(pairs$degradation),
      accuracy = table(pairs$accuracy)
    )
  )
}

full_recognition_smoke_cohort <- function(cohort, participant_n = 4L) {
  participant_n <- as.integer(participant_n)
  if (participant_n < 4L || participant_n > length(cohort$participants)) {
    stop("Smoke cohort participant_n must be between four and the cohort size.")
  }
  participants <- cohort$participants[seq_len(participant_n)]
  smoke <- cohort
  smoke$participants <- participants
  smoke$pairs <- cohort$pairs[
    cohort$pairs$participant %in% participants, , drop = FALSE
  ]
  smoke$balance <- list(
    probe_type = table(smoke$pairs$probe_type),
    degradation = table(smoke$pairs$degradation),
    accuracy = table(smoke$pairs$accuracy)
  )
  smoke
}

full_recognition_participant_map <- function(cohort,
                                             seed = full_recognition_seed) {
  counts <- as.integer(table(factor(
    cohort$pairs$participant, levels = cohort$participants
  )))
  tab <- data.frame(
    participant = cohort$participants,
    eligible_trials = counts,
    stringsAsFactors = FALSE
  )
  set.seed(seed)
  tab$tie_order <- stats::runif(nrow(tab))
  order_index <- order(-tab$eligible_trials, tab$tie_order, tab$participant)
  fold_total <- c(0L, 0L)
  fold_n <- c(0L, 0L)
  tab$participant_fold <- NA_integer_
  for (index in order_index) {
    eligible_fold <- which(fold_total == min(fold_total))
    if (length(eligible_fold) > 1L) {
      eligible_fold <- eligible_fold[
        fold_n[eligible_fold] == min(fold_n[eligible_fold])
      ]
    }
    chosen <- if (length(eligible_fold) == 1L) {
      eligible_fold
    } else {
      sample(eligible_fold, 1L)
    }
    tab$participant_fold[[index]] <- chosen
    fold_total[[chosen]] <- fold_total[[chosen]] + tab$eligible_trials[[index]]
    fold_n[[chosen]] <- fold_n[[chosen]] + 1L
  }
  tab$tie_order <- NULL
  tab[order(tab$participant), , drop = FALSE]
}

full_recognition_candidate_plan <- function(
    cohort, seed = full_recognition_seed,
    candidate_n = full_recognition_candidate_n) {
  groups <- split(
    seq_len(nrow(cohort$pairs)),
    interaction(
      cohort$pairs$participant, cohort$pairs$item_fold,
      drop = TRUE, lex.order = TRUE
    )
  )
  rows <- lapply(groups, function(index) {
    part <- cohort$pairs[index, , drop = FALSE]
    participant <- part$participant[[1L]]
    item_fold <- part$item_fold[[1L]]
    items <- sort(unique(part$item))
    if (length(items) < candidate_n) {
      stop("A retained participant-item fold has too few candidates.")
    }
    set.seed(seed + as.integer(participant) * 1009L + item_fold * 101L)
    ring <- sample(items)
    dplyr::bind_rows(lapply(items, function(target) {
      target_position <- match(target, ring)
      positions <- ((target_position - 1L + 0:(candidate_n - 1L)) %%
                      length(ring)) + 1L
      candidates <- ring[positions]
      shift <- (as.integer(participant) + target + item_fold) %% candidate_n
      if (shift > 0L) {
        candidates <- candidates[c(
          (shift + 1L):candidate_n, seq_len(shift)
        )]
      }
      data.frame(
        participant = participant,
        target_item = target,
        item_fold = item_fold,
        candidate_set_id = paste(participant, target, sep = ":"),
        candidate_position = seq_len(candidate_n),
        candidate_item = candidates,
        is_true = candidates == target,
        stringsAsFactors = FALSE
      )
    }))
  })
  plan <- dplyr::bind_rows(rows)
  plan <- plan[order(
    plan$participant, plan$target_item, plan$candidate_position
  ), ]
  rownames(plan) <- NULL
  plan
}

full_recognition_fold_plan <- function(cohort,
                                       seed = full_recognition_seed) {
  participant_map <- full_recognition_participant_map(cohort, seed)
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
        eval_items = cohort$item_map$item[
          cohort$item_map$item_fold == item_fold
        ]
      )
      index <- index + 1L
    }
  }
  list(
    folds = folds,
    participant_map = participant_map,
    item_map = cohort$item_map
  )
}

full_recognition_validate_design <- function(cohort, candidate_plan,
                                             fold_plan,
                                             strict = TRUE) {
  expected <- full_recognition_expected
  observed <- c(
    raw_participants = cohort$raw_participant_n,
    eligible_participants = cohort$eligible_participant_n,
    eligible_pairs = cohort$eligible_pair_n,
    retained_participants = length(cohort$participants),
    retained_pairs = nrow(cohort$pairs),
    items = nrow(cohort$item_map)
  )
  if (strict && !identical(as.integer(observed), as.integer(expected))) {
    stop("The maximal-cohort support differs from the frozen protocol.")
  }
  set_sizes <- table(candidate_plan$candidate_set_id)
  truth_sizes <- tapply(
    candidate_plan$is_true, candidate_plan$candidate_set_id, sum
  )
  if (any(set_sizes != full_recognition_candidate_n) ||
      any(truth_sizes != 1L) ||
      nrow(candidate_plan) != nrow(cohort$pairs) *
        full_recognition_candidate_n) {
    stop("Candidate sets violate the fixed five-candidate truth contract.")
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
    stop("Outer folds do not evaluate every retained target exactly once.")
  }
  list(observed = observed, fold_audit = dplyr::bind_rows(audits))
}

full_recognition_task_tables <- function(raw, cohort) {
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
  reference <- study[
    study$presentation == 4L,
    c(
      "participant", "age", "item", "item_fold", "study_image_version",
      "nfix", "fixgroup"
    )
  ]
  reference <- merge(
    reference, metadata, by = c("participant", "age", "item"),
    all.x = TRUE, sort = FALSE
  )
  signal <- retrieval
  names(signal)[names(signal) == "Accuracy"] <- "accuracy"
  reference$task <- signal$task <- "combined"
  reference$condition <- signal$condition <- "signal"
  reference <- reference[order(reference$participant, reference$item), ]
  signal <- signal[order(signal$participant, signal$item), ]
  rownames(reference) <- rownames(signal) <- NULL
  if (nrow(reference) != nrow(cohort$pairs) ||
      nrow(signal) != nrow(reference)) {
    stop("Combined task tables do not match the maximal cohort.")
  }
  list(
    reference = tibble::as_tibble(reference),
    signal = tibble::as_tibble(signal)
  )
}

full_recognition_fold_subsets <- function(tables, candidate_plan, fold) {
  eval_target <- with(
    tables$signal,
    participant %in% fold$eval_participants & item %in% fold$eval_items
  )
  train_target <- with(
    tables$signal,
    !participant %in% fold$eval_participants & !item %in% fold$eval_items
  )
  source_train <- tables$signal[train_target, , drop = FALSE]
  ref_train <- tables$reference[train_target, , drop = FALSE]
  source_eval <- tables$signal[eval_target, , drop = FALSE]
  plan <- candidate_plan[
    candidate_plan$participant %in% fold$eval_participants &
      candidate_plan$target_item %in% fold$eval_items, , drop = FALSE
  ]
  source_eval$candidate_set_id <- paste(
    source_eval$participant, source_eval$item, sep = ":"
  )
  reference_lookup <- tables$reference
  names(reference_lookup)[names(reference_lookup) == "item"] <- "candidate_item"
  ref_eval <- merge(
    plan, reference_lookup,
    by = c("participant", "candidate_item"), all.x = TRUE, sort = FALSE
  )
  ref_eval$item <- ref_eval$candidate_item
  ref_eval <- ref_eval[order(
    ref_eval$candidate_set_id, ref_eval$candidate_position
  ), , drop = FALSE]
  source_eval <- source_eval[order(source_eval$candidate_set_id), ]
  if (nrow(ref_eval) != nrow(source_eval) * full_recognition_candidate_n ||
      anyNA(ref_eval$fixgroup)) {
    stop("Could not expand the frozen evaluation candidate sets.")
  }
  list(
    ref_train = ref_train,
    source_train = source_train,
    ref_eval = tibble::as_tibble(ref_eval),
    source_eval = tibble::as_tibble(source_eval)
  )
}

full_recognition_specs <- function(smoke = FALSE) {
  specs <- real_specs(smoke)
  specs$baseline <- gaze_baseline_spec(
    screen = specs$screen,
    density_sigmas = c(30, 60, 80, 120, 160),
    density_grid = if (smoke) 12L else 24L,
    methods = c("multimatch", "density", "elastic"),
    warp = specs$warp,
    lambda_grid = c(0.01, 0.1, 1, 10),
    inner_folds = 2L,
    elastic_radii = c(consensus = 80, rigidity = 400, matching = 40),
    elastic_maxit = if (smoke) 10L else 25L,
    elastic_tolerance = 1e-4
  )
  specs
}

full_recognition_slim_candidates <- function(candidates) {
  columns <- intersect(
    c(
      "candidate_key", "is_true", "log_score", "prior", "posterior",
      "status", "reason"
    ),
    names(candidates)
  )
  candidates[columns]
}

full_recognition_baseline_row <- function(source_row, evidence, outer_fold) {
  row <- real_baseline_row(source_row, evidence, outer_fold)
  row$candidates <- list(full_recognition_slim_candidates(evidence$candidates))
  row$compatibility_probability_true <- if (!is.null(
    evidence$compatibility_probability_true
  )) evidence$compatibility_probability_true else NA_real_
  row$compatibility_log_loss <- if (!is.null(
    evidence$compatibility_log_loss
  )) evidence$compatibility_log_loss else NA_real_
  row
}

full_recognition_score_fold_method <- function(
    tables, candidate_plan, fold, specs, method,
    workers = 1L, seed = full_recognition_seed, smoke = FALSE,
    subset_function = full_recognition_fold_subsets) {
  subset <- subset_function(tables, candidate_plan, fold)
  match_on <- c("participant", "item")
  train_contrast <- "participant"
  eval_contrast <- "candidate_set_id"
  profile <- real_profile(method, "combined", fold$id, {
    if (identical(method, "transport_v2")) {
      model <- real_ns("fit_gaze_transport_v2_model")(
        subset$ref_train, subset$source_train,
        match_on, train_contrast, "fixgroup", "fixgroup",
        specs$transport, workers = workers
      )
      scored <- lapply(seq_len(nrow(subset$source_eval)), function(i) {
        real_ns("score_gaze_transport_v2_row")(
          subset$source_eval[i, , drop = FALSE], subset$ref_eval,
          match_on, eval_contrast, "fixgroup", "fixgroup", model,
          workers = workers
        )
      })
      list(model = model, scored = scored)
    } else if (identical(method, "replay")) {
      model <- fit_gaze_replay_model(
        subset$ref_train, subset$source_train, match_on,
        contrast_on = train_contrast,
        spec = specs$replay
      )
      scored <- lapply(seq_len(nrow(subset$source_eval)), function(i) {
        real_ns("score_gaze_replay_row")(
          subset$source_eval[i, , drop = FALSE], subset$ref_eval,
          match_on, eval_contrast, "fixgroup", "fixgroup", model
        )
      })
      list(model = model, scored = scored)
    } else if (identical(method, "baselines")) {
      model <- real_ns("fit_gaze_baseline_model")(
        subset$ref_train, subset$source_train,
        match_on, train_contrast, c("participant", "item"),
        "fixgroup", "fixgroup", specs$baseline, seed + fold$id
      )
      if (!all(model$availability$available)) {
        missing <- model$availability$method[!model$availability$available]
        stop("Required full-cohort comparators unavailable: ",
             paste(missing, collapse = ", "))
      }
      sets <- real_ns("baseline_build_feature_sets")(
        subset$ref_eval, subset$source_eval,
        match_on, eval_contrast, "fixgroup", "fixgroup",
        specs$baseline, model$warp, model$availability
      )
      list(
        model = model,
        scored = real_ns("score_gaze_baseline_sets")(sets, model)
      )
    } else {
      stop("Unknown full-cohort scoring method: ", method)
    }
  })
  model <- profile$value$model
  rows <- list()
  index <- 1L
  if (method %in% full_recognition_engines) {
    for (i in seq_along(profile$value$scored)) {
      rows[[index]] <- real_engine_row(
        subset$source_eval[i, , drop = FALSE],
        profile$value$scored[[i]], method, fold$id
      )
      index <- index + 1L
    }
  } else {
    for (i in seq_along(profile$value$scored)) {
      available <- intersect(
        c(
          full_recognition_native_methods,
          full_recognition_calibrated_comparators
        ),
        names(profile$value$scored[[i]])
      )
      for (baseline_method in available) {
        rows[[index]] <- full_recognition_baseline_row(
          subset$source_eval[i, , drop = FALSE],
          profile$value$scored[[i]][[baseline_method]], fold$id
        )
        index <- index + 1L
      }
    }
  }
  warp <- real_warp_audit(
    model$warp$info, "combined", fold$id, method
  )
  audit <- data.frame(
    outer_fold = fold$id,
    method_group = method,
    train_trials = nrow(subset$source_train),
    eval_trials = nrow(subset$source_eval),
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
  list(
    scored = dplyr::bind_rows(rows),
    resources = profile$resource,
    warp = warp,
    audit = audit,
    model_audit = if (identical(method, "baselines")) {
      lapply(model$composites, function(value) {
        if (!identical(value$status, "scored")) return(value)
        list(
          status = value$status,
          lambda = value$selection$lambda,
          coefficients = value$model$coefficients,
          inner_log_loss = value$selection$mean_log_loss
        )
      })
    } else if (identical(method, "transport_v2")) {
      list(
        coverage_penalty = model$coverage_policy$coverage_penalty,
        temperature = model$coverage_policy$temperature
      )
    } else {
      list(
        parameters = model$parameters,
        calibration = model$calibration
      )
    },
    protocol_version = full_recognition_protocol_version,
    smoke = smoke,
    seed = seed,
    method_group = method,
    outer_fold = fold$id
  )
}

full_recognition_checkpoint <- function(output_dir, method, fold_id) {
  file.path(
    output_dir,
    sprintf("checkpoint-%s-fold-%02d.rds", method, fold_id)
  )
}

full_recognition_bootstrap_p <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(NA_real_)
  min(1, 2 * min(mean(values <= 0), mean(values >= 0)))
}

full_recognition_analyse <- function(scored, draws, seed) {
  calibrated <- scored[
    scored$calibrated & scored$status == "scored" &
      scored$method %in% c(
        full_recognition_engines,
        full_recognition_calibrated_comparators
      ), , drop = FALSE
  ]
  reference <- calibrated[calibrated$method == "transport_v2", ]
  reference <- recognition_prepare_rows(reference)
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws, seed)
  effects <- list()
  audits <- list()
  lopo_rows <- list()
  index <- 1L
  methods <- c(
    full_recognition_engines,
    full_recognition_calibrated_comparators
  )
  for (method in methods) {
    message("Full-cohort conditional model: ", method)
    tab <- recognition_prepare_rows(calibrated[calibrated$method == method, ])
    tab <- tab[order(tab$participant, tab$item), ]
    if (!identical(tab$trial_key, reference$trial_key)) {
      stop("Calibrated methods do not share full-cohort trial support.")
    }
    fit <- recognition_fit_method_phase(tab, plan)
    effect <- fit$effects
    effect$method <- method
    effect$phase <- "combined"
    effect$bootstrap_p <- vapply(
      seq_len(nrow(effect)),
      function(i) full_recognition_bootstrap_p(
        fit$bootstrap[effect$contrast[[i]], ]
      ),
      numeric(1)
    )
    if (!method %in% full_recognition_engines) {
      effect$lower_family <- NA_real_
      effect$upper_family <- NA_real_
    }
    effects[[index]] <- effect
    audit <- fit$audit
    audit$method <- method
    audit$phase <- "combined"
    audits[[index]] <- audit
    lopo <- as.data.frame(as.table(fit$lopo), stringsAsFactors = FALSE)
    names(lopo) <- c("contrast", "participant", "estimate")
    lopo$method <- method
    lopo$phase <- "combined"
    lopo_rows[[index]] <- lopo
    index <- index + 1L
  }
  effects <- dplyr::bind_rows(effects)
  comparator <- !effects$method %in% full_recognition_engines
  effects$comparator_bh_q <- NA_real_
  effects$comparator_bh_q[comparator] <- stats::p.adjust(
    effects$bootstrap_p[comparator], method = "BH"
  )
  audit <- dplyr::bind_rows(audits)
  lopo <- dplyr::bind_rows(lopo_rows)
  primary <- effects[effects$method %in% full_recognition_engines, ]
  primary_audit <- audit[audit$method %in% full_recognition_engines, ]
  lopo_required <- ceiling(
    0.8 * length(unique(reference$participant))
  )
  convergence_rate <- vapply(
    full_recognition_engines,
    function(method) mean(scored$converged[scored$method == method]),
    numeric(1)
  )
  convergence_ok <- all(is.finite(convergence_rate)) &&
    all(convergence_rate >= 0.95)
  integrity <- all(primary$valid_fraction >= 0.95) &&
    all(primary$design_kappa <= 30) &&
    all(primary_audit$lmer_status == "scored") &&
    all(primary_audit$lmer_converged) && convergence_ok
  ordinary <- primary$lower_95 > 0 | primary$upper_95 < 0
  family <- primary$lower_family > 0 | primary$upper_family < 0
  stable <- primary$lopo_sign_n >= lopo_required &
    primary$lmer_sign_agrees
  supported <- family & stable & integrity
  label <- if (!integrity) {
    "not_supported"
  } else if (any(supported)) {
    "supported"
  } else if (any(ordinary)) {
    "suggestive"
  } else {
    "not_supported"
  }
  verdict <- data.frame(
    verdict = label,
    supported_effects = paste(
      paste(primary$method[supported], primary$contrast[supported], sep = ":"),
      collapse = ";"
    ),
    suggestive_effects = paste(
      paste(primary$method[ordinary], primary$contrast[ordinary], sep = ":"),
      collapse = ";"
    ),
    design_integrity = integrity,
    transport_convergence_rate = convergence_rate[["transport_v2"]],
    replay_convergence_rate = convergence_rate[["replay"]],
    lopo_required = lopo_required,
    supersedes_gw12_for_sensitivity = TRUE,
    independent_replication = FALSE,
    closes_wang_gate = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    effects = effects,
    audit = audit,
    lopo = lopo,
    verdict = verdict,
    plan = plan
  )
}

full_recognition_candidate_scores <- function(scored) {
  rows <- lapply(seq_len(nrow(scored)), function(i) {
    candidates <- scored$candidates[[i]]
    if (is.null(candidates) || nrow(candidates) == 0L) return(NULL)
    candidates <- full_recognition_slim_candidates(candidates)
    candidates$participant <- scored$participant[[i]]
    candidates$target_item <- scored$item[[i]]
    candidates$method <- scored$method[[i]]
    candidates$outer_fold <- scored$outer_fold[[i]]
    candidates
  })
  dplyr::bind_rows(rows)
}

full_recognition_configuration <- function(result) {
  cohort <- result$config$cohort
  data.frame(
    protocol_version = full_recognition_protocol_version,
    seed = result$config$seed,
    item_seed = full_recognition_item_seed,
    bootstrap_draws = result$config$draws,
    raw_participants = cohort$raw_participant_n,
    eligible_participants = cohort$eligible_participant_n,
    eligible_pairs = cohort$eligible_pair_n,
    retained_participants = length(cohort$participants),
    retained_pairs = nrow(cohort$pairs),
    candidate_n = full_recognition_candidate_n,
    density_sigmas_px = "80;160",
    multimatch_metrics = paste(full_recognition_mm_metrics, collapse = ";"),
    verdict = result$analysis$verdict$verdict,
    local_only = TRUE,
    stringsAsFactors = FALSE
  )
}

full_recognition_write_results <- function(result, output_dir) {
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  score_columns <- setdiff(names(result$scored), "candidates")
  write(result$scored[score_columns], "trial-scores.csv")
  write(full_recognition_candidate_scores(result$scored), "candidate-scores.csv")
  write(real_method_summary(result$scored), "method-summary.csv")
  write(real_frozen_summary(result$scored), "native-method-summary.csv")
  write(result$analysis$effects, "conditional-effects.csv")
  write(result$analysis$audit, "conditional-model-audit.csv")
  write(result$analysis$lopo, "leave-one-participant-out.csv")
  write(result$analysis$verdict, "sensitivity-verdict.csv")
  write(result$config$design$fold_audit, "design-fold-audit.csv")
  write(result$config$fold_plan$participant_map, "participant-folds.csv")
  write(result$config$fold_plan$item_map, "item-folds.csv")
  write(result$config$candidate_plan, "candidate-plan.csv")
  write(result$resources, "resources.csv")
  write(result$warps, "warp-audit.csv")
  write(result$fold_audit, "scoring-fold-audit.csv")
  write(full_recognition_configuration(result), "configuration.csv")
  saveRDS(
    result[c(
      "analysis", "resources", "warps", "fold_audit", "config"
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
  write(manifest, "manifest-md5.csv")
  invisible(result)
}

run_gaze_weave_recognition_full_cohort <- function(
    output_dir = file.path(
      "inst", "validation",
      "gaze-weave-recognition-full-cohort-results"
    ),
    seed = full_recognition_seed,
    draws = full_recognition_bootstrap_draws,
    smoke = FALSE,
    method_groups = c("replay", "baselines", "transport_v2"),
    fold_ids = NULL,
    workers = NULL,
    resume = TRUE,
    finalize = !smoke) {
  if (!smoke && (seed != full_recognition_seed ||
                 draws != full_recognition_bootstrap_draws)) {
    stop("Non-smoke full-cohort runs use the frozen seed and bootstrap draws.")
  }
  method_groups <- match.arg(
    method_groups,
    c("transport_v2", "replay", "baselines"),
    several.ok = TRUE
  )
  raw <- full_recognition_read_inputs(verify = !smoke)
  cohort <- full_recognition_select_cohort(raw)
  if (smoke) cohort <- full_recognition_smoke_cohort(cohort)
  candidate_plan <- full_recognition_candidate_plan(cohort, seed)
  fold_plan <- full_recognition_fold_plan(cohort, seed)
  design <- full_recognition_validate_design(
    cohort, candidate_plan, fold_plan, strict = !smoke
  )
  tables <- full_recognition_task_tables(raw, cohort)
  specs <- full_recognition_specs(smoke)
  if (is.null(workers)) {
    detected <- parallel::detectCores(logical = FALSE)
    if (!is.finite(detected) || detected < 1L) detected <- 1L
    workers <- if (smoke) min(2L, detected) else min(4L, detected)
  }
  workers <- as.integer(workers)
  folds <- fold_plan$folds
  if (!is.null(fold_ids)) {
    folds <- folds[vapply(
      folds, function(fold) fold$id %in% fold_ids, logical(1)
    )]
  }
  if (length(folds) == 0L) stop("fold_ids selected no outer folds.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (fold in folds) {
    for (method in method_groups) {
      path <- full_recognition_checkpoint(output_dir, method, fold$id)
      if (resume && file.exists(path)) {
        checkpoint <- readRDS(path)
        if (!identical(
          checkpoint$protocol_version,
          full_recognition_protocol_version
        ) || !identical(checkpoint$smoke, smoke) ||
            !identical(checkpoint$seed, seed)) {
          stop("A checkpoint uses another protocol version: ", path)
        }
        message("Using checkpoint: ", basename(path))
        next
      }
      message(
        "Full recognition scoring: fold ", fold$id,
        ", method ", method
      )
      checkpoint <- full_recognition_score_fold_method(
        tables, candidate_plan, fold, specs, method,
        workers = workers, seed = seed, smoke = smoke
      )
      saveRDS(checkpoint, path, version = 3)
    }
  }
  if (!finalize) {
    return(invisible(list(
      cohort = cohort,
      candidate_plan = candidate_plan,
      fold_plan = fold_plan,
      design = design
    )))
  }
  required <- expand.grid(
    method = c("replay", "baselines", "transport_v2"),
    fold = 1:4,
    stringsAsFactors = FALSE
  )
  paths <- mapply(
    full_recognition_checkpoint,
    MoreArgs = list(output_dir = output_dir),
    method = required$method,
    fold_id = required$fold,
    USE.NAMES = FALSE
  )
  if (!all(file.exists(paths))) {
    stop("Full-cohort finalization requires all 12 method-fold checkpoints.")
  }
  checkpoints <- lapply(paths, readRDS)
  scored <- dplyr::bind_rows(lapply(checkpoints, `[[`, "scored"))
  support <- table(scored$method)
  if (!all(full_recognition_all_methods %in% names(support)) ||
      any(support[full_recognition_all_methods] != nrow(cohort$pairs))) {
    stop("Full-cohort methods do not share the frozen trial support.")
  }
  calibrated <- scored[
    scored$method %in% c(
      full_recognition_engines,
      full_recognition_calibrated_comparators
    ), ]
  if (any(calibrated$status != "scored") ||
      any(!calibrated$calibrated) ||
      any(!is.finite(calibrated$gaze_info_bits)) ||
      any(calibrated$candidate_count != full_recognition_candidate_n)) {
    stop("A calibrated full-cohort method failed its score contract.")
  }
  analysis <- full_recognition_analyse(scored, draws, seed)
  result <- list(
    scored = scored,
    analysis = analysis,
    resources = dplyr::bind_rows(lapply(checkpoints, `[[`, "resources")),
    warps = dplyr::bind_rows(lapply(checkpoints, `[[`, "warp")),
    fold_audit = dplyr::bind_rows(lapply(checkpoints, `[[`, "audit")),
    model_audit = lapply(checkpoints, `[[`, "model_audit"),
    config = list(
      seed = seed,
      draws = draws,
      smoke = smoke,
      workers = workers,
      cohort = cohort,
      candidate_plan = candidate_plan,
      fold_plan = fold_plan,
      design = design
    )
  )
  full_recognition_write_results(result, output_dir)
  result
}
