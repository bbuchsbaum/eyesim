# Frozen probe-delay persistence court for GazeWeave.
#
# Protocol: inst/validation/GAZEWEAVE-PROBE-DELAY-CONTROLS.md
# Run from a source checkout after devtools::load_all():
#   source("inst/validation/gaze-weave-probe-delay.R")
#   run_gaze_weave_probe_delay_court(
#     "inst/validation/gaze-weave-probe-delay-results"
#   )

probe_delay_protocol_version <- 2L
probe_delay_protocol_seed <- 20260820L
probe_delay_expected_items <- list(
  `1015` = c(6, 17, 25, 48, 62, 68, 71, 95, 98, 116),
  `1018` = c(5, 27, 51, 62, 66, 68, 69, 98, 106, 117),
  `1019` = c(9, 39, 43, 53, 72, 77, 78, 97, 119, 120),
  `1026` = c(13, 17, 18, 24, 27, 40, 59, 67, 100, 108),
  `1038` = c(29, 35, 79, 83, 88, 89, 94, 96, 105, 112),
  `1039` = c(13, 20, 27, 35, 69, 75, 84, 90, 99, 112),
  `1042` = c(6, 22, 29, 41, 62, 73, 93, 100, 106, 111),
  `1044` = c(7, 8, 29, 65, 85, 93, 95, 101, 112, 114)
)
probe_delay_expected_fold1_items <- list(
  `1015` = c(6, 48, 62, 95, 98),
  `1018` = c(5, 62, 66, 69, 98),
  `1019` = c(53, 72, 77, 119, 120),
  `1026` = c(13, 24, 40, 67, 100),
  `1038` = c(29, 79, 94, 96, 105),
  `1039` = c(13, 20, 69, 90, 99),
  `1042` = c(6, 29, 62, 73, 100),
  `1044` = c(7, 29, 95, 101, 114)
)
probe_delay_study_window <- c(0, 2500)
probe_delay_windows <- list(
  probe = c(0, 500),
  delay = c(500, 3000),
  combined = c(0, 3000),
  early_delay = c(500, 1500),
  late_delay = c(1500, 3000)
)
probe_delay_screen_bounds <- c(xmin = 112, xmax = 912, ymin = 84, ymax = 684)
probe_delay_min_duration <- 80

probe_delay_helper_file <- system.file(
  "validation", "gaze-weave-real-controls.R", package = "eyesim"
)
if (!nzchar(probe_delay_helper_file)) {
  probe_delay_helper_file <- file.path(
    "inst", "validation", "gaze-weave-real-controls.R"
  )
}
if (!file.exists(probe_delay_helper_file)) {
  stop("Run the probe-delay court from the eyesim source root.")
}
source(probe_delay_helper_file, local = TRUE)

probe_delay_data_files <- function(
    data_dir = file.path("test_data", "wynn_probe_delay")) {
  c(
    study = file.path(data_dir, "study_fix_input_new.csv"),
    retrieval = file.path(data_dir, "testdelay_fix_input_matched.csv")
  )
}

probe_delay_sha256 <- function(path) {
  executable <- Sys.which("shasum")
  if (!nzchar(executable)) stop("shasum is required to verify frozen inputs.")
  result <- system2(
    executable, c("-a", "256", shQuote(normalizePath(path))),
    stdout = TRUE, stderr = TRUE, env = "LC_ALL=C"
  )
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop("Could not compute SHA-256 for ", path, ": ", paste(result, collapse = " "))
  }
  digest <- grep("^[[:xdigit:]]{64}[[:space:]]", result, value = TRUE)
  if (length(digest) != 1L) {
    stop("Could not identify one SHA-256 digest for ", path, ".")
  }
  tolower(sub("[[:space:]].*$", "", digest[[1L]]))
}

probe_delay_verify_inputs <- function(
    data_dir = file.path("test_data", "wynn_probe_delay")) {
  manifest_path <- file.path(data_dir, "manifest.csv")
  if (!file.exists(manifest_path)) stop("Probe-delay input manifest is missing.")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  required <- c("file", "bytes", "rows", "columns", "sha256")
  if (!all(required %in% names(manifest)) || nrow(manifest) != 2L) {
    stop("Probe-delay input manifest has an unexpected schema.")
  }
  paths <- file.path(data_dir, manifest$file)
  if (!all(file.exists(paths))) stop("One or more frozen probe-delay inputs are missing.")
  info <- file.info(paths)
  observed_sha <- vapply(paths, probe_delay_sha256, character(1))
  if (!identical(as.numeric(info$size), as.numeric(manifest$bytes)) ||
      !identical(unname(observed_sha), manifest$sha256)) {
    stop("Probe-delay input integrity check failed.")
  }
  manifest$path <- paths
  manifest
}

probe_delay_required_columns <- function() {
  list(
    study = c(
      "Subject", "Trial", "TrialTotal", "Repetition", "Saliency",
      "ImageVersion", "Run", "ImageNumber", "ImageSet", "FixX", "FixY",
      "FixDuration", "FixStartTime", "FixEndTime", "FixTrialOnset"
    ),
    retrieval = c(
      "Subject", "Trial", "TrialTotal", "Repetition", "Saliency",
      "ImageVersion", "Run", "ImageNumber", "ImageSet", "Response",
      "Accuracy", "FixX", "FixY", "DelayOnset", "FixStartTime",
      "FixDuration", "FixEndTime", "FixTrialOnset", "EpochSrc",
      "FixProbeOnset", "Epoch", "EpochStraddle"
    )
  )
}

probe_delay_read_csv <- function(path, kind) {
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (length(names(tab)) > 0L && names(tab)[[1L]] %in% c("", "X")) {
    tab <- tab[-1L]
  }
  required <- probe_delay_required_columns()[[kind]]
  missing <- setdiff(required, names(tab))
  if (length(missing) > 0L) {
    stop(kind, " input is missing columns: ", paste(missing, collapse = ", "))
  }
  tab
}

probe_delay_spatially_valid <- function(tab) {
  bounds <- probe_delay_screen_bounds
  is.finite(tab$FixX) & is.finite(tab$FixY) &
    tab$FixX >= bounds[["xmin"]] & tab$FixX <= bounds[["xmax"]] &
    tab$FixY >= bounds[["ymin"]] & tab$FixY <= bounds[["ymax"]]
}

# Intersect fixation intervals with one half-open analysis window. This pure
# helper is the oracle target for boundary and refinement tests.
probe_delay_clip_rows <- function(tab, start, duration, window,
                                  minimum_duration = probe_delay_min_duration) {
  stopifnot(length(window) == 2L, window[[1L]] < window[[2L]])
  interval_start <- as.numeric(start)
  interval_end <- interval_start + as.numeric(duration)
  clipped_start <- pmax(interval_start, window[[1L]])
  clipped_end <- pmin(interval_end, window[[2L]])
  clipped_duration <- clipped_end - clipped_start
  keep <- is.finite(interval_start) & is.finite(interval_end) &
    probe_delay_spatially_valid(tab) & is.finite(clipped_duration) &
    clipped_duration > minimum_duration
  result <- tab[keep, , drop = FALSE]
  result$path_onset <- clipped_start[keep] - window[[1L]]
  result$path_duration <- clipped_duration[keep]
  result
}

probe_delay_assign_presentations <- function(trials) {
  required <- c("Subject", "ImageNumber", "Run", "TrialTotal", "Trial")
  if (!all(required %in% names(trials))) stop("Study trial keys are incomplete.")
  key <- interaction(
    trials$Subject, trials$ImageNumber, drop = TRUE, lex.order = TRUE
  )
  trials$presentation <- NA_integer_
  for (index in split(seq_len(nrow(trials)), key)) {
    order_index <- order(
      trials$Run[index], trials$TrialTotal[index], trials$Trial[index],
      na.last = TRUE
    )
    trials$presentation[index[order_index]] <- seq_along(index)
  }
  trials
}

probe_delay_path <- function(rows) {
  rows <- rows[order(rows$path_onset, rows$FixEndTime, seq_len(nrow(rows))),
               , drop = FALSE]
  onset <- as.numeric(rows$path_onset)
  if (length(onset) > 1L && any(diff(onset) <= 0)) {
    onset <- onset + seq_along(onset) * 1e-7
  }
  fixation_group(
    x = rows$FixX - probe_delay_screen_bounds[["xmin"]],
    y = rows$FixY - probe_delay_screen_bounds[["ymin"]],
    duration = rows$path_duration,
    onset = onset
  )
}

probe_delay_split_paths <- function(tab, keys) {
  key <- do.call(interaction, c(tab[keys], list(drop = TRUE, lex.order = TRUE)))
  rows <- lapply(split(seq_len(nrow(tab)), key), function(index) {
    part <- tab[index, , drop = FALSE]
    metadata <- part[1L, keys, drop = FALSE]
    metadata$fixgroup <- list(probe_delay_path(part))
    metadata$nfix <- nrow(metadata$fixgroup[[1L]])
    metadata
  })
  dplyr::bind_rows(rows)
}

probe_delay_study_paths <- function(study) {
  trial_keys <- c(
    "Subject", "Trial", "TrialTotal", "Repetition", "Saliency",
    "ImageVersion", "Run", "ImageNumber", "ImageSet"
  )
  trials <- unique(study[trial_keys])
  trials <- probe_delay_assign_presentations(trials)
  join_key <- interaction(
    trials$Subject, trials$Trial, trials$Run, trials$ImageNumber,
    drop = TRUE, lex.order = TRUE
  )
  row_key <- interaction(
    study$Subject, study$Trial, study$Run, study$ImageNumber,
    drop = TRUE, lex.order = TRUE
  )
  matched <- match(row_key, join_key)
  if (anyNA(matched)) stop("Could not map study rows to presentation order.")
  study$presentation <- trials$presentation[matched]
  study <- probe_delay_clip_rows(
    study, study$FixTrialOnset, study$FixDuration, probe_delay_study_window
  )
  paths <- probe_delay_split_paths(
    study,
    c(
      "Subject", "Trial", "TrialTotal", "Repetition", "Saliency",
      "ImageVersion", "Run", "ImageNumber", "ImageSet", "presentation"
    )
  )
  names(paths)[names(paths) == "Subject"] <- "participant"
  names(paths)[names(paths) == "ImageNumber"] <- "item"
  names(paths)[names(paths) == "ImageVersion"] <- "study_image_version"
  names(paths)[names(paths) == "Repetition"] <- "probe_type"
  names(paths)[names(paths) == "Saliency"] <- "degradation"
  paths$participant <- as.character(paths$participant)
  paths$item <- as.integer(paths$item)
  paths$age <- "all"
  paths
}

probe_delay_retrieval_paths <- function(retrieval, phase) {
  window <- probe_delay_windows[[phase]]
  if (is.null(window)) stop("Unknown probe-delay phase: ", phase)
  clipped <- probe_delay_clip_rows(
    retrieval, retrieval$FixProbeOnset, retrieval$FixDuration, window
  )
  paths <- probe_delay_split_paths(
    clipped,
    c(
      "Subject", "Trial", "TrialTotal", "Repetition", "Saliency",
      "ImageVersion", "Run", "ImageNumber", "ImageSet", "Response",
      "Accuracy", "EpochSrc"
    )
  )
  names(paths)[names(paths) == "Subject"] <- "participant"
  names(paths)[names(paths) == "ImageNumber"] <- "item"
  names(paths)[names(paths) == "ImageVersion"] <- "retrieval_image_version"
  names(paths)[names(paths) == "Repetition"] <- "probe_type"
  names(paths)[names(paths) == "Saliency"] <- "degradation"
  paths$participant <- as.character(paths$participant)
  paths$item <- as.integer(paths$item)
  paths$age <- "all"
  paths$phase <- phase
  paths$cue_duration <- window[[2L]] - window[[1L]]
  paths
}

probe_delay_trial_tables <- function(study, retrieval) {
  study_keys <- c(
    "Subject", "Trial", "TrialTotal", "Repetition", "Saliency",
    "ImageVersion", "Run", "ImageNumber", "ImageSet"
  )
  retrieval_keys <- c(
    study_keys, "Response", "Accuracy"
  )
  list(
    study = unique(study[study_keys]),
    retrieval = unique(retrieval[retrieval_keys])
  )
}

probe_delay_read_inputs <- function(
    data_dir = file.path("test_data", "wynn_probe_delay"),
    verify = TRUE) {
  if (verify) probe_delay_verify_inputs(data_dir)
  files <- probe_delay_data_files(data_dir)
  study_raw <- probe_delay_read_csv(files[["study"]], "study")
  retrieval_raw <- probe_delay_read_csv(files[["retrieval"]], "retrieval")
  if (nrow(study_raw) != 101991L || nrow(retrieval_raw) != 29699L) {
    stop("Probe-delay row counts differ from the frozen inputs.")
  }
  if (length(unique(study_raw$Subject)) != 46L ||
      length(unique(retrieval_raw$Subject)) != 46L) {
    stop("Probe-delay participant counts differ from the frozen inputs.")
  }
  phases <- lapply(names(probe_delay_windows), function(phase) {
    probe_delay_retrieval_paths(retrieval_raw, phase)
  })
  names(phases) <- names(probe_delay_windows)
  list(
    study = probe_delay_study_paths(study_raw),
    retrieval = phases,
    trials = probe_delay_trial_tables(study_raw, retrieval_raw)
  )
}

probe_delay_complete_pairs <- function(raw) {
  study <- raw$study
  study_groups <- split(
    seq_len(nrow(study)),
    interaction(study$participant, study$item, drop = TRUE)
  )
  study_pairs <- dplyr::bind_rows(lapply(study_groups, function(index) {
    part <- study[index, , drop = FALSE]
    data.frame(
      participant = part$participant[[1L]], item = part$item[[1L]],
      four_presentations = identical(sort(unique(part$presentation)), 1:4),
      study_min_nfix = min(part$nfix), stringsAsFactors = FALSE
    )
  }))
  phase_requirement <- c(
    probe = 1L, delay = 3L, combined = 3L,
    early_delay = 1L, late_delay = 1L
  )
  result <- study_pairs
  for (phase in names(phase_requirement)) {
    paths <- raw$retrieval[[phase]]
    paths <- paths[paths$probe_type %in% c("old", "lure"), , drop = FALSE]
    keep <- paths$nfix >= phase_requirement[[phase]]
    phase_pairs <- unique(paths[keep, c("participant", "item")])
    phase_pairs[[paste0(phase, "_ok")]] <- TRUE
    result <- merge(result, phase_pairs, by = c("participant", "item"), all = FALSE)
  }
  result[
    result$four_presentations & result$study_min_nfix >= 3L,
    , drop = FALSE
  ]
}

probe_delay_select_cohort <- function(raw, seed = probe_delay_protocol_seed,
                                      participants_per_design = 2L,
                                      items_per_fold = 5L) {
  study_trials <- raw$trials$study
  retrieval_trials <- raw$trials$retrieval
  study_n <- table(study_trials$Subject)
  retrieval_n <- table(retrieval_trials$Subject)
  complete <- intersect(
    names(study_n)[study_n == 360L],
    names(retrieval_n)[retrieval_n == 90L]
  )
  signal_trials <- retrieval_trials[
    retrieval_trials$Repetition %in% c("old", "lure") &
      as.character(retrieval_trials$Subject) %in% complete,
    , drop = FALSE
  ]
  signature <- tapply(
    signal_trials$ImageNumber, signal_trials$Subject,
    function(value) paste(sort(unique(value)), collapse = ";")
  )
  signature_levels <- sort(unique(signature))
  signature_id <- stats::setNames(
    match(signature, signature_levels), names(signature)
  )
  pairs <- unique(probe_delay_complete_pairs(raw)[c("participant", "item")])

  # Assign item identities globally before selecting trials. Every participant
  # then contributes exactly five eligible items to each held-out item fold,
  # without requiring an infeasible common item rectangle.
  set.seed(seed)
  item_order <- sample(sort(unique(pairs$item)))
  item_map <- data.frame(
    item = item_order,
    item_fold = rep(1:2, length.out = length(item_order)),
    stringsAsFactors = FALSE
  )
  pair_folds <- merge(pairs, item_map, by = "item", sort = FALSE)
  has_fold_balance <- vapply(complete, function(participant) {
    counts <- table(factor(
      pair_folds$item_fold[pair_folds$participant == participant],
      levels = 1:2
    ))
    all(counts >= items_per_fold)
  }, logical(1))
  eligible_participants <- complete[has_fold_balance]

  set.seed(seed + 1L)
  participants <- unlist(lapply(seq_along(signature_levels), function(group) {
    candidates <- sort(intersect(
      names(signature_id)[signature_id == group], eligible_participants
    ))
    if (length(candidates) < participants_per_design) {
      stop("A counterbalancing group lacks enough fixation-eligible participants.")
    }
    sample(candidates, participants_per_design)
  }), use.names = FALSE)
  participants <- sort(participants)

  set.seed(seed + 2L)
  selected_pairs <- dplyr::bind_rows(lapply(participants, function(participant) {
    dplyr::bind_rows(lapply(1:2, function(item_fold) {
      candidates <- sort(pair_folds$item[
        pair_folds$participant == participant &
          pair_folds$item_fold == item_fold
      ])
      data.frame(
        participant = participant,
        item = sort(sample(candidates, items_per_fold)),
        item_fold = item_fold,
        stringsAsFactors = FALSE
      )
    }))
  }))
  selected_pairs$design_signature <- unname(
    signature_id[selected_pairs$participant]
  )
  selected_meta <- unique(raw$retrieval$combined[c(
    "participant", "item", "probe_type", "degradation"
  )])
  selected_meta <- merge(
    selected_pairs, selected_meta,
    by = c("participant", "item"), all.x = TRUE, sort = FALSE
  )
  list(
    participants = participants,
    pairs = selected_pairs,
    item_map = item_map,
    signature_id = signature_id,
    complete_participant_n = length(complete),
    eligible_participant_n = length(eligible_participants),
    design_group_n = length(signature_levels),
    balance = list(
      probe_type = table(selected_meta$probe_type),
      degradation = table(selected_meta$degradation)
    )
  )
}

probe_delay_fold_plan <- function(cohort, seed = probe_delay_protocol_seed) {
  participant_signature <- unique(cohort$pairs[c(
    "participant", "design_signature"
  )])
  set.seed(seed + 3L)
  participant_map <- dplyr::bind_rows(lapply(
    sort(unique(participant_signature$design_signature)), function(group) {
      participants <- sample(participant_signature$participant[
        participant_signature$design_signature == group
      ])
      data.frame(
        participant = participants,
        participant_fold = rep(1:2, length.out = length(participants)),
        design_signature = group,
        stringsAsFactors = FALSE
      )
    }
  ))
  item_map <- cohort$item_map
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

probe_delay_task_tables <- function(raw, cohort, task,
                                    seed = probe_delay_protocol_seed) {
  pair_keys <- cohort$pairs[c("participant", "item", "item_fold")]
  study <- merge(
    raw$study, pair_keys, by = c("participant", "item"),
    all = FALSE, sort = FALSE
  )
  retrieval_phase <- if (identical(task, "repeated_viewing")) {
    raw$retrieval$combined
  } else {
    raw$retrieval[[task]]
  }
  if (is.null(retrieval_phase)) stop("Unknown probe-delay task: ", task)
  retrieval_phase <- merge(
    retrieval_phase, pair_keys, by = c("participant", "item"),
    all = FALSE, sort = FALSE
  )
  retrieval_phase <- retrieval_phase[
    retrieval_phase$probe_type %in% c("old", "lure"), , drop = FALSE
  ]
  metadata <- retrieval_phase[
    !duplicated(retrieval_phase[c("participant", "item")]),
    c(
      "participant", "age", "item", "probe_type", "degradation",
      "cue_duration", "Accuracy", "retrieval_image_version", "phase"
    )
  ]
  names(metadata)[names(metadata) == "Accuracy"] <- "accuracy"
  reference_rep <- if (identical(task, "repeated_viewing")) 1L else 4L
  reference <- study[
    study$presentation == reference_rep,
    c(
      "participant", "age", "item", "item_fold", "study_image_version",
      "nfix", "fixgroup"
    )
  ]
  reference <- merge(
    reference, metadata, by = c("participant", "age", "item"),
    all.x = TRUE, sort = FALSE
  )
  if (identical(task, "repeated_viewing")) {
    signal <- study[
      study$presentation == 4L,
      c(
        "participant", "age", "item", "item_fold", "study_image_version",
        "nfix", "fixgroup"
      )
    ]
    signal <- merge(
      signal, metadata, by = c("participant", "age", "item"),
      all.x = TRUE, sort = FALSE
    )
    signal$phase <- "repeated_viewing"
  } else {
    signal <- retrieval_phase
    names(signal)[names(signal) == "Accuracy"] <- "accuracy"
  }
  reference$task <- task
  signal$task <- task
  reference <- reference[order(reference$participant, reference$item), ]
  signal <- signal[order(signal$participant, signal$item), ]
  rownames(reference) <- rownames(signal) <- NULL
  if (nrow(reference) != nrow(cohort$pairs) || nrow(signal) != nrow(reference)) {
    stop("Task tables do not match the frozen participant-item pairs for task ", task, ".")
  }

  signal_lookup <- split(
    seq_len(nrow(signal)),
    paste(signal$participant, signal$item, sep = "\r")
  )
  generic <- lapply(split(signal$fixgroup, signal$participant), real_generic_path)
  controls <- lapply(seq_len(nrow(signal)), function(i) {
    row <- signal[i, , drop = FALSE]
    participant <- row$participant[[1L]]
    item <- row$item[[1L]]
    item_fold <- row$item_fold[[1L]]
    participant_items <- sort(signal$item[
      signal$participant == participant & signal$item_fold == item_fold
    ])
    item_position <- match(item, participant_items)
    next_item <- participant_items[(item_position %% length(participant_items)) + 1L]
    wrong_key <- paste(participant, next_item, sep = "\r")
    wrong_row <- signal_lookup[[wrong_key]]
    if (is.null(wrong_row) || length(wrong_row) != 1L) {
      stop("Could not construct the frozen wrong-item control.")
    }
    paths <- list(
      signal = row$fixgroup[[1L]],
      shuffled_order = real_shuffle_path(
        row$fixgroup[[1L]],
        seed + as.integer(participant) * 1009L + item * 101L +
          match(task, c(
            "repeated_viewing", "probe", "delay", "combined",
            "early_delay", "late_delay"
          ))
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
  list(
    reference = tibble::as_tibble(reference),
    signal = tibble::as_tibble(signal),
    source = dplyr::bind_rows(controls)
  )
}

probe_delay_specs <- function(task, smoke = FALSE) {
  specs <- real_specs(smoke)
  if (task %in% c("probe", "early_delay", "late_delay")) {
    specs$baseline <- gaze_baseline_spec(
      screen = specs$screen,
      density_sigmas = c(30, 60, 120),
      density_grid = if (smoke) 12L else 24L,
      methods = c("density", "elastic"),
      warp = specs$warp,
      lambda_grid = c(0.01, 0.1, 1, 10),
      inner_folds = 2L,
      elastic_radii = c(consensus = 80, rigidity = 400, matching = 40),
      elastic_maxit = if (smoke) 10L else 25L,
      elastic_tolerance = 1e-4
    )
  }
  specs
}

probe_delay_persistence_verdict <- function(summary, bootstrap, operating,
                                            fold_audit, scored) {
  engines <- c("transport_v2", "replay")
  signal <- summary[
    summary$condition == "signal" & summary$method %in% engines,
    , drop = FALSE
  ]
  intervals <- bootstrap[
    bootstrap$condition == "signal" & bootstrap$method %in% engines &
      bootstrap$metric == "gaze_info_bits",
    , drop = FALSE
  ]
  merged <- merge(
    signal[c("task", "method", "mean_info_bits")],
    intervals[c("task", "method", "lower_95", "upper_95")],
    by = c("task", "method"), all = TRUE
  )
  delay <- merged[merged$task == "delay", ]
  late <- merged[merged$task == "late_delay", ]
  joint <- merge(delay, late, by = "method", suffixes = c("_delay", "_late"))
  null_ok <- operating[
    operating$task == "delay" & operating$method %in% engines,
    , drop = FALSE
  ]
  safe <- null_ok$method[null_ok$null_error <= 0.075]
  supported <- joint$method[
    joint$lower_95_delay > 0 & joint$lower_95_late > 0 &
      joint$method %in% safe
  ]
  suggestive <- joint$method[
    joint$mean_info_bits_delay > 0 & joint$mean_info_bits_late > 0 &
      joint$method %in% safe
  ]
  leakage_ok <- all(fold_audit$participant_overlap == 0L) &&
    all(fold_audit$item_overlap == 0L)
  convergence_ok <- all(vapply(engines, function(method) {
    values <- scored$converged[scored$method == method]
    length(values) > 0L && mean(values) >= 0.99
  }, logical(1)))
  label <- if (!leakage_ok || !convergence_ok) {
    "not_supported"
  } else if (length(supported) > 0L) {
    "supported"
  } else if (length(suggestive) > 0L) {
    "suggestive"
  } else {
    "not_supported"
  }
  data.frame(
    verdict = label,
    supported_engines = paste(supported, collapse = ";"),
    suggestive_engines = paste(suggestive, collapse = ";"),
    leakage_safe = leakage_ok,
    engine_convergence = convergence_ok,
    independent_replication = FALSE,
    closes_wang_gate = FALSE,
    stringsAsFactors = FALSE
  )
}

probe_delay_path_coverage <- function(tables) {
  dplyr::bind_rows(lapply(names(tables), function(task) {
    reference <- tables[[task]]$reference[c("participant", "item", "nfix")]
    signal <- tables[[task]]$signal[c("participant", "item", "nfix")]
    paired <- merge(
      reference, signal, by = c("participant", "item"),
      suffixes = c("_reference", "_source"), sort = FALSE
    )
    data.frame(
      task = task,
      n = nrow(paired),
      source_min_nfix = min(paired$nfix_source),
      source_median_nfix = stats::median(paired$nfix_source),
      source_max_nfix = max(paired$nfix_source),
      source_ge_3_n = sum(paired$nfix_source >= 3L),
      multimatch_pair_coverage = mean(
        paired$nfix_reference >= 3L & paired$nfix_source >= 3L
      ),
      stringsAsFactors = FALSE
    )
  }))
}

probe_delay_configuration <- function(result) {
  cohort <- result$config$cohort
  pair_receipt <- paste(
    paste(cohort$pairs$participant, cohort$pairs$item,
          cohort$pairs$item_fold, sep = ":"),
    collapse = ";"
  )
  data.frame(
    protocol_version = probe_delay_protocol_version,
    seed = result$config$seed,
    smoke = result$config$smoke,
    participants = paste(cohort$participants, collapse = ";"),
    participant_item_fold_pairs = pair_receipt,
    complete_participant_n = cohort$complete_participant_n,
    eligible_participant_n = cohort$eligible_participant_n,
    design_group_n = cohort$design_group_n,
    probe_type_balance = paste(
      names(cohort$balance$probe_type), as.integer(cohort$balance$probe_type),
      sep = ":", collapse = ";"
    ),
    saliency_balance = paste(
      names(cohort$balance$degradation),
      as.integer(cohort$balance$degradation), sep = ":", collapse = ";"
    ),
    verdict = result$verdict$verdict,
    independent_replication = FALSE,
    closes_wang_gate = FALSE,
    stringsAsFactors = FALSE
  )
}

probe_delay_write_results <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  write(probe_delay_configuration(result), "configuration.csv")
  write(result$verdict, "persistence-verdict.csv")
  write(result$method_summary, "method-summary.csv")
  write(result$frozen_summary, "frozen-summary.csv")
  write(result$bootstrap, "bootstrap-intervals.csv")
  write(result$operating, "operating-characteristics.csv")
  write(result$order_checks, "order-controls.csv")
  write(result$path_coverage, "path-coverage.csv")
  write(result$diagnostics, "alignment-diagnostics.csv")
  write(result$warps, "warp-estimates.csv")
  write(result$fold_audit, "fold-audit.csv")
  write(result$covariates, "covariate-sensitivity.csv")
  write(result$resources, "resources.csv")
  score_columns <- setdiff(names(result$scored), c("alignment", "candidates"))
  write(result$scored[score_columns], "trial-scores.csv")
  saveRDS(
    result[c(
      "method_summary", "frozen_summary", "bootstrap", "operating",
      "order_checks", "path_coverage", "diagnostics", "warps",
      "fold_audit", "covariates", "verdict", "resources", "config"
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

run_gaze_weave_probe_delay_court <- function(
    output_dir = NULL, seed = probe_delay_protocol_seed, smoke = FALSE,
    bootstrap_draws = if (smoke) 100L else 1000L,
    fold_ids = NULL,
    tasks = c(
      "repeated_viewing", "probe", "delay", "combined",
      "early_delay", "late_delay"
    )) {
  if (!smoke && seed != probe_delay_protocol_seed) {
    stop("The probe-delay court is frozen to seed 20260820; use smoke = TRUE for debugging.")
  }
  if (!smoke && !is.null(fold_ids)) {
    stop("Non-smoke runs must evaluate all four frozen outer folds.")
  }
  raw <- probe_delay_read_inputs(verify = !smoke)
  cohort <- probe_delay_select_cohort(raw, seed)
  observed_items <- lapply(
    split(cohort$pairs$item, cohort$pairs$participant),
    function(value) as.numeric(sort(value))
  )
  observed_fold1 <- lapply(split(
    cohort$pairs$item[cohort$pairs$item_fold == 1L],
    cohort$pairs$participant[cohort$pairs$item_fold == 1L]
  ), function(value) as.numeric(sort(value)))
  if (!smoke &&
      (!identical(cohort$participants, names(probe_delay_expected_items)) ||
       !identical(observed_items, probe_delay_expected_items) ||
       !identical(observed_fold1, probe_delay_expected_fold1_items))) {
    stop("The data-only cohort no longer matches the frozen protocol.")
  }
  tables <- lapply(tasks, function(task) {
    probe_delay_task_tables(raw, cohort, task, seed)
  })
  names(tables) <- tasks
  plan <- probe_delay_fold_plan(cohort, seed)
  folds_to_run <- if (is.null(fold_ids)) {
    plan$folds
  } else {
    plan$folds[vapply(plan$folds, function(fold) fold$id %in% fold_ids, logical(1))]
  }
  if (length(folds_to_run) == 0L) stop("fold_ids selected no outer folds.")
  detected <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected) || detected < 1L) detected <- 2L
  workers <- if (smoke) min(2L, detected) else min(4L, detected)
  fold_results <- list()
  index <- 1L
  for (task in tasks) {
    specs <- probe_delay_specs(task, smoke)
    for (fold in folds_to_run) {
      message("Probe-delay court: ", task, ", outer fold ", fold$id, "/4")
      conditions <- if (task %in% c("early_delay", "late_delay")) "signal" else
        c("signal", "shuffled_order", "wrong_item", "generic_gaze")
      fold_results[[index]] <- real_score_fold(
        tables[[task]], fold, specs, workers = workers,
        conditions = conditions, seed = seed
      )
      index <- index + 1L
    }
  }
  scored <- dplyr::bind_rows(lapply(fold_results, `[[`, "scored"))
  resources <- dplyr::bind_rows(lapply(fold_results, `[[`, "resources"))
  warps <- dplyr::bind_rows(lapply(fold_results, `[[`, "warps"))
  fold_audit <- real_fold_audit_table(lapply(fold_results, `[[`, "audit"))
  method_summary <- real_method_summary(scored)
  frozen_summary <- real_frozen_summary(scored)
  bootstrap <- real_bootstrap_summary(scored, bootstrap_draws, seed)
  operating <- real_operating_characteristics(
    scored[!scored$task %in% c("early_delay", "late_delay"), , drop = FALSE]
  )
  order_checks <- real_order_checks(scored)
  path_coverage <- probe_delay_path_coverage(tables)
  diagnostics <- real_diagnostic_summary(scored)
  covariates <- real_covariate_summary(scored)
  verdict <- probe_delay_persistence_verdict(
    method_summary, bootstrap, operating, fold_audit, scored
  )
  result <- list(
    scored = scored, method_summary = method_summary,
    frozen_summary = frozen_summary, bootstrap = bootstrap,
    operating = operating, order_checks = order_checks,
    path_coverage = path_coverage, diagnostics = diagnostics,
    warps = warps, fold_audit = fold_audit,
    covariates = covariates, verdict = verdict, resources = resources,
    config = list(
      protocol_version = probe_delay_protocol_version, seed = seed,
      smoke = smoke, bootstrap_draws = bootstrap_draws, cohort = cohort,
      fold_plan = plan,
      evaluated_fold_ids = vapply(folds_to_run, `[[`, integer(1), "id"),
      workers = workers, tasks = tasks
    )
  )
  if (!is.null(output_dir)) probe_delay_write_results(result, output_dir)
  result
}
