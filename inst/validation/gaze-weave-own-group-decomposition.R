# Own-versus-other-participant Transport decomposition for recognition gaze.

own_group_seed <- 20260830L
own_group_draws <- 2000L
own_group_result_dir <- file.path(
  "inst", "validation", "gaze-weave-transport-v3-retrieval-results",
  "own-group"
)

own_group_dependency <- system.file(
  "validation", "gaze-weave-transport-v3-retrieval.R", package = "eyesim"
)
if (!nzchar(own_group_dependency)) {
  own_group_dependency <- file.path(
    "inst", "validation", "gaze-weave-transport-v3-retrieval.R"
  )
}
if (!file.exists(own_group_dependency)) {
  stop("The canonical Transport retrieval helpers are unavailable.")
}
source(own_group_dependency, local = TRUE)

own_group_item_dependency <- system.file(
  "validation", "gaze-weave-item-effects.R", package = "eyesim"
)
if (!nzchar(own_group_item_dependency)) {
  own_group_item_dependency <- file.path(
    "inst", "validation", "gaze-weave-item-effects.R"
  )
}
if (!file.exists(own_group_item_dependency)) {
  stop("The crossed-bootstrap helpers are unavailable.")
}
source(own_group_item_dependency, local = TRUE)

own_group_trial_key <- function(participant, item) {
  paste(as.character(participant), as.integer(item), sep = ":")
}

own_group_candidate_margin <- function(candidates) {
  if (!is.data.frame(candidates) ||
      !all(c("log_score", "is_true") %in% names(candidates)) ||
      sum(candidates$is_true) != 1L || sum(!candidates$is_true) < 1L ||
      any(!is.finite(candidates$log_score))) {
    stop("Candidate evidence must contain one finite truth and nonmatches.")
  }
  mean(candidates$log_score[candidates$is_true]) -
    mean(candidates$log_score[!candidates$is_true])
}

own_group_add_margin <- function(scored) {
  scored$candidate_margin <- vapply(
    scored$candidates, own_group_candidate_margin, numeric(1)
  )
  scored$trial_key <- own_group_trial_key(scored$participant, scored$item)
  scored
}

own_group_warp_from_info <- function(info) {
  if (!is.list(info) || !identical(info$type, "contraction") ||
      length(info$groups) != 1L) {
    stop("Only a one-group frozen contraction warp can be reconstructed.")
  }
  group <- info$groups[[1L]]
  group_key <- as.character(group$group)
  group_model <- list(
    A = diag(as.numeric(group$scale), 2L),
    translation = as.numeric(group$translation),
    center = as.numeric(info$center),
    shift = as.numeric(group$shift),
    scale = as.numeric(group$scale),
    matched_pairs = as.integer(group$matched_pairs)
  )
  structure(
    list(
      type = "contraction", fit_by = info$fit_by,
      group_models = stats::setNames(list(group_model), group_key),
      info = info
    ),
    class = c("gaze_warp_model", "list")
  )
}

own_group_spec_from_frozen <- function(spec) {
  if (!inherits(spec, "gaze_transport_v3_spec") ||
      !identical(spec$version, 3L) ||
      !identical(spec$estimand, "edge_normalized_episode_transport")) {
    stop("The frozen specification is not the canonical Transport-v3 estimand.")
  }
  class(spec) <- c("gaze_transport_spec", "list")
  spec
}

own_group_complete_study <- function(study, minimum_fixations = 3L) {
  key <- interaction(
    study$participant, study$item, study$study_image_version,
    drop = TRUE, lex.order = TRUE
  )
  keep <- unlist(lapply(split(seq_len(nrow(study)), key), function(index) {
    part <- study[index, , drop = FALSE]
    complete <- identical(sort(unique(part$presentation)), 1:4) &&
      all(part$nfix >= minimum_fixations)
    if (complete) index else integer(0)
  }), use.names = FALSE)
  result <- study[sort(keep), , drop = FALSE]
  if (!nrow(result) || any(table(interaction(
    result$participant, result$item, result$study_image_version,
    drop = TRUE, lex.order = TRUE
  )) != 4L)) {
    stop("Complete study donors must contribute exactly four presentations.")
  }
  result
}

own_group_rotate <- function(value, offset) {
  if (!length(value)) return(value)
  start <- as.integer(offset %% length(value)) + 1L
  value[c(start:length(value), if (start > 1L) seq_len(start - 1L))]
}

own_group_select_donors <- function(
    study, source_participant, candidate_item, image_version,
    seed = own_group_seed) {
  part <- study[
    study$item == candidate_item &
      as.character(study$study_image_version) == as.character(image_version) &
      as.character(study$participant) != as.character(source_participant),
    , drop = FALSE
  ]
  donors <- sort(unique(as.character(part$participant)))
  donors <- donors[vapply(donors, function(donor) {
    identical(sort(unique(part$presentation[
      as.character(part$participant) == donor
    ])), 1:4)
  }, logical(1))]
  if (!length(donors)) {
    stop("No other-participant donor is available.")
  }
  offset <- as.integer(seed) + sum(utf8ToInt(as.character(source_participant))) +
    17L * as.integer(candidate_item) +
    sum(utf8ToInt(as.character(image_version)))
  donor_order <- own_group_rotate(donors, offset)
  selected <- donor_order[(seq_len(4L) - 1L) %% length(donor_order) + 1L]
  rows <- lapply(1:4, function(episode) {
    row <- part[
      as.character(part$participant) == selected[[episode]] &
        part$presentation == episode,
      , drop = FALSE
    ]
    if (nrow(row) != 1L) stop("A selected donor episode is not unique.")
    row$donor_participant <- as.character(row$participant)
    row$episode_id <- episode
    row
  })
  result <- dplyr::bind_rows(rows)
  if (length(unique(result$donor_participant)) != min(4L, length(donors)) ||
      as.character(source_participant) %in% result$donor_participant ||
      !identical(result$episode_id, 1:4)) {
    stop("Other-participant donor selection violated independence.")
  }
  result
}

own_group_build_pool <- function(
    study, source_row, candidate_plan, version_lookup,
    seed = own_group_seed) {
  source_participant <- as.character(source_row$participant[[1L]])
  source_item <- as.integer(source_row$item[[1L]])
  candidate_set_id <- own_group_trial_key(source_participant, source_item)
  plan <- candidate_plan[
    candidate_plan$candidate_set_id == candidate_set_id, , drop = FALSE
  ]
  if (nrow(plan) != 5L || sum(plan$is_true) != 1L) {
    stop("Every candidate plan must contain one truth among five candidates.")
  }
  rows <- lapply(seq_len(nrow(plan)), function(index) {
    candidate_item <- as.integer(plan$candidate_item[[index]])
    image_version <- version_lookup(candidate_item)
    donor <- own_group_select_donors(
      study, source_participant, candidate_item, image_version, seed
    )
    donor$participant <- source_participant
    donor$item <- candidate_item
    donor$candidate_item <- candidate_item
    donor$candidate_set_id <- candidate_set_id
    donor$candidate_position <- plan$candidate_position[[index]]
    donor$is_true <- plan$is_true[[index]]
    donor$prior_weight <- 1
    donor
  })
  result <- dplyr::bind_rows(rows)
  result <- result[order(result$candidate_position, result$episode_id), ]
  if (nrow(result) != 20L ||
      any(table(result$candidate_item) != 4L) ||
      sum(result$is_true) != 4L || anyNA(result$fixgroup) ||
      source_participant %in% result$donor_participant) {
    stop("Other-participant candidate pool violates K=5 by four episodes.")
  }
  tibble::as_tibble(result)
}

own_group_existing_plan <- function(candidate_plan, source_row) {
  participant <- as.character(source_row$participant[[1L]])
  item <- as.integer(source_row$item[[1L]])
  candidate_plan[
    as.character(candidate_plan$participant) == participant &
      candidate_plan$target_item == item,
    , drop = FALSE
  ]
}

own_group_newtest_plan <- function(source_row, item_map) {
  participant <- as.character(source_row$participant[[1L]])
  item <- as.integer(source_row$item[[1L]])
  fold <- item_map$item_fold[match(item, item_map$item)]
  candidates <- sort(as.integer(item_map$item[item_map$item_fold == fold]))
  position <- match(item, candidates)
  rotated <- candidates[c(
    position:length(candidates),
    if (position > 1L) seq_len(position - 1L)
  )]
  selected <- c(item, head(rotated[rotated != item], 4L))
  data.frame(
    participant = participant, target_item = item, item_fold = fold,
    candidate_set_id = own_group_trial_key(participant, item),
    candidate_position = seq_along(selected), candidate_item = selected,
    is_true = selected == item, stringsAsFactors = FALSE
  )
}

own_group_version_lookup <- function(study, participant) {
  own <- unique(study[
    as.character(study$participant) == as.character(participant),
    c("item", "study_image_version")
  ])
  if (anyDuplicated(own$item)) {
    stop("A participant has multiple study versions for one base item.")
  }
  function(candidate_item) {
    value <- own$study_image_version[match(candidate_item, own$item)]
    if (length(value) != 1L || is.na(value)) {
      stop("A planned old/lure candidate lacks the participant study version.")
    }
    as.character(value)
  }
}

own_group_newtest_version_lookup <- function(source_row) {
  suffix <- sub("^[0-9]+", "", as.character(
    source_row$retrieval_image_version[[1L]]
  ))
  if (!suffix %in% c("A", "B")) stop("Unknown retrieval image version suffix.")
  function(candidate_item) paste0(as.integer(candidate_item), suffix)
}

own_group_score_one <- function(
    source_row, plan, study, spec, warp, calibration, version_lookup,
    seed = own_group_seed) {
  source_row$candidate_set_id <- own_group_trial_key(
    source_row$participant, source_row$item
  )
  pool <- own_group_build_pool(
    study, source_row, plan, version_lookup, seed
  )
  scored <- v3_retrieval_score_row(
    source_row, pool, spec, warp, calibration
  )
  scored$template_source <- "other_participants"
  scored$donor_participants <- I(list(as.character(pool$donor_participant)))
  scored$candidate_donor_counts <- I(list(vapply(
    split(pool$donor_participant, pool$candidate_position),
    function(value) length(unique(value)), integer(1)
  )))
  scored
}

own_group_checkpoint <- function(output_dir, family, fold_id) {
  file.path(output_dir, sprintf(
    "checkpoint-other-participants-%s-fold-%02d.rds", family, fold_id
  ))
}

own_group_verify_canonical_checkpoints <- function(
    output_dir = v3_retrieval_output_dir) {
  manifest <- utils::read.csv(
    file.path(output_dir, "checkpoint-manifest.csv"),
    stringsAsFactors = FALSE
  )
  paths <- file.path(output_dir, manifest$file)
  if (!all(file.exists(paths)) ||
      !identical(unname(tools::md5sum(paths)), manifest$md5)) {
    stop("Canonical Transport-v3 checkpoint hashes do not verify.")
  }
  invisible(TRUE)
}

own_group_context <- function() {
  own_group_verify_canonical_checkpoints()
  raw <- full_recognition_read_inputs(verify = TRUE)
  study <- own_group_complete_study(raw$study)
  private <- file.path(v3_retrieval_output_dir, "freeze-private")
  list(
    raw = raw, study = study,
    tables = readRDS(file.path(private, "score-blind-tables.rds")),
    candidate_plan = utils::read.csv(
      file.path(private, "candidate-plan.csv"), stringsAsFactors = FALSE
    ),
    participant_map = utils::read.csv(
      file.path(private, "participant-folds.csv"), stringsAsFactors = FALSE
    ),
    item_map = utils::read.csv(
      file.path(private, "item-folds.csv"), stringsAsFactors = FALSE
    ),
    spec = own_group_spec_from_frozen(
      readRDS(file.path(private, "spec.rds"))
    )
  )
}

own_group_newtest_source <- function(context, fold) {
  source <- context$raw$retrieval
  source <- source[
    source$probe_type == "newtest" & source$nfix >= 3L &
      as.character(source$participant) %in%
        as.character(context$participant_map$participant) &
      source$item %in% context$item_map$item,
    , drop = FALSE
  ]
  source$item_fold <- context$item_map$item_fold[
    match(source$item, context$item_map$item)
  ]
  source <- source[
    as.character(source$participant) %in%
      as.character(fold$eval_participants) &
      source$item %in% fold$eval_items,
    , drop = FALSE
  ]
  source <- source[order(source$participant, source$item), ]
  source$task <- "retrieval_0_3000"
  source$condition <- "newtest_other_participants"
  tibble::as_tibble(source)
}

own_group_fold_plan <- function(context) {
  cohort <- full_recognition_select_cohort(context$raw)
  full_recognition_fold_plan(cohort, v3_retrieval_seed)$folds
}

own_group_score_fold <- function(
    family = c("old_lure", "newtest"), fold_id,
    output_dir = own_group_result_dir, resume = TRUE,
    seed = own_group_seed) {
  family <- match.arg(family)
  path <- own_group_checkpoint(output_dir, family, fold_id)
  if (resume && file.exists(path)) {
    checkpoint <- readRDS(path)
    if (identical(checkpoint$protocol, "own-group/1.1.0")) {
      message("Using checkpoint: ", basename(path))
      return(invisible(checkpoint))
    }
    message("Replacing checkpoint from an earlier own-group protocol: ",
            basename(path))
  }
  context <- own_group_context()
  folds <- own_group_fold_plan(context)
  fold_position <- match(fold_id, vapply(folds, `[[`, integer(1), "id"))
  if (is.na(fold_position)) stop("Unknown outer fold.")
  fold <- folds[[fold_position]]
  canonical <- readRDS(v3_retrieval_checkpoint(
    v3_retrieval_output_dir, fold_id
  ))
  warp <- own_group_warp_from_info(canonical$warp)
  if (family == "old_lure") {
    source <- context$tables$source
    source <- source[
      as.character(source$participant) %in%
        as.character(fold$eval_participants) &
        source$item %in% fold$eval_items,
      , drop = FALSE
    ]
  } else {
    source <- own_group_newtest_source(context, fold)
  }
  started <- proc.time()[["elapsed"]]
  rows <- lapply(seq_len(nrow(source)), function(index) {
    if (index %% 25L == 0L) {
      message(family, " fold ", fold_id, ": ", index, "/", nrow(source))
    }
    source_row <- source[index, , drop = FALSE]
    if (family == "old_lure") {
      plan <- own_group_existing_plan(context$candidate_plan, source_row)
      version_lookup <- own_group_version_lookup(
        context$study, source_row$participant[[1L]]
      )
    } else {
      plan <- own_group_newtest_plan(source_row, context$item_map)
      version_lookup <- own_group_newtest_version_lookup(source_row)
    }
    own_group_score_one(
      source_row, plan, context$study, context$spec, warp,
      canonical$calibration, version_lookup, seed
    )
  })
  scored <- dplyr::bind_rows(rows)
  scored$outer_fold <- fold_id
  scored$probe_family <- family
  elapsed <- proc.time()[["elapsed"]] - started
  donor_overlap <- vapply(seq_len(nrow(scored)), function(index) {
    as.character(scored$participant[[index]]) %in%
      scored$donor_participants[[index]]
  }, logical(1))
  audit <- data.frame(
    family = family, outer_fold = fold_id, trials = nrow(scored),
    donor_overlap = sum(donor_overlap),
    min_donors = min(vapply(scored$donor_participants, length, integer(1))),
    max_donors = max(vapply(scored$donor_participants, length, integer(1))),
    min_unique_donors_per_candidate = min(unlist(
      scored$candidate_donor_counts, use.names = FALSE
    )),
    max_unique_donors_per_candidate = max(unlist(
      scored$candidate_donor_counts, use.names = FALSE
    )),
    candidate_count_min = min(scored$candidate_count),
    candidate_count_max = max(scored$candidate_count),
    episode_count_min = min(scored$common_episode_count),
    episode_count_max = max(scored$common_episode_count),
    convergence = sum(scored$converged_candidates) /
      sum(scored$candidate_alignments),
    elapsed_seconds = elapsed, stringsAsFactors = FALSE
  )
  if (audit$donor_overlap != 0L || audit$min_donors != 20L ||
      audit$max_donors != 20L || audit$candidate_count_min != 5L ||
      audit$candidate_count_max != 5L || audit$episode_count_min != 4L ||
      audit$episode_count_max != 4L) {
    stop("Other-participant scoring violated its support contract.")
  }
  result <- list(
    scored = scored, audit = audit, seed = seed,
    protocol = "own-group/1.1.0"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, path, version = 3)
  invisible(result)
}

run_gaze_weave_own_group_measurement <- function(
    output_dir = own_group_result_dir, families = c("old_lure", "newtest"),
    fold_ids = 1:4, resume = TRUE) {
  for (family in families) {
    for (fold_id in fold_ids) {
      message("Other-participant Transport: ", family, " fold ", fold_id)
      own_group_score_fold(family, fold_id, output_dir, resume)
    }
  }
  invisible(TRUE)
}

own_group_read_scores <- function(
    output_dir = own_group_result_dir, family = c("old_lure", "newtest")) {
  family <- match.arg(family)
  paths <- vapply(1:4, function(fold) {
    own_group_checkpoint(output_dir, family, fold)
  }, character(1))
  if (!all(file.exists(paths))) stop("Other-participant checkpoints incomplete.")
  values <- lapply(paths, readRDS)
  if (any(vapply(values, `[[`, integer(1), "seed") != own_group_seed) ||
      any(vapply(values, `[[`, character(1), "protocol") != "own-group/1.1.0")) {
    stop("Other-participant checkpoints belong to another protocol.")
  }
  list(
    scored = own_group_add_margin(dplyr::bind_rows(lapply(values, `[[`, "scored"))),
    audit = dplyr::bind_rows(lapply(values, `[[`, "audit")),
    paths = paths
  )
}

own_group_interval <- function(value, reference, draws, seed) {
  plan <- transport_exploratory_bootstrap_plan(reference, draws, seed)
  weights <- transport_exploratory_align_weights(reference, plan)
  bootstrap <- vapply(seq_len(ncol(weights)), function(draw) {
    transport_exploratory_weighted_mean(value, weights[, draw])
  }, numeric(1))
  interval <- transport_exploratory_interval(bootstrap, 0.95)
  c(
    estimate = mean(value), lower_95 = interval[["lower"]],
    upper_95 = interval[["upper"]], probability_le_zero = mean(bootstrap <= 0)
  )
}

own_group_condition_contrast <- function(
    first, second, label, draws, seed) {
  first$..condition <- "first"
  second$..condition <- "second"
  combined <- dplyr::bind_rows(first, second)
  plan <- transport_exploratory_bootstrap_plan(combined, draws, seed)
  weights <- transport_exploratory_align_weights(combined, plan)
  first_index <- combined$..condition == "first"
  bootstrap <- vapply(seq_len(ncol(weights)), function(draw) {
    transport_exploratory_weighted_mean(
      combined$candidate_margin[first_index], weights[first_index, draw]
    ) - transport_exploratory_weighted_mean(
      combined$candidate_margin[!first_index], weights[!first_index, draw]
    )
  }, numeric(1))
  interval <- transport_exploratory_interval(bootstrap, 0.95)
  data.frame(
    contrast = label, trials = nrow(combined),
    estimate = mean(first$candidate_margin) - mean(second$candidate_margin),
    lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
    probability_le_zero = mean(bootstrap <= 0), stringsAsFactors = FALSE
  )
}

own_group_measurement_summary <- function(
    own, group, newtest, metadata, draws = own_group_draws,
    seed = own_group_seed) {
  own <- own_group_add_margin(own)
  key_group <- match(own$trial_key, group$trial_key)
  if (anyNA(key_group)) stop("Own and group old/lure support do not align.")
  group <- group[key_group, , drop = FALSE]
  meta_key <- match(own$trial_key, metadata$trial_key)
  if (anyNA(meta_key)) stop("Old/lure metadata support does not align.")
  own$probe_type <- metadata$probe_type[meta_key]
  rows <- list()
  add <- function(label, value, reference = own) {
    summary <- own_group_interval(value, reference, draws, seed)
    data.frame(
      contrast = label, trials = length(value),
      estimate = summary[["estimate"]],
      lower_95 = summary[["lower_95"]],
      upper_95 = summary[["upper_95"]],
      probability_le_zero = summary[["probability_le_zero"]],
      stringsAsFactors = FALSE
    )
  }
  rows[[1L]] <- add("own_margin", own$candidate_margin)
  rows[[2L]] <- add("own_top1_minus_chance", own$top1_credit - 0.2)
  rows[[3L]] <- add("group_margin", group$candidate_margin)
  rows[[4L]] <- add("group_top1_minus_chance", group$top1_credit - 0.2)
  rows[[5L]] <- add(
    "own_minus_group_margin", own$candidate_margin - group$candidate_margin
  )
  rows[[6L]] <- add(
    "own_minus_group_top1", own$top1_credit - group$top1_credit
  )
  index <- 7L
  for (probe in c("old", "lure")) {
    keep <- own$probe_type == probe
    rows[[index]] <- add(
      paste0("own_minus_group_margin_", probe),
      own$candidate_margin[keep] - group$candidate_margin[keep],
      own[keep, , drop = FALSE]
    )
    index <- index + 1L
  }
  rows[[index]] <- add(
    "newtest_group_margin", newtest$candidate_margin, newtest
  )
  rows[[index + 1L]] <- add(
    "newtest_group_top1_minus_chance", newtest$top1_credit - 0.2, newtest
  )
  index <- index + 2L
  for (probe in c("old", "lure")) {
    keep <- own$probe_type == probe
    rows[[index]] <- own_group_condition_contrast(
      group[keep, , drop = FALSE], newtest,
      paste0("group_margin_", probe, "_minus_newtest"), draws, seed
    )
    index <- index + 1L
  }
  dplyr::bind_rows(rows)
}

own_group_behavior_fold <- function(tab, fold_id, probe_type) {
  evaluate <- tab[
    tab$outer_fold == fold_id & tab$probe_type == probe_type &
      !is.na(tab$oldness), , drop = FALSE
  ]
  train <- tab[
    !tab$participant %in% evaluate$participant &
      !tab$item %in% evaluate$item & tab$probe_type == probe_type &
      !is.na(tab$oldness), , drop = FALSE
  ]
  if (!nrow(evaluate) || !nrow(train)) return(NULL)
  train$log_effective_fixations <- log(train$effective_fixations)
  evaluate$log_effective_fixations <- log(evaluate$effective_fixations)
  mappings <- list(
    c("log_effective_fixations", "z_quality"),
    c("total_duration", "z_duration"),
    c("group_margin", "z_group"), c("own_margin", "z_own")
  )
  for (mapping in mappings) {
    value <- item_effect_safe_standardize(
      train, evaluate, mapping[[1L]], mapping[[2L]]
    )
    train <- value$train
    evaluate <- value$evaluate
  }
  base_formula <- oldness ~ saliency_z + z_quality + z_duration
  group_formula <- stats::update.formula(base_formula, ". ~ . + z_group")
  own_formula <- stats::update.formula(group_formula, ". ~ . + z_own")
  fits <- lapply(
    list(base = base_formula, group = group_formula, own_beyond_group = own_formula),
    transport_density_ordinal_fit, data = train
  )
  losses <- lapply(fits, transport_density_ordinal_loss, evaluate = evaluate)
  data.frame(
    participant = evaluate$participant, item = evaluate$item,
    trial_key = evaluate$trial_key, probe_type = evaluate$probe_type,
    oldness = as.integer(evaluate$oldness), outer_fold = fold_id,
    loss_base = losses$base, loss_group = losses$group,
    loss_own_beyond_group = losses$own_beyond_group,
    stringsAsFactors = FALSE
  )
}

own_group_behavior <- function(tab, draws = own_group_draws,
                               seed = own_group_seed) {
  predictions <- dplyr::bind_rows(lapply(c("old", "lure"), function(probe) {
    dplyr::bind_rows(lapply(1:4, function(fold) {
      own_group_behavior_fold(tab, fold, probe)
    }))
  }))
  rows <- lapply(c("old", "lure"), function(probe) {
    part <- predictions[predictions$probe_type == probe, , drop = FALSE]
    plan <- transport_exploratory_bootstrap_plan(part, draws, seed)
    weights <- transport_exploratory_align_weights(part, plan)
    gains <- list(
      group_over_base = (part$loss_base - part$loss_group) / log(2),
      own_over_group = (part$loss_group - part$loss_own_beyond_group) / log(2)
    )
    dplyr::bind_rows(lapply(names(gains), function(contrast) {
      value <- gains[[contrast]]
      bootstrap <- vapply(seq_len(ncol(weights)), function(draw) {
        transport_exploratory_weighted_mean(value, weights[, draw])
      }, numeric(1))
      interval <- transport_exploratory_interval(bootstrap, 0.95)
      data.frame(
        probe_type = probe, contrast = contrast, trials = nrow(part),
        information_gain_bits = mean(value),
        lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
        stringsAsFactors = FALSE
      )
    }))
  })
  list(predictions = predictions, summary = dplyr::bind_rows(rows))
}

own_group_behavior_panel <- function(own, context = own_group_context()) {
  metadata <- unique(context$raw$retrieval[c(
    "participant", "item", "probe_type", "degradation", "Response"
  )])
  metadata <- metadata[metadata$probe_type %in% c("old", "lure"), ]
  if (anyDuplicated(metadata[c("participant", "item")])) {
    stop("Behavior metadata are not unique by old/lure trial.")
  }
  quality <- own[c(
    "participant", "item", "outer_fold", "effective_fixations",
    "total_duration"
  )]
  panel <- merge(
    quality, metadata, by = c("participant", "item"),
    all.x = TRUE, sort = FALSE
  )
  if (nrow(panel) != nrow(own) || anyNA(panel$probe_type)) {
    stop("Behavior metadata do not align with frozen Transport support.")
  }
  panel$saliency_z <- (as.numeric(panel$degradation) - 60) / 20
  panel$oldness <- transport_density_oldness(panel$Response)
  panel$trial_key <- own_group_trial_key(panel$participant, panel$item)
  panel[order(panel$participant, panel$item), , drop = FALSE]
}

own_group_finalize <- function(
    output_dir = own_group_result_dir, draws = own_group_draws,
    seed = own_group_seed) {
  old_group <- own_group_read_scores(output_dir, "old_lure")
  newtest <- own_group_read_scores(output_dir, "newtest")
  own_paths <- vapply(1:4, function(fold) {
    v3_retrieval_checkpoint(v3_retrieval_output_dir, fold)
  }, character(1))
  own <- own_group_add_margin(dplyr::bind_rows(lapply(own_paths, function(path) {
    readRDS(path)$scored
  })))
  panel <- own_group_behavior_panel(own)
  summary <- own_group_measurement_summary(
    own, old_group$scored, newtest$scored, panel, draws, seed
  )
  group_key <- match(own$trial_key, old_group$scored$trial_key)
  behavior_tab <- panel[match(own$trial_key, panel$trial_key), ]
  behavior_tab$own_margin <- own$candidate_margin
  behavior_tab$group_margin <- old_group$scored$candidate_margin[group_key]
  behavior <- own_group_behavior(behavior_tab, draws, seed)
  audit <- rbind(old_group$audit, newtest$audit)
  configuration <- data.frame(
    protocol = "own-group/1.1.0", seed = seed, draws = draws,
    old_lure_trials = nrow(old_group$scored),
    newtest_trials = nrow(newtest$scored), episodes_per_candidate = 4L,
    exact_masks_used = FALSE, response_blind_measurement = TRUE,
    stringsAsFactors = FALSE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, file.path(output_dir, "measurement-summary.csv"),
                   row.names = FALSE)
  utils::write.csv(behavior$summary, file.path(output_dir, "behavior-summary.csv"),
                   row.names = FALSE)
  utils::write.csv(audit, file.path(output_dir, "fold-audit.csv"),
                   row.names = FALSE)
  utils::write.csv(configuration, file.path(output_dir, "configuration.csv"),
                   row.names = FALSE)
  saveRDS(
    list(
      measurement = summary, behavior = behavior$summary,
      audit = audit, configuration = configuration
    ),
    file.path(output_dir, "scientific-results.rds"), version = 3
  )
  saveRDS(
    list(
      own = own, group = old_group$scored, newtest = newtest$scored,
      behavior = behavior$predictions
    ),
    file.path(output_dir, "local-predictions.rds"), version = 3
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[basename(files) != "manifest-md5.csv"]
  utils::write.csv(
    data.frame(
      file = basename(files), md5 = unname(tools::md5sum(files)),
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "manifest-md5.csv"), row.names = FALSE
  )
  invisible(list(
    measurement = summary, behavior = behavior$summary,
    audit = audit, configuration = configuration
  ))
}
