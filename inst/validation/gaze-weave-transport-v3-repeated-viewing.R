# Frozen response-blind repeated-viewing court for GazeWeave Transport v3.
#
# This script deliberately reads study fixations only. Trial-linked checkpoints
# are written below the Git-ignored results directory; only aggregate summaries
# may be copied into a tracked report.

v3_repeated_protocol_version <- 1L
v3_repeated_seed <- 20260823L
v3_repeated_shuffle_seed <- 20260824L
v3_repeated_bootstrap_draws <- 2000L
v3_repeated_output_dir <- file.path(
  "inst", "validation",
  "gaze-weave-transport-v3-repeated-viewing-results"
)
v3_repeated_contrasts <- data.frame(
  task = c(
    "study_p1_p2", "study_p1_p3", "study_p1_p4",
    "study_p2_p3", "study_p3_p4"
  ),
  reference_presentation = c(1L, 1L, 1L, 2L, 3L),
  source_presentation = c(2L, 3L, 4L, 3L, 4L),
  stringsAsFactors = FALSE
)

v3_repeated_followup_file <- system.file(
  "validation", "gaze-weave-followup-full-cohort.R", package = "eyesim"
)
if (!nzchar(v3_repeated_followup_file)) {
  v3_repeated_followup_file <- file.path(
    "inst", "validation", "gaze-weave-followup-full-cohort.R"
  )
}
if (!file.exists(v3_repeated_followup_file)) {
  stop("Run the repeated-viewing court from the eyesim source root.")
}
source(v3_repeated_followup_file, local = TRUE)

v3_repeated_spec <- function(smoke = FALSE) {
  inherited <- followup_specs("study_p1_p4", smoke)
  gaze_transport_spec(
    spatial = inherited$transport$spatial,
    chronology = inherited$transport$chronology,
    coverage_nodes = if (smoke) 2L else 12L,
    temporal_weight = 2,
    entropy_schedule = if (smoke) 0.03 else c(0.05, 0.015),
    warp = inherited$warp,
    screen = inherited$screen,
    maxit = 1000L,
    tolerance = 5e-5,
    projection_maxit = 1000L,
    projection_tolerance = 1e-8,
    backend = "auto",
    reliability = "effective_fixations",
    calibration_folds = 2L,
    calibration_seed = 20260822L
  )
}

v3_repeated_task_tables <- function(raw, cohort, task) {
  row <- v3_repeated_contrasts[
    v3_repeated_contrasts$task == task, , drop = FALSE
  ]
  if (nrow(row) != 1L) stop("Unknown repeated-viewing task: ", task)
  followup_study_tables(
    raw, cohort, row$reference_presentation, row$source_presentation
  )
}

v3_repeated_path_order <- function(path, order) {
  order <- as.integer(order)
  if (!inherits(path, "fixation_group") ||
      length(order) != nrow(path) || !setequal(order, seq_len(nrow(path)))) {
    stop("A path perturbation requires one complete fixation permutation.")
  }
  duration <- path$duration[order]
  fixation_group(
    x = path$x[order], y = path$y[order], duration = duration,
    onset = c(0, head(cumsum(duration), -1L))
  )
}

v3_repeated_perturb_path <- function(path, control, key = "row",
                                     seed = v3_repeated_shuffle_seed) {
  control <- match.arg(control, c("intact", "reversed", "shuffled"))
  if (identical(control, "intact")) return(path)
  order <- rev(seq_len(nrow(path)))
  if (identical(control, "shuffled")) {
    key_code <- sum(utf8ToInt(as.character(key))) %% 100000L
    set.seed(as.integer(seed + key_code))
    order <- sample(seq_len(nrow(path)), nrow(path), replace = FALSE)
    if (nrow(path) > 2L && identical(order, seq_len(nrow(path)))) {
      order <- c(order[-1L], order[[1L]])
    }
  }
  v3_repeated_path_order(path, order)
}

v3_repeated_expand_pool <- function(reference, source, candidate_plan) {
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
  pool <- pool[order(pool$candidate_set_id, pool$candidate_position), ]
  source <- source[order(source$candidate_set_id), ]
  if (nrow(pool) != 5L * nrow(source) || anyNA(pool$fixgroup) ||
      any(table(pool$candidate_set_id) != 5L) ||
      any(tapply(pool$is_true, pool$candidate_set_id, sum) != 1L)) {
    stop("Repeated-viewing candidate expansion violates the frozen K=5 plan.")
  }
  list(reference = tibble::as_tibble(pool), source = tibble::as_tibble(source))
}

v3_repeated_fold_data <- function(tables, candidate_plan, fold) {
  subset <- followup_fold_subsets(tables, candidate_plan, fold)
  train <- v3_repeated_expand_pool(
    tables$reference, subset$source_train, candidate_plan
  )
  eval <- v3_repeated_expand_pool(
    tables$reference, subset$source_eval, candidate_plan
  )
  list(
    ref_true_train = subset$ref_train,
    source_train = train$source,
    ref_train_pool = train$reference,
    source_eval = eval$source,
    ref_eval_pool = eval$reference
  )
}

v3_repeated_score_row <- function(source_row, reference_pool, spec, warp,
                                  calibration, control = "intact") {
  key <- paste(source_row$participant, source_row$item, sep = ":")
  source_row$fixgroup[[1L]] <- v3_repeated_perturb_path(
    source_row$fixgroup[[1L]], control, key
  )
  scored <- eyesim:::score_transport_v3_cv_row(
    source_row, reference_pool,
    c("participant", "item"), "candidate_set_id",
    "fixgroup", "fixgroup", NULL, "prior_weight", spec, warp,
    temperature = calibration$temperature,
    kappa = calibration$kappa
  )
  evidence <- scored$evidence
  candidate_converged <- vapply(scored$candidates, function(candidate) {
    isTRUE(candidate$convergence$converged)
  }, logical(1))
  data.frame(
    task = source_row$task[[1L]], control = control,
    participant = source_row$participant[[1L]], item = source_row$item[[1L]],
    gaze_info_bits = evidence$gaze_info_bits,
    log_loss = evidence$log_loss,
    brier_score = evidence$brier_score,
    posterior_true = evidence$posterior_true,
    prior_true = evidence$prior_true,
    template_rank = evidence$template_rank,
    top1_credit = evidence$top1_credit,
    candidate_count = evidence$candidate_count,
    temperature = evidence$temperature,
    reliability = evidence$reliability,
    effective_fixations = scored$quality$effective_fixations,
    converged = scored$all_converged,
    converged_candidates = sum(candidate_converged),
    candidate_alignments = length(candidate_converged),
    common_episode_count = length(scored$common_episode_ids),
    candidates = I(list(evidence$candidates)),
    stringsAsFactors = FALSE
  )
}

v3_repeated_score_fold <- function(tables, candidate_plan, fold, spec,
                                   controls = "intact") {
  total_started <- proc.time()[["elapsed"]]
  data <- v3_repeated_fold_data(tables, candidate_plan, fold)
  fit_started <- proc.time()[["elapsed"]]
  warp <- eyesim:::fit_transport_v3_warp(
    data$ref_true_train, data$source_train,
    c("participant", "item"), "fixgroup", "fixgroup", NULL, spec
  )
  data$source_train$..gaze_row_id <- seq_len(nrow(data$source_train))
  inner <- eyesim:::fit_transport_v3_inner_calibration(
    data$ref_train_pool, data$source_train,
    c("participant", "item"), "candidate_set_id",
    "fixgroup", "fixgroup", NULL, "prior_weight", spec,
    fold_contrast_on = NULL
  )
  fit_elapsed <- proc.time()[["elapsed"]] - fit_started
  scoring_started <- proc.time()[["elapsed"]]
  rows <- unlist(lapply(controls, function(control) {
    lapply(seq_len(nrow(data$source_eval)), function(index) {
      v3_repeated_score_row(
        data$source_eval[index, , drop = FALSE], data$ref_eval_pool,
        spec, warp, inner$calibration, control
      )
    })
  }), recursive = FALSE)
  scoring_elapsed <- proc.time()[["elapsed"]] - scoring_started
  total_elapsed <- proc.time()[["elapsed"]] - total_started
  scored <- dplyr::bind_rows(rows)
  scored$outer_fold <- fold$id
  list(
    scored = scored,
    resources = data.frame(
      task = unique(scored$task), outer_fold = fold$id,
      controls = paste(controls, collapse = ","),
      eval_trials = nrow(data$source_eval),
      scored_rows = nrow(scored), elapsed_seconds = total_elapsed,
      fit_calibration_seconds = fit_elapsed,
      scoring_seconds = scoring_elapsed,
      seconds_per_score = scoring_elapsed / nrow(scored),
      stringsAsFactors = FALSE
    ),
    audit = data.frame(
      task = unique(scored$task), outer_fold = fold$id,
      train_trials = nrow(data$source_train),
      eval_trials = nrow(data$source_eval),
      participant_overlap = length(intersect(
        unique(data$source_train$participant),
        unique(data$source_eval$participant)
      )),
      item_overlap = length(intersect(
        unique(data$source_train$item), unique(data$source_eval$item)
      )),
      candidate_count_min = min(scored$candidate_count),
      candidate_count_max = max(scored$candidate_count),
      row_convergence = mean(scored$converged),
      alignment_convergence = sum(scored$converged_candidates) /
        sum(scored$candidate_alignments),
      stringsAsFactors = FALSE
    ),
    warp = warp$info,
    calibration = inner$calibration,
    protocol_version = v3_repeated_protocol_version,
    seed = v3_repeated_seed,
    shuffle_seed = v3_repeated_shuffle_seed
  )
}

v3_repeated_checkpoint <- function(output_dir, task, fold_id) {
  file.path(output_dir, sprintf("checkpoint-%s-fold-%02d.rds", task, fold_id))
}

run_gaze_weave_transport_v3_repeated_viewing <- function(
    output_dir = v3_repeated_output_dir, tasks = v3_repeated_contrasts$task,
    fold_ids = NULL, smoke = FALSE, resume = TRUE) {
  tasks <- match.arg(tasks, v3_repeated_contrasts$task, several.ok = TRUE)
  raw <- full_recognition_read_inputs(verify = !smoke)
  cohort <- full_recognition_select_cohort(raw)
  if (smoke) cohort <- full_recognition_smoke_cohort(cohort)
  candidate_plan <- full_recognition_candidate_plan(cohort, followup_seed)
  fold_plan <- full_recognition_fold_plan(cohort, followup_seed)
  folds <- fold_plan$folds
  if (!is.null(fold_ids)) {
    folds <- folds[vapply(folds, function(x) x$id %in% fold_ids, logical(1))]
  }
  if (length(folds) == 0L) stop("fold_ids selected no outer folds.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  spec <- v3_repeated_spec(smoke)
  for (task in tasks) {
    tables <- v3_repeated_task_tables(raw, cohort, task)
    plan <- followup_available_candidate_plan(candidate_plan, tables$signal)
    controls <- if (identical(task, "study_p1_p4")) {
      c("intact", "reversed", "shuffled")
    } else {
      "intact"
    }
    for (fold in folds) {
      path <- v3_repeated_checkpoint(output_dir, task, fold$id)
      if (resume && file.exists(path)) {
        old <- readRDS(path)
        if (!identical(old$protocol_version, v3_repeated_protocol_version) ||
            !identical(old$seed, v3_repeated_seed) ||
            !identical(old$shuffle_seed, v3_repeated_shuffle_seed)) {
          stop("A checkpoint belongs to another repeated-viewing protocol: ", path)
        }
        message("Using checkpoint: ", basename(path))
        next
      }
      message("Transport v3 repeated viewing: ", task, ", fold ", fold$id)
      value <- v3_repeated_score_fold(tables, plan, fold, spec, controls)
      saveRDS(value, path, version = 3)
    }
  }
  invisible(list(cohort = cohort, candidate_plan = candidate_plan,
                 fold_plan = fold_plan, tasks = tasks))
}

v3_repeated_read <- function(
    output_dir = v3_repeated_output_dir,
    tasks = v3_repeated_contrasts$task, fold_ids = 1:4,
    strict = TRUE) {
  paths <- unlist(lapply(tasks, function(task) {
    vapply(fold_ids, function(fold) {
      v3_repeated_checkpoint(output_dir, task, fold)
    }, character(1))
  }), use.names = FALSE)
  if (strict && !all(file.exists(paths))) {
    stop("One or more repeated-viewing checkpoints are missing.")
  }
  paths <- paths[file.exists(paths)]
  values <- lapply(paths, readRDS)
  list(
    scored = dplyr::bind_rows(lapply(values, `[[`, "scored")),
    resources = dplyr::bind_rows(lapply(values, `[[`, "resources")),
    audit = dplyr::bind_rows(lapply(values, `[[`, "audit")),
    checkpoints = values, paths = paths
  )
}

v3_repeated_interval <- function(value, weight) {
  draws <- vapply(seq_len(ncol(weight)), function(index) {
    followup_weighted_mean(value, weight[, index])
  }, numeric(1))
  recognition_interval(draws, 0.95)
}

v3_repeated_summary <- function(scored, cohort,
                                draws = v3_repeated_bootstrap_draws,
                                seed = v3_repeated_seed) {
  intact <- scored[scored$control == "intact", , drop = FALSE]
  plan <- followup_bootstrap_plan(cohort, draws, seed)
  groups <- split(seq_len(nrow(intact)), intact$task)
  dplyr::bind_rows(lapply(groups, function(index) {
    tab <- intact[index, , drop = FALSE]
    weights <- followup_score_weights(tab, plan)
    interval <- v3_repeated_interval(tab$gaze_info_bits, weights)
    data.frame(
      task = tab$task[[1L]], trials = nrow(tab),
      mean_info_bits = mean(tab$gaze_info_bits),
      lower_95 = interval[["lower"]], upper_95 = interval[["upper"]],
      mean_log_loss = mean(tab$log_loss),
      mean_rank = mean(tab$template_rank),
      top1_credit = mean(tab$top1_credit),
      row_convergence = mean(tab$converged),
      alignment_convergence = sum(tab$converged_candidates) /
        sum(tab$candidate_alignments),
      stringsAsFactors = FALSE
    )
  }))
}

v3_repeated_advance <- function(scored, cohort,
                                draws = v3_repeated_bootstrap_draws,
                                seed = v3_repeated_seed) {
  tab <- scored[scored$task == "study_p1_p4", , drop = FALSE]
  wide <- reshape(
    tab[c("participant", "item", "control", "gaze_info_bits")],
    idvar = c("participant", "item"), timevar = "control", direction = "wide"
  )
  required <- paste0("gaze_info_bits.", c("intact", "reversed", "shuffled"))
  if (!all(required %in% names(wide)) || anyNA(wide[required])) {
    stop("The P1-to-P4 intact/control panel is incomplete.")
  }
  wide$intact_control <- wide$gaze_info_bits.intact -
    (wide$gaze_info_bits.reversed + wide$gaze_info_bits.shuffled) / 2
  plan <- followup_bootstrap_plan(cohort, draws, seed)
  weights <- followup_score_weights(wide, plan)
  intact_ci <- v3_repeated_interval(wide$gaze_info_bits.intact, weights)
  control_ci <- v3_repeated_interval(wide$intact_control, weights)
  result <- data.frame(
    endpoint = c("intact_information", "intact_minus_mean_control"),
    estimate = c(mean(wide$gaze_info_bits.intact), mean(wide$intact_control)),
    lower_95 = c(intact_ci[["lower"]], control_ci[["lower"]]),
    upper_95 = c(intact_ci[["upper"]], control_ci[["upper"]]),
    stringsAsFactors = FALSE
  )
  result$positive <- result$estimate > 0
  result$lower_above_zero <- result$lower_95 > 0
  attr(result, "advance") <- all(result$positive) && any(result$lower_above_zero)
  result
}

v3_repeated_calibration <- function(scored) {
  intact <- scored[scored$control == "intact", , drop = FALSE]
  groups <- split(seq_len(nrow(intact)), intact$task)
  dplyr::bind_rows(lapply(groups, function(index) {
    tab <- intact[index, , drop = FALSE]
    evidence <- lapply(seq_len(nrow(tab)), function(i) {
      structure(
        list(candidates = tab$candidates[[i]]),
        class = c("gaze_candidate_evidence", "list")
      )
    })
    curve <- eyesim:::gaze_reliability_table(evidence)
    data.frame(
      task = tab$task[[1L]],
      heldout_log_loss = mean(tab$log_loss),
      heldout_ece = eyesim:::transport_v3_expected_calibration_error(curve),
      stringsAsFactors = FALSE
    )
  }))
}

v3_repeated_comparator_summary <- function(
    v3_scored, cohort, draws = v3_repeated_bootstrap_draws,
    seed = v3_repeated_seed,
    comparator_dir = file.path(
      "inst", "validation", "gaze-weave-followup-full-cohort-results"
    )) {
  prior <- followup_read_checkpoints(
    output_dir = comparator_dir, tasks = "study_p1_p4",
    method_groups = c("replay", "baselines", "transport_v2")
  )$scored
  prior <- prior[
    prior$status == "scored" &
      (prior$calibrated | prior$method %in% c("transport_v2", "replay")),
    , drop = FALSE
  ]
  v3 <- v3_scored[
    v3_scored$task == "study_p1_p4" & v3_scored$control == "intact",
    , drop = FALSE
  ]
  v3$method <- "transport_v3"
  v3$calibrated <- TRUE
  v3$status <- "scored"
  combined <- dplyr::bind_rows(prior, v3)
  followup_study_summary(combined, cohort, draws, seed)
}
