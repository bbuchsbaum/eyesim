# Prospectively frozen participant- and item-disjoint GazeWeave replication.
#
# Protocol: inst/validation/GAZEWEAVE-REPLICATION-CONTROLS.md
# Run from a source checkout after devtools::load_all():
#   source("inst/validation/gaze-weave-real-replication.R")
#   run_gaze_weave_real_replication(
#     "inst/validation/gaze-weave-replication-results"
#   )

source("inst/validation/gaze-weave-real-controls.R")

replication_exclusions <- function() {
  list(
    participants = c("104", "121", "124", "128", "315", "317", "319", "327"),
    items = c(5L, 18L, 69L, 73L, 81L, 82L, 88L, 105L)
  )
}

real_select_replication_cohort <- function(
    raw, seed = 20260818L, n_participants_per_age = 4L, n_items = 8L) {
  excluded <- replication_exclusions()
  study <- raw$study
  recall <- raw$recall
  study$nfix <- vapply(study$fixgroup, nrow, integer(1))
  recall$nfix <- vapply(recall$fixgroup, nrow, integer(1))
  first <- study[
    study$repetition == 1L & study$nfix >= 3L,
    c("participant", "age", "item")
  ]
  fourth <- study[
    study$repetition == 4L & study$nfix >= 3L,
    c("participant", "age", "item")
  ]
  retrieval <- recall[
    recall$nfix >= 3L, c("participant", "age", "item")
  ]
  complete <- merge(
    merge(first, fourth, by = c("participant", "age", "item")),
    retrieval,
    by = c("participant", "age", "item")
  )
  counts <- aggregate(item ~ participant + age, complete, function(x) {
    length(unique(x))
  })
  eligible <- counts[
    counts$item >= 100L &
      !counts$participant %in% excluded$participants,
    , drop = FALSE
  ]
  if (length(unique(eligible$age)) < 2L) {
    stop("The age-stratified replication cohort cannot be formed.")
  }

  set.seed(seed)
  chosen_participants <- unlist(
    lapply(sort(unique(eligible$age)), function(group) {
      ids <- sort(eligible$participant[eligible$age == group])
      if (length(ids) < n_participants_per_age) {
        stop("Too few unused eligible participants in age group ", group, ".")
      }
      sample(ids, n_participants_per_age)
    }),
    use.names = FALSE
  )
  complete_selected <- complete[
    complete$participant %in% chosen_participants &
      !complete$item %in% excluded$items,
    , drop = FALSE
  ]
  item_counts <- table(complete_selected$item)
  eligible_items <- as.integer(names(item_counts)[
    item_counts == length(chosen_participants)
  ])
  if (length(eligible_items) < n_items) {
    stop("Too few unused items are complete for the replication participants.")
  }
  chosen_items <- sample(sort(eligible_items), n_items)
  list(
    participants = sort(as.character(chosen_participants)),
    items = sort(as.integer(chosen_items)),
    eligibility = list(
      participants_per_age = table(eligible$age),
      common_unused_item_n = length(eligible_items),
      excluded_participants = excluded$participants,
      excluded_items = excluded$items
    )
  )
}

replication_write_results <- function(result, output_dir) {
  real_write_results(result, output_dir)
  configuration <- data.frame(
    protocol_version = 2L,
    seed = result$config$seed,
    smoke = result$config$smoke,
    participants = paste(result$config$cohort$participants, collapse = ";"),
    items = paste(result$config$cohort$items, collapse = ";"),
    excluded_participants = paste(
      result$config$cohort$eligibility$excluded_participants,
      collapse = ";"
    ),
    excluded_items = paste(
      result$config$cohort$eligibility$excluded_items,
      collapse = ";"
    ),
    replay_score_scale = "mean_log_likelihood_per_duration_bin",
    advance = result$verdict$advance,
    provisional_default = result$verdict$provisional_default,
    strongest_baseline = result$verdict$strongest_baseline,
    superiority_supported = result$verdict$superiority_supported,
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    configuration,
    file.path(output_dir, "configuration.csv"),
    row.names = FALSE
  )
  saveRDS(
    result[c(
      "method_summary", "frozen_summary", "bootstrap", "operating",
      "order_checks", "diagnostics", "warps", "fold_audit", "covariates",
      "comparison", "verdict", "resources", "sensitivity", "config"
    )],
    file.path(output_dir, "scientific-results.rds"),
    version = 3
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[basename(files) != "manifest-md5.csv"]
  manifest <- data.frame(
    file = basename(files),
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    manifest,
    file.path(output_dir, "manifest-md5.csv"),
    row.names = FALSE
  )
  invisible(result)
}

run_gaze_weave_real_replication <- function(
    output_dir = NULL, seed = 20260818L, smoke = FALSE,
    bootstrap_draws = if (smoke) 100L else 1000L) {
  if (!smoke && seed != 20260818L) {
    stop("The final replication is frozen to seed 20260818; use smoke = TRUE for debugging.")
  }
  raw <- real_read_wynn(c(0, 3000))
  cohort <- real_select_replication_cohort(raw, seed)
  expected_participants <- c("111", "130", "18", "300", "302", "316", "325", "9")
  expected_items <- c(4L, 11L, 24L, 62L, 63L, 80L, 92L, 101L)
  if (!smoke && (!identical(cohort$participants, expected_participants) ||
                 !identical(cohort$items, expected_items))) {
    stop("The data-only replication cohort no longer matches the frozen protocol.")
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
      message("Real replication: ", task, ", outer fold ", fold$id, "/4")
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
  full_task <- real_task_tables(
    raw_full, cohort, "blank_screen_imagery", seed
  )
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
    scored = scored,
    method_summary = method_summary,
    frozen_summary = frozen_summary,
    bootstrap = bootstrap,
    operating = operating,
    order_checks = order_checks,
    diagnostics = diagnostics,
    warps = warps,
    fold_audit = fold_audit,
    covariates = covariates,
    comparison = comparison,
    verdict = verdict,
    resources = resources,
    sensitivity = sensitivity,
    config = list(
      protocol_version = 2L,
      seed = seed,
      smoke = smoke,
      bootstrap_draws = bootstrap_draws,
      cohort = cohort,
      fold_plan = plan,
      workers = workers,
      replay_score_scale = "mean_log_likelihood_per_duration_bin"
    ),
    audits = audits
  )
  if (!is.null(output_dir)) replication_write_results(result, output_dir)
  result
}
