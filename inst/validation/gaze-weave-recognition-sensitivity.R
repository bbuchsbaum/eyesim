# Secondary saliency/correctness sensitivity court for frozen GW-11 scores.
#
# Protocol: inst/validation/GAZEWEAVE-RECOGNITION-SENSITIVITY.md
# Run from the eyesim source root after devtools::load_all().

recognition_protocol_version <- 1L
recognition_protocol_seed <- 20260824L
recognition_bootstrap_draws <- 2000L
recognition_family_size <- 6L
recognition_methods <- c(
  "transport_v2", "replay", "density_ridge_registered",
  "multimatch_ridge_registered", "elastic_ridge_registered"
)
recognition_engines <- c("transport_v2", "replay")
recognition_phases <- c("combined", "probe", "delay")
recognition_contrast_names <- c(
  "saliency_20_to_100", "correct_at_60", "interaction_per_20"
)

recognition_result_dir <- function() {
  file.path("inst", "validation", "gaze-weave-probe-delay-results")
}

recognition_verify_source <- function(
    result_dir = recognition_result_dir()) {
  manifest_path <- file.path(result_dir, "manifest-md5.csv")
  if (!file.exists(manifest_path)) {
    stop("The local-only GW-11 result manifest is missing.")
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  if (!all(c("file", "md5") %in% names(manifest)) || nrow(manifest) < 1L) {
    stop("The GW-11 result manifest has an unexpected schema.")
  }
  paths <- file.path(result_dir, manifest$file)
  if (!all(file.exists(paths))) stop("A manifested GW-11 result is missing.")
  observed <- unname(tools::md5sum(paths))
  if (!identical(observed, manifest$md5)) {
    stop("The GW-11 result manifest integrity check failed.")
  }
  manifest$path <- paths
  manifest
}

recognition_methods_for_phase <- function(phase) {
  if (identical(phase, "probe")) {
    setdiff(recognition_methods, "multimatch_ridge_registered")
  } else {
    recognition_methods
  }
}

recognition_prepare_rows <- function(tab) {
  tab$participant <- as.character(tab$participant)
  tab$item <- as.integer(tab$item)
  tab$saliency <- as.numeric(tab$degradation)
  tab$accuracy <- as.integer(tab$accuracy)
  tab$saliency_z <- (tab$saliency - 60) / 20
  tab$correct_ec <- tab$accuracy - 0.5
  tab$probe_type_ec <- ifelse(tab$probe_type == "old", 0.5, -0.5)
  tab$trial_key <- paste(tab$participant, tab$item, sep = ":")
  tab[order(tab$participant, tab$item, tab$method), , drop = FALSE]
}

recognition_read_scores <- function(
    result_dir = recognition_result_dir(), verify = TRUE) {
  if (verify) recognition_verify_source(result_dir)
  path <- file.path(result_dir, "trial-scores.csv")
  scores <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c(
    "task", "condition", "participant", "item", "method", "calibrated",
    "status", "gaze_info_bits", "candidate_count", "probe_type",
    "degradation", "accuracy"
  )
  missing <- setdiff(required, names(scores))
  if (length(missing) > 0L) {
    stop("GW-11 scores are missing columns: ", paste(missing, collapse = ", "))
  }
  keep <- scores$task %in% recognition_phases &
    scores$condition == "signal" & scores$calibrated &
    scores$status == "scored" & scores$method %in% recognition_methods
  scores <- recognition_prepare_rows(scores[keep, , drop = FALSE])
  if (any(!is.finite(scores$gaze_info_bits)) ||
      any(!scores$accuracy %in% 0:1) ||
      any(!scores$saliency %in% c(20, 40, 60, 80, 100)) ||
      any(!scores$probe_type %in% c("old", "lure")) ||
      any(scores$candidate_count != 5L)) {
    stop("GW-11 recognition scores violate the frozen value contract.")
  }
  for (phase in recognition_phases) {
    phase_methods <- recognition_methods_for_phase(phase)
    expected_keys <- NULL
    for (method in phase_methods) {
      part <- scores[scores$task == phase & scores$method == method, ]
      if (nrow(part) != 80L || anyDuplicated(part$trial_key)) {
        stop("Unexpected trial support for ", phase, "/", method, ".")
      }
      keys <- sort(part$trial_key)
      if (is.null(expected_keys)) expected_keys <- keys
      if (!identical(keys, expected_keys)) {
        stop("Methods do not share trial support within phase ", phase, ".")
      }
    }
  }
  reference <- scores[
    scores$task == "combined" & scores$method == "transport_v2", ]
  if (length(unique(reference$participant)) != 8L ||
      length(unique(reference$item)) != 61L ||
      !identical(as.integer(table(reference$accuracy)), c(21L, 59L))) {
    stop("Recognition support differs from the frozen protocol.")
  }
  scores
}

recognition_design_matrix <- function(tab) {
  stats::model.matrix(
    ~ saliency_z * correct_ec + probe_type_ec,
    data = tab
  )
}

recognition_fit_fixed <- function(tab, weights = NULL) {
  design <- recognition_design_matrix(tab)
  response <- as.numeric(tab$gaze_info_bits)
  if (is.null(weights)) weights <- rep(1, nrow(tab))
  keep <- is.finite(weights) & weights > 0
  if (sum(keep) < ncol(design)) {
    return(list(status = "insufficient", coefficients = NULL, rank = 0L))
  }
  fit <- stats::lm.wfit(design[keep, , drop = FALSE], response[keep], weights[keep])
  coefficients <- stats::setNames(as.numeric(fit$coefficients), colnames(design))
  valid <- fit$rank == ncol(design) && all(is.finite(coefficients))
  list(
    status = if (valid) "scored" else "rank_deficient",
    coefficients = if (valid) coefficients else NULL,
    rank = fit$rank
  )
}

recognition_contrasts <- function(coefficients) {
  required <- c("saliency_z", "correct_ec", "saliency_z:correct_ec")
  if (is.null(coefficients) || !all(required %in% names(coefficients))) {
    return(stats::setNames(rep(NA_real_, 3L), recognition_contrast_names))
  }
  c(
    saliency_20_to_100 = 4 * coefficients[["saliency_z"]],
    correct_at_60 = coefficients[["correct_ec"]],
    interaction_per_20 = coefficients[["saliency_z:correct_ec"]]
  )
}

recognition_bootstrap_plan <- function(tab, draws = recognition_bootstrap_draws,
                                       seed = recognition_protocol_seed) {
  participants <- sort(unique(as.character(tab$participant)))
  items <- sort(unique(as.integer(tab$item)))
  set.seed(seed)
  weights <- vapply(seq_len(draws), function(draw) {
    participant_frequency <- table(sample(
      participants, length(participants), replace = TRUE
    ))
    item_frequency <- table(sample(items, length(items), replace = TRUE))
    participant_weight <- as.numeric(
      participant_frequency[as.character(tab$participant)]
    )
    item_weight <- as.numeric(item_frequency[as.character(tab$item)])
    participant_weight[is.na(participant_weight)] <- 0
    item_weight[is.na(item_weight)] <- 0
    participant_weight * item_weight
  }, numeric(nrow(tab)))
  if (draws == 1L) weights <- matrix(weights, ncol = 1L)
  list(
    weights = weights, participants = participants, items = items,
    draws = draws, seed = seed, trial_keys = tab$trial_key
  )
}

recognition_align_plan <- function(tab, plan) {
  position <- match(tab$trial_key, plan$trial_keys)
  if (anyNA(position) || anyDuplicated(position) || length(position) != nrow(tab)) {
    stop("The bootstrap plan does not match this phase/method trial support.")
  }
  plan$weights[position, , drop = FALSE]
}

recognition_bootstrap_contrasts <- function(tab, plan) {
  weights <- recognition_align_plan(tab, plan)
  values <- vapply(seq_len(ncol(weights)), function(draw) {
    recognition_contrasts(
      recognition_fit_fixed(tab, weights[, draw])$coefficients
    )
  }, numeric(length(recognition_contrast_names)))
  if (ncol(weights) == 1L) values <- matrix(values, ncol = 1L)
  rownames(values) <- recognition_contrast_names
  values
}

recognition_lopo <- function(tab) {
  participants <- sort(unique(tab$participant))
  values <- vapply(participants, function(participant) {
    recognition_contrasts(recognition_fit_fixed(
      tab[tab$participant != participant, , drop = FALSE]
    )$coefficients)
  }, numeric(length(recognition_contrast_names)))
  if (length(participants) == 1L) values <- matrix(values, ncol = 1L)
  rownames(values) <- recognition_contrast_names
  colnames(values) <- participants
  values
}

recognition_fit_lmer <- function(tab) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    return(list(status = "missing_lme4", fit = NULL))
  }
  warnings <- character()
  messages <- character()
  fit <- tryCatch(
    withCallingHandlers(
      lme4::lmer(
        gaze_info_bits ~ saliency_z * correct_ec + probe_type_ec +
          (1 | participant) + (1 | item),
        data = tab, REML = FALSE,
        control = lme4::lmerControl(
          optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)
        )
      ),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      },
      message = function(message) {
        messages <<- c(messages, conditionMessage(message))
        invokeRestart("muffleMessage")
      }
    ),
    error = function(error) error
  )
  if (inherits(fit, "error")) {
    return(list(
      status = "error", fit = NULL, error = conditionMessage(fit),
      warnings = warnings, messages = messages
    ))
  }
  convergence_messages <- fit@optinfo$conv$lme4$messages
  if (is.null(convergence_messages)) convergence_messages <- character()
  variances <- as.data.frame(lme4::VarCorr(fit))
  optimizer_code <- fit@optinfo$conv$opt
  list(
    status = "scored", fit = fit,
    coefficients = lme4::fixef(fit),
    contrasts = recognition_contrasts(lme4::fixef(fit)),
    singular = lme4::isSingular(fit, tol = 1e-4),
    converged = is.numeric(optimizer_code) && length(optimizer_code) == 1L &&
      optimizer_code == 0L,
    optimizer_code = optimizer_code,
    convergence_messages = convergence_messages,
    warnings = warnings, messages = messages,
    participant_variance = variances$vcov[variances$grp == "participant"][[1L]],
    item_variance = variances$vcov[variances$grp == "item"][[1L]],
    residual_variance = variances$vcov[variances$grp == "Residual"][[1L]]
  )
}

recognition_interval <- function(values, probability) {
  finite <- values[is.finite(values)]
  if (length(finite) == 0L) return(c(lower = NA_real_, upper = NA_real_))
  alpha <- 1 - probability
  stats::setNames(
    unname(stats::quantile(finite, c(alpha / 2, 1 - alpha / 2), type = 8)),
    c("lower", "upper")
  )
}

recognition_fit_method_phase <- function(tab, plan) {
  fixed <- recognition_fit_fixed(tab)
  if (!identical(fixed$status, "scored")) {
    stop("The full recognition design is rank deficient.")
  }
  estimate <- recognition_contrasts(fixed$coefficients)
  bootstrap <- recognition_bootstrap_contrasts(tab, plan)
  lopo <- recognition_lopo(tab)
  lmer <- recognition_fit_lmer(tab)
  kappa_value <- kappa(recognition_design_matrix(tab), exact = TRUE)
  family_probability <- 1 - 0.05 / recognition_family_size
  effects <- dplyr::bind_rows(lapply(recognition_contrast_names, function(name) {
    values <- bootstrap[name, ]
    ordinary <- recognition_interval(values, 0.95)
    family <- recognition_interval(values, family_probability)
    lmer_value <- if (identical(lmer$status, "scored")) {
      lmer$contrasts[[name]]
    } else {
      NA_real_
    }
    data.frame(
      contrast = name,
      estimate = estimate[[name]],
      lower_95 = ordinary[["lower"]], upper_95 = ordinary[["upper"]],
      lower_family = family[["lower"]], upper_family = family[["upper"]],
      valid_draws = sum(is.finite(values)),
      valid_fraction = mean(is.finite(values)),
      design_kappa = kappa_value,
      lmer_estimate = lmer_value,
      lmer_sign_agrees = is.finite(lmer_value) &&
        sign(lmer_value) == sign(estimate[[name]]),
      lopo_sign_n = sum(sign(lopo[name, ]) == sign(estimate[[name]]), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  audit <- data.frame(
    fixed_rank = fixed$rank,
    fixed_columns = ncol(recognition_design_matrix(tab)),
    design_kappa = kappa_value,
    lmer_status = lmer$status,
    lmer_singular = if (identical(lmer$status, "scored")) lmer$singular else NA,
    lmer_converged = if (identical(lmer$status, "scored")) lmer$converged else FALSE,
    participant_variance = if (identical(lmer$status, "scored")) {
      lmer$participant_variance
    } else NA_real_,
    item_variance = if (identical(lmer$status, "scored")) {
      lmer$item_variance
    } else NA_real_,
    residual_variance = if (identical(lmer$status, "scored")) {
      lmer$residual_variance
    } else NA_real_,
    lmer_optimizer_code = if (identical(lmer$status, "scored")) {
      lmer$optimizer_code
    } else NA_integer_,
    lmer_warnings = paste(
      if (is.null(lmer$warnings)) character() else lmer$warnings,
      collapse = " | "
    ),
    lmer_messages = paste(
      if (is.null(lmer$messages)) character() else lmer$messages,
      collapse = " | "
    ),
    lmer_convergence_messages = paste(
      if (is.null(lmer$convergence_messages)) {
        character()
      } else {
        lmer$convergence_messages
      },
      collapse = " | "
    ),
    stringsAsFactors = FALSE
  )
  list(
    fixed = fixed, effects = effects, bootstrap = bootstrap,
    lopo = lopo, lmer = lmer, audit = audit
  )
}

recognition_cell_predictions <- function(coefficients, method, phase) {
  grid <- expand.grid(
    saliency = c(20, 40, 60, 80, 100),
    accuracy = 0:1,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$saliency_z <- (grid$saliency - 60) / 20
  grid$correct_ec <- grid$accuracy - 0.5
  grid$probe_type_ec <- 0
  design <- recognition_design_matrix(grid)
  grid$predicted_info_bits <- as.numeric(design %*% coefficients[colnames(design)])
  grid$method <- method
  grid$phase <- phase
  grid[c("phase", "method", "saliency", "accuracy", "predicted_info_bits")]
}

recognition_support_table <- function(scores) {
  reference <- scores[
    scores$task == "combined" & scores$method == "transport_v2", ]
  counts <- as.data.frame(table(
    saliency = reference$saliency,
    accuracy = reference$accuracy,
    probe_type = reference$probe_type
  ), stringsAsFactors = FALSE)
  counts[counts$Freq > 0L, , drop = FALSE]
}

recognition_verdict <- function(effects, audit) {
  primary <- effects[
    effects$phase == "combined" & effects$method %in% recognition_engines,
    , drop = FALSE
  ]
  audit_primary <- audit[
    audit$phase == "combined" & audit$method %in% recognition_engines,
    , drop = FALSE
  ]
  integrity <- all(primary$valid_fraction >= 0.95) &&
    all(primary$design_kappa <= 30) &&
    all(audit_primary$lmer_status == "scored") &&
    all(audit_primary$lmer_converged)
  primary$ordinary_excludes_zero <- primary$lower_95 > 0 | primary$upper_95 < 0
  primary$family_excludes_zero <- primary$lower_family > 0 |
    primary$upper_family < 0
  primary$stable <- primary$lopo_sign_n >= 7L & primary$lmer_sign_agrees
  supported <- primary$family_excludes_zero & primary$stable & integrity
  suggestive <- primary$ordinary_excludes_zero
  label <- if (!integrity) {
    "not_supported"
  } else if (any(supported)) {
    "supported"
  } else if (any(suggestive)) {
    "suggestive"
  } else {
    "not_supported"
  }
  format_effect <- function(index) {
    paste(primary$method[index], primary$contrast[index], sep = ":")
  }
  data.frame(
    verdict = label,
    supported_effects = paste(format_effect(which(supported)), collapse = ";"),
    suggestive_effects = paste(format_effect(which(suggestive)), collapse = ";"),
    design_integrity = integrity,
    post_primary_secondary = TRUE,
    rescues_persistence_gate = FALSE,
    independent_replication = FALSE,
    closes_wang_gate = FALSE,
    stringsAsFactors = FALSE
  )
}

recognition_configuration <- function(result) {
  data.frame(
    protocol_version = recognition_protocol_version,
    seed = result$config$seed,
    bootstrap_draws = result$config$draws,
    family_size = recognition_family_size,
    participants = paste(result$config$plan$participants, collapse = ";"),
    item_n = length(result$config$plan$items),
    trial_n = length(result$config$plan$trial_keys),
    methods = paste(recognition_methods, collapse = ";"),
    phases = paste(recognition_phases, collapse = ";"),
    verdict = result$verdict$verdict,
    local_only = TRUE,
    stringsAsFactors = FALSE
  )
}

recognition_write_results <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  write(recognition_configuration(result), "configuration.csv")
  write(result$verdict, "sensitivity-verdict.csv")
  write(result$effects, "conditional-effects.csv")
  write(result$audit, "model-audit.csv")
  write(result$lopo, "leave-one-participant-out.csv")
  write(result$predictions, "cell-predictions.csv")
  write(result$support, "support-counts.csv")
  saveRDS(
    result[c("verdict", "effects", "audit", "lopo", "predictions", "support", "config")],
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

run_gaze_weave_recognition_sensitivity <- function(
    output_dir = NULL,
    seed = recognition_protocol_seed,
    draws = recognition_bootstrap_draws,
    smoke = FALSE) {
  if (!smoke && (seed != recognition_protocol_seed ||
                 draws != recognition_bootstrap_draws)) {
    stop("Non-smoke recognition sensitivity runs use the frozen seed and draws.")
  }
  scores <- recognition_read_scores(verify = !smoke)
  reference <- scores[
    scores$task == "combined" & scores$method == "transport_v2", ]
  reference <- reference[order(reference$participant, reference$item), ]
  plan <- recognition_bootstrap_plan(reference, draws, seed)
  effects <- list()
  audits <- list()
  lopo_rows <- list()
  predictions <- list()
  index <- 1L
  for (phase in recognition_phases) {
    for (method in recognition_methods_for_phase(phase)) {
      message("Recognition sensitivity: ", phase, ", ", method)
      tab <- scores[scores$task == phase & scores$method == method, ]
      tab <- tab[order(tab$participant, tab$item), ]
      fit <- recognition_fit_method_phase(tab, plan)
      effect <- fit$effects
      effect$phase <- phase
      effect$method <- method
      effects[[index]] <- effect
      audit <- fit$audit
      audit$phase <- phase
      audit$method <- method
      audits[[index]] <- audit
      lopo <- as.data.frame(as.table(fit$lopo), stringsAsFactors = FALSE)
      names(lopo) <- c("contrast", "participant", "estimate")
      lopo$phase <- phase
      lopo$method <- method
      lopo_rows[[index]] <- lopo
      predictions[[index]] <- recognition_cell_predictions(
        fit$fixed$coefficients, method, phase
      )
      index <- index + 1L
    }
  }
  effects <- dplyr::bind_rows(effects)
  effects <- effects[c(
    "phase", "method", "contrast", "estimate", "lower_95", "upper_95",
    "lower_family", "upper_family", "valid_draws", "valid_fraction",
    "design_kappa", "lmer_estimate", "lmer_sign_agrees", "lopo_sign_n"
  )]
  audit <- dplyr::bind_rows(audits)
  audit <- audit[c(
    "phase", "method", setdiff(names(audit), c("phase", "method"))
  )]
  lopo <- dplyr::bind_rows(lopo_rows)
  predictions <- dplyr::bind_rows(predictions)
  verdict <- recognition_verdict(effects, audit)
  result <- list(
    effects = effects, audit = audit, lopo = lopo,
    predictions = predictions,
    support = recognition_support_table(scores),
    verdict = verdict,
    config = list(seed = seed, draws = draws, smoke = smoke, plan = plan)
  )
  if (!is.null(output_dir)) recognition_write_results(result, output_dir)
  result
}
