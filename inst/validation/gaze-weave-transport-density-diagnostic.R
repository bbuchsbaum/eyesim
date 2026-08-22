# Post-court Transport-versus-density saliency and confidence diagnostics.

transport_density_seed <- 20260827L
transport_density_draws <- 2000L
transport_density_methods <- c(
  transport = "transport",
  density_sigma_80 = "density_sigma_80_raw_calibrated",
  density_sigma_160 = "density_sigma_160_raw_calibrated",
  density_registered = "density_ridge_registered"
)
transport_density_result_dir <- file.path(
  "inst", "validation", "gaze-weave-transport-v3-retrieval-results",
  "transport-density-diagnostic"
)

transport_density_dependency <- system.file(
  "validation", "gaze-weave-transport-exploratory.R", package = "eyesim"
)
if (!nzchar(transport_density_dependency)) {
  transport_density_dependency <- file.path(
    "inst", "validation", "gaze-weave-transport-exploratory.R"
  )
}
if (!file.exists(transport_density_dependency)) {
  stop("The Transport exploratory helpers are unavailable.")
}
source(transport_density_dependency, local = TRUE)

transport_density_oldness <- function(response) {
  response <- suppressWarnings(as.integer(response))
  value <- ifelse(response %in% 1:4, 5L - response, NA_integer_)
  ordered(value, levels = 1:4)
}

transport_density_verify_manifest <- function(directory) {
  manifest_path <- file.path(directory, "manifest-md5.csv")
  if (!file.exists(manifest_path)) stop("A required result manifest is missing.")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  paths <- file.path(directory, manifest$file)
  if (!all(file.exists(paths)) ||
      !identical(unname(tools::md5sum(paths)), manifest$md5)) {
    stop("A required result manifest does not verify.")
  }
  invisible(manifest)
}

transport_density_read_panel <- function(
    retrieval_dir = dirname(transport_density_result_dir),
    comparator_dir = file.path(
      "inst", "validation", "gaze-weave-recognition-full-cohort-results"
    )) {
  v3 <- transport_exploratory_verify_scores(retrieval_dir)
  transport_density_verify_manifest(comparator_dir)
  comparator <- utils::read.csv(
    file.path(comparator_dir, "trial-scores.csv"), stringsAsFactors = FALSE
  )
  comparator <- comparator[
    comparator$task == "combined" & comparator$condition == "signal" &
      comparator$method %in% unname(transport_density_methods[-1L]),
    c(
      "participant", "item", "outer_fold", "method", "gaze_info_bits",
      "probe_type", "degradation", "accuracy"
    ),
    drop = FALSE
  ]

  full_file <- system.file(
    "validation", "gaze-weave-recognition-full-cohort.R", package = "eyesim"
  )
  if (!nzchar(full_file)) {
    full_file <- file.path(
      "inst", "validation", "gaze-weave-recognition-full-cohort.R"
    )
  }
  full_file <- normalizePath(full_file, mustWork = TRUE)
  source_root <- normalizePath(
    file.path(dirname(full_file), "..", ".."), mustWork = TRUE
  )
  previous_directory <- setwd(source_root)
  on.exit(setwd(previous_directory), add = TRUE)
  environment <- new.env(parent = globalenv())
  sys.source(full_file, envir = environment)
  raw <- environment$full_recognition_read_inputs(verify = TRUE)
  metadata <- unique(raw$retrieval[c(
    "participant", "item", "probe_type", "degradation", "Response",
    "Accuracy"
  )])
  metadata <- metadata[metadata$probe_type %in% c("old", "lure"), ]
  if (anyDuplicated(metadata[c("participant", "item")])) {
    stop("Behavior metadata are not unique by trial.")
  }

  quality <- v3[c(
    "participant", "item", "effective_fixations", "total_duration"
  )]
  comparator <- merge(
    comparator, quality, by = c("participant", "item"),
    all.x = TRUE, sort = FALSE
  )
  comparator <- merge(
    comparator,
    metadata[c("participant", "item", "Response", "Accuracy")],
    by = c("participant", "item"), all.x = TRUE, sort = FALSE
  )
  names(comparator)[names(comparator) == "Accuracy"] <- "accuracy_source"
  if (any(comparator$accuracy != comparator$accuracy_source)) {
    stop("Comparator accuracy differs from the verified raw metadata.")
  }
  comparator$accuracy_source <- NULL

  v3 <- merge(
    v3, metadata, by = c("participant", "item"),
    all.x = TRUE, sort = FALSE
  )
  v3$method <- "transport"
  names(v3)[names(v3) == "Accuracy"] <- "accuracy"
  columns <- c(
    "participant", "item", "outer_fold", "method", "gaze_info_bits",
    "probe_type", "degradation", "accuracy", "Response",
    "effective_fixations", "total_duration"
  )
  panel <- rbind(v3[columns], comparator[columns])
  panel$method <- names(transport_density_methods)[match(
    panel$method, unname(transport_density_methods)
  )]
  panel$saliency_z <- (as.numeric(panel$degradation) - 60) / 20
  panel$z_quality <- as.numeric(scale(log(panel$effective_fixations)))
  panel$z_duration <- as.numeric(scale(panel$total_duration))
  panel$trial_key <- paste(panel$participant, panel$item, sep = ":")
  panel$oldness <- transport_density_oldness(panel$Response)

  support <- table(panel$method)
  if (!identical(sort(names(support)), sort(names(transport_density_methods))) ||
      any(support != 1295L)) {
    stop("Transport and density methods do not share 1,295-trial support.")
  }
  key_sets <- lapply(names(transport_density_methods), function(method) {
    sort(panel$trial_key[panel$method == method])
  })
  if (!all(vapply(key_sets[-1L], identical, logical(1), key_sets[[1L]]))) {
    stop("Transport and density trial keys differ.")
  }
  panel[order(panel$method, panel$participant, panel$item), , drop = FALSE]
}

transport_density_fit_saliency <- function(tab, weights = NULL,
                                           adjusted = FALSE) {
  if (is.null(weights)) weights <- rep(1, nrow(tab))
  formula <- if (adjusted) {
    gaze_info_bits ~ saliency_z + z_quality + z_duration
  } else {
    gaze_info_bits ~ saliency_z
  }
  design <- stats::model.matrix(formula, data = tab)
  keep <- is.finite(weights) & weights > 0
  fit <- stats::lm.wfit(
    design[keep, , drop = FALSE], tab$gaze_info_bits[keep], weights[keep]
  )
  if (fit$rank != ncol(design) ||
      !is.finite(fit$coefficients[["saliency_z"]])) {
    return(NA_real_)
  }
  4 * fit$coefficients[["saliency_z"]]
}

transport_density_saliency_probe <- function(tab, probe, plan) {
  methods <- names(transport_density_methods)
  subset_probe <- function(method) {
    part <- tab[tab$method == method, , drop = FALSE]
    if (!identical(probe, "all")) {
      part <- part[part$probe_type == probe, , drop = FALSE]
    }
    part[order(part$participant, part$item), , drop = FALSE]
  }
  weights <- transport_exploratory_align_weights(
    subset_probe("transport"), plan
  )
  estimates <- bootstrap <- adjusted <- vector("list", length(methods))
  names(estimates) <- names(bootstrap) <- names(adjusted) <- methods
  for (method in methods) {
    part <- subset_probe(method)
    estimates[[method]] <- transport_density_fit_saliency(part)
    adjusted[[method]] <- transport_density_fit_saliency(
      part, adjusted = TRUE
    )
    bootstrap[[method]] <- vapply(seq_len(ncol(weights)), function(draw) {
      transport_density_fit_saliency(part, weights[, draw])
    }, numeric(1))
  }
  method_rows <- do.call(rbind, lapply(methods, function(method) {
    interval <- transport_exploratory_interval(bootstrap[[method]], 0.95)
    data.frame(
      probe = probe,
      method = method,
      trials = nrow(subset_probe(method)),
      saliency_20_to_100_bits = estimates[[method]],
      adjusted_saliency_bits = adjusted[[method]],
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  densities <- setdiff(methods, "transport")
  family_probability <- 1 - 0.05 / (length(densities) * 3)
  contrast_rows <- do.call(rbind, lapply(densities, function(method) {
    values <- bootstrap[["transport"]] - bootstrap[[method]]
    interval <- transport_exploratory_interval(values, 0.95)
    family <- transport_exploratory_interval(values, family_probability)
    data.frame(
      probe = probe,
      contrast = paste("transport_minus", method, sep = "_"),
      estimate = estimates[["transport"]] - estimates[[method]],
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      lower_family = family[["lower"]],
      upper_family = family[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  list(methods = method_rows, contrasts = contrast_rows)
}

transport_density_saliency <- function(tab, draws, seed) {
  results <- lapply(c("all", "old", "lure"), function(probe) {
    reference <- tab[tab$method == "transport", , drop = FALSE]
    if (!identical(probe, "all")) {
      reference <- reference[reference$probe_type == probe, , drop = FALSE]
    }
    reference <- reference[order(reference$participant, reference$item), ]
    plan <- transport_exploratory_bootstrap_plan(
      reference, draws = draws, seed = seed
    )
    transport_density_saliency_probe(tab, probe, plan)
  })
  list(
    methods = do.call(rbind, lapply(results, `[[`, "methods")),
    contrasts = do.call(rbind, lapply(results, `[[`, "contrasts"))
  )
}

transport_density_wide <- function(panel) {
  methods <- names(transport_density_methods)
  reference <- panel[panel$method == "transport", c(
    "participant", "item", "outer_fold", "probe_type", "degradation",
    "Response", "accuracy", "effective_fixations", "total_duration",
    "saliency_z", "trial_key", "oldness"
  )]
  for (method in methods) {
    part <- panel[panel$method == method, c("trial_key", "gaze_info_bits")]
    names(part)[[2L]] <- paste0("score_", method)
    reference <- merge(reference, part, by = "trial_key", all.x = TRUE,
                       sort = FALSE)
  }
  reference[order(reference$participant, reference$item), , drop = FALSE]
}

transport_density_ordinal_probabilities <- function(parameter, design) {
  coefficient_n <- ncol(design)
  coefficient <- parameter[seq_len(coefficient_n)]
  threshold_parameter <- parameter[coefficient_n + 1:3]
  threshold <- c(
    threshold_parameter[[1L]],
    threshold_parameter[[1L]] + exp(threshold_parameter[[2L]]),
    threshold_parameter[[1L]] + exp(threshold_parameter[[2L]]) +
      exp(threshold_parameter[[3L]])
  )
  eta <- as.numeric(design %*% coefficient)
  cumulative <- vapply(threshold, function(value) {
    stats::plogis(value - eta)
  }, numeric(nrow(design)))
  probability <- cbind(
    cumulative[, 1L],
    cumulative[, 2L] - cumulative[, 1L],
    cumulative[, 3L] - cumulative[, 2L],
    1 - cumulative[, 3L]
  )
  probability <- pmax(probability, 1e-12)
  probability <- probability / rowSums(probability)
  colnames(probability) <- as.character(1:4)
  probability
}

transport_density_ordinal_fit <- function(formula, data, penalty = 1) {
  terms <- stats::terms(formula, data = data)
  design <- stats::model.matrix(terms, data = data)
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  response <- as.integer(stats::model.response(
    stats::model.frame(terms, data = data)
  ))
  proportions <- tabulate(response, nbins = 4L) / length(response)
  cumulative <- pmin(pmax(cumsum(proportions)[1:3], 1e-4), 1 - 1e-4)
  threshold <- stats::qlogis(cumulative)
  difference <- pmax(diff(threshold), 1e-3)
  initial <- c(
    rep(0, ncol(design)), threshold[[1L]], log(difference[[1L]]),
    log(difference[[2L]])
  )
  objective <- function(parameter) {
    probability <- transport_density_ordinal_probabilities(parameter, design)
    loss <- -sum(log(probability[cbind(seq_along(response), response)]))
    coefficient <- parameter[seq_len(ncol(design))]
    loss + penalty * sum(coefficient^2) / 2
  }
  fit <- stats::optim(initial, objective, method = "BFGS",
                      control = list(maxit = 2000, reltol = 1e-10))
  if (fit$convergence != 0L || !is.finite(fit$value)) {
    stop("Penalized ordinal model did not converge.")
  }
  structure(list(
    parameter = fit$par,
    terms = terms,
    design_columns = colnames(design),
    penalty = penalty,
    convergence = fit$convergence
  ), class = "transport_density_ordinal_fit")
}

transport_density_ordinal_predict <- function(model, newdata) {
  design <- stats::model.matrix(
    stats::delete.response(model$terms), data = newdata
  )
  design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  design <- design[, model$design_columns, drop = FALSE]
  transport_density_ordinal_probabilities(model$parameter, design)
}

transport_density_ordinal_log_loss <- function(probability, oldness) {
  category <- match(as.character(oldness), colnames(probability))
  if (anyNA(category)) stop("Ordinal prediction omitted an observed category.")
  selected <- probability[cbind(seq_len(nrow(probability)), category)]
  -log(pmax(selected, 1e-12))
}

transport_density_ordinal_loss <- function(model, evaluate) {
  probability <- transport_density_ordinal_predict(model, evaluate)
  transport_density_ordinal_log_loss(probability, evaluate$oldness)
}

transport_density_ordinal_fold <- function(train, evaluate) {
  if (length(intersect(train$participant, evaluate$participant)) ||
      length(intersect(train$item, evaluate$item))) {
    stop("Ordinal training and evaluation support overlaps.")
  }
  train$log_effective_fixations <- log(train$effective_fixations)
  evaluate$log_effective_fixations <- log(evaluate$effective_fixations)
  mappings <- list(
    c("log_effective_fixations", "z_quality_fold"),
    c("total_duration", "z_duration_fold")
  )
  for (method in names(transport_density_methods)) {
    mappings[[length(mappings) + 1L]] <- c(
      paste0("score_", method), paste0("z_", method)
    )
  }
  for (mapping in mappings) {
    value <- transport_exploratory_standardize(
      train, evaluate, mapping[[1L]], mapping[[2L]]
    )
    train <- value$train
    evaluate <- value$evaluate
  }
  base_formula <- oldness ~ saliency_z + z_quality_fold + z_duration_fold
  formulas <- list(base = base_formula)
  for (method in names(transport_density_methods)) {
    formulas[[method]] <- stats::update.formula(
      base_formula, paste(". ~ . +", paste0("z_", method))
    )
  }
  for (density in setdiff(names(transport_density_methods), "transport")) {
    formulas[[paste0("transport_plus_", density)]] <- stats::update.formula(
      base_formula,
      paste(". ~ . + z_transport +", paste0("z_", density))
    )
  }
  fits <- lapply(formulas, transport_density_ordinal_fit, data = train)
  losses <- lapply(fits, transport_density_ordinal_loss, evaluate = evaluate)
  result <- data.frame(
    participant = evaluate$participant,
    item = evaluate$item,
    trial_key = evaluate$trial_key,
    outer_fold = evaluate$outer_fold,
    probe_type = evaluate$probe_type,
    response = evaluate$Response,
    oldness = as.integer(evaluate$oldness),
    stringsAsFactors = FALSE
  )
  for (name in names(losses)) result[[paste0("loss_", name)]] <- losses[[name]]
  result
}

transport_density_ordinal_crossfit <- function(tab, probe) {
  part <- tab[tab$probe_type == probe & !is.na(tab$oldness), , drop = FALSE]
  rows <- lapply(sort(unique(part$outer_fold)), function(fold_id) {
    evaluate <- part[part$outer_fold == fold_id, , drop = FALSE]
    train <- part[
      !part$participant %in% evaluate$participant &
        !part$item %in% evaluate$item, , drop = FALSE
    ]
    transport_density_ordinal_fold(train, evaluate)
  })
  result <- do.call(rbind, rows)
  if (nrow(result) != nrow(part) || anyDuplicated(result$trial_key)) {
    stop("Ordinal cross-fitting did not return one row per response.")
  }
  result[order(result$participant, result$item), , drop = FALSE]
}

transport_density_confidence_summary <- function(predictions, plan) {
  weights <- transport_exploratory_align_weights(predictions, plan)
  methods <- names(transport_density_methods)
  contrasts <- list()
  for (method in methods) {
    contrasts[[paste0(method, "_vs_base")]] <- c("base", method)
  }
  for (density in setdiff(methods, "transport")) {
    joint <- paste0("transport_plus_", density)
    contrasts[[paste0("transport_vs_", density)]] <- c(density, "transport")
    contrasts[[paste0("transport_beyond_", density)]] <- c(density, joint)
    contrasts[[paste0(density, "_beyond_transport")]] <- c("transport", joint)
  }
  family_probability <- 1 - 0.05 / length(contrasts)
  do.call(rbind, lapply(names(contrasts), function(name) {
    pair <- contrasts[[name]]
    gain <- (
      predictions[[paste0("loss_", pair[[1L]])]] -
        predictions[[paste0("loss_", pair[[2L]])]]
    ) / log(2)
    bootstrap <- vapply(seq_len(ncol(weights)), function(draw) {
      transport_exploratory_weighted_mean(gain, weights[, draw])
    }, numeric(1))
    interval <- transport_exploratory_interval(bootstrap, 0.95)
    family <- transport_exploratory_interval(bootstrap, family_probability)
    data.frame(
      probe_type = predictions$probe_type[[1L]],
      contrast = name,
      trials = nrow(predictions),
      information_gain_bits = mean(gain),
      lower_95 = interval[["lower"]],
      upper_95 = interval[["upper"]],
      lower_family = family[["lower"]],
      upper_family = family[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
}

transport_density_response_counts <- function(tab) {
  reference <- tab[tab$method == "transport", ]
  counts <- as.data.frame(table(
    probe_type = reference$probe_type,
    response = reference$Response,
    useNA = "ifany"
  ), stringsAsFactors = FALSE)
  counts[counts$Freq > 0L, , drop = FALSE]
}

transport_density_saliency_means <- function(tab) {
  stats::aggregate(
    gaze_info_bits ~ method + probe_type + degradation,
    data = tab, FUN = mean
  )
}

transport_density_fold_audit <- function(predictions) {
  do.call(rbind, lapply(names(predictions), function(probe) {
    tab <- predictions[[probe]]
    do.call(rbind, lapply(names(transport_density_methods), function(method) {
      gain <- (
        tab$loss_base - tab[[paste0("loss_", method)]]
      ) / log(2)
      data.frame(
        probe_type = probe,
        method = method,
        outer_fold = sort(unique(tab$outer_fold)),
        trials = as.integer(table(tab$outer_fold)),
        information_gain_bits = as.numeric(tapply(
          gain, tab$outer_fold, mean
        )),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

transport_density_write <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(object, name) {
    utils::write.csv(object, file.path(output_dir, name), row.names = FALSE)
  }
  write(result$saliency$methods, "saliency-slopes.csv")
  write(result$saliency$contrasts, "saliency-contrasts.csv")
  write(result$saliency_means, "saliency-means.csv")
  write(result$confidence, "confidence-prediction.csv")
  write(result$confidence_fold_audit, "confidence-fold-audit.csv")
  write(result$response_counts, "response-counts.csv")
  write(result$configuration, "configuration.csv")
  saveRDS(
    result$predictions,
    file.path(output_dir, "confidence-predictions.rds"), version = 3
  )
  saveRDS(
    result[c(
      "saliency", "saliency_means", "confidence", "confidence_fold_audit",
      "response_counts", "configuration"
    )],
    file.path(output_dir, "scientific-results.rds"), version = 3
  )
  files <- sort(list.files(output_dir, full.names = TRUE))
  files <- files[basename(files) != "manifest-md5.csv"]
  write(data.frame(
    file = basename(files), md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  ), "manifest-md5.csv")
  invisible(result)
}

run_gaze_weave_transport_density_diagnostic <- function(
    output_dir = transport_density_result_dir,
    draws = transport_density_draws,
    seed = transport_density_seed) {
  panel <- transport_density_read_panel()
  saliency <- transport_density_saliency(panel, draws, seed)
  wide <- transport_density_wide(panel)
  predictions <- lapply(c("old", "lure"), function(probe) {
    message("Ordinal response cross-fit: ", probe)
    transport_density_ordinal_crossfit(wide, probe)
  })
  names(predictions) <- c("old", "lure")
  confidence <- do.call(rbind, lapply(names(predictions), function(probe) {
    plan <- transport_exploratory_bootstrap_plan(
      predictions[[probe]], draws = draws, seed = seed
    )
    transport_density_confidence_summary(predictions[[probe]], plan)
  }))
  result <- list(
    saliency = saliency,
    saliency_means = transport_density_saliency_means(panel),
    confidence = confidence,
    confidence_fold_audit = transport_density_fold_audit(predictions),
    predictions = predictions,
    response_counts = transport_density_response_counts(panel),
    configuration = data.frame(
      status = "post_court_exploratory_diagnostic",
      seed = seed,
      bootstrap_draws = draws,
      trials = 1295L,
      old_responses = nrow(predictions$old),
      lure_responses = nrow(predictions$lure),
      response_scale = "newness_1_to_4",
      frozen_primary_unchanged = TRUE,
      local_only = TRUE,
      stringsAsFactors = FALSE
    )
  )
  transport_density_write(result, output_dir)
}
