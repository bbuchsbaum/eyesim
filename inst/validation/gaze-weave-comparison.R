# Reproducible comparative validation for GazeWeave.
#
# Run from a source checkout after loading eyesim:
#   devtools::load_all()
#   source("inst/validation/gaze-weave-comparison.R")
#   run_gaze_weave_validation(output_dir = "inst/validation/results")

validation_fixation_group <- function(coords, duration = NULL, time_scale = 1) {
  if (is.null(duration)) {
    duration <- rep(1, nrow(coords))
  }
  onset <- c(0, head(cumsum(duration), -1)) * time_scale
  fixation_group(
    x = coords[, 1],
    y = coords[, 2],
    duration = duration * time_scale,
    onset = onset
  )
}

validation_rotate <- function(points, angle) {
  rotation <- matrix(
    c(cos(angle), -sin(angle), sin(angle), cos(angle)),
    nrow = 2,
    byrow = TRUE
  )
  points %*% t(rotation)
}

validation_unique_permutations <- function(n_items, n_points) {
  permutations <- list(seq_len(n_points))
  while (length(permutations) < n_items) {
    candidate <- sample.int(n_points)
    duplicate <- any(vapply(
      permutations,
      function(existing) identical(existing, candidate),
      logical(1)
    ))
    if (!duplicate) {
      permutations[[length(permutations) + 1L]] <- candidate
    }
  }
  permutations
}

validation_spec <- function(warp = gaze_warp_none()) {
  gaze_weave_spec(
    spatial = gaze_gaussian_mixture(
      sigmas = c(2.5, 5, 10),
      weights = c(0.45, 0.35, 0.20),
      unit = "px"
    ),
    chronology = gaze_local_order(horizon = 0.3),
    transport = gaze_partial_transport(
      temporal_weight = 4,
      unmatched = c(2, 2),
      entropy = 0.03
    ),
    warp = warp,
    screen = gaze_screen(100, 100, unit = "px"),
    control = gaze_weave_control(maxit = 250, multistart = 2)
  )
}

validation_spatial_templates <- function(n_items, phase_shift = 0,
                                         motif_scale = 1) {
  motif <- rbind(
    c(-8, -4), c(-3, 7), c(4, 3), c(8, -6), c(1, -9), c(-6, 2)
  )
  angles <- seq(0, 2 * pi, length.out = n_items + 1L)[-1L] + phase_shift
  lapply(seq_len(n_items), function(item) {
    anchor <- c(50, 50) + 24 * c(cos(angles[[item]]), sin(angles[[item]]))
    shape <- validation_rotate(motif * motif_scale, angles[[item]] * 0.7)
    shape[, 1] <- shape[, 1] * (0.85 + item / (4 * n_items))
    validation_fixation_group(
      sweep(shape, 2, anchor, FUN = "+"),
      duration = c(1, 1.8, 0.8, 1.4, 1.1, 0.9)
    )
  })
}

validation_fit_warp <- function(reference, source, spec) {
  ref_tab <- tibble::tibble(
    calibration_id = seq_along(reference),
    fixgroup = reference
  )
  source_tab <- tibble::tibble(
    calibration_id = seq_along(source),
    fixgroup = source
  )
  eyesim:::fit_gaze_warp_model(
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = "calibration_id",
    refvar = "fixgroup",
    sourcevar = "fixgroup",
    spec = spec
  )
}

validation_register_path <- function(path, warp_model) {
  parameters <- eyesim:::warp_parameters(warp_model)
  coords <- cbind(path$x, path$y) %*% t(parameters$A)
  coords <- sweep(coords, 2, parameters$translation, FUN = "+")
  fixation_group(
    x = coords[, 1],
    y = coords[, 2],
    duration = path$duration,
    onset = path$onset
  )
}

validation_geometry_case <- function(participant, n_items) {
  scale_true <- stats::rnorm(1, mean = 1.25, sd = 0.04)
  translation_true <- stats::rnorm(2, mean = c(-6, 4), sd = 0.8)
  spec <- validation_spec(
    gaze_warp_contraction(center = c(50, 50), translation = TRUE)
  )

  calibration_reference <- validation_spatial_templates(
    n_items,
    phase_shift = pi / n_items,
    motif_scale = 0.9
  )
  calibration_source <- lapply(calibration_reference, function(path) {
    coords <- sweep(cbind(path$x, path$y), 2, translation_true, FUN = "-") /
      scale_true
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 0.25), ncol = 2)
    validation_fixation_group(coords, path$duration, time_scale = 1.6)
  })
  warp_model <- validation_fit_warp(
    calibration_reference,
    calibration_source,
    spec
  )

  reference <- validation_spatial_templates(n_items)
  source <- lapply(seq_along(reference), function(item) {
    coords <- sweep(cbind(reference[[item]]$x, reference[[item]]$y),
                    2, translation_true, FUN = "-") / scale_true
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 0.7), ncol = 2)
    validation_fixation_group(
      coords,
      reference[[item]]$duration,
      time_scale = stats::runif(1, 1.4, 2.2)
    )
  })

  list(
    participant = participant,
    scenario = "registered_geometry",
    reference = reference,
    source = source,
    mappings = lapply(reference, function(path) seq_len(nrow(path))),
    omitted = lapply(reference, function(path) integer()),
    spec = spec,
    warp_model = warp_model,
    true_scale = scale_true,
    fitted_scale = eyesim:::warp_parameters(warp_model)$scale
  )
}

validation_order_case <- function(participant, n_items, partial = FALSE) {
  common_points <- rbind(
    c(24, 24), c(39, 71), c(53, 31),
    c(72, 66), c(81, 19), c(18, 53)
  )
  common_points <- common_points +
    matrix(stats::rnorm(length(common_points), sd = 1.2), ncol = 2)
  permutations <- validation_unique_permutations(n_items, nrow(common_points))
  reference <- lapply(permutations, function(permutation) {
    validation_fixation_group(
      common_points[permutation, , drop = FALSE],
      duration = rep(1, nrow(common_points))
    )
  })

  source_records <- lapply(seq_len(n_items), function(item) {
    coords <- cbind(reference[[item]]$x, reference[[item]]$y)
    mapping <- seq_len(nrow(coords))
    missing <- integer()
    if (partial) {
      missing <- sort(sample(2:(nrow(coords) - 1L), 2L))
      keep <- setdiff(seq_len(nrow(coords)), missing)
      coords <- coords[keep, , drop = FALSE]
      mapping <- keep
      insertion <- (coords[2, ] + coords[3, ]) / 2 + stats::rnorm(2, sd = 2)
      coords <- rbind(coords[1:2, , drop = FALSE], insertion,
                      coords[3:nrow(coords), , drop = FALSE])
      mapping <- append(mapping, NA_integer_, after = 2L)
    }
    coords <- coords + matrix(stats::rnorm(length(coords), sd = 0.8), ncol = 2)
    list(
      path = validation_fixation_group(
        coords,
        duration = rep(1, nrow(coords)),
        time_scale = stats::runif(1, 1.4, 2.2)
      ),
      mapping = mapping,
      omitted = missing
    )
  })
  source <- lapply(source_records, `[[`, "path")
  mappings <- lapply(source_records, `[[`, "mapping")
  omitted <- lapply(source_records, `[[`, "omitted")

  list(
    participant = participant,
    scenario = if (partial) "partial_replay" else "order_at_fixed_density",
    reference = reference,
    source = source,
    mappings = mappings,
    omitted = omitted,
    spec = validation_spec(),
    warp_model = eyesim:::identity_gaze_warp_model(),
    true_scale = 1,
    fitted_scale = 1
  )
}

validation_density_map <- function(path, sigma, bounds = c(0, 100), grid_n = 35L) {
  grid <- seq(bounds[[1]], bounds[[2]], length.out = grid_n)
  locations <- as.matrix(expand.grid(x = grid, y = grid))
  coords <- cbind(path$x, path$y)
  mass <- path$duration / sum(path$duration)
  squared_distance <- outer(
    locations[, 1], coords[, 1], FUN = "-"
  )^2 + outer(locations[, 2], coords[, 2], FUN = "-")^2
  density <- as.numeric(exp(-squared_distance / (2 * sigma^2)) %*% mass)
  total <- sum(density)
  if (!is.finite(total) || total <= 0) {
    stop("Density validation baseline produced invalid mass.")
  }
  density / total
}

validation_density_signature <- function(path, sigmas = c(2.5, 5, 10)) {
  lapply(sigmas, function(sigma) validation_density_map(path, sigma))
}

validation_cosine <- function(x, y) {
  denominator <- sqrt(sum(x^2) * sum(y^2))
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  sum(x * y) / denominator
}

validation_density_similarity <- function(x, y) {
  mean(mapply(validation_cosine, x, y))
}

validation_score_density <- function(reference, source) {
  ref_density <- lapply(reference, validation_density_signature)
  source_density <- lapply(source, validation_density_signature)
  matrix(
    vapply(seq_along(source_density), function(source_index) {
      vapply(seq_along(ref_density), function(ref_index) {
        validation_density_similarity(
          source_density[[source_index]], ref_density[[ref_index]]
        )
      }, numeric(1))
    }, numeric(length(reference))),
    nrow = length(source),
    byrow = TRUE
  )
}

validation_score_multimatch <- function(reference, source) {
  metric_names <- c(
    "mm_vector", "mm_direction", "mm_length", "mm_position",
    "mm_duration", "mm_position_emd"
  )
  arrays <- lapply(metric_names, function(name) {
    matrix(NA_real_, nrow = length(source), ncol = length(reference))
  })
  names(arrays) <- metric_names
  for (source_index in seq_along(source)) {
    source_scanpath <- scanpath(source[[source_index]])
    for (ref_index in seq_along(reference)) {
      scores <- multi_match(
        scanpath(reference[[ref_index]]),
        source_scanpath,
        screensize = c(100, 100)
      )
      for (name in metric_names) {
        arrays[[name]][source_index, ref_index] <- scores[[name]]
      }
    }
  }
  arrays
}

validation_gaze_diagnostics <- function(alignment, mapping, omitted) {
  genuine <- which(!is.na(mapping))
  linked_mass <- if (length(genuine) == 0L) {
    NA_real_
  } else {
    sum(alignment$coupling[cbind(mapping[genuine], genuine)]) /
      sum(alignment$coupling[, genuine, drop = FALSE])
  }

  omission_auc <- NA_real_
  if (length(omitted) > 0L) {
    retained <- setdiff(seq_along(alignment$reference$mass), omitted)
    missing_mass <- pmax(
      alignment$reference$mass - alignment$reference_marginal,
      0
    )
    comparisons <- outer(
      missing_mass[omitted], missing_mass[retained], FUN = "-"
    )
    omission_auc <- mean((comparisons > 0) + 0.5 * (comparisons == 0))
  }

  c(
    coupling_link_mass = linked_mass,
    omission_auc = omission_auc,
    transported_mass = alignment$transported_mass,
    converged = as.numeric(isTRUE(alignment$convergence$converged))
  )
}

validation_score_case <- function(case) {
  registered_source <- lapply(
    case$source,
    validation_register_path,
    warp_model = case$warp_model
  )

  gaze_scores <- matrix(
    NA_real_,
    nrow = length(case$source),
    ncol = length(case$reference)
  )
  diagnostics <- vector("list", length(case$source))
  for (source_index in seq_along(case$source)) {
    for (ref_index in seq_along(case$reference)) {
      alignment <- gaze_align(
        case$reference[[ref_index]],
        case$source[[source_index]],
        spec = case$spec,
        warp_model = case$warp_model
      )
      gaze_scores[source_index, ref_index] <- -alignment$energy
      if (source_index == ref_index) {
        diagnostics[[source_index]] <- validation_gaze_diagnostics(
          alignment,
          mapping = case$mappings[[source_index]],
          omitted = case$omitted[[source_index]]
        )
      }
    }
  }

  mm_raw <- validation_score_multimatch(case$reference, case$source)
  mm_registered <- validation_score_multimatch(
    case$reference,
    registered_source
  )
  scores <- list(
    gaze_weave = gaze_scores,
    density_raw = validation_score_density(case$reference, case$source),
    density_registered = validation_score_density(
      case$reference,
      registered_source
    )
  )
  for (name in names(mm_raw)) {
    scores[[paste0(name, "_raw")]] <- mm_raw[[name]]
    scores[[paste0(name, "_registered")]] <- mm_registered[[name]]
  }

  list(
    participant = case$participant,
    scenario = case$scenario,
    scores = scores,
    target = seq_along(case$source),
    diagnostics = diagnostics,
    true_scale = case$true_scale,
    fitted_scale = case$fitted_scale
  )
}

validation_log_mean_exp <- function(x) {
  maximum <- max(x)
  maximum + log(mean(exp(x - maximum)))
}

validation_score_tolerance <- function(x) {
  1e-10 * max(1, max(abs(x), na.rm = TRUE))
}

validation_row_statistics <- function(scores, target, gaze = FALSE) {
  vapply(seq_len(nrow(scores)), function(row) {
    true_score <- scores[row, target[[row]]]
    false_scores <- scores[row, -target[[row]]]
    value <- if (gaze) {
      (true_score - validation_log_mean_exp(false_scores)) / log(2)
    } else {
      true_score - mean(false_scores)
    }
    if (abs(value) <= validation_score_tolerance(scores[row, ])) 0 else value
  }, numeric(1))
}

validation_retrieval_statistics <- function(scores, target) {
  row_values <- lapply(seq_len(nrow(scores)), function(row) {
    true_score <- scores[row, target[[row]]]
    false_scores <- scores[row, -target[[row]]]
    delta <- true_score - false_scores
    tolerance <- validation_score_tolerance(scores[row, ])
    better <- sum(delta < -tolerance)
    tied <- sum(abs(delta) <= tolerance)
    c(
      pairwise_auc = mean((delta > tolerance) +
                            0.5 * (abs(delta) <= tolerance)),
      top1_credit = if (better == 0L) 1 / (tied + 1L) else 0,
      reciprocal_rank = 1 / (1 + better + tied / 2)
    )
  })
  colMeans(do.call(rbind, row_values))
}

validation_summarize_case <- function(scored_case) {
  rows <- lapply(names(scored_case$scores), function(method) {
    scores <- scored_case$scores[[method]]
    retrieval <- validation_retrieval_statistics(scores, scored_case$target)
    statistic <- mean(validation_row_statistics(
      scores,
      scored_case$target,
      gaze = identical(method, "gaze_weave")
    ))
    data.frame(
      participant = scored_case$participant,
      scenario = scored_case$scenario,
      method = method,
      pairwise_auc = retrieval[["pairwise_auc"]],
      top1_credit = retrieval[["top1_credit"]],
      reciprocal_rank = retrieval[["reciprocal_rank"]],
      participant_statistic = statistic,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validation_candidate_statistics <- function(scores, gaze = FALSE) {
  result <- matrix(NA_real_, nrow = nrow(scores), ncol = ncol(scores))
  for (candidate in seq_len(ncol(scores))) {
    result[, candidate] <- validation_row_statistics(
      scores,
      rep(candidate, nrow(scores)),
      gaze = gaze
    )
  }
  result
}

validation_case_statistic_matrices <- function(scored_case) {
  result <- lapply(names(scored_case$scores), function(method) {
    validation_candidate_statistics(
      scored_case$scores[[method]],
      gaze = identical(method, "gaze_weave")
    )
  })
  stats::setNames(result, names(scored_case$scores))
}

validation_null_study_distribution <- function(
    statistic_matrices, sampled, methods, draws) {
  distribution <- matrix(
    0,
    nrow = draws,
    ncol = length(methods),
    dimnames = list(NULL, methods)
  )
  for (position in seq_along(sampled)) {
    case_index <- sampled[[position]]
    case_matrices <- statistic_matrices[[case_index]]
    n_rows <- nrow(case_matrices[[1]])
    n_candidates <- ncol(case_matrices[[1]])
    target_draws <- matrix(
      sample.int(n_candidates, draws * n_rows, replace = TRUE),
      nrow = draws,
      ncol = n_rows
    )
    for (method_index in seq_along(methods)) {
      statistic_matrix <- case_matrices[[methods[[method_index]]]]
      selected <- vapply(seq_len(n_rows), function(row) {
        statistic_matrix[row, target_draws[, row]]
      }, numeric(draws))
      distribution[, method_index] <- distribution[, method_index] +
        rowMeans(selected) / length(sampled)
    }
  }
  distribution
}

validation_observed_study <- function(
    statistic_matrices, scored_cases, sampled, methods,
    signal_probability) {
  observed <- stats::setNames(numeric(length(methods)), methods)
  for (position in seq_along(sampled)) {
    case_index <- sampled[[position]]
    case <- scored_cases[[case_index]]
    case_matrices <- statistic_matrices[[case_index]]
    n_rows <- nrow(case_matrices[[1]])
    retain_signal <- stats::runif(n_rows) < signal_probability
    random_target <- sample.int(ncol(case_matrices[[1]]), n_rows, replace = TRUE)
    target <- ifelse(retain_signal, case$target, random_target)
    for (method in methods) {
      statistic_matrix <- case_matrices[[method]]
      observed[[method]] <- observed[[method]] + mean(
        statistic_matrix[cbind(seq_len(n_rows), target)]
      ) / length(sampled)
    }
  }
  observed
}

validation_randomization_operating_characteristics <- function(
    scored_cases, sample_sizes = c(8L, 16L, 32L),
    repetitions = 500L, randomization_draws = 199L,
    signal_probabilities = c(0.25, 0.5, 0.75, 1),
    alpha = 0.05, seed = 20260815L) {
  if (randomization_draws < 39L) {
    stop("randomization_draws must be at least 39 for an alpha-0.05 test.")
  }
  if (!is.numeric(signal_probabilities) ||
      any(!is.finite(signal_probabilities)) ||
      any(signal_probabilities < 0 | signal_probabilities > 1)) {
    stop("signal_probabilities must be finite values between zero and one.")
  }
  set.seed(seed)
  methods <- names(scored_cases[[1]]$scores)
  statistic_matrices <- lapply(scored_cases, validation_case_statistic_matrices)
  family_methods <- list(
    multimatch_raw_holm_any = paste0(
      c("mm_vector", "mm_direction", "mm_length", "mm_position",
        "mm_duration", "mm_position_emd"),
      "_raw"
    ),
    multimatch_registered_holm_any = paste0(
      c("mm_vector", "mm_direction", "mm_length", "mm_position",
        "mm_duration", "mm_position_emd"),
      "_registered"
    )
  )
  family_methods <- Filter(
    function(family) all(family %in% methods),
    family_methods
  )
  rows <- list()
  row_index <- 1L
  for (signal_probability in signal_probabilities) {
    for (sample_size in sample_sizes) {
      alternative_rejections <- matrix(
        FALSE, nrow = repetitions, ncol = length(methods),
        dimnames = list(NULL, methods)
      )
      null_rejections <- alternative_rejections
      null_observed <- matrix(
        NA_real_, nrow = repetitions, ncol = length(methods),
        dimnames = list(NULL, methods)
      )
      family_alternative <- matrix(
        FALSE, nrow = repetitions, ncol = length(family_methods),
        dimnames = list(NULL, names(family_methods))
      )
      family_null <- family_alternative

      for (iteration in seq_len(repetitions)) {
        sampled <- sample.int(length(scored_cases), sample_size, replace = TRUE)
        observed_alternative <- validation_observed_study(
          statistic_matrices,
          scored_cases = scored_cases,
          sampled = sampled,
          methods = methods,
          signal_probability = signal_probability
        )
        null_distribution <- validation_null_study_distribution(
          statistic_matrices,
          sampled = sampled,
          methods = methods,
          draws = randomization_draws
        )
        observed_null <- null_distribution[1, ]
        reference_null <- null_distribution[-1, , drop = FALSE]
        alternative_p <- (1 + colSums(sweep(
          null_distribution,
          2,
          observed_alternative,
          FUN = ">="
        ))) / (randomization_draws + 1)
        null_p <- (1 + colSums(sweep(
          reference_null,
          2,
          observed_null,
          FUN = ">="
        ))) / randomization_draws

        alternative_rejections[iteration, ] <- alternative_p <= alpha
        null_rejections[iteration, ] <- null_p <= alpha
        null_observed[iteration, ] <- observed_null
        for (family in names(family_methods)) {
          family_alternative[iteration, family] <- any(
            stats::p.adjust(
              alternative_p[family_methods[[family]]],
              method = "holm"
            ) <= alpha
          )
          family_null[iteration, family] <- any(
            stats::p.adjust(
              null_p[family_methods[[family]]],
              method = "holm"
            ) <= alpha
          )
        }
      }

      for (method in methods) {
        rows[[row_index]] <- data.frame(
          method = method,
          signal_probability = signal_probability,
          sample_size = sample_size,
          repetitions = repetitions,
          randomization_draws = randomization_draws,
          estimated_power = mean(alternative_rejections[, method]),
          estimated_type1 = mean(null_rejections[, method]),
          mean_null_statistic = mean(null_observed[, method]),
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
      for (family in names(family_methods)) {
        rows[[row_index]] <- data.frame(
          method = family,
          signal_probability = signal_probability,
          sample_size = sample_size,
          repetitions = repetitions,
          randomization_draws = randomization_draws,
          estimated_power = mean(family_alternative[, family]),
          estimated_type1 = mean(family_null[, family]),
          mean_null_statistic = NA_real_,
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }
  }
  do.call(rbind, rows)
}

validation_diagnostic_table <- function(scored_cases) {
  do.call(rbind, lapply(scored_cases, function(case) {
    diagnostics <- do.call(rbind, case$diagnostics)
    data.frame(
      participant = case$participant,
      scenario = case$scenario,
      image_id = seq_len(nrow(diagnostics)),
      coupling_link_mass = diagnostics[, "coupling_link_mass"],
      omission_auc = diagnostics[, "omission_auc"],
      transported_mass = diagnostics[, "transported_mass"],
      converged = diagnostics[, "converged"],
      true_scale = case$true_scale,
      fitted_scale = case$fitted_scale,
      stringsAsFactors = FALSE
    )
  }))
}

validation_aggregate_retrieval <- function(participant_results) {
  aggregate(
    cbind(pairwise_auc, top1_credit, reciprocal_rank, participant_statistic) ~
      scenario + method,
    data = participant_results,
    FUN = mean
  )
}

validation_write_results <- function(result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    result$participant_results,
    file.path(output_dir, "participant-results.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result$retrieval,
    file.path(output_dir, "retrieval-summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result$power,
    file.path(output_dir, "power-summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    result$diagnostics,
    file.path(output_dir, "diagnostic-recovery.csv"),
    row.names = FALSE
  )
  config <- data.frame(
    parameter = names(result$config),
    value = vapply(result$config, function(value) {
      paste(value, collapse = ",")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    config,
    file.path(output_dir, "run-config.csv"),
    row.names = FALSE
  )
  utils::capture.output(
    utils::sessionInfo(),
    file = file.path(output_dir, "session-info.txt")
  )
  invisible(result)
}

run_gaze_weave_validation <- function(
    n_participants = 24L, n_items = 6L,
    sample_sizes = c(8L, 16L, 32L), bootstrap_repetitions = 500L,
    randomization_draws = 199L,
    signal_probabilities = c(0.25, 0.5, 0.75, 1),
    seed = 20260815L, output_dir = NULL, verbose = interactive()) {
  if (n_items < 4L || n_items > factorial(6)) {
    stop("n_items must be between four and 720.")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Comparative validation requires the suggested package 'igraph'.")
  }
  emd_available <- vapply(
    c("emdist", "T4transport", "transport"),
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
  if (!any(emd_available)) {
    stop("Comparative validation requires emdist, T4transport, or transport for MultiMatch EMD.")
  }
  set.seed(seed)
  scenario_generators <- list(
    registered_geometry = function(participant) {
      validation_geometry_case(participant, n_items)
    },
    order_at_fixed_density = function(participant) {
      validation_order_case(participant, n_items, partial = FALSE)
    },
    partial_replay = function(participant) {
      validation_order_case(participant, n_items, partial = TRUE)
    }
  )

  scored_by_scenario <- list()
  for (scenario in names(scenario_generators)) {
    if (isTRUE(verbose)) {
      message("Scoring ", scenario, " ...")
    }
    scored_by_scenario[[scenario]] <- lapply(seq_len(n_participants), function(participant) {
      scored <- validation_score_case(scenario_generators[[scenario]](participant))
      if (isTRUE(verbose) &&
          (participant %% 4L == 0L || participant == n_participants)) {
        message("  participants scored: ", participant, "/", n_participants)
      }
      scored
    })
  }

  participant_results <- do.call(rbind, lapply(scored_by_scenario, function(cases) {
    do.call(rbind, lapply(cases, validation_summarize_case))
  }))
  retrieval <- validation_aggregate_retrieval(participant_results)
  diagnostics <- validation_diagnostic_table(unlist(scored_by_scenario, recursive = FALSE))
  power <- do.call(rbind, lapply(names(scored_by_scenario), function(scenario) {
    if (isTRUE(verbose)) {
      message("Estimating operating characteristics for ", scenario, " ...")
    }
    cases <- scored_by_scenario[[scenario]]
    combined <- validation_randomization_operating_characteristics(
      cases,
      sample_sizes = sample_sizes,
      repetitions = bootstrap_repetitions,
      randomization_draws = randomization_draws,
      signal_probabilities = signal_probabilities,
      seed = seed + match(scenario, names(scored_by_scenario))
    )
    combined$scenario <- scenario
    combined
  }))

  result <- list(
    config = list(
      seed = seed,
      n_participants = n_participants,
      n_items = n_items,
      sample_sizes = sample_sizes,
      bootstrap_repetitions = bootstrap_repetitions,
      randomization_draws = randomization_draws,
      signal_probabilities = signal_probabilities,
      density_sigmas_px = c(2.5, 5, 10),
      chronology_horizon = 0.3,
      temporal_weight = 4,
      unmatched = c(2, 2),
      entropy = 0.03
    ),
    participant_results = participant_results,
    retrieval = retrieval,
    power = power,
    diagnostics = diagnostics
  )
  if (!is.null(output_dir)) {
    validation_write_results(result, output_dir)
  }
  result
}
