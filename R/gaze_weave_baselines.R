# Fair predictive baselines for GazeWeave -----------------------------------

gaze_namespace_available <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

#' Specify the fair GazeWeave comparator court
#'
#' The comparator court preserves the native MultiMatch dimensions, uses a
#' fixed multiscale density score, and includes order-free elastic-consensus
#' matching. Ridge composites for registered MultiMatch and density features
#' are tuned in inner folds of each outer training split.
#'
#' @param screen Screen geometry created by [gaze_screen()].
#' @param density_sigmas Positive density bandwidths in `screen$unit`.
#' @param density_grid Number of grid points per spatial axis.
#' @param methods Any of `"multimatch"`, `"density"`, and `"elastic"`.
#' @param warp A candidate-invariant warp specification.
#' @param lambda_grid Non-negative ridge penalties.
#' @param inner_folds Number of inner folds used to select the ridge penalty.
#' @param elastic_radii Named consensus, rigidity, and matching radii in
#'   `screen$unit`. Defaults follow the 2, 10, and 1 degree values proposed by
#'   Wang, Holmqvist, and Alexa when the declared unit is degrees.
#' @param elastic_maxit Maximum elastic-consensus iterations.
#' @param elastic_tolerance Convergence tolerance in screen coordinates.
#'
#' @return A frozen `gaze_baseline_spec`.
#' @export
gaze_baseline_spec <- function(
    screen,
    density_sigmas,
    density_grid = 32L,
    methods = c("multimatch", "density", "elastic"),
    warp = gaze_warp_none(),
    lambda_grid = c(0, 0.01, 0.1, 1, 10),
    inner_folds = 2L,
    elastic_radii = c(consensus = 2, rigidity = 10, matching = 1),
    elastic_maxit = 25L,
    elastic_tolerance = 1e-5) {
  if (!inherits(screen, "gaze_screen")) {
    stop("screen must be created by gaze_screen().")
  }
  if (!is.numeric(density_sigmas) || length(density_sigmas) == 0L ||
      any(!is.finite(density_sigmas)) || any(density_sigmas <= 0)) {
    stop("density_sigmas must contain finite positive values.")
  }
  density_grid <- as.integer(density_grid)
  if (length(density_grid) != 1L || is.na(density_grid) || density_grid < 8L) {
    stop("density_grid must be an integer of at least eight.")
  }
  methods <- unique(match.arg(
    methods, c("multimatch", "density", "elastic"), several.ok = TRUE
  ))
  if (!inherits(warp, "gaze_warp_spec")) {
    stop("warp must be created by gaze_warp_none() or gaze_warp_contraction().")
  }
  if (identical(warp$type, "contraction") &&
      is.character(warp$center) && is.null(screen)) {
    stop("A screen-centred contraction requires screen geometry.")
  }
  if (!is.numeric(lambda_grid) || length(lambda_grid) == 0L ||
      any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) {
    stop("lambda_grid must contain finite non-negative values.")
  }
  inner_folds <- as.integer(inner_folds)
  if (length(inner_folds) != 1L || is.na(inner_folds) || inner_folds < 2L) {
    stop("inner_folds must be an integer of at least two.")
  }
  required_radii <- c("consensus", "rigidity", "matching")
  if (!is.numeric(elastic_radii) ||
      !all(required_radii %in% names(elastic_radii)) ||
      any(!is.finite(elastic_radii[required_radii])) ||
      any(elastic_radii[required_radii] <= 0)) {
    stop("elastic_radii must contain positive consensus, rigidity, and matching values.")
  }
  elastic_maxit <- as.integer(elastic_maxit)
  if (length(elastic_maxit) != 1L || is.na(elastic_maxit) ||
      elastic_maxit < 1L || !is.numeric(elastic_tolerance) ||
      length(elastic_tolerance) != 1L || !is.finite(elastic_tolerance) ||
      elastic_tolerance <= 0) {
    stop("Invalid elastic solver controls.")
  }

  structure(
    list(
      screen = screen,
      density_sigmas = sort(unique(as.numeric(density_sigmas))),
      density_grid = density_grid,
      methods = methods,
      warp = warp,
      chronology = gaze_order_neighbours(),
      lambda_grid = sort(unique(as.numeric(lambda_grid))),
      inner_folds = inner_folds,
      elastic = list(
        radii = stats::setNames(
          as.numeric(elastic_radii[required_radii]), required_radii
        ),
        maxit = elastic_maxit,
        tolerance = as.numeric(elastic_tolerance)
      ),
      version = 1L
    ),
    class = c("gaze_baseline_spec", "list")
  )
}

#' Report comparator availability
#'
#' Optional MultiMatch dependencies are checked explicitly. Density and the
#' built-in elastic-consensus reference implementation have no optional
#' dependency.
#'
#' @param spec A [gaze_baseline_spec()].
#'
#' @return A data frame with one row per requested comparator.
#' @export
gaze_baseline_availability <- function(spec) {
  if (!inherits(spec, "gaze_baseline_spec")) {
    stop("spec must be created by gaze_baseline_spec().")
  }
  rows <- lapply(spec$methods, function(method) {
    available <- TRUE
    reason <- NA_character_
    if (identical(method, "multimatch")) {
      graph_available <- gaze_namespace_available("igraph")
      emd_available <- any(vapply(
        c("emdist", "T4transport", "transport"),
        gaze_namespace_available,
        logical(1)
      ))
      available <- graph_available && emd_available
      if (!available) {
        missing <- c(
          if (!graph_available) "igraph" else character(),
          if (!emd_available) "one of emdist, T4transport, or transport" else character()
        )
        reason <- paste("missing optional dependency:", paste(missing, collapse = "; "))
      }
    }
    data.frame(
      method = method,
      available = available,
      reason = reason,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Density features ----------------------------------------------------------

baseline_density_signature <- function(path, spec) {
  valid <- is.finite(path$x) & is.finite(path$y) & is.finite(path$duration) &
    path$duration > 0
  if (!any(valid)) stop("Density baselines require positive finite fixation mass.")
  coords <- cbind(path$x[valid], path$y[valid])
  mass <- path$duration[valid] / sum(path$duration[valid])
  x_grid <- seq(spec$screen$xlim[[1]], spec$screen$xlim[[2]],
                length.out = spec$density_grid)
  y_grid <- seq(spec$screen$ylim[[1]], spec$screen$ylim[[2]],
                length.out = spec$density_grid)
  locations <- as.matrix(expand.grid(x = x_grid, y = y_grid))
  squared_distance <- outer(locations[, 1], coords[, 1], FUN = "-")^2 +
    outer(locations[, 2], coords[, 2], FUN = "-")^2
  lapply(spec$density_sigmas, function(sigma) {
    value <- as.numeric(exp(-squared_distance / (2 * sigma^2)) %*% mass)
    value / sum(value)
  })
}

baseline_cosine <- function(first, second) {
  denominator <- sqrt(sum(first^2) * sum(second^2))
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sum(first * second) / denominator
}

baseline_density_features <- function(reference, source, spec) {
  reference_density <- baseline_density_signature(reference, spec)
  source_density <- baseline_density_signature(source, spec)
  value <- mapply(baseline_cosine, reference_density, source_density)
  stats::setNames(
    as.numeric(value),
    paste0("density_sigma_", format(spec$density_sigmas, trim = TRUE))
  )
}

# Elastic consensus reference implementation -------------------------------

baseline_weighted_rigid_map <- function(control, target, weight, query) {
  weight <- pmax(as.numeric(weight), 0)
  if (sum(weight) <= 0) weight[] <- 1
  weight <- weight / sum(weight)
  control_center <- colSums(control * weight)
  target_center <- colSums(target * weight)
  centered_control <- sweep(control, 2L, control_center, FUN = "-")
  centered_target <- sweep(target, 2L, target_center, FUN = "-")
  covariance <- t(centered_control * weight) %*% centered_target
  decomposition <- svd(covariance)
  rotation <- decomposition$u %*% t(decomposition$v)
  if (det(rotation) < 0) {
    decomposition$u[, ncol(decomposition$u)] <-
      -decomposition$u[, ncol(decomposition$u)]
    rotation <- decomposition$u %*% t(decomposition$v)
  }
  as.numeric((query - control_center) %*% rotation + target_center)
}

baseline_consensus_targets <- function(encoding, relocated, radius) {
  dx <- outer(relocated[, 1], encoding[, 1], FUN = "-")
  dy <- outer(relocated[, 2], encoding[, 2], FUN = "-")
  log_weight <- -(dx^2 + dy^2) / radius^2
  # Equation 4 normalizes only over encoding fixations. Subtracting each
  # row maximum leaves the normalized Gaussian weights unchanged and avoids
  # underflow when the two configurations begin far apart.
  log_weight <- sweep(log_weight, 1L, apply(log_weight, 1L, max), FUN = "-")
  weight <- exp(log_weight)
  (weight %*% encoding) / rowSums(weight)
}

baseline_mls_relocate <- function(recall, targets, rigidity_radius) {
  relocated <- matrix(NA_real_, nrow(recall), 2L)
  for (i in seq_len(nrow(recall))) {
    squared_distance <- rowSums((sweep(recall, 2L, recall[i, ], FUN = "-"))^2)
    log_weight <- -squared_distance / rigidity_radius^2
    weight <- exp(log_weight - max(log_weight))
    relocated[i, ] <- baseline_weighted_rigid_map(
      recall, targets, weight, recall[i, ]
    )
  }
  relocated
}

#' Elastic consensus matching for encoding and recall fixations
#'
#' This is an inspectable reference implementation of the order-free
#' consensus and moving-least-squares relocation described by Wang,
#' Holmqvist, and Alexa (2021). It is a comparator, not a GazeWeave engine.
#' The iterative relocation follows their Gaussian weighting equations. Their
#' method returns relocated recall and matched encoding fixations rather than a
#' pairwise scalar; `score` is therefore an eyesim-specific, duration-weighted
#' soft match used only as a frozen compatibility comparator.
#'
#' @param reference Encoding fixation path.
#' @param source Recall fixation path.
#' @param spec A [gaze_baseline_spec()].
#'
#' @return An `elastic_consensus_alignment` with relocated points, soft spatial
#'   similarity, matching coverage, deformation, and convergence diagnostics.
#'
#' @references Wang, X., Holmqvist, K., & Alexa, M. (2021). The recorded
#'   trajectories of visually guided eye movements can be used to recover the
#'   spatial location of previously attended objects. *Attention, Perception,
#'   & Psychophysics*, 83, 1700-1716.
#' @export
elastic_consensus_align <- function(reference, source, spec) {
  if (!inherits(spec, "gaze_baseline_spec")) {
    stop("spec must be created by gaze_baseline_spec().")
  }
  required <- c("x", "y", "duration")
  if (!all(required %in% names(reference)) || !all(required %in% names(source))) {
    stop("Elastic consensus paths require x, y, and duration columns.")
  }
  encoding <- cbind(reference$x, reference$y)
  recall <- cbind(source$x, source$y)
  if (nrow(encoding) < 1L || nrow(recall) < 1L ||
      any(!is.finite(encoding)) || any(!is.finite(recall))) {
    stop("Elastic consensus paths must contain finite fixations.")
  }
  if (any(!is.finite(reference$duration)) ||
      any(!is.finite(source$duration)) ||
      any(reference$duration < 0) || any(source$duration < 0) ||
      sum(reference$duration) <= 0 || sum(source$duration) <= 0) {
    stop("Elastic consensus durations must be finite, non-negative, and have positive total mass.")
  }
  radius <- spec$elastic$radii
  relocated <- recall
  converged <- FALSE
  change <- Inf
  iteration <- 0L
  for (iteration in seq_len(spec$elastic$maxit)) {
    targets <- baseline_consensus_targets(
      encoding, relocated, radius[["consensus"]]
    )
    proposal <- baseline_mls_relocate(
      recall, targets, radius[["rigidity"]]
    )
    change <- max(sqrt(rowSums((proposal - relocated)^2)))
    relocated <- proposal
    if (change <= spec$elastic$tolerance) {
      converged <- TRUE
      break
    }
  }

  dx <- outer(encoding[, 1], relocated[, 1], FUN = "-")
  dy <- outer(encoding[, 2], relocated[, 2], FUN = "-")
  distance <- sqrt(dx^2 + dy^2)
  encoding_mass <- reference$duration / sum(reference$duration)
  recall_mass <- source$duration / sum(source$duration)
  nearest_encoding <- apply(distance, 1L, min)
  nearest_recall <- apply(distance, 2L, min)
  matching_radius <- radius[["matching"]]
  soft_similarity <- 0.5 * (
    sum(encoding_mass * exp(-nearest_encoding^2 / matching_radius^2)) +
      sum(recall_mass * exp(-nearest_recall^2 / matching_radius^2))
  )
  global_target <- t(vapply(seq_len(nrow(recall)), function(i) {
    baseline_weighted_rigid_map(recall, relocated, rep(1, nrow(recall)), recall[i, ])
  }, numeric(2)))
  local_deformation <- sqrt(mean(rowSums((relocated - global_target)^2)))

  structure(
    list(
      score = log(pmax(soft_similarity, .Machine$double.xmin)),
      soft_similarity = soft_similarity,
      encoding_coverage = sum(encoding_mass[nearest_encoding <= matching_radius]),
      recall_coverage = sum(recall_mass[nearest_recall <= matching_radius]),
      encoding_matched = nearest_encoding <= matching_radius,
      recall_matched = nearest_recall <= matching_radius,
      local_deformation = local_deformation,
      relocated = relocated,
      original_source = recall,
      reference = encoding,
      convergence = list(
        converged = converged,
        iterations = iteration,
        final_change = change
      ),
      parameters = list(
        consensus_radius = radius[["consensus"]],
        rigidity_radius = radius[["rigidity"]],
        matching_radius = matching_radius,
        unit = spec$screen$unit
      )
    ),
    class = c("elastic_consensus_alignment", "list")
  )
}

# Pair features -------------------------------------------------------------

baseline_register_path <- function(path, warp_model) {
  parameters <- warp_parameters(warp_model)
  coords <- cbind(path$x, path$y) %*% t(parameters$A)
  coords <- sweep(coords, 2L, parameters$translation, FUN = "+")
  registered <- path
  registered$x <- coords[, 1]
  registered$y <- coords[, 2]
  registered
}

baseline_pair_features <- function(reference, source, registered_source,
                                   spec, availability) {
  features <- numeric()
  diagnostics <- list()
  if ("multimatch" %in% spec$methods &&
      availability$available[availability$method == "multimatch"]) {
    raw <- suppressWarnings(multi_match(
      scanpath(reference), scanpath(source),
      screensize = c(spec$screen$width, spec$screen$height)
    ))
    registered <- suppressWarnings(multi_match(
      scanpath(reference), scanpath(registered_source),
      screensize = c(spec$screen$width, spec$screen$height)
    ))
    features <- c(
      features,
      stats::setNames(raw, paste0("raw_", names(raw))),
      stats::setNames(registered, paste0("registered_", names(registered)))
    )
    diagnostics$multimatch_raw <- raw
    diagnostics$multimatch_registered <- registered
  }
  if ("density" %in% spec$methods) {
    raw <- baseline_density_features(reference, source, spec)
    registered <- baseline_density_features(reference, registered_source, spec)
    features <- c(
      features,
      stats::setNames(raw, paste0("raw_", names(raw))),
      stats::setNames(registered, paste0("registered_", names(registered)))
    )
    diagnostics$density_raw <- raw
    diagnostics$density_registered <- registered
  }
  if ("elastic" %in% spec$methods) {
    elastic <- elastic_consensus_align(reference, registered_source, spec)
    elastic_features <- c(
      registered_elastic_score = elastic$score,
      registered_elastic_similarity = elastic$soft_similarity,
      registered_elastic_encoding_coverage = elastic$encoding_coverage,
      registered_elastic_recall_coverage = elastic$recall_coverage,
      registered_elastic_deformation = elastic$local_deformation
    )
    features <- c(features, elastic_features)
    diagnostics$elastic <- elastic
  }
  list(features = features, diagnostics = diagnostics)
}

baseline_group_warp <- function(warp_model, spec, source_row) {
  if (identical(warp_model$type, "none")) return(warp_model)
  group <- if (is.null(spec$warp$fit_by)) {
    names(warp_model$group_models)[[1]]
  } else {
    gaze_key(source_row, spec$warp$fit_by, "warp fit_by")
  }
  subset_gaze_warp_model(warp_model, group)
}

baseline_build_feature_sets <- function(ref_tab, source_tab, match_on,
                                        contrast_on, refvar, sourcevar,
                                        spec, warp_model, availability) {
  ref_key <- gaze_key(ref_tab, match_on, "match_on")
  source_key <- gaze_key(source_tab, match_on, "match_on")
  ref_contrast <- if (is.null(contrast_on)) rep("all", nrow(ref_tab)) else
    gaze_key(ref_tab, contrast_on, "contrast_on")
  source_contrast <- if (is.null(contrast_on)) rep("all", nrow(source_tab)) else
    gaze_key(source_tab, contrast_on, "contrast_on")

  lapply(seq_len(nrow(source_tab)), function(i) {
    candidate_rows <- which(ref_contrast == source_contrast[[i]])
    true_index <- match(source_key[[i]], ref_key[candidate_rows])
    if (is.na(true_index) || length(candidate_rows) < 2L) {
      stop("Every baseline row needs one true candidate and at least one nonmatch.")
    }
    source <- source_tab[[sourcevar]][[i]]
    group_warp <- baseline_group_warp(
      warp_model, spec, source_tab[i, , drop = FALSE]
    )
    registered_source <- baseline_register_path(source, group_warp)
    pairs <- lapply(candidate_rows, function(row) {
      baseline_pair_features(
        ref_tab[[refvar]][[row]], source, registered_source,
        spec, availability
      )
    })
    feature_names <- unique(unlist(lapply(pairs, function(pair) {
      names(pair$features)
    })))
    features <- matrix(
      NA_real_, nrow = length(pairs), ncol = length(feature_names),
      dimnames = list(NULL, feature_names)
    )
    for (candidate in seq_along(pairs)) {
      features[candidate, names(pairs[[candidate]]$features)] <-
        pairs[[candidate]]$features
    }
    list(
      features = features,
      diagnostics = lapply(pairs, `[[`, "diagnostics"),
      candidate_key = ref_key[candidate_rows],
      candidate_rows = candidate_rows,
      true_index = true_index,
      true_key = source_key[[i]],
      contrast_key = source_contrast[[i]],
      source_row = i,
      warp = warp_parameters(group_warp)
    )
  })
}

# Ridge ranking composites --------------------------------------------------

baseline_validate_feature_sets <- function(sets, columns) {
  length(sets) >= 1L && length(columns) >= 1L &&
    all(vapply(sets, function(set) {
      all(columns %in% colnames(set$features)) &&
        all(is.finite(set$features[, columns, drop = FALSE]))
    }, logical(1)))
}

fit_baseline_ridge <- function(sets, columns, lambda) {
  if (!baseline_validate_feature_sets(sets, columns)) {
    stop("Ridge training features are missing or non-finite.")
  }
  stacked <- do.call(rbind, lapply(sets, function(set) {
    set$features[, columns, drop = FALSE]
  }))
  center <- colMeans(stacked)
  scale <- apply(stacked, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= sqrt(.Machine$double.eps)] <- 1
  standardized <- lapply(sets, function(set) {
    sweep(
      sweep(set$features[, columns, drop = FALSE], 2L, center, FUN = "-"),
      2L, scale, FUN = "/"
    )
  })
  objective <- function(beta) {
    losses <- vapply(seq_along(sets), function(i) {
      score <- as.numeric(standardized[[i]] %*% beta)
      gaze_log_sum_exp(score) - score[[sets[[i]]$true_index]]
    }, numeric(1))
    mean(losses) + lambda * sum(beta^2) / 2
  }
  gradient <- function(beta) {
    values <- lapply(seq_along(sets), function(i) {
      design <- standardized[[i]]
      score <- as.numeric(design %*% beta)
      probability <- exp(score - gaze_log_sum_exp(score))
      truth <- numeric(length(score))
      truth[[sets[[i]]$true_index]] <- 1
      as.numeric(crossprod(design, probability - truth))
    })
    Reduce(`+`, values) / length(values) + lambda * beta
  }
  fit <- stats::optim(
    rep(0, length(columns)), objective, gr = gradient,
    method = "BFGS", control = list(reltol = 1e-10, maxit = 500)
  )
  structure(
    list(
      coefficients = stats::setNames(fit$par, columns),
      columns = columns,
      center = center,
      scale = scale,
      lambda = lambda,
      log_loss = fit$value - lambda * sum(fit$par^2) / 2,
      convergence = fit$convergence,
      counts = fit$counts
    ),
    class = c("gaze_baseline_ridge", "list")
  )
}

predict_baseline_ridge <- function(model, set) {
  design <- set$features[, model$columns, drop = FALSE]
  design <- sweep(sweep(design, 2L, model$center, FUN = "-"),
                  2L, model$scale, FUN = "/")
  as.numeric(design %*% model$coefficients)
}

select_baseline_lambda <- function(inner_splits, columns, lambda_grid) {
  if (length(inner_splits) < 2L) {
    stop("At least two valid inner splits are required for ridge selection.")
  }
  fold_loss <- matrix(
    NA_real_, nrow = length(inner_splits), ncol = length(lambda_grid)
  )
  for (fold in seq_along(inner_splits)) {
    split <- inner_splits[[fold]]
    if (!baseline_validate_feature_sets(split$train, columns) ||
        !baseline_validate_feature_sets(split$eval, columns)) {
      stop("Inner-fold ridge features are missing or non-finite.")
    }
    for (index in seq_along(lambda_grid)) {
      model <- fit_baseline_ridge(
        split$train, columns, lambda_grid[[index]]
      )
      fold_loss[fold, index] <- mean(vapply(split$eval, function(set) {
        score <- predict_baseline_ridge(model, set)
        probability <- gaze_candidate_probabilities(score)
        -log(probability[[set$true_index]])
      }, numeric(1)))
    }
  }
  mean_loss <- colMeans(fold_loss)
  minimum <- min(mean_loss)
  tolerance <- 1e-10 * max(1, abs(minimum))
  eligible <- which(mean_loss <= minimum + tolerance)
  selected <- eligible[[which.max(lambda_grid[eligible])]]
  list(
    lambda = lambda_grid[[selected]],
    selected = selected,
    lambda_grid = lambda_grid,
    fold_log_loss = fold_loss,
    mean_log_loss = mean_loss
  )
}

baseline_composite_columns <- function(spec) {
  result <- list()
  if ("multimatch" %in% spec$methods) {
    metric_names <- c(
      "mm_vector", "mm_direction", "mm_length", "mm_position",
      "mm_duration", "mm_position_emd"
    )
    result$multimatch_ridge_registered <- paste0("registered_", metric_names)
    for (metric in metric_names) {
      name <- paste("multimatch", metric, "raw_calibrated", sep = "_")
      result[[name]] <- paste0("raw_", metric)
    }
  }
  if ("density" %in% spec$methods) {
    result$density_ridge_registered <- paste0(
      "registered_density_sigma_",
      format(spec$density_sigmas, trim = TRUE)
    )
    for (sigma in spec$density_sigmas) {
      sigma_name <- format(sigma, trim = TRUE)
      name <- paste0("density_sigma_", sigma_name, "_raw_calibrated")
      result[[name]] <- paste0("raw_density_sigma_", sigma_name)
    }
  }
  if ("elastic" %in% spec$methods) {
    # A one-feature conditional ridge model gives elastic consensus the same
    # training-only calibration contract as the multivariate comparators.
    # The frozen native score remains available as an audit diagnostic.
    result$elastic_ridge_registered <- "registered_elastic_score"
  }
  result
}

baseline_composite_court <- function(method) {
  if (grepl("_raw_calibrated$", method)) {
    "supervised_raw"
  } else {
    "supervised_registered"
  }
}

baseline_frozen_methods <- function(spec, availability) {
  result <- list()
  if ("multimatch" %in% spec$methods) {
    metric_names <- c(
      "mm_vector", "mm_direction", "mm_length", "mm_position",
      "mm_duration", "mm_position_emd"
    )
    mm_available <- availability$available[availability$method == "multimatch"]
    mm_reason <- availability$reason[availability$method == "multimatch"]
    for (registration in c("raw", "registered")) {
      for (metric in metric_names) {
        name <- paste("multimatch", metric, registration, sep = "_")
        result[[name]] <- list(
          columns = paste0(registration, "_", metric),
          aggregate = "identity",
          court = paste0("frozen_", registration),
          available = mm_available,
          reason = mm_reason
        )
      }
    }
  }
  if ("density" %in% spec$methods) {
    for (registration in c("raw", "registered")) {
      name <- paste("density", registration, sep = "_")
      result[[name]] <- list(
        columns = paste0(
          registration, "_density_sigma_",
          format(spec$density_sigmas, trim = TRUE)
        ),
        aggregate = "mean",
        court = paste0("frozen_", registration),
        available = TRUE,
        reason = NA_character_
      )
      for (sigma in spec$density_sigmas) {
        sigma_name <- format(sigma, trim = TRUE)
        scale_name <- paste0(
          "density_sigma_", sigma_name, "_", registration
        )
        result[[scale_name]] <- list(
          columns = paste0(
            registration, "_density_sigma_", sigma_name
          ),
          aggregate = "identity",
          court = paste0("frozen_", registration),
          available = TRUE,
          reason = NA_character_
        )
      }
    }
  }
  if ("elastic" %in% spec$methods) {
    result$elastic_consensus_registered <- list(
      columns = "registered_elastic_score",
      aggregate = "identity",
      court = "frozen_registered",
      available = TRUE,
      reason = NA_character_
    )
  }
  result
}

baseline_method_scores <- function(set, definition, model = NULL) {
  if (!is.null(model)) return(predict_baseline_ridge(model, set))
  value <- set$features[, definition$columns, drop = FALSE]
  if (identical(definition$aggregate, "mean")) rowMeans(value) else value[, 1]
}

# Nested fitting ------------------------------------------------------------

fit_baseline_warp <- function(ref_tab, source_tab, match_on,
                              refvar, sourcevar, spec) {
  fit_gaze_warp_model(
    ref_tab, source_tab, match_on, refvar, sourcevar, spec
  )
}

baseline_make_inner_splits <- function(ref_tab, source_tab, match_on,
                                       contrast_on, split_on,
                                       refvar, sourcevar, spec,
                                       availability, seed) {
  folds <- make_gaze_weave_folds(
    source_tab, split_on, contrast_on, spec$inner_folds, seed
  )
  ref_key <- gaze_key(ref_tab, match_on, "match_on")
  source_key <- gaze_key(source_tab, match_on, "match_on")
  splits <- vector("list", folds$n_folds)
  for (fold in seq_len(folds$n_folds)) {
    eval_rows <- which(folds$fold_id == fold)
    eval_keys <- unique(source_key[eval_rows])
    train_rows <- which(folds$fold_id != fold & !source_key %in% eval_keys)
    train_keys <- unique(source_key[train_rows])
    if (length(train_rows) < 2L || length(eval_rows) < 2L) {
      stop("Inner folds need at least two training and two evaluation rows.")
    }
    ref_train <- ref_tab[ref_key %in% train_keys, , drop = FALSE]
    ref_eval <- ref_tab[ref_key %in% eval_keys, , drop = FALSE]
    source_train <- source_tab[train_rows, , drop = FALSE]
    source_eval <- source_tab[eval_rows, , drop = FALSE]
    warp <- fit_baseline_warp(
      ref_train, source_train, match_on, refvar, sourcevar, spec
    )
    splits[[fold]] <- list(
      train = baseline_build_feature_sets(
        ref_train, source_train, match_on, contrast_on,
        refvar, sourcevar, spec, warp, availability
      ),
      eval = baseline_build_feature_sets(
        ref_eval, source_eval, match_on, contrast_on,
        refvar, sourcevar, spec, warp, availability
      ),
      train_match_keys = train_keys,
      eval_match_keys = eval_keys,
      overlap_match_n = length(intersect(train_keys, eval_keys)),
      warp = warp$info
    )
  }
  splits
}

fit_gaze_baseline_model <- function(ref_tab, source_tab, match_on,
                                    contrast_on, split_on,
                                    refvar, sourcevar, spec, seed) {
  availability <- gaze_baseline_availability(spec)
  warp <- fit_baseline_warp(
    ref_tab, source_tab, match_on, refvar, sourcevar, spec
  )
  training_sets <- baseline_build_feature_sets(
    ref_tab, source_tab, match_on, contrast_on,
    refvar, sourcevar, spec, warp, availability
  )
  inner_error <- NULL
  inner_splits <- tryCatch(
    baseline_make_inner_splits(
      ref_tab, source_tab, match_on, contrast_on, split_on,
      refvar, sourcevar, spec, availability, seed
    ),
    error = function(error) {
      inner_error <<- conditionMessage(error)
      list()
    }
  )
  composites <- list()
  for (method in names(baseline_composite_columns(spec))) {
    columns <- baseline_composite_columns(spec)[[method]]
    available <- baseline_validate_feature_sets(training_sets, columns) &&
      length(inner_splits) >= 2L
    if (!available) {
      reason <- if (!is.null(inner_error)) inner_error else
        "composite features are unavailable or non-finite"
      composites[[method]] <- list(status = "skipped", reason = reason)
      next
    }
    selection <- tryCatch(
      select_baseline_lambda(inner_splits, columns, spec$lambda_grid),
      error = function(error) error
    )
    if (inherits(selection, "error")) {
      composites[[method]] <- list(
        status = "skipped", reason = conditionMessage(selection)
      )
      next
    }
    fitted <- fit_baseline_ridge(
      training_sets, columns, selection$lambda
    )
    composites[[method]] <- list(
      status = "scored",
      reason = NA_character_,
      selection = selection,
      model = fitted
    )
  }
  structure(
    list(
      spec = spec,
      availability = availability,
      warp = warp,
      composites = composites,
      inner_splits = inner_splits,
      inner_error = inner_error,
      training_n = nrow(source_tab)
    ),
    class = c("gaze_baseline_model", "list")
  )
}

baseline_skipped_evidence <- function(set, method, court, reason) {
  list(
    method = method,
    court = court,
    status = "skipped",
    reason = reason,
    calibrated = FALSE,
    compatibility_bits = NA_real_,
    compatibility_odds_bits = NA_real_,
    compatibility_probability_true = NA_real_,
    compatibility_log_loss = NA_real_,
    gaze_info_bits = NA_real_,
    odds_bits = NA_real_,
    posterior_true = NA_real_,
    log_loss = NA_real_,
    brier_score = NA_real_,
    template_rank = NA_real_,
    top1_credit = NA_real_,
    candidate_count = length(set$candidate_key),
    candidates = data.frame(
      candidate_key = set$candidate_key,
      is_true = seq_along(set$candidate_key) == set$true_index,
      status = "skipped",
      reason = reason,
      stringsAsFactors = FALSE
    )
  )
}

baseline_scored_evidence <- function(set, method, court, score,
                                     calibrated, definition = NULL) {
  if (any(!is.finite(score))) {
    return(baseline_skipped_evidence(
      set, method, court, "candidate scores are non-finite"
    ))
  }
  evidence <- score_gaze_candidates(
    score,
    true_index = set$true_index,
    candidate_key = set$candidate_key,
    prior = rep(1 / length(score), length(score)),
    temperature = 1,
    candidate_pool_id = paste0("held-out:", set$contrast_key)
  )
  candidates <- evidence$candidates
  candidates <- cbind(
    candidates,
    as.data.frame(set$features, check.names = FALSE),
    stringsAsFactors = FALSE
  )
  candidates$native_diagnostics <- I(set$diagnostics)
  list(
    method = method,
    court = court,
    status = "scored",
    reason = NA_character_,
    calibrated = calibrated,
    compatibility_bits = evidence$gaze_info_bits,
    compatibility_odds_bits = evidence$odds_bits,
    compatibility_probability_true = evidence$posterior_true,
    compatibility_log_loss = evidence$log_loss,
    gaze_info_bits = if (calibrated) evidence$gaze_info_bits else NA_real_,
    odds_bits = if (calibrated) evidence$odds_bits else NA_real_,
    posterior_true = if (calibrated) evidence$posterior_true else NA_real_,
    log_loss = if (calibrated) evidence$log_loss else NA_real_,
    brier_score = if (calibrated) evidence$brier_score else NA_real_,
    template_rank = evidence$template_rank,
    top1_credit = evidence$top1_credit,
    candidate_count = evidence$candidate_count,
    candidates = candidates,
    definition = definition
  )
}

score_gaze_baseline_sets <- function(sets, model) {
  frozen <- baseline_frozen_methods(model$spec, model$availability)
  lapply(sets, function(set) {
    results <- list()
    for (method in names(frozen)) {
      definition <- frozen[[method]]
      if (!isTRUE(definition$available)) {
        results[[method]] <- baseline_skipped_evidence(
          set, method, definition$court, definition$reason
        )
      } else {
        score <- baseline_method_scores(set, definition)
        results[[method]] <- baseline_scored_evidence(
          set, method, definition$court, score,
          calibrated = FALSE, definition = definition
        )
      }
    }
    for (method in names(model$composites)) {
      composite <- model$composites[[method]]
      court <- baseline_composite_court(method)
      if (!identical(composite$status, "scored")) {
        results[[method]] <- baseline_skipped_evidence(
          set, method, court, composite$reason
        )
      } else {
        score <- predict_baseline_ridge(composite$model, set)
        results[[method]] <- baseline_scored_evidence(
          set, method, court, score,
          calibrated = TRUE,
          definition = list(
            lambda = composite$selection$lambda,
            coefficients = composite$model$coefficients,
            inner_log_loss = composite$selection$mean_log_loss
          )
        )
      }
    }
    results
  })
}

#' Nested cross-validated comparator court for GazeWeave
#'
#' Every method receives the same outer folds, held-out candidate sets,
#' candidate priors, and cross-fitted warp. Registered ridge composites select
#' their penalty using inner folds only. Frozen scores remain explicitly
#' uncalibrated compatibility scores.
#'
#' @inheritParams gaze_weave_cv
#' @param spec A [gaze_baseline_spec()].
#' @param fit_source_filter,eval_source_filter Optional logical vectors or
#'   functions evaluated on `source_tab`. They permit a model to be trained on
#'   positive-control rows and evaluated on both signal and independently
#'   generated null rows without admitting the null rows to preprocessing or
#'   hyperparameter selection.
#' @param inner_split_on Columns defining inner-fold groups. Defaults to
#'   `split_on`; specify a finer training-only unit when an outer fold holds out
#'   an entire session or replication.
#'
#' @return A `gaze_baseline_fit` with one long-format row per held-out trial and
#'   comparator.
#' @export
gaze_baseline_cv <- function(ref_tab, source_tab, match_on,
                             contrast_on = NULL,
                             refvar = "fixgroup", sourcevar = "fixgroup",
                             spec, split_on = match_on, n_folds = NULL,
                             seed = 1, fit_source_filter = NULL,
                             eval_source_filter = NULL,
                             inner_split_on = split_on) {
  if (!inherits(spec, "gaze_baseline_spec")) {
    stop("spec must be created by gaze_baseline_spec().")
  }
  required_ref <- unique(c(match_on, contrast_on, spec$warp$fit_by, refvar))
  required_source <- unique(c(
    match_on, contrast_on, split_on, inner_split_on,
    spec$warp$fit_by, sourcevar
  ))
  if (!all(required_ref %in% names(ref_tab)) ||
      !all(required_source %in% names(source_tab))) {
    stop("Baseline matching, split, warp, and fixation columns must exist.")
  }
  ref_key <- gaze_key(ref_tab, match_on, "match_on")
  if (anyDuplicated(ref_key)) stop("Reference match_on keys must be unique.")
  source_tab <- dplyr::ungroup(source_tab)
  source_tab[["..gaze_row_id"]] <- seq_len(nrow(source_tab))
  source_key <- gaze_key(source_tab, match_on, "match_on")
  if (any(!source_key %in% ref_key)) {
    stop("Every baseline source key must exist in ref_tab.")
  }
  fit_mask <- resolve_gaze_weave_filter(
    source_tab, fit_source_filter, "fit_source_filter"
  )
  eval_mask <- resolve_gaze_weave_filter(
    source_tab, eval_source_filter, "eval_source_filter"
  )
  folds <- make_gaze_weave_folds(
    source_tab, split_on, contrast_on, n_folds, seed
  )
  id_columns <- unique(c(match_on, contrast_on, split_on))
  result_rows <- list()
  fold_info <- vector("list", folds$n_folds)
  result_index <- 1L
  for (fold in seq_len(folds$n_folds)) {
    eval_rows <- which(folds$fold_id == fold & eval_mask)
    if (length(eval_rows) == 0L) next
    eval_source <- source_tab[eval_rows, , drop = FALSE]
    eval_keys <- unique(gaze_key(eval_source, match_on, "match_on"))
    train_rows <- which(
      folds$fold_id != fold & fit_mask & !source_key %in% eval_keys
    )
    train_source <- source_tab[train_rows, , drop = FALSE]
    train_keys <- unique(gaze_key(train_source, match_on, "match_on"))
    ref_train <- ref_tab[ref_key %in% train_keys, , drop = FALSE]
    ref_eval <- ref_tab[ref_key %in% eval_keys, , drop = FALSE]
    model <- fit_gaze_baseline_model(
      ref_train, train_source, match_on, contrast_on, inner_split_on,
      refvar, sourcevar, spec, seed + fold
    )
    eval_sets <- baseline_build_feature_sets(
      ref_eval, eval_source, match_on, contrast_on,
      refvar, sourcevar, spec, model$warp, model$availability
    )
    scored <- score_gaze_baseline_sets(eval_sets, model)
    for (i in seq_along(scored)) {
      retained_columns <- setdiff(
        names(eval_source), c(sourcevar, "..gaze_row_id")
      )
      source_identifiers <- eval_source[i, retained_columns, drop = FALSE]
      for (method in names(scored[[i]])) {
        evidence <- scored[[i]][[method]]
        row <- source_identifiers
        row$.cv_fold <- fold
        row$method <- method
        row$court <- evidence$court
        row$status <- evidence$status
        row$reason <- evidence$reason
        row$calibrated <- evidence$calibrated
        for (field in c(
          "compatibility_bits", "compatibility_odds_bits",
          "compatibility_probability_true", "compatibility_log_loss",
          "gaze_info_bits", "odds_bits",
          "posterior_true", "log_loss", "brier_score", "template_rank",
          "top1_credit", "candidate_count"
        )) {
          row[[field]] <- evidence[[field]]
        }
        row$candidates <- list(evidence$candidates)
        row$definition <- list(
          if (is.null(evidence$definition)) NULL else evidence$definition
        )
        row$source_row_id <- eval_source[["..gaze_row_id"]][[i]]
        result_rows[[result_index]] <- row
        result_index <- result_index + 1L
      }
    }
    fold_info[[fold]] <- list(
      fold = fold,
      train_rows = train_rows,
      eval_rows = eval_rows,
      train_match_keys = train_keys,
      eval_match_keys = eval_keys,
      overlap_match_n = length(intersect(train_keys, eval_keys)),
      inner_folds = lapply(model$inner_splits, function(split) {
        split[c("train_match_keys", "eval_match_keys", "overlap_match_n", "warp")]
      }),
      availability = model$availability,
      composites = lapply(model$composites, function(composite) {
        if (!identical(composite$status, "scored")) return(composite)
        list(
          status = composite$status,
          lambda = composite$selection$lambda,
          inner_log_loss = composite$selection$mean_log_loss,
          coefficients = composite$model$coefficients,
          convergence = composite$model$convergence
        )
      }),
      warp = model$warp$info
    )
  }
  results <- dplyr::bind_rows(result_rows)
  results <- results[order(results$source_row_id, results$court, results$method), , drop = FALSE]
  results$source_row_id <- NULL
  structure(
    list(
      results = tibble::as_tibble(results),
      spec = spec,
      folds = fold_info,
      keys = list(
        match_on = match_on,
        contrast_on = contrast_on,
        split_on = split_on,
        id_columns = id_columns
      ),
      provenance = list(
        engine = "fair_baseline_court",
        outer_folds = folds$n_folds,
        inner_folds = spec$inner_folds,
        candidate_prior = "uniform within held-out candidate set",
        frozen_score_semantics = "uncalibrated compatibility",
        supervised_score_semantics = "ridge conditional-logit probability",
        training_filter = !is.null(fit_source_filter),
        evaluation_filter = !is.null(eval_source_filter),
        seed = seed
      )
    ),
    class = c("gaze_baseline_fit", "list")
  )
}

#' @export
print.gaze_baseline_fit <- function(x, ...) {
  cat("Nested GazeWeave comparator court\n")
  cat("  held-out result rows:", nrow(x$results), "\n")
  cat("  methods:", length(unique(x$results$method)), "\n")
  cat("  scored:", sum(x$results$status == "scored"), "\n")
  cat("  skipped:", sum(x$results$status == "skipped"), "\n")
  invisible(x)
}

#' Tidy GazeWeave comparator results
#'
#' @param x A `gaze_baseline_fit`.
#' @param diagnostics Include calibration, rank, status, and fold diagnostics.
#' @param ... Unused.
#'
#' @return A tibble of held-out comparator results.
#' @export
tidy.gaze_baseline_fit <- function(x, diagnostics = FALSE, ...) {
  columns <- unique(c(
    x$keys$id_columns, "method", "court", "gaze_info_bits"
  ))
  if (isTRUE(diagnostics)) {
    columns <- unique(c(
      columns, "status", "reason", "calibrated", "odds_bits",
      "compatibility_bits", "compatibility_odds_bits",
      "compatibility_probability_true", "compatibility_log_loss",
      "posterior_true", "log_loss",
      "brier_score", "template_rank", "top1_credit", "candidate_count",
      ".cv_fold"
    ))
  }
  x$results[, columns, drop = FALSE]
}
