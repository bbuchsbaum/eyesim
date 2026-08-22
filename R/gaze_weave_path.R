# Internal gaze representation -----------------------------------------------

as_gaze_measure <- function(x, chronology) {
  if (!inherits(x, "fixation_group") && !is.data.frame(x)) {
    stop("GazeWeave inputs must be fixation_group objects or data frames.")
  }
  required <- c("x", "y", "onset", "duration")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0L) {
    stop("GazeWeave input is missing: ", paste(missing_cols, collapse = ", "), ".")
  }
  if (nrow(x) < 1L) {
    stop("GazeWeave requires at least one fixation.")
  }

  vals <- x[required]
  valid_numeric <- vapply(vals, is.numeric, logical(1))
  if (!all(valid_numeric)) {
    stop("x, y, onset, and duration must be numeric.")
  }
  if (any(!is.finite(as.matrix(vals)))) {
    stop("GazeWeave inputs must contain only finite values.")
  }
  if (any(x$duration < 0)) {
    stop("Fixation durations must be non-negative.")
  }
  if (sum(x$duration) <= 0) {
    stop("Total fixation duration must be positive.")
  }

  row_order <- order(x$onset, seq_len(nrow(x)))
  sorted <- x[row_order, , drop = FALSE]
  if (nrow(sorted) > 1L && any(diff(sorted$onset) <= 0)) {
    stop("Fixation onsets must be unique; simultaneous onsets are ambiguous.")
  }

  keep <- sorted$duration > 0
  filtered <- sorted[keep, , drop = FALSE]
  kept_index <- row_order[keep]
  if (identical(chronology$type, "order_neighbours")) {
    coalesced <- coalesce_adjacent_gaze_fixations(
      filtered, chronology$coalesce_distance, kept_index
    )
    filtered <- coalesced$fixations
    kept_index <- coalesced$kept_index
  }
  mass <- filtered$duration / sum(filtered$duration)

  trial_start <- min(filtered$onset)
  trial_end <- max(filtered$onset + filtered$duration)
  trial_span <- trial_end - trial_start
  if (!is.finite(trial_span) || trial_span <= 0) {
    stop("The fixation sequence must span positive time.")
  }
  if (identical(chronology$type, "order_neighbours")) {
    normalized_time <- (seq_len(nrow(filtered)) - 0.5) / nrow(filtered)
    relation <- gaze_order_relation(nrow(filtered), chronology$neighbours)
  } else {
    midpoint <- filtered$onset + filtered$duration / 2
    normalized_time <- (midpoint - trial_start) / trial_span
    relation <- gaze_local_relation(normalized_time, chronology$horizon)
  }

  structure(
    list(
      coords = cbind(x = filtered$x, y = filtered$y),
      mass = as.numeric(mass),
      time = as.numeric(normalized_time),
      relation = relation,
      fixations = filtered,
      original = x,
      kept_index = kept_index,
      trial = c(start = trial_start, end = trial_end, span = trial_span)
    ),
    class = c("gaze_measure", "list")
  )
}

coalesce_adjacent_gaze_fixations <- function(fixations, distance, kept_index) {
  if (nrow(fixations) <= 1L) {
    return(list(fixations = fixations, kept_index = kept_index))
  }
  rows <- list(fixations[1, , drop = FALSE])
  indices <- list(kept_index[[1]])
  for (i in 2:nrow(fixations)) {
    previous <- rows[[length(rows)]]
    spatial_distance <- sqrt(
      (fixations$x[[i]] - previous$x[[1]])^2 +
        (fixations$y[[i]] - previous$y[[1]])^2
    )
    if (spatial_distance <= distance) {
      total_duration <- previous$duration[[1]] + fixations$duration[[i]]
      previous$x[[1]] <- (
        previous$x[[1]] * previous$duration[[1]] +
          fixations$x[[i]] * fixations$duration[[i]]
      ) / total_duration
      previous$y[[1]] <- (
        previous$y[[1]] * previous$duration[[1]] +
          fixations$y[[i]] * fixations$duration[[i]]
      ) / total_duration
      previous$duration[[1]] <- total_duration
      rows[[length(rows)]] <- previous
      indices[[length(indices)]] <- c(indices[[length(indices)]], kept_index[[i]])
    } else {
      rows[[length(rows) + 1L]] <- fixations[i, , drop = FALSE]
      indices[[length(indices) + 1L]] <- kept_index[[i]]
    }
  }
  list(
    fixations = do.call(rbind, rows),
    kept_index = indices
  )
}

gaze_order_relation <- function(n_fixations, neighbours) {
  relation <- matrix(0, n_fixations, n_fixations)
  if (n_fixations <= 1L) return(relation)
  for (i in seq_len(n_fixations - 1L)) {
    successors <- seq.int(i + 1L, min(n_fixations, i + neighbours))
    relation[i, successors] <- 1 / (successors - i)
  }
  relation
}

gaze_local_relation <- function(time, horizon) {
  delta <- outer(time, time, FUN = function(left, right) right - left)
  out <- matrix(0, nrow = length(time), ncol = length(time))
  forward <- delta > 0
  out[forward] <- exp(-delta[forward] / horizon)
  out
}

gaze_measure_moments <- function(measures) {
  if (length(measures) == 0L) {
    stop("Cannot compute gaze moments from an empty set of paths.")
  }
  moments <- lapply(measures, function(measure) {
    coords <- measure$coords
    mass <- measure$mass
    mean_vec <- colSums(coords * mass)
    centered <- sweep(coords, 2, mean_vec, FUN = "-")
    cov_mat <- t(centered * mass) %*% centered
    list(mean = mean_vec, cov = cov_mat)
  })

  means <- do.call(rbind, lapply(moments, `[[`, "mean"))
  grand_mean <- colMeans(means)
  grand_cov <- Reduce(`+`, lapply(moments, function(moment) {
    centered_mean <- moment$mean - grand_mean
    moment$cov + tcrossprod(centered_mean)
  })) / length(moments)

  list(mean = as.numeric(grand_mean), cov = grand_cov)
}

transform_gaze_measure <- function(measure, A, translation) {
  if (!inherits(measure, "gaze_measure")) {
    stop("measure must be a gaze_measure.")
  }
  if (!is.matrix(A) || any(dim(A) != c(2L, 2L)) || any(!is.finite(A))) {
    stop("A must be a finite 2 by 2 matrix.")
  }
  if (!is.numeric(translation) || length(translation) != 2L ||
      any(!is.finite(translation))) {
    stop("translation must be a finite vector of length two.")
  }

  transformed <- measure
  transformed$raw_coords <- measure$coords
  transformed$coords <- sweep(measure$coords %*% t(A), 2, translation, FUN = "+")
  transformed$fixations$x <- transformed$coords[, 1]
  transformed$fixations$y <- transformed$coords[, 2]
  transformed
}
