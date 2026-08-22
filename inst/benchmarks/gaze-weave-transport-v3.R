# Benchmark Transport-v3 scaling across episodes and exhaustive candidates.

benchmark_gaze_transport_episodes <- function(
    source, reference_candidates, chronology, spec, true_key, aligner,
    episode_counts = 1:4, candidate_counts = length(reference_candidates),
    repetitions = 3L) {
  repetitions <- as.integer(repetitions)
  if (length(repetitions) != 1L || is.na(repetitions) || repetitions < 1L) {
    stop("repetitions must be a positive integer.")
  }
  candidate_keys <- names(reference_candidates)
  if (is.null(candidate_keys)) stop("reference_candidates must be named.")
  episode_scorer <- getFromNamespace(
    "score_transport_v3_episode_candidates", "eyesim"
  )
  rows <- list()
  row_index <- 0L
  for (episode_count in as.integer(episode_counts)) {
    for (candidate_count in as.integer(candidate_counts)) {
      if (candidate_count > length(reference_candidates) || candidate_count < 2L) {
        stop("candidate_counts must be between two and the available pool size.")
      }
      selected_keys <- candidate_keys[seq_len(candidate_count)]
      if (!true_key %in% selected_keys) {
        selected_keys[[candidate_count]] <- true_key
      }
      selected <- reference_candidates[selected_keys]
      selected <- lapply(selected, function(episodes) {
        episodes[seq_len(min(episode_count, length(episodes)))]
      })
      elapsed <- numeric(repetitions)
      result <- NULL
      for (iteration in seq_len(repetitions)) {
        timing <- system.time({
          result <- episode_scorer(
            source = source,
            reference_candidates = selected,
            chronology = chronology,
            spec = spec,
            true_key = true_key,
            candidate_pool_id = paste0("benchmark-", candidate_count),
            aligner = aligner
          )
        })
        elapsed[[iteration]] <- unname(timing[["elapsed"]])
      }
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        episode_count = episode_count,
        candidate_count = candidate_count,
        pair_evaluations = episode_count * candidate_count,
        repetitions = repetitions,
        median_seconds = stats::median(elapsed),
        max_seconds = max(elapsed),
        result_bytes = as.numeric(utils::object.size(result)),
        all_converged = all(vapply(
          result$candidates,
          function(candidate) candidate$convergence$converged,
          logical(1)
        )),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

benchmark_gaze_transport_backend <- function(
    repetitions = 9L, fold_candidates = 6L) {
  path <- function(n, offset = c(0, 0), phase = 0) {
    time <- seq(0, 3 * pi, length.out = n) + phase
    duration <- 0.5 + (seq_len(n) %% 5) / 5
    fixation_group(
      x = cumsum(cos(time)) + offset[[1]],
      y = cumsum(sin(time)) + offset[[2]],
      duration = duration,
      onset = cumsum(c(0, head(duration, -1)))
    )
  }
  source <- path(15, c(0.05, 0))
  reference <- path(15)
  make_spec <- function(backend) {
    gaze_transport_spec(backend = backend, multistart = 1)
  }
  reference_spec <- make_spec("reference")
  optimized_spec <- make_spec("optimized")
  invisible(gaze_transport_align(reference, source, optimized_spec))
  pair_reference <- replicate(repetitions, system.time(
    gaze_transport_align(reference, source, reference_spec)
  )[["elapsed"]])
  pair_optimized <- replicate(repetitions, system.time(
    gaze_transport_align(reference, source, optimized_spec)
  )[["elapsed"]])

  references <- lapply(seq_len(fold_candidates), function(index) {
    path(15, c(index / 5, -index / 10), phase = index / 30)
  })
  names(references) <- paste0("item_", seq_len(fold_candidates))
  fold_time <- function(spec) {
    system.time(gaze_transport_align_batch(
      references, source, spec, batch_size = fold_candidates
    ))[["elapsed"]]
  }
  fold_reference <- replicate(3L, fold_time(reference_spec))
  fold_optimized <- replicate(3L, fold_time(optimized_spec))
  data.frame(
    scope = c("pair_15x15", "representative_fold"),
    evaluations = c(1L, fold_candidates),
    reference_median_seconds = c(
      stats::median(pair_reference), stats::median(fold_reference)
    ),
    optimized_median_seconds = c(
      stats::median(pair_optimized), stats::median(fold_optimized)
    ),
    speedup = c(
      stats::median(pair_reference) / stats::median(pair_optimized),
      stats::median(fold_reference) / stats::median(fold_optimized)
    ),
    stringsAsFactors = FALSE
  )
}
