# Reproducible directional Replay performance smoke benchmark.
# Run after loading eyesim:
# source(system.file("benchmarks/gaze-weave-replay.R", package = "eyesim"))

benchmark_replay_path <- function(n, offset = 0) {
  time <- seq(0, 3 * pi, length.out = n)
  fixation_group(
    x = cumsum(cos(time)) + offset,
    y = cumsum(sin(time)),
    duration = 0.5 + (seq_len(n) %% 5) / 5,
    onset = cumsum(c(0, head(0.5 + (seq_len(n) %% 5) / 5, -1)))
  )
}

benchmark_gaze_replay <- function(sizes = c(20L, 100L),
                                  grid_size = 64L,
                                  repetitions = 5L) {
  results <- lapply(sizes, function(n) {
    reference <- benchmark_replay_path(n)
    source <- benchmark_replay_path(n, offset = 0.05)
    spec <- gaze_replay_spec(
      grid_size = grid_size,
      max_skip = 2,
      transition_grid = list(
        background = 0.03,
        restart = 0.03,
        advance = 0.35,
        background_stay = 0.9
      )
    )
    model <- fit_gaze_replay_model(
      tibble::tibble(item = 1L, fixgroup = list(reference)),
      tibble::tibble(item = 1L, fixgroup = list(source)),
      match_on = "item",
      spec = spec
    )
    elapsed <- replicate(repetitions, {
      unname(system.time(
        gaze_replay_align(reference, source, model)
      )[["elapsed"]])
    })
    data.frame(
      encoding_fixations = n,
      recall_fixations = n,
      grid_size = grid_size,
      max_skip = spec$max_skip,
      repetitions = repetitions,
      median_seconds = stats::median(elapsed),
      max_seconds = max(elapsed)
    )
  })
  do.call(rbind, results)
}

if (interactive()) {
  print(benchmark_gaze_replay())
}
