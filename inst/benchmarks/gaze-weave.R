# Reproducible GazeWeave performance smoke benchmark.
# Run after installing or loading eyesim:
# source(system.file("benchmarks/gaze-weave.R", package = "eyesim"))

benchmark_gaze_path <- function(n, offset = 0) {
  t <- seq(0, 2 * pi, length.out = n)
  fixation_group(
    x = cumsum(cos(t)) + offset,
    y = cumsum(sin(t)),
    duration = rep(1, n),
    onset = seq_len(n) - 1
  )
}

benchmark_gaze_weave <- function(sizes = c(20L, 50L), repetitions = 3L) {
  spec <- gaze_weave_spec(
    spatial = gaze_gaussian_mixture(c(0.75, 1.5, 3)),
    chronology = gaze_local_order(0.2),
    transport = gaze_partial_transport(
      temporal_weight = 1,
      unmatched = 1,
      entropy = 0.05
    ),
    control = gaze_weave_control(maxit = 200, multistart = 1)
  )

  results <- lapply(sizes, function(n) {
    reference <- benchmark_gaze_path(n)
    source <- benchmark_gaze_path(n, offset = 0.1)
    elapsed <- replicate(repetitions, {
      unname(system.time(gaze_align(reference, source, spec))[["elapsed"]])
    })
    data.frame(
      reference_fixations = n,
      source_fixations = n,
      repetitions = repetitions,
      median_seconds = stats::median(elapsed),
      max_seconds = max(elapsed)
    )
  })
  do.call(rbind, results)
}

if (interactive()) {
  print(benchmark_gaze_weave())
}
