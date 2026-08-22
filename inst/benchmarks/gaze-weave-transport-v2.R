# Reproducible coverage-conditioned Transport v2 performance benchmark.
# Run after loading eyesim:
# source(system.file("benchmarks/gaze-weave-transport-v2.R", package = "eyesim"))

benchmark_transport_v2_path <- function(n, offset = 0) {
  time <- seq(0, 3 * pi, length.out = n)
  duration <- 0.5 + (seq_len(n) %% 5) / 5
  fixation_group(
    x = cumsum(cos(time)) + offset,
    y = cumsum(sin(time)),
    duration = duration,
    onset = cumsum(c(0, head(duration, -1)))
  )
}

benchmark_gaze_transport_v2 <- function(
    sizes = c(10L, 20L), repetitions = 3L,
    projection_methods = c("standard", "log")) {
  results <- lapply(sizes, function(n) {
    reference <- benchmark_transport_v2_path(n)
    source <- benchmark_transport_v2_path(n, offset = 0.05)
    method_results <- lapply(projection_methods, function(projection_method) {
      spec <- gaze_transport_v2_spec(
        spatial = gaze_gaussian_mixture(c(0.35, 0.8)),
        coverage_grid = c(0.5, 1),
        coverage_penalty_grid = 1,
        entropy_schedule = c(0.03, 0.01),
        maxit = 150,
        tolerance = 1e-4,
        projection_maxit = 300,
        projection_tolerance = 1e-7,
        projection_method = projection_method,
        multistart = 1
      )
      fits <- replicate(repetitions, {
        elapsed <- system.time({
          result <- gaze_transport_v2_align(reference, source, spec)
        })[["elapsed"]]
        c(
          elapsed = unname(elapsed),
          converged = result$convergence$converged,
          log_score = result$log_score
        )
      })
      data.frame(
        reference_fixations = n,
        source_fixations = n,
        projection_method = projection_method,
        coverage_values = length(spec$coverage_grid),
        continuation_values = length(spec$entropy_schedule),
        repetitions = repetitions,
        median_seconds = stats::median(fits["elapsed", ]),
        max_seconds = max(fits["elapsed", ]),
        convergence_rate = mean(fits["converged", ]),
        log_score = stats::median(fits["log_score", ])
      )
    })
    do.call(rbind, method_results)
  })
  result <- do.call(rbind, results)
  reference_time <- stats::setNames(
    result$median_seconds[result$projection_method == "log"],
    result$reference_fixations[result$projection_method == "log"]
  )[as.character(result$reference_fixations)]
  result$speedup_vs_log <- reference_time / result$median_seconds
  result
}

if (interactive()) {
  print(benchmark_gaze_transport_v2())
}
