# Reproducible fair-comparator smoke court for GazeWeave.
#
# Run from a source checkout after loading eyesim:
#   devtools::load_all()
#   source("inst/validation/gaze-weave-baseline-court.R")
#   run_gaze_weave_baseline_court("inst/validation/baseline-results")

baseline_court_path <- function(coords, duration = NULL) {
  if (is.null(duration)) duration <- rep(1, nrow(coords))
  fixation_group(
    x = coords[, 1],
    y = coords[, 2],
    duration = duration,
    onset = c(0, head(cumsum(duration), -1))
  )
}

baseline_court_data <- function(n_participants = 4L, n_items = 6L,
                                seed = 20260815L) {
  set.seed(seed)
  rows <- expand.grid(
    participant = paste0("p", seq_len(n_participants)),
    scenario = c("registered_geometry", "order_at_fixed_density"),
    image_id = seq_len(n_items),
    stringsAsFactors = FALSE
  )
  order_points <- rbind(
    c(20, 20), c(35, 72), c(51, 29),
    c(70, 67), c(82, 18), c(17, 52)
  )
  permutations <- lapply(seq_len(n_items), function(item) {
    c(seq.int(item, n_items), if (item > 1L) seq_len(item - 1L))
  })
  reference <- vector("list", nrow(rows))
  source <- vector("list", nrow(rows))
  for (row in seq_len(nrow(rows))) {
    item <- rows$image_id[[row]]
    participant <- match(rows$participant[[row]], unique(rows$participant))
    if (rows$scenario[[row]] == "order_at_fixed_density") {
      coords <- order_points[permutations[[item]], , drop = FALSE]
    } else {
      angle <- 2 * pi * item / n_items
      anchor <- c(50, 50) + 25 * c(cos(angle), sin(angle))
      motif <- rbind(
        c(-7, -4), c(-2, 7), c(5, 3),
        c(8, -5), c(1, -8), c(-6, 2)
      )
      coords <- sweep(motif, 2, anchor, FUN = "+")
    }
    coords <- coords + participant * c(0.15, -0.1)
    duration <- if (rows$scenario[[row]] == "order_at_fixed_density") {
      rep(1, 6)
    } else {
      c(1, 1.8, 0.8, 1.4, 1.1, 0.9)
    }
    reference[[row]] <- baseline_court_path(coords, duration = duration)
    if (rows$scenario[[row]] == "registered_geometry") {
      source_coords <- sweep(coords, 2, c(-5, 4), FUN = "-") / 1.22
    } else {
      source_coords <- coords
    }
    source_coords <- source_coords + matrix(
      stats::rnorm(length(source_coords), sd = 0.25), ncol = 2
    )
    source[[row]] <- baseline_court_path(source_coords, duration = duration)
  }
  list(
    reference = tibble::tibble(
      participant = rows$participant,
      scenario = rows$scenario,
      image_id = rows$image_id,
      fixgroup = reference
    ),
    source = tibble::tibble(
      participant = rows$participant,
      scenario = rows$scenario,
      image_id = rows$image_id,
      fixgroup = source
    )
  )
}

run_gaze_weave_baseline_court <- function(
    output_dir = NULL,
    n_participants = 4L,
    n_items = 6L,
    seed = 20260815L) {
  data <- baseline_court_data(n_participants, n_items, seed)
  spec <- gaze_baseline_spec(
    screen = gaze_screen(100, 100, unit = "px"),
    density_sigmas = c(3, 6, 12),
    density_grid = 24,
    warp = gaze_warp_contraction(
      center = "screen",
      translation = TRUE,
      fit_by = c("participant", "scenario")
    ),
    lambda_grid = c(0.01, 0.1, 1, 10),
    inner_folds = 2,
    elastic_radii = c(consensus = 5, rigidity = 20, matching = 3),
    elastic_maxit = 20,
    elastic_tolerance = 1e-3
  )
  fit <- gaze_baseline_cv(
    data$reference,
    data$source,
    match_on = c("participant", "scenario", "image_id"),
    contrast_on = c("participant", "scenario"),
    split_on = c("participant", "scenario", "image_id"),
    n_folds = 3,
    seed = seed,
    spec = spec
  )
  scored <- fit$results[fit$results$status == "scored", , drop = FALSE]
  scored$reported_bits <- ifelse(
    scored$calibrated, scored$gaze_info_bits, scored$compatibility_bits
  )
  scored$reported_log_loss <- ifelse(
    scored$calibrated, scored$log_loss, scored$compatibility_log_loss
  )
  summary <- aggregate(
    cbind(
      reported_bits, reported_log_loss, template_rank, top1_credit
    ) ~ scenario + court + method + calibrated,
    data = scored,
    FUN = mean
  )
  fold_audit <- do.call(rbind, lapply(fit$folds, function(fold) {
    data.frame(
      outer_fold = fold$fold,
      outer_overlap = fold$overlap_match_n,
      inner_overlap = max(vapply(
        fold$inner_folds, `[[`, integer(1), "overlap_match_n"
      )),
      multimatch_lambda = fold$composites$multimatch_ridge_registered$lambda,
      density_lambda = fold$composites$density_ridge_registered$lambda
    )
  }))
  result <- list(
    fit = fit,
    summary = summary,
    fold_audit = fold_audit,
    availability = gaze_baseline_availability(spec)
  )
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      fit$results[, setdiff(names(fit$results), c("candidates", "definition"))],
      file.path(output_dir, "baseline-results.csv"), row.names = FALSE
    )
    utils::write.csv(
      summary, file.path(output_dir, "baseline-summary.csv"), row.names = FALSE
    )
    utils::write.csv(
      fold_audit, file.path(output_dir, "fold-audit.csv"), row.names = FALSE
    )
    utils::write.csv(
      result$availability, file.path(output_dir, "availability.csv"),
      row.names = FALSE
    )
    utils::capture.output(
      utils::sessionInfo(), file = file.path(output_dir, "session-info.txt")
    )
  }
  result
}
