# Generate the public Transport-v2 estimator-lock fixture used by Transport v3.
#
# From the repository root:
# Rscript inst/validation/gaze-weave-transport-v3-golden.R

load_eyesim_for_v3_golden <- function() {
  if ("package:eyesim" %in% search()) return(invisible(TRUE))
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(eyesim))
  }
  invisible(TRUE)
}

v3_golden_path <- function(coords, duration) {
  stopifnot(nrow(coords) == length(duration), all(duration > 0))
  fixation_group(
    x = coords[, 1],
    y = coords[, 2],
    duration = duration,
    onset = cumsum(c(0, head(duration, -1)))
  )
}

v3_golden_inputs <- function() {
  base_coords <- rbind(
    c(-4.0, -2.0),
    c(-2.7, 0.6),
    c(-0.4, -0.8),
    c(1.2, 1.9),
    c(3.8, 0.3),
    c(5.1, 2.6)
  )
  base_duration <- c(1, 2, 1, 1.5, 0.5, 2)

  cases <- list(
    identical = list(
      reference = v3_golden_path(base_coords, base_duration),
      source = v3_golden_path(base_coords, base_duration)
    ),
    reversed = list(
      reference = v3_golden_path(base_coords, base_duration),
      source = v3_golden_path(
        base_coords[6:1, , drop = FALSE],
        base_duration[6:1]
      )
    ),
    local_swap = list(
      reference = v3_golden_path(base_coords, base_duration),
      source = v3_golden_path(
        base_coords[c(1, 2, 4, 3, 5, 6), , drop = FALSE],
        base_duration[c(1, 2, 4, 3, 5, 6)]
      )
    ),
    offset = list(
      reference = v3_golden_path(base_coords, base_duration),
      source = v3_golden_path(
        sweep(base_coords, 2, c(0.45, -0.30), FUN = "+"),
        base_duration
      )
    ),
    partial_with_intrusions = list(
      reference = v3_golden_path(base_coords, base_duration),
      source = v3_golden_path(
        rbind(base_coords[1:4, , drop = FALSE], c(8.5, -6), c(-8, 6.5)),
        c(base_duration[1:4], 1.25, 1.75)
      )
    ),
    split_adjacent = list(
      reference = v3_golden_path(base_coords, base_duration),
      source = v3_golden_path(
        rbind(base_coords[1, ], base_coords[1, ], base_coords[-1, , drop = FALSE]),
        c(0.4, 0.6, base_duration[-1])
      )
    ),
    two_fixation_reversal = list(
      reference = v3_golden_path(base_coords[1:2, , drop = FALSE], base_duration[1:2]),
      source = v3_golden_path(base_coords[2:1, , drop = FALSE], base_duration[2:1])
    )
  )
  cases
}

v3_golden_spec <- function() {
  gaze_transport_v2_spec(
    spatial = gaze_gaussian_mixture(c(0.75, 1.5), weights = c(0.7, 0.3)),
    chronology = gaze_order_neighbours(neighbours = 2),
    coverage_grid = c(0.5, 0.75, 1),
    coverage_penalty_grid = c(0.25, 0.75, 1.5),
    temporal_weight = 2,
    entropy_schedule = c(0.05, 0.015),
    maxit = 80,
    tolerance = 5e-4,
    projection_maxit = 500,
    projection_tolerance = 1e-8,
    projection_method = "auto",
    multistart = 1
  )
}

v3_golden_input_table <- function(cases) {
  rows <- lapply(names(cases), function(fixture_id) {
    roles <- lapply(c("reference", "source"), function(role) {
      path <- cases[[fixture_id]][[role]]
      data.frame(
        fixture_id = fixture_id,
        role = role,
        fixation_index = seq_len(nrow(path)),
        x = path$x,
        y = path$y,
        onset = path$onset,
        duration = path$duration,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, roles)
  })
  do.call(rbind, rows)
}

v3_golden_barycentric <- function(fit) {
  coupling <- fit$alignment$coupling
  reference_coords <- fit$alignment$reference$coords
  selected <- colSums(coupling)
  result <- matrix(NA_real_, nrow = ncol(coupling), ncol = 2)
  positive <- selected > 0
  result[positive, ] <- t(coupling[, positive, drop = FALSE]) %*%
    reference_coords / selected[positive]
  result
}

generate_transport_v3_golden <- function(output_dir = file.path(
    "inst", "validation", "gaze-weave-transport-v3-golden")) {
  load_eyesim_for_v3_golden()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  old_options <- options(digits = 17, scipen = 999)
  on.exit(options(old_options), add = TRUE)
  set.seed(20260820)

  cases <- v3_golden_inputs()
  spec <- v3_golden_spec()
  fits <- lapply(names(cases), function(fixture_id) {
    pair <- cases[[fixture_id]]
    gaze_transport_v2_align(
      pair$reference,
      pair$source,
      spec,
      candidate_key = fixture_id
    )
  })
  names(fits) <- names(cases)

  inputs <- v3_golden_input_table(cases)
  spec_table <- data.frame(
    key = c(
      "engine_version", "seed", "coordinate_unit", "spatial_scales",
      "spatial_weights", "ordinal_neighbours", "coverage_grid",
      "coverage_penalty_grid", "temporal_weight", "entropy_schedule",
      "maxit", "tolerance", "projection_maxit", "projection_tolerance",
      "projection_method", "multistart"
    ),
    value = c(
      "2", "20260820", "degrees_visual_angle", "0.75;1.5", "0.7;0.3",
      "2", "0.5;0.75;1", "0.25;0.75;1.5", "2", "0.05;0.015",
      "80", "0.0005", "500", "0.00000001", "auto", "1"
    ),
    stringsAsFactors = FALSE
  )
  summary_rows <- lapply(names(fits), function(fixture_id) {
    fit <- fits[[fixture_id]]
    data.frame(
      fixture_id = fixture_id,
      log_score = fit$log_score,
      replay_coverage = fit$diagnostics$replay_coverage,
      map_coverage = fit$diagnostics$map_coverage,
      spatial_rmse = fit$diagnostics$spatial_rmse,
      conditional_spatial = fit$diagnostics$conditional_spatial,
      local_order_error = fit$diagnostics$local_order_error,
      correspondence_information = fit$diagnostics$correspondence_information,
      stationarity = fit$diagnostics$stationarity,
      coupling_change = fit$diagnostics$coupling_change,
      converged = fit$convergence$converged,
      stringsAsFactors = FALSE
    )
  })
  summaries <- do.call(rbind, summary_rows)
  profile_rows <- lapply(names(fits), function(fixture_id) {
    profile <- fits[[fixture_id]]$alignment$profile
    data.frame(
      fixture_id = fixture_id,
      coverage = profile$coverage,
      scientific_energy = profile$scientific_energy,
      regularized_energy = profile$regularized_energy,
      conditional_spatial = profile$conditional_spatial,
      conditional_temporal = profile$conditional_temporal,
      mutual_information = profile$mutual_information,
      stringsAsFactors = FALSE
    )
  })
  profiles <- do.call(rbind, profile_rows)
  coupling_rows <- lapply(names(fits), function(fixture_id) {
    coupling <- fits[[fixture_id]]$alignment$coupling
    grid <- expand.grid(
      reference_index = seq_len(nrow(coupling)),
      source_index = seq_len(ncol(coupling))
    )
    grid$fixture_id <- fixture_id
    grid$coupling <- as.vector(coupling)
    grid[, c("fixture_id", "reference_index", "source_index", "coupling")]
  })
  couplings <- do.call(rbind, coupling_rows)
  barycentric_rows <- lapply(names(fits), function(fixture_id) {
    coords <- v3_golden_barycentric(fits[[fixture_id]])
    data.frame(
      fixture_id = fixture_id,
      source_index = seq_len(nrow(coords)),
      reference_x = coords[, 1],
      reference_y = coords[, 2],
      stringsAsFactors = FALSE
    )
  })
  barycentric <- do.call(rbind, barycentric_rows)

  tables <- list(
    "inputs.csv" = inputs,
    "specification.csv" = spec_table,
    "expected-summary.csv" = summaries,
    "expected-profile.csv" = profiles,
    "expected-coupling.csv" = couplings,
    "expected-barycentric.csv" = barycentric
  )
  for (filename in names(tables)) {
    utils::write.csv(
      tables[[filename]],
      file.path(output_dir, filename),
      row.names = FALSE,
      na = "NA"
    )
  }
  readme <- c(
    "# Transport v2 estimator-lock fixture",
    "",
    "These public deterministic paths and Transport-v2 outputs were frozen",
    "before Transport-v3 algorithm edits. Later backends regenerate the tables",
    "and compare scientific values within 1e-8. The coupling is an optimized",
    "correspondence, not a posterior probability.",
    "",
    "Regenerate from the repository root with:",
    "",
    "```sh",
    "Rscript inst/validation/gaze-weave-transport-v3-golden.R",
    "```"
  )
  writeLines(readme, file.path(output_dir, "README.md"), useBytes = TRUE)

  files <- sort(setdiff(list.files(output_dir), "manifest-md5.csv"))
  hashes <- unname(tools::md5sum(file.path(output_dir, files)))
  manifest <- data.frame(file = files, md5 = hashes, stringsAsFactors = FALSE)
  utils::write.csv(
    manifest,
    file.path(output_dir, "manifest-md5.csv"),
    row.names = FALSE
  )
  invisible(list(
    output_dir = normalizePath(output_dir),
    summary = summaries,
    profile = profiles,
    manifest = manifest
  ))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_arg <- grep("^--output-dir=", args, value = TRUE)
  output_dir <- if (length(output_arg)) {
    sub("^--output-dir=", "", output_arg[[1]])
  } else {
    file.path("inst", "validation", "gaze-weave-transport-v3-golden")
  }
  result <- generate_transport_v3_golden(output_dir)
  print(result$manifest)
}
