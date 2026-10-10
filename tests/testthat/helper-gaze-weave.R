make_gaze_fixations <- function(coords,
                                duration = rep(1, nrow(coords)),
                                onset = seq(0, length(duration) - 1)) {
  fixation_group(
    x = coords[, 1],
    y = coords[, 2],
    duration = duration,
    onset = onset
  )
}

# inst/validation (research scripts, reports and frozen results) is kept out
# of the built package, so R CMD check finds it only through
# EYESIM_SOURCE_DIR, a source checkout that CI sets. devtools::test() run in a
# checkout finds it through the source tree.
gaze_weave_test_inst_path <- function(...) {
  source_dir <- Sys.getenv("EYESIM_SOURCE_DIR")
  candidates <- c(
    system.file(..., package = "eyesim"),
    if (nzchar(source_dir)) file.path(source_dir, "inst", ...),
    testthat::test_path("..", "..", "inst", ...)
  )
  found <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(found) > 0L) return(found[[1L]])
  testthat::test_path("..", "..", "inst", ...)
}

# Skip a test file whose checks need inst/validation when it is unavailable,
# as on CRAN.
skip_if_no_validation <- function() {
  testthat::skip_if_not(
    dir.exists(gaze_weave_test_inst_path("validation")),
    "inst/validation is not in the built package (set EYESIM_SOURCE_DIR)"
  )
}

make_small_gaze_replay_spec <- function() {
  gaze_replay_spec(
    grid_size = 8,
    transition_grid = list(
      background = 0.05,
      restart = 0.05,
      advance = 0.3,
      background_stay = 0.9
    )
  )
}

make_gaze_weave_cv_tables <- function() {
  scale_true <- 0.75
  translation_true <- c(0.4, -0.25)
  participants <- c("p1", "p2")
  items <- seq_len(6)

  rows <- expand.grid(
    participant = participants,
    image_id = items,
    stringsAsFactors = FALSE
  )
  reference_paths <- lapply(seq_len(nrow(rows)), function(i) {
    item <- rows$image_id[[i]]
    participant_shift <- if (rows$participant[[i]] == "p1") c(0, 0) else c(0.2, 0.15)
    anchor <- c(item * 2.4, (item %% 3) * 2.7) + participant_shift
    coords <- rbind(
      anchor + c(0, 0),
      anchor + c(0.8, 1.1 + item * 0.03),
      anchor + c(1.7, -0.2),
      anchor + c(2.1, 0.65)
    )
    make_gaze_fixations(coords, duration = c(1, 2, 1, 1), onset = c(0, 1, 3, 4))
  })
  source_paths <- lapply(reference_paths, function(path) {
    ref_coords <- cbind(path$x, path$y)
    source_coords <- sweep(ref_coords, 2, translation_true, FUN = "-") / scale_true
    make_gaze_fixations(
      source_coords,
      duration = path$duration * 1.5,
      onset = path$onset * 1.5
    )
  })

  list(
    ref_tab = tibble::tibble(
      participant = rows$participant,
      image_id = rows$image_id,
      fixgroup = reference_paths
    ),
    source_tab = tibble::tibble(
      row_id = seq_len(nrow(rows)),
      participant = rows$participant,
      image_id = rows$image_id,
      phase = "recall",
      fixgroup = source_paths
    ),
    scale = scale_true,
    translation = translation_true
  )
}
