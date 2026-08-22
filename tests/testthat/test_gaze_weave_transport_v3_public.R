make_transport_v3_public_fit <- function(episode_count = 2L) {
  make_path <- function(offset = c(0, 0), episode = 1L) {
    xy <- rbind(c(0, 0), c(1, 1.2), c(2.2, 0.2), c(3.1, 1.4))
    xy <- xy + matrix(
      offset + c(episode / 50, -episode / 60),
      nrow(xy), 2, byrow = TRUE
    )
    fixation_group(
      x = xy[, 1], y = xy[, 2],
      duration = c(1, 2, 1, 2), onset = c(0, 1, 3, 4)
    )
  }
  ids <- paste0("P", seq_len(episode_count))
  references <- list(
    item_a = stats::setNames(lapply(seq_along(ids), function(index) {
      make_path(c(0, 0), index)
    }), ids),
    item_b = stats::setNames(lapply(seq_along(ids), function(index) {
      make_path(c(3, -2), index)
    }), ids)
  )
  source <- make_path(c(0.08, -0.04), 1L)
  spec <- gaze_transport_spec(
    spatial = gaze_gaussian_mixture(0.8, unit = "deg"),
    coverage_nodes = 2,
    entropy_schedule = 0.03,
    maxit = 40,
    tolerance = 1e-3,
    projection_maxit = 300,
    projection_tolerance = 1e-7,
    backend = "optimized",
    reliability = "none"
  )
  scored <- eyesim:::score_transport_v3_episode_candidates(
    source = source,
    reference_candidates = references,
    chronology = gaze_order_neighbours(neighbours = 2),
    spec = spec,
    true_key = "item_a",
    temperature = 1,
    reliability = 1,
    candidate_pool_id = "public-synthetic",
    aligner = gaze_transport_align
  )
  results <- tibble::tibble(
    participant = "synthetic-1",
    item = "item_a",
    gaze_info_bits = scored$evidence$gaze_info_bits,
    posterior_true = scored$evidence$posterior_true,
    prior_true = scored$evidence$prior_true,
    log_loss = scored$evidence$log_loss,
    template_rank = scored$evidence$template_rank,
    top1_credit = scored$evidence$top1_credit,
    candidate_count = scored$evidence$candidate_count,
    common_episode_count = length(scored$common_episode_ids),
    temperature = 1,
    reliability = 1,
    .cv_fold = 1L,
    all_converged = all(vapply(
      scored$candidates, function(candidate) candidate$convergence$converged,
      logical(1)
    )),
    candidates = list(scored$evidence$candidates),
    alignments = list(scored$candidates)
  )
  structure(
    list(
      results = results,
      spec = spec,
      folds = list(),
      calibration = list(),
      keys = list(
        match_on = c("participant", "item"),
        contrast_on = "participant",
        split_on = c("participant", "item"),
        episode_on = "presentation",
        priorvar = NULL,
        id_columns = c("participant", "item")
      ),
      provenance = list(
        engine = "transport",
        primary_score = "gaze_info_bits=log2(p_true/prior_true)",
        candidate_pool = "public_synthetic",
        candidate_prior = "declared_uniform_design",
        calibration = "test_identity",
        reliability = "none",
        generic_gaze_quality = "separate_response_blind_diagnostic_channel",
        seed = 1L,
        n_folds = 1L
      )
    ),
    class = c("gaze_transport_fit", "list")
  )
}

test_that("Transport defaults to one score and nests diagnostics", {
  fit <- make_transport_v3_public_fit(2L)
  primary <- broom::tidy(fit)
  audited <- broom::tidy(fit, diagnostics = TRUE)
  result <- gaze_transport_result(fit, row = 1L)

  expect_named(primary, c("participant", "item", "gaze_info_bits"))
  expect_true("diagnostics" %in% names(audited))
  expect_named(
    result$diagnostics,
    c(
      "matched_coverage", "spatial_rmse", "spatial_rmse_unit",
      "local_order_preservation", "contraction_scale", "template_rank",
      "candidate_count", "episode_count", "episode_normalization",
      "episodes", "solver_stability"
    )
  )
  expect_identical(result$diagnostics$spatial_rmse_unit, "deg")
  expect_equal(
    result$diagnostics$episodes$equal_episode_weight,
    rep(0.5, 2), tolerance = 0
  )
  expect_true(result$diagnostics$solver_stability$all_converged)
  expect_identical(
    result$provenance$correspondence_semantics,
    "optimized_alignment_not_posterior_probability"
  )
})

test_that("one episode is the explicit identity case", {
  fit <- make_transport_v3_public_fit(1L)
  result <- gaze_transport_result(fit)
  episode <- result$diagnostics$episodes

  expect_identical(result$diagnostics$episode_count, 1L)
  expect_identical(
    result$diagnostics$episode_normalization,
    "fixed_equal_prior_likelihood_mixture"
  )
  expect_equal(episode$equal_episode_weight, 1, tolerance = 0)
})

test_that("Transport plots name alignment and evidence semantics", {
  fit <- make_transport_v3_public_fit(2L)
  alignment <- fit$results$alignments[[1L]]$item_a$alignment$episodes$P1$alignment

  overlay <- ggplot2::autoplot(alignment, type = "overlay")
  braid <- ggplot2::autoplot(alignment, type = "braid")
  raw_registered <- ggplot2::autoplot(alignment, type = "raw_registered")
  diagnostics <- ggplot2::autoplot(alignment, type = "diagnostics")
  evidence <- ggplot2::autoplot(fit, type = "evidence")
  fit_braid <- ggplot2::autoplot(fit, type = "braid", episode = "P2")

  for (plot in list(
    overlay, braid, raw_registered, diagnostics, evidence, fit_braid
  )) {
    expect_s3_class(plot, "ggplot")
  }
  expect_match(overlay$labels$subtitle, "not posterior")
  expect_match(braid$labels$subtitle, "not a posterior")
  expect_match(diagnostics$labels$subtitle, "not a posterior")
  expect_match(evidence$labels$title, "One-score evidence ledger")
  expect_match(evidence$labels$subtitle, "equal-prior episode")
  expect_match(fit_braid$labels$caption, "episode P2 of 2")
  expect_error(
    ggplot2::autoplot(fit, type = "braid", episode = "P3"),
    "one retained study presentation"
  )
})

test_that("private Transport result directories are excluded", {
  old <- setwd(testthat::test_path("..", ".."))
  on.exit(setwd(old), add = TRUE)
  skip_if_not(
    dir.exists(".git") && file.exists(".Rbuildignore"),
    "Source ignore policy is checked only in a source checkout"
  )
  build_ignore <- readLines(".Rbuildignore", warn = FALSE)
  git_status <- system2(
    "git",
    c(
      "check-ignore", "-q",
      "inst/validation/gaze-weave-transport-v3-retrieval-results/.probe"
    ),
    stdout = FALSE, stderr = FALSE
  )

  expect_true(any(grepl(
    "gaze-weave-transport-v3-retrieval-results", build_ignore,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "gaze-weave-transport-v3-repeated-viewing-results", build_ignore,
    fixed = TRUE
  )))
  expect_identical(git_status, 0L)
})
