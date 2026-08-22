recognition_sensitivity_source <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    gaze_weave_test_inst_path(
      "validation", "gaze-weave-recognition-sensitivity.R"
    ),
    envir = environment
  )
  environment
}

recognition_sensitivity_source_dir <- function() {
  gaze_weave_test_inst_path(
    "validation",
    "gaze-weave-probe-delay-results"
  )
}

recognition_sensitivity_result_path <- function(...) {
  gaze_weave_test_inst_path(
    "validation",
    "gaze-weave-recognition-sensitivity-results", ...
  )
}

recognition_skip_without_source <- function() {
  testthat::skip_if_not(
    file.exists(file.path(
      recognition_sensitivity_source_dir(), "manifest-md5.csv"
    )),
    "Local-only GW-11 results are not present."
  )
}

recognition_fixture <- function() {
  grid <- expand.grid(
    saliency_z = c(-2, -1, 0, 1, 2),
    correct_ec = c(-0.5, 0.5),
    probe_type_ec = c(-0.5, 0.5),
    replicate = 1:4,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$participant <- paste0("p", grid$replicate)
  grid$item <- rep(1:10, length.out = nrow(grid))
  grid$trial_key <- paste(grid$participant, grid$item, seq_len(nrow(grid)), sep = ":")
  design <- stats::model.matrix(
    ~ saliency_z * correct_ec + probe_type_ec, data = grid
  )
  beta <- c(
    `(Intercept)` = 0.1, saliency_z = 0.2, correct_ec = 0.3,
    probe_type_ec = -0.1, `saliency_z:correct_ec` = 0.4
  )
  grid$gaze_info_bits <- as.numeric(design %*% beta[colnames(design)])
  grid
}

test_that("fixed recognition model recovers an analytic linear fixture", {
  court <- recognition_sensitivity_source()
  fixture <- recognition_fixture()
  fit <- court$recognition_fit_fixed(fixture)

  expect_identical(fit$status, "scored")
  expect_equal(
    fit$coefficients,
    c(
      `(Intercept)` = 0.1, saliency_z = 0.2, correct_ec = 0.3,
      probe_type_ec = -0.1, `saliency_z:correct_ec` = 0.4
    ),
    tolerance = 1e-12
  )
  expect_equal(
    court$recognition_contrasts(fit$coefficients),
    c(
      saliency_20_to_100 = 0.8,
      correct_at_60 = 0.3,
      interaction_per_20 = 0.4
    ),
    tolerance = 1e-12
  )
})

test_that("rank-deficient conditional designs fail explicitly", {
  court <- recognition_sensitivity_source()
  fixture <- recognition_fixture()
  fixture$saliency_z <- 0
  fixture$correct_ec <- 0.5
  fixture$probe_type_ec <- 0.5

  expect_identical(
    court$recognition_fit_fixed(fixture)$status, "rank_deficient"
  )
})

test_that("crossed bootstrap is deterministic and row-permutation invariant", {
  court <- recognition_sensitivity_source()
  fixture <- recognition_fixture()
  set.seed(101)
  shuffled <- fixture[sample(seq_len(nrow(fixture))), ]
  plan <- court$recognition_bootstrap_plan(fixture, draws = 50L, seed = 99L)
  shuffled_plan <- court$recognition_bootstrap_plan(
    shuffled, draws = 50L, seed = 99L
  )
  original <- court$recognition_bootstrap_contrasts(fixture, plan)
  permuted <- court$recognition_bootstrap_contrasts(shuffled, shuffled_plan)

  expect_equal(original, permuted, tolerance = 1e-12)
  expect_true(all(plan$weights >= 0))
  expect_true(all(plan$weights == floor(plan$weights)))
  expect_true(any(plan$weights == 0))
  expect_equal(dim(plan$weights), c(nrow(fixture), 50L))
})

test_that("local GW-11 scores satisfy the frozen sensitivity support", {
  recognition_skip_without_source()
  court <- recognition_sensitivity_source()
  scores <- court$recognition_read_scores(
    recognition_sensitivity_source_dir()
  )
  reference <- scores[
    scores$task == "combined" & scores$method == "transport_v2", ]

  expect_equal(nrow(scores), 1120L)
  expect_equal(nrow(reference), 80L)
  expect_equal(length(unique(reference$participant)), 8L)
  expect_equal(length(unique(reference$item)), 61L)
  expect_identical(as.integer(table(reference$accuracy)), c(21L, 59L))
  expect_equal(
    matrix(
      as.integer(table(reference$saliency, reference$accuracy)), nrow = 5
    ),
    matrix(c(4, 6, 3, 4, 4, 9, 15, 11, 8, 16), nrow = 5)
  )
  expect_lt(kappa(court$recognition_design_matrix(reference), exact = TRUE), 30)
})

test_that("crossed lmer audit reports rather than rewrites singular fits", {
  recognition_skip_without_source()
  court <- recognition_sensitivity_source()
  scores <- court$recognition_read_scores(
    recognition_sensitivity_source_dir()
  )
  reference <- scores[
    scores$task == "combined" & scores$method == "transport_v2", ]
  fit <- court$recognition_fit_lmer(reference)

  expect_identical(fit$status, "scored")
  expect_type(fit$singular, "logical")
  expect_type(fit$converged, "logical")
  expect_true(all(is.finite(fit$contrasts)))
})

test_that("local sensitivity artifacts preserve their manifest and verdict", {
  testthat::skip_if_not(
    file.exists(recognition_sensitivity_result_path("manifest-md5.csv")),
    "Local-only sensitivity results are not present."
  )
  manifest <- utils::read.csv(
    recognition_sensitivity_result_path("manifest-md5.csv"),
    stringsAsFactors = FALSE
  )
  paths <- vapply(
    manifest$file, recognition_sensitivity_result_path, character(1)
  )
  verdict <- utils::read.csv(
    recognition_sensitivity_result_path("sensitivity-verdict.csv"),
    stringsAsFactors = FALSE
  )
  effects <- utils::read.csv(
    recognition_sensitivity_result_path("conditional-effects.csv"),
    stringsAsFactors = FALSE
  )

  expect_true(all(file.exists(paths)))
  expect_identical(unname(tools::md5sum(paths)), manifest$md5)
  expect_true(all(c(
    "configuration.csv", "conditional-effects.csv", "model-audit.csv",
    "cell-predictions.csv", "sensitivity-verdict.csv",
    "scientific-results.rds"
  ) %in% manifest$file))
  expect_true(verdict$post_primary_secondary)
  expect_false(verdict$rescues_persistence_gate)
  expect_false(verdict$closes_wang_gate)
  expect_true(all(effects$valid_fraction >= 0.95))
})
