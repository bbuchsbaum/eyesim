transport_exploratory_source <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    gaze_weave_test_inst_path(
      "validation", "gaze-weave-transport-exploratory.R"
    ),
    envir = environment
  )
  environment
}

transport_exploratory_fixture <- function() {
  grid <- expand.grid(
    saliency_z = c(-2, -1, 0, 1, 2),
    correct_ec = c(-0.5, 0.5),
    probe_type_ec = c(-0.5, 0.5),
    replicate = 1:4,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$participant <- paste0("p", grid$replicate)
  grid$item <- rep(1:20, length.out = nrow(grid))
  grid$z_quality <- rep(c(-1, 1), length.out = nrow(grid))
  grid$z_duration <- rep(c(-1, -1, 1, 1), length.out = nrow(grid))
  grid$trial_key <- paste(grid$participant, grid$item, seq_len(nrow(grid)), sep = ":")
  design <- stats::model.matrix(
    ~ saliency_z * correct_ec * probe_type_ec + z_quality + z_duration,
    data = grid
  )
  beta <- stats::setNames(rep(0, ncol(design)), colnames(design))
  beta[c("(Intercept)", "saliency_z", "correct_ec", "probe_type_ec",
         "saliency_z:correct_ec:probe_type_ec")] <- c(0.1, 0.2, 0.3, -0.1, 0.4)
  grid$gaze_info_bits <- as.numeric(design %*% beta)
  grid
}

test_that("exploratory interaction fit recovers a three-way fixture", {
  court <- transport_exploratory_source()
  fit <- court$transport_exploratory_fit_association(
    transport_exploratory_fixture()
  )
  contrast <- court$transport_exploratory_contrasts(fit$coefficients)

  expect_identical(fit$status, "scored")
  expect_equal(contrast[["saliency_20_to_100"]], 0.8, tolerance = 1e-12)
  expect_equal(contrast[["correct_at_60"]], 0.3, tolerance = 1e-12)
  expect_equal(contrast[["old_minus_lure_at_60"]], -0.1, tolerance = 1e-12)
  expect_equal(contrast[["three_way_per_20"]], 0.4, tolerance = 1e-12)
})

test_that("crossed exploratory bootstrap is row-order invariant", {
  court <- transport_exploratory_source()
  fixture <- transport_exploratory_fixture()
  set.seed(44)
  shuffled <- fixture[sample(seq_len(nrow(fixture))), ]
  first <- court$transport_exploratory_bootstrap_plan(
    fixture, draws = 30L, seed = 73L
  )
  second <- court$transport_exploratory_bootstrap_plan(
    shuffled, draws = 30L, seed = 73L
  )

  expect_equal(first$weights, second$weights)
  expect_true(all(first$weights >= 0))
  expect_true(any(first$weights == 0))
})

test_that("reliability sensitivity has exact boundary behavior", {
  court <- transport_exploratory_source()
  tab <- data.frame(
    effective_fixations = c(2, 5),
    gaze_info_bits = c(NA_real_, NA_real_)
  )
  tab$candidates <- I(list(
    data.frame(base_posterior = c(0.6, 0.4), prior = c(0.5, 0.5),
               is_true = c(TRUE, FALSE)),
    data.frame(base_posterior = c(0.25, 0.75), prior = c(0.5, 0.5),
               is_true = c(TRUE, FALSE))
  ))

  unshrunk <- court$transport_exploratory_kappa_scores(tab, 0)
  shrunken <- court$transport_exploratory_kappa_scores(tab, 1e9)

  expect_equal(unshrunk$gaze_info_bits, log2(c(0.6, 0.25) / 0.5))
  expect_equal(shrunken$gaze_info_bits, c(0, 0), tolerance = 1e-7)
})

test_that("deeper alignment grid is complete and one-factor-at-a-time", {
  court <- transport_exploratory_source()
  grid <- court$transport_exploratory_alignment_grid()
  factors <- data.frame(
    temporal_weight = grid$temporal_weight,
    neighbours = grid$neighbours,
    spatial_scale = grid$spatial_scale,
    coverage = paste(grid$coverage_a, grid$coverage_b, sep = ":")
  )
  default <- factors[grid$spec_id == "default", , drop = FALSE]
  changed <- vapply(seq_len(nrow(grid)), function(index) {
    sum(as.character(factors[index, ]) != as.character(default[1L, ]))
  }, integer(1))

  expect_equal(nrow(grid), 10L)
  expect_equal(anyDuplicated(grid$spec_id), 0L)
  expect_identical(changed, c(0L, rep(1L, 9L)))
})

test_that("local exploratory inputs match frozen checkpoint support", {
  testthat::skip_if_not(
    file.exists(gaze_weave_test_inst_path(
      "validation", "gaze-weave-transport-v3-retrieval-results",
      "checkpoint-manifest.csv"
    )),
    "Local-only frozen Transport checkpoints are not present."
  )
  court <- transport_exploratory_source()
  tab <- court$transport_exploratory_read_data(
    gaze_weave_test_inst_path(
      "validation", "gaze-weave-transport-v3-retrieval-results"
    )
  )

  expect_equal(nrow(tab), 1295L)
  expect_equal(length(unique(tab$participant)), 36L)
  expect_equal(length(unique(tab$item)), 120L)
  expect_equal(nrow(court$transport_exploratory_support(tab)), 20L)
})
