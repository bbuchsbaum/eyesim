transport_v3_court_result_path <- function(file) {
  gaze_weave_test_inst_path(
    "validation", "gaze-weave-transport-v3-court-results", file
  )
}

test_that("the frozen Transport court ledger passes without moved gates", {
  ledger <- utils::read.csv(
    transport_v3_court_result_path("gate-ledger.csv"),
    stringsAsFactors = FALSE
  )
  expected_threshold <- c(
    representation_invariance = 1,
    null_mean_information = 0.02,
    null_top1_prior = 0.02,
    heldout_ece = 0.10,
    matched_false_positive_rate = 0.075,
    finite_converged_rate = 0.99,
    harmless_time_dilation = 1e-8,
    fixed_density_chronology_destruction = 0,
    partial_coherent_intermediate = 1,
    registration_wrong_candidate_guard = 1
  )

  expect_equal(ledger$gate, names(expected_threshold))
  expect_equal(ledger$threshold, unname(expected_threshold), tolerance = 0)
  expect_true(all(ledger$passed))
  expect_true(all(ledger$failure_action == "none"))
})

test_that("the court shares exhaustive support, folds, and priors", {
  receipt <- utils::read.csv(
    transport_v3_court_result_path("fold-receipts.csv"),
    stringsAsFactors = FALSE
  )
  manifest <- utils::read.csv(
    transport_v3_court_result_path("court-manifest.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(receipt$overlap_item_n, c(0L, 0L))
  expect_equal(receipt$v3_inner_overlap_max, c(0L, 0L))
  expect_equal(receipt$baseline_inner_overlap_max, c(0L, 0L))
  expect_equal(receipt$candidate_count, c(8L, 8L))
  expect_equal(receipt$candidate_prior_sum, c(1, 1), tolerance = 1e-14)
  expect_identical(manifest$candidate_policy, "shared_exhaustive")
  expect_identical(manifest$prior_policy, "shared_actual_uniform_design")
  expect_false(manifest$retrieval_fields_read)
  expect_false(manifest$private_data_read)
  expect_true(manifest$passed)
})

test_that("the full comparator panel is present on identical rows", {
  summary <- utils::read.csv(
    transport_v3_court_result_path("predictive-summary.csv"),
    stringsAsFactors = FALSE
  )
  required <- c(
    "transport_v3", "transport_v2", "replay",
    "density_ridge_registered",
    "density_sigma_80_raw_calibrated",
    "density_sigma_160_raw_calibrated",
    "multimatch_ridge_registered",
    "multimatch_mm_vector_raw_calibrated",
    "multimatch_mm_direction_raw_calibrated",
    "multimatch_mm_length_raw_calibrated",
    "multimatch_mm_position_raw_calibrated",
    "multimatch_mm_duration_raw_calibrated",
    "multimatch_mm_position_emd_raw_calibrated"
  )
  expect_true(all(required %in% summary$method))
  counts <- table(summary$method)
  expect_true(all(counts[required] == 2L))
  expect_true(all(summary$n[summary$method %in% required] == 8L))
  expect_true(all(summary$candidate_count[summary$method %in% required] == 8L))
})

test_that("invariance and targeted perturbation relationships are frozen", {
  invariance <- utils::read.csv(
    transport_v3_court_result_path("invariance.csv"),
    stringsAsFactors = FALSE
  )
  perturbation <- utils::read.csv(
    transport_v3_court_result_path("perturbations.csv"),
    stringsAsFactors = FALSE
  )
  factors <- c(
    "contraction", "translation", "spatial_noise",
    "duration_heterogeneity", "global_time_dilation", "local_swap",
    "block_reorder", "reversal", "deletion", "insertion",
    "central_bias", "fixation_count", "candidate_difficulty",
    "coherent_short_subsequence"
  )
  value <- function(factor, column = "transport_v3_info_bits") {
    perturbation[[column]][perturbation$factor == factor]
  }

  expect_true(all(invariance$passed))
  expect_true(all(factors %in% perturbation$factor))
  expect_equal(
    value("global_time_dilation"), value("intact"), tolerance = 1e-8
  )
  expect_lt(value("reversal"), value("intact"))
  expect_gt(value("coherent_short_subsequence"), value("reversal"))
  expect_lt(value("coherent_short_subsequence"), value("intact"))
  expect_equal(
    value("reversal", "density_info_bits"),
    value("intact", "density_info_bits"), tolerance = 1e-12
  )
})

test_that("the court artifact manifest verifies byte for byte", {
  result_dir <- dirname(transport_v3_court_result_path("manifest-md5.csv"))
  manifest <- utils::read.csv(
    file.path(result_dir, "manifest-md5.csv"), stringsAsFactors = FALSE
  )
  paths <- file.path(result_dir, manifest$file)

  expect_true(all(file.exists(paths)))
  expect_equal(
    unname(tools::md5sum(paths)), manifest$md5, tolerance = 0
  )
})
