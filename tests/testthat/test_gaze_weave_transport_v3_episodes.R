make_v3_episode_set <- function(offset = c(0, 0), invalid = character()) {
  ids <- paste0("P", 1:4)
  episodes <- lapply(seq_along(ids), function(index) {
    coords <- rbind(
      c(0, 0), c(1, 1.5), c(2.5, 0.5), c(3.5, 2)
    ) + matrix(offset + c(index / 20, -index / 30), 4, 2, byrow = TRUE)
    duration <- if (ids[[index]] %in% invalid) rep(0, 4) else c(1, 2, 1, 2)
    make_gaze_fixations(coords, duration = duration, onset = c(0, 1, 3, 4))
  })
  stats::setNames(episodes, ids)
}

v3_mock_aligner <- function(reference, source, spec, candidate_key) {
  difference <- reference$coords - source$coords
  score <- -mean(rowSums(difference^2))
  eyesim:::new_gaze_engine_result(
    engine = "transport",
    candidate_key = candidate_key,
    log_score = score,
    diagnostics = list(
      matched_coverage = 1,
      local_order_residual = 0
    ),
    alignment = list(reference_n = nrow(reference$coords)),
    convergence = list(converged = TRUE),
    provenance = list(
      engine_version = 3L,
      directionality = "symmetric",
      score_semantics = "test_pair_energy",
      duration_semantics = "unit_duration_mass",
      candidate_invariant = TRUE
    )
  )
}

v3_episode_test_spec <- function() {
  list(chronology = gaze_order_neighbours(neighbours = 2))
}

test_that("study episodes are prepared independently without cross edges", {
  episodes <- make_v3_episode_set()
  prepared <- eyesim:::prepare_transport_v3_episodes(
    episodes, gaze_order_neighbours(neighbours = 2)
  )

  expect_s3_class(prepared, "gaze_transport_episodes")
  expect_identical(prepared$valid_ids, paste0("P", 1:4))
  expect_equal(prepared$valid_count, 4L)
  expect_false(prepared$provenance$concatenated)
  expect_false(prepared$provenance$cross_episode_edges)
  expect_true(all(vapply(prepared$episodes, function(episode) {
    nrow(episode$relation) == 4L && sum(episode$mass) == 1
  }, logical(1))))
  expect_true(all(vapply(prepared$episodes, function(episode) {
    all(episode$relation[lower.tri(episode$relation, diag = TRUE)] == 0)
  }, logical(1))))
})

test_that("equal-prior episode mixture is reorder invariant and reduces at one", {
  source <- make_v3_episode_set()[["P1"]]
  spec <- v3_episode_test_spec()
  original <- eyesim:::prepare_transport_v3_episodes(
    make_v3_episode_set(), spec$chronology
  )
  reordered <- eyesim:::prepare_transport_v3_episodes(
    make_v3_episode_set()[c("P4", "P2", "P1", "P3")], spec$chronology
  )
  first <- eyesim:::score_transport_v3_episode_candidate(
    source, original, spec, "item-a",
    common_episode_ids = paste0("P", 1:4),
    aligner = v3_mock_aligner
  )
  second <- eyesim:::score_transport_v3_episode_candidate(
    source, reordered, spec, "item-a",
    common_episode_ids = paste0("P", 1:4),
    aligner = v3_mock_aligner
  )
  single <- eyesim:::score_transport_v3_episode_candidate(
    source, original, spec, "item-a",
    common_episode_ids = "P1",
    aligner = v3_mock_aligner
  )

  expect_equal(first$log_score, second$log_score, tolerance = 1e-14)
  expect_equal(
    single$log_score,
    single$alignment$episodes[["P1"]]$log_score,
    tolerance = 1e-14
  )
  expect_equal(unname(first$diagnostics$equal_episode_weights), rep(0.25, 4))
  expect_identical(first$provenance$episode_weighting, "fixed_equal_prior")
})

test_that("missing episodes use the common valid intersection", {
  source <- make_v3_episode_set()[["P1"]]
  candidates <- list(
    item_a = make_v3_episode_set(),
    item_b = make_v3_episode_set(offset = c(2, 0), invalid = "P2")
  )
  scored <- eyesim:::score_transport_v3_episode_candidates(
    source = source,
    reference_candidates = candidates,
    chronology = gaze_order_neighbours(neighbours = 2),
    spec = v3_episode_test_spec(),
    true_key = "item_a",
    aligner = v3_mock_aligner
  )

  expect_identical(scored$common_episode_ids, c("P1", "P3", "P4"))
  expect_identical(scored$omitted_by_candidate$item_b, "P2")
  expect_equal(
    unname(scored$candidates$item_a$diagnostics$equal_episode_weights),
    rep(1 / 3, 3),
    tolerance = 0
  )
  expect_identical(
    scored$provenance$common_valid_rule,
    "intersection_then_equal_renormalization"
  )
})

test_that("episode candidates return deterministic ranking and diagnostics", {
  source <- make_v3_episode_set()[["P1"]]
  candidates <- list(
    item_a = make_v3_episode_set(),
    item_b = make_v3_episode_set(offset = c(2.5, -1)),
    item_c = make_v3_episode_set(offset = c(-3, 2))
  )
  score_once <- function() {
    eyesim:::score_transport_v3_episode_candidates(
      source = source,
      reference_candidates = candidates,
      chronology = gaze_order_neighbours(neighbours = 2),
      spec = v3_episode_test_spec(),
      true_key = "item_a",
      prior = c(0.5, 0.3, 0.2),
      candidate_pool_id = "episode-test",
      aligner = v3_mock_aligner
    )
  }
  first <- score_once()
  second <- score_once()

  expect_equal(first$evidence$gaze_info_bits, second$evidence$gaze_info_bits)
  expect_gt(first$evidence$gaze_info_bits, 0)
  expect_equal(first$evidence$template_rank, 1)
  expect_identical(names(first$episode_evidence), paste0("P", 1:4))
  expect_true(all(vapply(first$episode_evidence, function(evidence) {
    evidence$template_rank == 1 && is.finite(evidence$gaze_info_bits)
  }, logical(1))))
  expect_true(all(vapply(first$candidates, function(candidate) {
    length(candidate$alignment$episodes) == 4L
  }, logical(1))))
})

test_that("episode benchmark covers one through four presentations", {
  benchmark_path <- gaze_weave_test_inst_path(
    "benchmarks", "gaze-weave-transport-v3.R"
  )
  environment <- new.env(parent = globalenv())
  sys.source(benchmark_path, envir = environment)
  candidates <- list(
    item_a = make_v3_episode_set(),
    item_b = make_v3_episode_set(offset = c(2, 0))
  )
  real_spec <- gaze_transport_spec(
    coverage_nodes = 2,
    entropy_schedule = 0.03,
    maxit = 30,
    tolerance = 1e-3,
    projection_maxit = 300,
    projection_tolerance = 1e-7,
    backend = "reference"
  )
  benchmark <- environment$benchmark_gaze_transport_episodes(
    source = candidates$item_a[["P1"]],
    reference_candidates = candidates,
    chronology = gaze_order_neighbours(neighbours = 2),
    spec = real_spec,
    true_key = "item_a",
    aligner = gaze_transport_align,
    episode_counts = 1:4,
    candidate_counts = 2,
    repetitions = 1
  )

  expect_identical(benchmark$episode_count, 1:4)
  expect_equal(benchmark$pair_evaluations, 2 * (1:4))
  expect_true(all(benchmark$all_converged))
  expect_true(all(is.finite(benchmark$median_seconds)))
  expect_true(all(benchmark$result_bytes > 0))
})
