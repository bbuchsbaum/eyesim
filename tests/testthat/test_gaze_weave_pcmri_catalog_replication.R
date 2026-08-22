source(
  gaze_weave_test_inst_path(
    "validation", "gaze-weave-pcmri-catalog-replication.R"
  ),
  local = TRUE
)

test_that("catalog response coding follows the published report", {
  expect_equal(
    pcmri_catalog_response(c(0, 1, 2, 3, 4, NA)),
    c(NA_integer_, 1L, 1L, 0L, 0L, NA_integer_)
  )
})

test_that("catalog Pearson scorer uses exhaustive within-participant nonmatches", {
  reference <- tibble::tibble(
    participant = rep("p1", 3),
    item = 1:3
  )
  source <- tibble::tibble(
    participant = rep("p1", 3),
    item = 1:3
  )
  reference_vectors <- diag(3)
  source_vectors <- diag(3)
  result <- pcmri_catalog_pearson(
    reference, source, reference_vectors, source_vectors
  )
  expect_equal(result$sim_self, rep(1, 3))
  expect_equal(result$sim_self_perm, rep(0, 3))
  expect_equal(result$sim_self_diff, rep(1, 3))
  expect_equal(result$n_perm, rep(2L, 3))
})

test_that("catalog study pooling preserves four paths without transitions", {
  path <- function(x) {
    structure(
      data.frame(
        x = x, y = x, duration = 100,
        onset = seq_along(x) - 1
      ),
      class = c("fixation_group", "data.frame")
    )
  }
  study <- tibble::tibble(
    participant = rep("p1", 4),
    item = rep(7L, 4),
    presentation = 1:4,
    fixgroup = list(path(1:2), path(3:4), path(5:6), path(7:8))
  )
  pooled <- pcmri_catalog_pool_study(study)
  expect_equal(nrow(pooled), 1L)
  expect_equal(pooled$presentations, 4L)
  expect_equal(pooled$nfix, 8L)
  expect_equal(pooled$fixgroup[[1L]]$x, 1:8)
  expect_equal(pooled$fixgroup[[1L]]$onset, 0:7)
})
