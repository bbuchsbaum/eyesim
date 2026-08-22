v3_permutation_correspondence <- function(permutation) {
  n <- length(permutation)
  result <- matrix(0, n, n)
  result[cbind(seq_len(n), permutation)] <- 1 / n
  result
}

v3_ordinal_relation <- function(n, neighbours = 2L) {
  result <- matrix(0, n, n)
  for (index in seq_len(n)) {
    if (index < n) {
      next_index <- seq.int(index + 1L, min(n, index + neighbours))
      result[index, next_index] <- 1
    }
  }
  result
}

v3_edge_residual_slow <- function(correspondence, reference_relation,
                                  source_relation) {
  reference_selected <- rowSums(correspondence)
  source_selected <- colSums(correspondence)
  a <- sum(reference_relation * outer(reference_selected, reference_selected))
  b <- sum(source_relation * outer(source_selected, source_selected))
  agreement <- 0
  for (i in seq_len(nrow(correspondence))) {
    for (ip in seq_len(nrow(correspondence))) {
      for (j in seq_len(ncol(correspondence))) {
        for (jp in seq_len(ncol(correspondence))) {
          agreement <- agreement + correspondence[i, j] * correspondence[ip, jp] *
            reference_relation[i, ip] * source_relation[j, jp]
        }
      }
    }
  }
  if (a + b <= 1e-15) 0 else 1 - 2 * agreement / (a + b)
}

test_that("edge-normalized residual matches the analytic slow form", {
  set.seed(20260820)
  correspondence <- matrix(stats::runif(20), 4, 5)
  correspondence <- correspondence / sum(correspondence)
  reference_relation <- v3_ordinal_relation(4)
  source_relation <- v3_ordinal_relation(5)

  expect_equal(
    eyesim:::transport_v3_edge_residual(
      correspondence, reference_relation, source_relation
    ),
    v3_edge_residual_slow(
      correspondence, reference_relation, source_relation
    ),
    tolerance = 1e-14
  )
})

test_that("edge-normalized analytic gradient matches finite differences", {
  set.seed(20260820)
  correspondence <- matrix(stats::runif(20, 0.2, 1), 4, 5)
  correspondence <- correspondence / sum(correspondence)
  direction <- matrix(stats::rnorm(20), 4, 5)
  direction <- direction - mean(direction)
  direction <- direction / max(abs(direction))
  reference_relation <- v3_ordinal_relation(4)
  source_relation <- v3_ordinal_relation(5)
  gradient <- eyesim:::transport_v3_edge_gradient(
    correspondence, reference_relation, source_relation
  )
  objective <- function(value) {
    eyesim:::transport_v3_edge_residual(
      value, reference_relation, source_relation
    )
  }
  step <- 1e-7
  numerical <- (
    objective(correspondence + step * direction) -
      objective(correspondence - step * direction)
  ) / (2 * step)

  expect_equal(numerical, sum(gradient * direction), tolerance = 1e-7)
})

test_that("complete reversal sensitivity is stable from 4 to 64 fixations", {
  sizes <- c(4L, 8L, 16L, 32L, 64L)
  residuals <- vapply(sizes, function(n) {
    relation <- v3_ordinal_relation(n)
    eyesim:::transport_v3_edge_residual(
      v3_permutation_correspondence(n:1L), relation, relation
    )
  }, numeric(1))

  expect_true(all(residuals >= 0.95))
  expect_lte(diff(range(residuals)), 0.05)
  expect_equal(residuals, rep(1, length(sizes)), tolerance = 1e-14)
})

test_that("identical order is zero and ordinal chronology ignores dilation", {
  n <- 12L
  relation <- v3_ordinal_relation(n)
  correspondence <- diag(1 / n, n)
  expect_lte(
    eyesim:::transport_v3_edge_residual(
      correspondence, relation, relation
    ),
    1e-10
  )

  coords <- cbind(seq_len(n), sin(seq_len(n)))
  original <- make_gaze_fixations(
    coords,
    duration = seq_len(n),
    onset = cumsum(c(0, head(seq_len(n), -1)))
  )
  dilated <- original
  dilated$duration <- original$duration * 17
  dilated$onset <- original$onset * 17
  chronology <- gaze_order_neighbours(neighbours = 2)
  original_measure <- eyesim:::as_gaze_measure(original, chronology)
  dilated_measure <- eyesim:::as_gaze_measure(dilated, chronology)
  expect_equal(original_measure$relation, dilated_measure$relation, tolerance = 0)
  expect_equal(original_measure$mass, dilated_measure$mass, tolerance = 1e-14)
})

test_that("declared order perturbations degrade monotonically", {
  n <- 16L
  relation <- v3_ordinal_relation(n)
  permutations <- list(
    intact = seq_len(n),
    local_swap = c(1:5, 7, 6, 8:16),
    block_reorder = c(1:4, 9:12, 5:8, 13:16),
    complete_shuffle = c(1, 9, 3, 12, 6, 15, 2, 10, 5, 14, 8, 16, 4, 13, 7, 11)
  )
  residuals <- vapply(permutations, function(permutation) {
    eyesim:::transport_v3_edge_residual(
      v3_permutation_correspondence(permutation), relation, relation
    )
  }, numeric(1))

  expect_equal(residuals[["intact"]], 0, tolerance = 1e-14)
  expect_lt(residuals[["local_swap"]], residuals[["block_reorder"]])
  expect_lt(residuals[["block_reorder"]], residuals[["complete_shuffle"]])
})

test_that("short and empty-edge chronology cases are finite and declared", {
  one <- matrix(1, 1, 1)
  no_edge <- matrix(0, 1, 1)
  one_terms <- eyesim:::transport_v3_edge_terms(one, no_edge, no_edge)
  expect_equal(one_terms$residual, 0)
  expect_equal(one_terms$gradient, NULL)

  two_relation <- v3_ordinal_relation(2)
  intact <- eyesim:::transport_v3_edge_residual(
    diag(0.5, 2), two_relation, two_relation
  )
  reversed <- eyesim:::transport_v3_edge_residual(
    v3_permutation_correspondence(2:1), two_relation, two_relation
  )
  expect_true(is.finite(intact))
  expect_true(is.finite(reversed))
  expect_equal(intact, 0, tolerance = 1e-14)
  expect_equal(reversed, 1, tolerance = 1e-14)

  one_sided <- eyesim:::transport_v3_edge_residual(
    matrix(c(0.5, 0.5), 1, 2), no_edge, two_relation
  )
  expect_equal(one_sided, 1, tolerance = 1e-14)
})

test_that("coalescing makes proportional adjacent split invariant", {
  coords <- rbind(c(0, 0), c(2, 1), c(4, -1), c(5, 2))
  original <- make_gaze_fixations(
    coords,
    duration = c(2, 1, 3, 2),
    onset = c(0, 2, 3, 6)
  )
  split <- make_gaze_fixations(
    rbind(coords[1, ], coords[1, ], coords[-1, , drop = FALSE]),
    duration = c(0.5, 1.5, 1, 3, 2),
    onset = c(0, 0.5, 2, 3, 6)
  )
  chronology <- gaze_order_neighbours(neighbours = 2)
  original_measure <- eyesim:::as_gaze_measure(original, chronology)
  split_measure <- eyesim:::as_gaze_measure(split, chronology)

  expect_equal(split_measure$coords, original_measure$coords, tolerance = 0)
  expect_equal(split_measure$mass, original_measure$mass, tolerance = 1e-14)
  expect_equal(split_measure$relation, original_measure$relation, tolerance = 0)
  correspondence <- diag(original_measure$mass)
  expect_equal(
    eyesim:::transport_v3_edge_residual(
      correspondence,
      original_measure$relation,
      original_measure$relation
    ),
    eyesim:::transport_v3_edge_residual(
      correspondence,
      split_measure$relation,
      split_measure$relation
    ),
    tolerance = 1e-14
  )
})

test_that("edge residual validates normalized matrix inputs", {
  relation <- v3_ordinal_relation(2)
  expect_error(
    eyesim:::transport_v3_edge_residual(
      diag(1, 2), relation, relation
    ),
    "unit mass"
  )
  expect_error(
    eyesim:::transport_v3_edge_residual(
      matrix(c(1, -1, 0, 1), 2), relation, relation
    ),
    "non-negative"
  )
})
