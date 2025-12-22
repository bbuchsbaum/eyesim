make_density_stub <- function(val) {
  structure(
    list(
      z = matrix(val, nrow = 2, ncol = 2),
      x = 1:2,
      y = 1:2,
      sigma = 50
    ),
    class = c("density", "eye_density")
  )
}

make_density_vec <- function(vec) {
  structure(
    list(
      z = matrix(vec, nrow = 2, ncol = 2, byrow = TRUE),
      x = 1:2,
      y = 1:2,
      sigma = 50
    ),
    class = c("density", "eye_density")
  )
}

test_that("latent_pca_transform produces numeric vectors", {
  ref_tab <- tibble::tibble(id = 1:2, density = list(make_density_stub(1), make_density_stub(2)))
  source_tab <- tibble::tibble(id = 1:2, density = list(make_density_stub(1.1), make_density_stub(1.9)))

  res <- latent_pca_transform(ref_tab, source_tab, match_on = "id", comps = 2)
  expect_true(is.numeric(res$ref_tab$density[[1]]))
  expect_length(res$ref_tab$density[[1]], 2)
  expect_length(res$source_tab$density[[1]], 2)
})

test_that("coral_transform adapts source scores without changing shape", {
  ref_tab <- tibble::tibble(id = 1:2, density = list(make_density_stub(1), make_density_stub(2)))
  source_tab <- tibble::tibble(id = 1:2, density = list(make_density_stub(1.1), make_density_stub(1.9)))

  res <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 2, shrink = 1e-2)
  expect_length(res$ref_tab$density[[1]], 2)
  expect_length(res$source_tab$density[[1]], 2)
})

test_that("template_similarity accepts a transform hook", {
  ref_tab <- tibble::tibble(id = 1:2, density = list(make_density_stub(1), make_density_stub(2)))
  source_tab <- tibble::tibble(id = 1:2, density = list(make_density_stub(1.1), make_density_stub(1.9)))

  sim <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = latent_pca_transform,
    similarity_transform_args = list(comps = 2)
  )

  expect_true("eye_sim" %in% names(sim))
  expect_equal(nrow(sim), nrow(source_tab))
})

test_that("latent transforms reject multiscale inputs with differing lengths", {
  ms1 <- structure(list(
    list(z = matrix(1, 2, 2), sigma = 10, x = 1:2, y = 1:2),
    list(z = matrix(1, 3, 3), sigma = 20, x = 1:3, y = 1:3)
  ), class = c("eye_density_multiscale", "list"))

  ms2 <- structure(list(
    list(z = matrix(2, 2, 2), sigma = 10, x = 1:2, y = 1:2),
    list(z = matrix(2, 3, 3), sigma = 20, x = 1:3, y = 1:3)
  ), class = c("eye_density_multiscale", "list"))

  ref_tab <- tibble::tibble(id = 1, density = list(ms1))
  source_tab <- tibble::tibble(id = 1, density = list(ms2))

  expect_error(
    latent_pca_transform(ref_tab, source_tab, match_on = "id"),
    "grid dimensions"
  )
})

test_that("coral_transform improves similarity under feature-wise scaling", {
  set.seed(123)
  n <- 50
  base_vecs <- replicate(n, rnorm(4, 0, 1), simplify = FALSE)
  scale_mat <- diag(c(2, 0.6, 1.7, 0.8))
  source_vecs <- lapply(base_vecs, function(v) as.vector(scale_mat %*% v))

  ref_tab <- tibble::tibble(id = seq_len(n), density = lapply(base_vecs, make_density_vec))
  source_tab <- tibble::tibble(id = seq_len(n), density = lapply(source_vecs, make_density_vec))

  transformed <- coral_transform(
    ref_tab,
    source_tab,
    match_on = "id",
    refvar = "density",
    sourcevar = "density",
    comps = 4,
    shrink = 1e-6
  )

  ref_mat <- do.call(rbind, base_vecs)
  raw_src_mat <- do.call(rbind, source_vecs)
  adapted_mat <- do.call(rbind, transformed$source_tab$density)

  cov_ref <- stats::cov(ref_mat)
  cov_raw <- stats::cov(raw_src_mat)
  cov_adapted <- stats::cov(adapted_mat)

  frob <- function(m) sqrt(sum(m^2))
  raw_gap <- frob(cov_raw - cov_ref)
  adapted_gap <- frob(cov_adapted - cov_ref)

  expect_lt(adapted_gap, raw_gap)
})

test_that("cca_transform recovers a linear mixing between domains", {
  set.seed(456)
  n <- 50
  base_vecs <- replicate(n, runif(4, 0.5, 1.5), simplify = FALSE)
  mix_mat <- matrix(
    c(0.4, 0.6, 0.1, 0.2,
      0.2, 0.5, 0.3, 0.1,
      0.1, 0.2, 0.7, 0.3,
      0.3, 0.2, 0.2, 0.6),
    nrow = 4, byrow = TRUE
  )
  source_vecs <- lapply(base_vecs, function(v) as.vector(mix_mat %*% v))

  ref_tab <- tibble::tibble(id = seq_len(n), density = lapply(base_vecs, make_density_vec))
  source_tab <- tibble::tibble(id = seq_len(n), density = lapply(source_vecs, make_density_vec))

  raw <- template_similarity(ref_tab, source_tab, match_on = "id", permutations = 0, method = "cosine")
  cca_res <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = cca_transform,
    similarity_transform_args = list(comps = 4, shrink = 1e-6)
  )

  raw_mean <- mean(raw$eye_sim, na.rm = TRUE)
  cca_mean <- mean(cca_res$eye_sim, na.rm = TRUE)

  expect_gt(cca_mean, raw_mean)
})
