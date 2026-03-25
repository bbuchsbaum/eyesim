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

make_gaussian_density <- function(mean = c(0, 0), cov = diag(c(0.2, 0.2)),
                                  x = seq(-2, 2, length.out = 25),
                                  y = seq(-2, 2, length.out = 25)) {
  coords <- as.matrix(expand.grid(x = x, y = y))
  inv_cov <- solve(cov)
  centered <- sweep(coords, 2, mean, FUN = "-")
  expo <- rowSums((centered %*% inv_cov) * centered)
  z <- exp(-0.5 * expo)
  z <- matrix(z, nrow = length(x), ncol = length(y))
  z <- z / sum(z)

  structure(
    list(z = z, x = x, y = y, sigma = 1),
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

test_that("latent_pca_projection truncates the PCA basis to requested comps", {
  set.seed(99)
  ref_tab <- tibble::tibble(
    id = 1:5,
    density = replicate(5, make_density_vec(rnorm(4)), simplify = FALSE)
  )
  source_tab <- tibble::tibble(
    id = 1:5,
    density = replicate(5, make_density_vec(rnorm(4)), simplify = FALSE)
  )

  proj <- eyesim:::latent_pca_projection(
    ref_tab,
    source_tab,
    refvar = "density",
    sourcevar = "density",
    comps = 2
  )

  expect_equal(proj$k, 2)
  expect_equal(ncol(proj$basis$x), 2)
  expect_equal(ncol(proj$basis$rotation), 2)
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

  raw_sim <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine"
  )
  coral_sim <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = coral_transform,
    similarity_transform_args = list(comps = 4, shrink = 1e-6)
  )

  expect_lt(adapted_gap, raw_gap)
  expect_gt(mean(coral_sim$eye_sim), mean(raw_sim$eye_sim))
})

test_that("coral_transform is invariant to source row order", {
  set.seed(321)
  n <- 12
  base_vecs <- replicate(n, rnorm(4), simplify = FALSE)
  scale_mat <- diag(c(1.8, 0.7, 1.3, 0.9))
  source_vecs <- lapply(base_vecs, function(v) as.vector(scale_mat %*% v))

  ref_tab <- tibble::tibble(id = seq_len(n), density = lapply(base_vecs, make_density_vec))
  source_tab <- tibble::tibble(
    row_id = seq_len(n),
    id = seq_len(n),
    density = lapply(source_vecs, make_density_vec)
  )
  shuffled_source_tab <- source_tab[sample(seq_len(n)), , drop = FALSE]

  ordered <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 4, shrink = 1e-6)
  shuffled <- coral_transform(ref_tab, shuffled_source_tab, match_on = "id", comps = 4, shrink = 1e-6)

  ordered_mat <- do.call(rbind, ordered$source_tab$density)
  shuffled_mat <- do.call(rbind, shuffled$source_tab$density[order(shuffled$source_tab$row_id)])

  expect_equal(ordered_mat, shuffled_mat, tolerance = 1e-8)
})

test_that("coral_transform remains finite for near-singular covariance", {
  ref_tab <- tibble::tibble(
    id = 1:6,
    density = lapply(
      list(
        c(1, 2, 3, 4),
        c(1, 2, 3, 4.1),
        c(1, 2, 3, 4.2),
        c(2, 3, 4, 5),
        c(2, 3, 4, 5.1),
        c(2, 3, 4, 5.2)
      ),
      make_density_vec
    )
  )
  source_tab <- tibble::tibble(
    id = 1:6,
    density = lapply(rep(list(c(10, 10, 10, 10)), 6), make_density_vec)
  )

  res <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 4, shrink = 1e-3)
  src_mat <- do.call(rbind, res$source_tab$density)

  expect_true(all(is.finite(src_mat)))
  expect_equal(dim(src_mat), c(6, 4))
})

test_that("coral_transform with a single fit_by stratum matches pooled CORAL", {
  set.seed(111)
  n <- 20
  base_vecs <- replicate(n, rnorm(4), simplify = FALSE)
  scale_mat <- diag(c(1.7, 0.9, 1.2, 0.8))
  source_vecs <- lapply(base_vecs, function(v) as.vector(scale_mat %*% v))

  ref_tab <- tibble::tibble(
    id = seq_len(n),
    pid = rep("p1", n),
    density = lapply(base_vecs, make_density_vec)
  )
  source_tab <- tibble::tibble(
    id = seq_len(n),
    pid = rep("p1", n),
    density = lapply(source_vecs, make_density_vec)
  )

  pooled <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 4, shrink = 1e-6)
  grouped <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 4, shrink = 1e-6, fit_by = "pid")

  expect_equal(do.call(rbind, grouped$ref_tab$density), do.call(rbind, pooled$ref_tab$density), tolerance = 1e-10)
  expect_equal(do.call(rbind, grouped$source_tab$density), do.call(rbind, pooled$source_tab$density), tolerance = 1e-8)
})

test_that("contract_transform improves similarity under isotropic contraction and shift", {
  set.seed(901)
  ref_means <- list(c(-0.6, -0.2), c(0.4, -0.1), c(-0.2, 0.5), c(0.7, 0.4))
  ref_cov <- matrix(c(0.18, 0.02, 0.02, 0.12), nrow = 2)
  scale_true <- 0.72
  shift_true <- c(0.18, -0.12)

  ref_tab <- tibble::tibble(
    id = seq_along(ref_means),
    density = lapply(ref_means, function(mu) make_gaussian_density(mean = mu, cov = ref_cov))
  )
  source_tab <- tibble::tibble(
    id = seq_along(ref_means),
    density = lapply(ref_means, function(mu_ref) {
      mu_src <- as.numeric((mu_ref - shift_true) / scale_true)
      cov_src <- ref_cov / (scale_true^2)
      make_gaussian_density(mean = mu_src, cov = cov_src)
    })
  )

  raw <- template_similarity(ref_tab, source_tab, match_on = "id", permutations = 0, method = "cosine")
  contracted <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = contract_transform,
    similarity_transform_args = list(shrink = 1e-6)
  )

  expect_gt(mean(contracted$eye_sim), mean(raw$eye_sim))
})

test_that("contract_transform approximately recovers known scale and translation", {
  set.seed(902)
  ref_means <- replicate(10, runif(2, -0.8, 0.8), simplify = FALSE)
  ref_cov <- matrix(c(0.18, 0.02, 0.02, 0.12), nrow = 2)
  scale_true <- 0.72
  shift_true <- c(0.18, -0.12)

  ref_tab <- tibble::tibble(
    id = seq_along(ref_means),
    density = lapply(ref_means, function(mu) make_gaussian_density(mean = mu, cov = ref_cov))
  )
  source_tab <- tibble::tibble(
    id = seq_along(ref_means),
    density = lapply(ref_means, function(mu_ref) {
      mu_src <- as.numeric((mu_ref - shift_true) / scale_true)
      cov_src <- ref_cov / (scale_true^2)
      make_gaussian_density(mean = mu_src, cov = cov_src)
    })
  )

  model <- eyesim:::fit_similarity_transform_model(
    similarity_transform = contract_transform,
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = "id",
    refvar = "density",
    sourcevar = "density",
    shrink = 1e-6
  )
  fit <- model$group_models$all

  expect_equal(fit$A, diag(scale_true, 2), tolerance = 0.03)
  expect_equal(unname(fit$scale), scale_true, tolerance = 0.03)
  expect_equal(unname(fit$t), shift_true, tolerance = 0.04)
})

test_that("affine_transform improves similarity under affine distortion", {
  ref_means <- list(c(-0.7, -0.1), c(0.5, 0.2), c(-0.1, 0.6), c(0.6, -0.5))
  ref_cov <- matrix(c(0.16, 0.05, 0.05, 0.11), nrow = 2)
  A_true <- matrix(c(0.82, 0.18, -0.12, 1.08), nrow = 2, byrow = TRUE)
  t_true <- c(0.14, -0.09)
  A_inv <- solve(A_true)

  ref_tab <- tibble::tibble(
    id = seq_along(ref_means),
    density = lapply(ref_means, function(mu) make_gaussian_density(mean = mu, cov = ref_cov))
  )
  source_tab <- tibble::tibble(
    id = seq_along(ref_means),
    density = lapply(ref_means, function(mu_ref) {
      mu_src <- as.numeric(A_inv %*% (mu_ref - t_true))
      cov_src <- A_inv %*% ref_cov %*% t(A_inv)
      make_gaussian_density(mean = mu_src, cov = cov_src)
    })
  )

  raw <- template_similarity(ref_tab, source_tab, match_on = "id", permutations = 0, method = "cosine")
  affine <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = affine_transform,
    similarity_transform_args = list(shrink = 1e-6)
  )

  expect_gt(mean(affine$eye_sim), mean(raw$eye_sim))
})

test_that("affine_transform approximately recovers known affine warp", {
  set.seed(903)
  ref_means <- replicate(10, runif(2, -0.8, 0.8), simplify = FALSE)
  ref_cov <- matrix(c(0.16, 0.05, 0.05, 0.11), nrow = 2)
  A_true <- matrix(c(0.82, 0.18, 0.18, 1.08), nrow = 2, byrow = TRUE)
  t_true <- c(0.14, -0.09)
  A_inv <- solve(A_true)

  ref_tab <- tibble::tibble(
    id = seq_along(ref_means),
    density = lapply(ref_means, function(mu) make_gaussian_density(mean = mu, cov = ref_cov))
  )
  source_tab <- tibble::tibble(
    id = seq_along(ref_means),
    density = lapply(ref_means, function(mu_ref) {
      mu_src <- as.numeric(A_inv %*% (mu_ref - t_true))
      cov_src <- A_inv %*% ref_cov %*% t(A_inv)
      make_gaussian_density(mean = mu_src, cov = cov_src)
    })
  )

  model <- eyesim:::fit_similarity_transform_model(
    similarity_transform = affine_transform,
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = "id",
    refvar = "density",
    sourcevar = "density",
    shrink = 1e-6
  )
  fit <- model$group_models$all

  expect_equal(fit$A, A_true, tolerance = 0.05)
  expect_equal(unname(fit$t), t_true, tolerance = 0.05)
})

test_that("coral_transform supports grouped covariance alignment via fit_by", {
  set.seed(222)
  n_per_group <- 30
  base_a <- replicate(n_per_group, rnorm(4), simplify = FALSE)
  base_b <- replicate(n_per_group, rnorm(4), simplify = FALSE)
  scale_a <- diag(c(2.2, 0.7, 1.0, 1.3))
  scale_b <- diag(c(0.8, 1.9, 1.4, 0.6))

  ref_tab <- tibble::tibble(
    id = seq_len(2 * n_per_group),
    pid = rep(c("p1", "p2"), each = n_per_group),
    density = c(lapply(base_a, make_density_vec), lapply(base_b, make_density_vec))
  )
  source_tab <- tibble::tibble(
    id = seq_len(2 * n_per_group),
    pid = rep(c("p1", "p2"), each = n_per_group),
    density = c(
      lapply(base_a, function(v) make_density_vec(as.vector(scale_a %*% v))),
      lapply(base_b, function(v) make_density_vec(as.vector(scale_b %*% v)))
    )
  )

  pooled <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 4, shrink = 1e-6)
  grouped <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 4, shrink = 1e-6, fit_by = "pid")

  cov_gap_sum <- function(transformed_source) {
    sum(vapply(c("p1", "p2"), function(pid) {
      ref_mat <- do.call(rbind, lapply(ref_tab$density[ref_tab$pid == pid], function(x) as.vector(x$z)))
      src_mat <- do.call(rbind, transformed_source$density[transformed_source$pid == pid])
      sqrt(sum((stats::cov(src_mat) - stats::cov(ref_mat))^2))
    }, numeric(1)))
  }

  grouped_notes <- stats::setNames(
    vapply(grouped$info$groups, `[[`, character(1), "note"),
    vapply(grouped$info$groups, `[[`, character(1), "group")
  )

  expect_lt(cov_gap_sum(grouped$source_tab), cov_gap_sum(pooled$source_tab))
  expect_equal(unname(grouped_notes[c("p1", "p2")]), c("ok", "ok"))
})

test_that("grouped coral_transform is invariant to source row order", {
  set.seed(223)
  n_per_group <- 10
  base_a <- replicate(n_per_group, rnorm(4), simplify = FALSE)
  base_b <- replicate(n_per_group, rnorm(4), simplify = FALSE)
  scale_a <- diag(c(2.1, 0.8, 1.1, 1.2))
  scale_b <- diag(c(0.9, 1.8, 1.5, 0.7))

  ref_tab <- tibble::tibble(
    id = seq_len(2 * n_per_group),
    pid = rep(c("p1", "p2"), each = n_per_group),
    density = c(lapply(base_a, make_density_vec), lapply(base_b, make_density_vec))
  )
  source_tab <- tibble::tibble(
    row_id = seq_len(2 * n_per_group),
    id = seq_len(2 * n_per_group),
    pid = rep(c("p1", "p2"), each = n_per_group),
    density = c(
      lapply(base_a, function(v) make_density_vec(as.vector(scale_a %*% v))),
      lapply(base_b, function(v) make_density_vec(as.vector(scale_b %*% v)))
    )
  )
  shuffled_source <- source_tab[sample(seq_len(nrow(source_tab))), , drop = FALSE]

  ordered <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 4, shrink = 1e-6, fit_by = "pid")
  shuffled <- coral_transform(ref_tab, shuffled_source, match_on = "id", comps = 4, shrink = 1e-6, fit_by = "pid")

  ordered_mat <- do.call(rbind, ordered$source_tab$density)
  shuffled_mat <- do.call(rbind, shuffled$source_tab$density[order(shuffled$source_tab$row_id)])

  expect_equal(ordered_mat, shuffled_mat, tolerance = 1e-8)
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

test_that("cca_transform is invariant to source row order when match_on is unchanged", {
  set.seed(789)
  n <- 40
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
  shuffled_source_tab <- source_tab[sample(seq_len(n)), , drop = FALSE]

  ordered <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = cca_transform,
    similarity_transform_args = list(comps = 4, shrink = 1e-6)
  )
  shuffled <- template_similarity(
    ref_tab,
    shuffled_source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = cca_transform,
    similarity_transform_args = list(comps = 4, shrink = 1e-6)
  )

  ordered <- ordered[order(ordered$id), , drop = FALSE]
  shuffled <- shuffled[order(shuffled$id), , drop = FALSE]

  expect_equal(ordered$id, shuffled$id)
  expect_equal(ordered$eye_sim, shuffled$eye_sim, tolerance = 1e-8)
})

test_that("cca_transform strict mode tracks grouped unique matched pairs", {
  ref_tab <- tibble::tibble(
    id = 1:4,
    phase = c("delay", "delay", "scene", "scene"),
    density = lapply(
      list(c(1, 2, 3, 4), c(2, 3, 4, 5), c(5, 4, 3, 2), c(4, 3, 2, 1)),
      make_density_vec
    )
  )
  source_tab <- tibble::tibble(
    id = c(1, 1, 2, 3, 3, 4),
    phase = c("delay", "delay", "delay", "scene", "scene", "scene"),
    density = lapply(
      list(c(1.1, 2.1, 3.1, 4.1), c(1.2, 2.2, 3.2, 4.2), c(2.1, 3.1, 4.1, 5.1),
           c(5.1, 4.1, 3.1, 2.1), c(5.2, 4.2, 3.2, 2.2), c(4.1, 3.1, 2.1, 1.1)),
      make_density_vec
    )
  )

  res <- cca_transform(
    ref_tab,
    source_tab,
    match_on = "id",
    comps = 2,
    fit_by = "phase",
    unique_match_only = TRUE,
    shrink = 1e-6
  )

  group_counts <- stats::setNames(
    vapply(res$info$groups, `[[`, integer(1), "matched_pairs"),
    vapply(res$info$groups, `[[`, character(1), "group")
  )

  expect_equal(res$info$matched_pairs, 4)
  expect_equal(unname(group_counts[c("delay", "scene")]), c(2, 2))
  expect_true(all(lengths(res$ref_tab$density) == 2))
  expect_true(all(lengths(res$source_tab$density) == 2))
})

test_that("cca_transform strict mode records insufficient groups without breaking output shape", {
  ref_tab <- tibble::tibble(
    id = 1:4,
    phase = c("delay", "delay", "scene", "scene"),
    density = lapply(
      list(c(1, 2, 3, 4), c(2, 3, 4, 5), c(5, 4, 3, 2), c(4, 3, 2, 1)),
      make_density_vec
    )
  )
  source_tab <- tibble::tibble(
    id = c(1, 1, 3, 4),
    phase = c("delay", "delay", "scene", "scene"),
    density = lapply(
      list(c(1.1, 2.1, 3.1, 4.1), c(1.2, 2.2, 3.2, 4.2), c(5.1, 4.1, 3.1, 2.1), c(4.1, 3.1, 2.1, 1.1)),
      make_density_vec
    )
  )

  res <- cca_transform(
    ref_tab,
    source_tab,
    match_on = "id",
    comps = 2,
    fit_by = "phase",
    unique_match_only = TRUE,
    shrink = 1e-6
  )

  group_notes <- stats::setNames(
    vapply(res$info$groups, `[[`, character(1), "note"),
    vapply(res$info$groups, `[[`, character(1), "group")
  )

  expect_equal(group_notes[["delay"]], "insufficient matched observations")
  expect_equal(group_notes[["scene"]], "ok")
  expect_true(all(lengths(res$ref_tab$density) == 2))
  expect_true(all(lengths(res$source_tab$density) == 2))
  expect_true(all(vapply(res$source_tab$density, function(x) all(is.finite(x)), logical(1))))
})

test_that("cca_transform strict mode is stable under source row reordering", {
  base_vecs <- list(
    c(1, 2, 3, 4),
    c(2, 3, 4, 5),
    c(5, 4, 3, 2),
    c(4, 3, 2, 1)
  )
  mix_mat <- matrix(
    c(0.5, 0.2, 0.1, 0.2,
      0.2, 0.6, 0.1, 0.1,
      0.1, 0.2, 0.6, 0.1,
      0.2, 0.1, 0.2, 0.5),
    nrow = 4,
    byrow = TRUE
  )

  ref_tab <- tibble::tibble(
    id = 1:4,
    phase = c("delay", "delay", "scene", "scene"),
    density = lapply(base_vecs, make_density_vec)
  )

  source_rows <- tibble::tibble(
    row_id = 1:8,
    id = c(1, 1, 2, 2, 3, 3, 4, 4),
    phase = c("delay", "delay", "delay", "delay", "scene", "scene", "scene", "scene"),
    density = lapply(rep(base_vecs, each = 2), function(v) make_density_vec(as.vector(mix_mat %*% v)))
  )
  shuffled_source <- source_rows[sample(seq_len(nrow(source_rows))), , drop = FALSE]

  ordered <- template_similarity(
    ref_tab,
    source_rows,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = cca_transform,
    similarity_transform_args = list(
      comps = 2,
      fit_by = "phase",
      unique_match_only = TRUE,
      shrink = 1e-6
    )
  )
  shuffled <- template_similarity(
    ref_tab,
    shuffled_source,
    match_on = "id",
    permutations = 0,
    method = "cosine",
    similarity_transform = cca_transform,
    similarity_transform_args = list(
      comps = 2,
      fit_by = "phase",
      unique_match_only = TRUE,
      shrink = 1e-6
    )
  )

  ordered <- ordered[order(ordered$row_id), , drop = FALSE]
  shuffled <- shuffled[order(shuffled$row_id), , drop = FALSE]

  expect_equal(ordered$row_id, shuffled$row_id)
  expect_equal(ordered$eye_sim, shuffled$eye_sim, tolerance = 1e-8)
})
