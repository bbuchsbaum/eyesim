make_cv_density_vec <- function(vec) {
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

make_cv_density_tables <- function() {
  ref_base <- list(
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
    participant = c("p1", "p1", "p2", "p2"),
    phase = c("scene", "scene", "scene", "scene"),
    density = lapply(ref_base, make_cv_density_vec)
  )

  source_tab <- tibble::tibble(
    row_id = 1:8,
    id = rep(1:4, each = 2),
    participant = rep(c("p1", "p1", "p2", "p2"), each = 2),
    phase = rep(c("scene", "delay"), times = 4),
    density = lapply(rep(ref_base, each = 2), function(v) make_cv_density_vec(as.vector(mix_mat %*% v)))
  )

  list(ref_tab = ref_tab, source_tab = source_tab)
}

make_cv_warp_density <- function(mean = c(0, 0), cov = diag(c(0.2, 0.2)),
                                 x = seq(-2, 2, length.out = 25),
                                 y = seq(-2, 2, length.out = 25)) {
  coords <- as.matrix(expand.grid(x = x, y = y))
  inv_cov <- solve(cov)
  centered <- sweep(coords, 2, mean, FUN = "-")
  expo <- rowSums((centered %*% inv_cov) * centered)
  z <- exp(-0.5 * expo)
  z <- matrix(z, nrow = length(x), ncol = length(y))
  z <- z / sum(z)
  structure(list(z = z, x = x, y = y, sigma = 1), class = c("density", "eye_density"))
}

make_contract_cv_tables <- function(n = 16, seed = 1) {
  set.seed(seed)
  ref_cov <- matrix(c(0.16, 0.03, 0.03, 0.12), nrow = 2)
  scale_true <- 0.74
  shift_true <- c(0.16, -0.11)
  ref_means <- replicate(n, runif(2, -0.8, 0.8), simplify = FALSE)

  ref_tab <- tibble::tibble(
    id = seq_len(n),
    density = lapply(ref_means, function(mu) make_cv_warp_density(mean = mu, cov = ref_cov))
  )
  source_tab <- tibble::tibble(
    row_id = seq_len(n),
    id = seq_len(n),
    density = lapply(ref_means, function(mu_ref) {
      mu_src <- as.numeric((mu_ref - shift_true) / scale_true)
      cov_src <- ref_cov / (scale_true^2)
      make_cv_warp_density(mean = mu_src, cov = cov_src)
    })
  )

  list(ref_tab = ref_tab, source_tab = source_tab)
}

make_affine_cv_tables <- function(n = 16, seed = 2) {
  set.seed(seed)
  ref_cov <- matrix(c(0.16, 0.05, 0.05, 0.11), nrow = 2)
  A_true <- matrix(c(0.82, 0.18, -0.12, 1.08), nrow = 2, byrow = TRUE)
  t_true <- c(0.14, -0.09)
  A_inv <- solve(A_true)
  ref_means <- replicate(n, runif(2, -0.8, 0.8), simplify = FALSE)

  ref_tab <- tibble::tibble(
    id = seq_len(n),
    density = lapply(ref_means, function(mu) make_cv_warp_density(mean = mu, cov = ref_cov))
  )
  source_tab <- tibble::tibble(
    row_id = seq_len(n),
    id = seq_len(n),
    density = lapply(ref_means, function(mu_ref) {
      mu_src <- as.numeric(A_inv %*% (mu_ref - t_true))
      cov_src <- A_inv %*% ref_cov %*% t(A_inv)
      make_cv_warp_density(mean = mu_src, cov = cov_src)
    })
  )

  list(ref_tab = ref_tab, source_tab = source_tab)
}

test_that("template_similarity_cv excludes held-out match keys from transform fitting", {
  tabs <- make_cv_density_tables()

  res <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    method = "cosine",
    similarity_transform = cca_transform,
    similarity_transform_args = list(comps = 2, fit_by = "phase", unique_match_only = TRUE, shrink = 1e-6),
    split_on = "id",
    n_folds = 2,
    permutations = 0
  )

  cv_info <- attr(res, "similarity_cv")
  overlap_counts <- vapply(cv_info$folds, `[[`, integer(1), "overlap_match_n")

  expect_equal(nrow(res), nrow(tabs$source_tab))
  expect_true(all(overlap_counts == 0L))
  expect_true(all(vapply(cv_info$folds, function(x) x$train_source_n > 0L, logical(1))))
  expect_true(all(vapply(cv_info$folds, function(x) x$eval_source_n > 0L, logical(1))))
})

test_that("template_similarity_cv matches direct similarity when no transform is fit", {
  tabs <- make_cv_density_tables()

  direct <- template_similarity(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    method = "cosine",
    permutations = 0
  )
  cv <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    method = "cosine",
    split_on = "id",
    n_folds = 2,
    permutations = 0
  )

  direct <- direct[order(direct$row_id), , drop = FALSE]
  cv <- cv[order(cv$row_id), , drop = FALSE]

  expect_equal(cv$row_id, direct$row_id)
  expect_equal(cv$eye_sim, direct$eye_sim, tolerance = 1e-10)
})

test_that("template_similarity_cv supports fit and evaluation source filters", {
  tabs <- make_cv_density_tables()

  res <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    method = "cosine",
    similarity_transform = cca_transform,
    similarity_transform_args = list(comps = 2, fit_by = "phase", unique_match_only = TRUE, shrink = 1e-6),
    split_on = "id",
    n_folds = 2,
    permutations = 0,
    fit_source_filter = function(tab) tab$phase == "scene",
    eval_source_filter = function(tab) tab$phase == "delay"
  )

  cv_info <- attr(res, "similarity_cv")

  expect_true(all(res$phase == "delay"))
  expect_equal(nrow(res), sum(tabs$source_tab$phase == "delay"))
  expect_true(all(vapply(cv_info$folds, function(x) x$eval_source_n <= 2L, logical(1))))
  expect_true(all(vapply(cv_info$folds, function(x) x$overlap_match_n == 0L, logical(1))))
})

test_that("template_similarity_cv matches manual held-out CCA computation for one fold", {
  tabs <- make_cv_density_tables()
  fold_spec <- eyesim:::make_similarity_cv_folds(tabs$source_tab, split_on = "id", n_folds = 2, seed = 1)
  fold <- 1L
  eval_rows <- which(fold_spec$fold_id == fold & tabs$source_tab$phase == "delay")
  eval_source <- tabs$source_tab[eval_rows, , drop = FALSE]
  eval_keys <- unique(eval_source$id)

  train_candidate_rows <- which(fold_spec$fold_id != fold & tabs$source_tab$phase == "scene")
  train_rows <- train_candidate_rows[!tabs$source_tab$id[train_candidate_rows] %in% eval_keys]
  train_source <- tabs$source_tab[train_rows, , drop = FALSE]
  train_keys <- unique(train_source$id)

  ref_train <- tabs$ref_tab[tabs$ref_tab$id %in% train_keys, , drop = FALSE]
  ref_eval <- tabs$ref_tab[tabs$ref_tab$id %in% eval_keys, , drop = FALSE]

  model <- eyesim:::fit_similarity_transform_model(
    similarity_transform = cca_transform,
    ref_tab = ref_train,
    source_tab = train_source,
    match_on = "id",
    refvar = "density",
    sourcevar = "density",
    comps = 2,
    shrink = 1e-6
  )
  transformed <- eyesim:::apply_similarity_transform_model(model, ref_tab = ref_eval, source_tab = eval_source)
  manual <- eyesim:::run_similarity_analysis(
    ref_tab = transformed$ref_tab,
    source_tab = transformed$source_tab,
    match_on = "id",
    permutations = 0,
    permute_on = NULL,
    method = "cosine",
    refvar = transformed$refvar,
    sourcevar = transformed$sourcevar
  )
  manual$.cv_fold <- fold
  manual <- manual[order(manual$row_id), , drop = FALSE]

  cv <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    method = "cosine",
    similarity_transform = cca_transform,
    similarity_transform_args = list(comps = 2, shrink = 1e-6),
    split_on = "id",
    n_folds = 2,
    permutations = 0,
    fit_source_filter = function(tab) tab$phase == "scene",
    eval_source_filter = function(tab) tab$phase == "delay",
    seed = 1
  )
  cv_fold <- cv[cv$.cv_fold == fold, , drop = FALSE]
  cv_fold <- cv_fold[order(cv_fold$row_id), , drop = FALSE]

  expect_equal(cv_fold$row_id, manual$row_id)
  expect_equal(cv_fold$eye_sim, manual$eye_sim, tolerance = 1e-10)
})

test_that("template_similarity_cv is invariant to source row order", {
  tabs <- make_cv_density_tables()
  shuffled_source <- tabs$source_tab[sample(seq_len(nrow(tabs$source_tab))), , drop = FALSE]

  ordered <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    method = "cosine",
    similarity_transform = coral_transform,
    similarity_transform_args = list(comps = 2, shrink = 1e-6),
    split_on = "id",
    n_folds = 2,
    permutations = 0
  )
  shuffled <- template_similarity_cv(
    tabs$ref_tab,
    shuffled_source,
    match_on = "id",
    method = "cosine",
    similarity_transform = coral_transform,
    similarity_transform_args = list(comps = 2, shrink = 1e-6),
    split_on = "id",
    n_folds = 2,
    permutations = 0
  )

  ordered <- ordered[order(ordered$row_id), , drop = FALSE]
  shuffled <- shuffled[order(shuffled$row_id), , drop = FALSE]

  expect_equal(ordered$row_id, shuffled$row_id)
  expect_equal(ordered$eye_sim, shuffled$eye_sim, tolerance = 1e-8)
})

test_that("template_similarity_cv supports grouped CORAL fits without held-out leakage", {
  tabs <- make_cv_density_tables()

  res <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    method = "cosine",
    similarity_transform = coral_transform,
    similarity_transform_args = list(comps = 2, shrink = 1e-6, fit_by = "participant"),
    split_on = "id",
    n_folds = 4,
    permutations = 0,
    fit_source_filter = function(tab) tab$phase == "scene",
    eval_source_filter = function(tab) tab$phase == "delay",
    seed = 1
  )

  cv_info <- attr(res, "similarity_cv")
  transform_infos <- lapply(cv_info$folds, `[[`, "transform_info")

  expect_true(all(vapply(cv_info$folds, function(x) x$overlap_match_n == 0L, logical(1))))
  expect_true(all(vapply(transform_infos, function(info) identical(info$fit_by, "participant"), logical(1))))
  expect_true(all(vapply(transform_infos, function(info) {
    group_notes <- stats::setNames(
      vapply(info$groups, `[[`, character(1), "note"),
      vapply(info$groups, `[[`, character(1), "group")
    )
    identical(unname(group_notes[c("p1", "p2")]), c("ok", "ok"))
  }, logical(1))))
  expect_true(all(res$phase == "delay"))
})

test_that("template_similarity_cv matches manual held-out grouped CORAL computation", {
  tabs <- make_cv_density_tables()
  fold_spec <- eyesim:::make_similarity_cv_folds(tabs$source_tab, split_on = "id", n_folds = 2, seed = 1)
  fold <- 1L
  eval_rows <- which(fold_spec$fold_id == fold & tabs$source_tab$phase == "delay")
  eval_source <- tabs$source_tab[eval_rows, , drop = FALSE]
  eval_keys <- unique(eval_source$id)

  train_candidate_rows <- which(fold_spec$fold_id != fold & tabs$source_tab$phase == "scene")
  train_rows <- train_candidate_rows[!tabs$source_tab$id[train_candidate_rows] %in% eval_keys]
  train_source <- tabs$source_tab[train_rows, , drop = FALSE]
  train_keys <- unique(train_source$id)

  ref_train <- tabs$ref_tab[tabs$ref_tab$id %in% train_keys, , drop = FALSE]
  ref_eval <- tabs$ref_tab[tabs$ref_tab$id %in% eval_keys, , drop = FALSE]

  model <- eyesim:::fit_similarity_transform_model(
    similarity_transform = coral_transform,
    ref_tab = ref_train,
    source_tab = train_source,
    match_on = "id",
    refvar = "density",
    sourcevar = "density",
    comps = 2,
    shrink = 1e-6,
    fit_by = "participant"
  )
  transformed <- eyesim:::apply_similarity_transform_model(model, ref_tab = ref_eval, source_tab = eval_source)
  manual <- eyesim:::run_similarity_analysis(
    ref_tab = transformed$ref_tab,
    source_tab = transformed$source_tab,
    match_on = "id",
    permutations = 0,
    permute_on = NULL,
    method = "cosine",
    refvar = transformed$refvar,
    sourcevar = transformed$sourcevar
  )
  manual$.cv_fold <- fold
  manual <- manual[order(manual$row_id), , drop = FALSE]

  cv <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    method = "cosine",
    similarity_transform = coral_transform,
    similarity_transform_args = list(comps = 2, shrink = 1e-6, fit_by = "participant"),
    split_on = "id",
    n_folds = 2,
    permutations = 0,
    fit_source_filter = function(tab) tab$phase == "scene",
    eval_source_filter = function(tab) tab$phase == "delay",
    seed = 1
  )
  cv_fold <- cv[cv$.cv_fold == fold, , drop = FALSE]
  cv_fold <- cv_fold[order(cv_fold$row_id), , drop = FALSE]

  expect_equal(cv_fold$row_id, manual$row_id)
  expect_equal(cv_fold$eye_sim, manual$eye_sim, tolerance = 1e-10)
})

test_that("template_similarity_cv matches manual held-out contract transform computation", {
  ref_means <- list(c(-0.5, -0.2), c(0.4, -0.1), c(-0.2, 0.4), c(0.6, 0.3))
  ref_cov <- matrix(c(0.16, 0.03, 0.03, 0.12), nrow = 2)
  scale_true <- 0.74
  shift_true <- c(0.16, -0.11)

  ref_tab <- tibble::tibble(
    id = 1:4,
    participant = c("p1", "p1", "p2", "p2"),
    phase = "scene",
    density = lapply(ref_means, function(mu) make_cv_warp_density(mean = mu, cov = ref_cov))
  )
  source_tab <- tibble::tibble(
    row_id = 1:8,
    id = rep(1:4, each = 2),
    participant = rep(c("p1", "p1", "p2", "p2"), each = 2),
    phase = rep(c("scene", "delay"), times = 4),
    density = lapply(rep(ref_means, each = 2), function(mu_ref) {
      mu_src <- as.numeric((mu_ref - shift_true) / scale_true)
      cov_src <- ref_cov / (scale_true^2)
      make_cv_warp_density(mean = mu_src, cov = cov_src)
    })
  )

  fold_spec <- eyesim:::make_similarity_cv_folds(source_tab, split_on = "id", n_folds = 2, seed = 1)
  fold <- 1L
  eval_rows <- which(fold_spec$fold_id == fold & source_tab$phase == "delay")
  eval_source <- source_tab[eval_rows, , drop = FALSE]
  eval_keys <- unique(eval_source$id)

  train_candidate_rows <- which(fold_spec$fold_id != fold & source_tab$phase == "scene")
  train_rows <- train_candidate_rows[!source_tab$id[train_candidate_rows] %in% eval_keys]
  train_source <- source_tab[train_rows, , drop = FALSE]
  train_keys <- unique(train_source$id)

  ref_train <- ref_tab[ref_tab$id %in% train_keys, , drop = FALSE]
  ref_eval <- ref_tab[ref_tab$id %in% eval_keys, , drop = FALSE]

  model <- eyesim:::fit_similarity_transform_model(
    similarity_transform = contract_transform,
    ref_tab = ref_train,
    source_tab = train_source,
    match_on = "id",
    refvar = "density",
    sourcevar = "density",
    shrink = 1e-6
  )
  transformed <- eyesim:::apply_similarity_transform_model(model, ref_tab = ref_eval, source_tab = eval_source)
  manual <- eyesim:::run_similarity_analysis(
    ref_tab = transformed$ref_tab,
    source_tab = transformed$source_tab,
    match_on = "id",
    permutations = 0,
    permute_on = NULL,
    method = "cosine",
    refvar = transformed$refvar,
    sourcevar = transformed$sourcevar
  )
  manual$.cv_fold <- fold
  manual <- manual[order(manual$row_id), , drop = FALSE]

  cv <- template_similarity_cv(
    ref_tab,
    source_tab,
    match_on = "id",
    method = "cosine",
    similarity_transform = contract_transform,
    similarity_transform_args = list(shrink = 1e-6),
    split_on = "id",
    n_folds = 2,
    permutations = 0,
    fit_source_filter = function(tab) tab$phase == "scene",
    eval_source_filter = function(tab) tab$phase == "delay",
    seed = 1
  )
  cv_fold <- cv[cv$.cv_fold == fold, , drop = FALSE]
  cv_fold <- cv_fold[order(cv_fold$row_id), , drop = FALSE]

  expect_equal(cv_fold$row_id, manual$row_id)
  expect_equal(cv_fold$eye_sim, manual$eye_sim, tolerance = 1e-10)
})

test_that("template_similarity_cv improves held-out similarity under contract distortion", {
  tabs <- make_contract_cv_tables(n = 16, seed = 11)

  raw <- template_similarity(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine"
  )
  cv <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    split_on = "id",
    n_folds = 4,
    permutations = 0,
    method = "cosine",
    similarity_transform = contract_transform,
    similarity_transform_args = list(shrink = 1e-6),
    seed = 1
  )

  expect_gt(mean(cv$eye_sim), mean(raw$eye_sim))
})

test_that("template_similarity_cv improves held-out similarity under affine distortion", {
  tabs <- make_affine_cv_tables(n = 16, seed = 12)

  raw <- template_similarity(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    permutations = 0,
    method = "cosine"
  )
  cv <- template_similarity_cv(
    tabs$ref_tab,
    tabs$source_tab,
    match_on = "id",
    split_on = "id",
    n_folds = 4,
    permutations = 0,
    method = "cosine",
    similarity_transform = affine_transform,
    similarity_transform_args = list(shrink = 1e-6),
    seed = 1
  )

  expect_gt(mean(cv$eye_sim), mean(raw$eye_sim))
})
