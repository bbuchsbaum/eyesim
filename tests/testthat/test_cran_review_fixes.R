# Regressions for the pre-CRAN correctness review.

test_that("calcangle returns 0, not NaN, for parallel vectors", {
  expect_equal(calcangle(c(.1, .7), c(.3, 2.1)), 0, tolerance = 1e-5)
  expect_equal(calcangle(c(1, 2), c(-2, -4)), 180)
  set.seed(1)
  angles <- replicate(200, {
    v <- runif(2)
    calcangle(v, runif(1, 0.1, 10) * v)
  })
  expect_false(anyNA(angles))
})

test_that("density_by works without groups", {
  fixations <- data.frame(
    subject = rep(c("s1", "s2"), each = 10),
    x = c(seq(100, 900, length.out = 10), seq(150, 850, length.out = 10)),
    y = c(seq(200, 800, length.out = 10), seq(250, 750, length.out = 10)),
    duration = 200,
    onset = rep(seq(0, 1800, by = 200), 2)
  )
  et <- eye_table("x", "y", "duration", "onset", groupvar = "subject", data = fixations)
  res <- density_by(et, sigma = 80, xbounds = c(0, 1000), ybounds = c(0, 1000))
  expect_equal(nrow(res), 1L)
  expect_s3_class(res$density[[1]], "eye_density")
})

test_that("scanpath orders fixations by onset", {
  sorted <- fixation_group(x = c(0, 100, 100), y = c(0, 0, 100),
                           onset = c(0, 200, 400), duration = rep(200, 3))
  shuffled <- sorted[c(3, 1, 2), ]
  expect_equal(scanpath(shuffled)$lenx, scanpath(sorted)$lenx)
  expect_equal(scanpath(shuffled)$theta, scanpath(sorted)$theta)
})

test_that("fixation_overlap ignores time before the first fixation", {
  fg <- fixation_group(x = c(100, 200, 300), y = c(100, 200, 300),
                       onset = c(1000, 2000, 3000), duration = rep(1000, 3))
  expect_equal(fixation_overlap(fg, fg)$perc, 1)
})

test_that("empty fixation groups do not crash", {
  fg <- fixation_group(x = numeric(0), y = numeric(0), onset = numeric(0),
                       duration = numeric(0))
  expect_equal(nrow(fg), 0L)
  expect_equal(nrow(rep_fixations(fg)), 0L)
})

test_that("suggest_sigma returns NA for coincident fixations without bounds", {
  fg <- fixation_group(x = rep(500, 5), y = rep(500, 5), onset = seq(0, 800, by = 200),
                       duration = rep(200, 5))
  expect_true(is.na(suggest_sigma(fg)))
  expect_gt(suggest_sigma(fg, xbounds = c(0, 1000), ybounds = c(0, 1000)), 0)
})

test_that("template_regression names a missing baseline key", {
  mk <- function(v) gen_density(x = 1:2, y = 1:2, z = matrix(v, 2))
  ref <- tibble::tibble(key = c("a", "b"), density = list(mk(c(1, 2, 3, 4)), mk(c(4, 3, 2, 1))))
  src <- tibble::tibble(key = c("a", "b"), base = c("x", "z"),
                        density = list(mk(c(1, 2, 3, 5)), mk(c(4, 3, 2, 2))))
  base <- tibble::tibble(base = "x", density = list(mk(c(1, 1, 2, 3))))
  expect_error(template_regression(ref, src, "key", base, "base"), "no row for base = z")
})

test_that("template_regression reads a custom density column", {
  mk <- function(v) gen_density(x = 1:2, y = 1:2, z = matrix(v, 2))
  ref <- tibble::tibble(key = c("a", "b"), dens = list(mk(c(1, 2, 3, 4)), mk(c(4, 3, 2, 1))))
  src <- tibble::tibble(key = c("a", "b"), base = "x",
                        dens = list(mk(c(1, 2, 3, 5)), mk(c(4, 3, 2, 2))))
  base <- tibble::tibble(base = "x", dens = list(mk(c(1, 1, 2, 3))))
  res <- template_regression(ref, src, "key", base, "base", density_var = "dens")
  expect_true(all(is.finite(res$beta_source)))
})

latent_density <- function(vec) {
  structure(list(z = matrix(vec, 3, 3), x = 1:3, y = 1:3, sigma = 1),
            class = c("density", "eye_density"))
}

test_that("coral_transform gives source scores the reference covariance", {
  set.seed(11)
  n <- 60
  ref_mat <- matrix(rnorm(n * 9), n, 9) %*% matrix(rnorm(81), 9, 9)
  src_mat <- matrix(rnorm(n * 9), n, 9) %*% matrix(rnorm(81), 9, 9)
  ref_tab <- tibble::tibble(id = seq_len(n), density = lapply(seq_len(n), function(i) latent_density(ref_mat[i, ])))
  source_tab <- tibble::tibble(id = seq_len(n), density = lapply(seq_len(n), function(i) latent_density(src_mat[i, ])))

  res <- coral_transform(ref_tab, source_tab, match_on = "id", comps = 3, shrink = 1e-10)
  ref_scores <- do.call(rbind, res$ref_tab$density)
  src_scores <- do.call(rbind, res$source_tab$density)
  expect_equal(cov(src_scores), cov(ref_scores), tolerance = 1e-6)
})

test_that("latent matrix square roots handle a single component", {
  m <- matrix(4, 1, 1)
  expect_equal(mat_sqrt(m), matrix(2, 1, 1))
  expect_equal(mat_inv_sqrt(m, 1e-8), matrix(0.5, 1, 1))
})

test_that("cca_transform outputs canonical variates only", {
  set.seed(12)
  n <- 6
  ref_mat <- matrix(rnorm(n * 9), n, 9)
  src_mat <- ref_mat %*% matrix(rnorm(81), 9, 9)
  ref_tab <- tibble::tibble(id = seq_len(n), density = lapply(seq_len(n), function(i) latent_density(ref_mat[i, ])))
  source_tab <- tibble::tibble(id = seq_len(n), density = lapply(seq_len(n), function(i) latent_density(src_mat[i, ])))

  res <- cca_transform(ref_tab, source_tab, match_on = "id", comps = 8)
  k <- res$info$groups[[1]]$comps
  width <- length(res$ref_tab$density[[1]])
  expect_lt(k, width)
  # Columns past the canonical variates are zero, not leftover PCA scores.
  tail_cols <- function(v) v[-seq_len(k)]
  expect_true(all(vapply(res$ref_tab$density, function(v) all(tail_cols(v) == 0), logical(1))))
  expect_true(all(vapply(res$source_tab$density, function(v) all(tail_cols(v) == 0), logical(1))))
})

test_that("Replay folds keep a true candidate and a nonmatch in every fold", {
  tab <- data.frame(item = sprintf("i%02d", 1:9), stratum = "s")
  all_rows <- rep(TRUE, nrow(tab))
  folds <- make_gaze_weave_candidate_folds(tab, "item", "stratum", NULL, 1, "item", all_rows)
  expect_equal(folds$n_folds, 4L)
  expect_true(all(table(folds$fold_id) >= 2L))
  expect_error(
    make_gaze_weave_candidate_folds(tab, "item", "stratum", 5, 1, "item", all_rows),
    "fewer than two"
  )
  # Holding out participants keeps every item in each fold: the default stays.
  tab2 <- expand.grid(participant = sprintf("p%d", 1:6), item = sprintf("i%d", 1:3),
                      stringsAsFactors = FALSE)
  folds2 <- make_gaze_weave_candidate_folds(tab2, "participant", NULL, NULL, 1, "item",
                                            rep(TRUE, nrow(tab2)))
  expect_equal(folds2$n_folds, 5L)
  # The underlying fold builder is unchanged.
  expect_equal(make_gaze_weave_folds(tab, "item", "stratum")$n_folds, 5L)
})
