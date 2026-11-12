# Regression tests for bugs found by differential testing against peyesim
# (GitHub issues #42-#51).

has_emd_backend <- function() {
  any(vapply(c("emdist", "transport", "T4transport"), requireNamespace,
             logical(1), quietly = TRUE))
}

ms_density <- function(...) {
  structure(list(...), class = c("eye_density_multiscale", "list"))
}

sigma_density <- function(v, sigma) {
  d <- gen_density(1:3, 1:2, matrix(v, 3, 2))
  d$sigma <- sigma
  d
}

test_that("match_scale forwards only the caller's extra arguments (#42)", {
  fg <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                       onset = c(0, 200, 400), duration = c(200, 200, 200))
  ref <- tibble::tibble(trial = 1, fixgroup = list(fg))
  src <- tibble::tibble(trial = 1, fixgroup = list(rescale(fg, 2, 2)))

  res <- eyesim:::match_scale(ref, src, match_on = "trial")
  expect_equal(res$scale_x, 0.5, tolerance = 1e-3)
  expect_equal(res$scale_y, 0.5, tolerance = 1e-3)
})

test_that("crqa wrapper gives a 1x1 plot when embedding leaves one point (#43)", {
  skip_if_not_installed("crqa")
  fg1 <- fixation_group(x = c(100, 200), y = c(100, 150), onset = c(0, 200), duration = c(200, 200))
  fg2 <- fixation_group(x = c(110, 205), y = c(105, 150), onset = c(0, 200), duration = c(200, 200))

  res <- eyesim:::crqa(fg1, fg2, radius = 60, embed = 2)
  expect_equal(dim(res$RP), c(1L, 1L))
  expect_equal(res$RR, 100)

  expect_error(eyesim:::crqa(fg1, fg2, radius = 60, embed = 3), "no embedded points")
})

test_that("c() skips an empty leading fixation group (#44)", {
  empty <- fixation_group(x = numeric(0), y = numeric(0), onset = numeric(0), duration = numeric(0))
  fg <- fixation_group(x = c(1, 2), y = c(3, 4), onset = c(0, 100), duration = c(100, 100))

  out <- expect_silent(c(empty, fg))
  expect_equal(out$onset, c(0, 100))
})

test_that("fast cosine path scores every all-zero map pair as 1 (#45)", {
  expect_equal(eyesim:::cosine_similarity_matrix(matrix(0, 3, 6), matrix(0, 3, 6)),
               matrix(1, 3, 3))

  zero <- gen_density(1:3, 1:2, matrix(0, 3, 2))
  tab <- tibble::tibble(k = 1:3, density = list(zero, zero, zero))
  res <- template_similarity(tab, tab, match_on = "k", method = "cosine", permutations = 0)
  expect_equal(res$eye_sim, c(1, 1, 1))
})

test_that("multiscale 'none' aggregation tolerates a NULL map (#46)", {
  a <- ms_density(sigma_density(c(1, 5, 2, 8, 3, 3), 20), sigma_density(c(2, 2, 3, 1, 0, 4), 10))
  b <- ms_density(sigma_density(c(2, 1, 4, 3, 9, 0), 20), sigma_density(c(1, 3, 3, 2, 2, 5), 10))
  ref <- tibble::tibble(k = 1:2, density = list(a, b))
  src <- tibble::tibble(k = 1:2, density = list(b, NULL))

  res <- expect_warnings(
    template_similarity(ref, src, match_on = "k", method = "pearson",
                        permutations = 0, multiscale_aggregation = "none"),
    regexp = "Invalid or NULL"
  )
  expect_equal(nrow(res), 2L)
  expect_length(res$eye_sim[[1]], 2L)
  expect_true(is.na(res$eye_sim[[2]]))
})

test_that("NULL templates are never drawn as permutation controls (#47)", {
  d1 <- gen_density(1:3, 1:2, matrix(c(1, 5, 2, 8, 3, 3), 3, 2))
  d2 <- gen_density(1:3, 1:2, matrix(c(2, 1, 4, 3, 9, 1), 3, 2))
  ref <- tibble::tibble(k = 1:2, density = list(d1, NULL))
  src <- tibble::tibble(k = 1:2, density = list(d2, d1))

  methods <- c("pearson", if (has_emd_backend()) "emd")
  for (m in methods) {
    res <- suppressWarnings(
      template_similarity(ref, src, match_on = "k", method = m, permutations = 5)
    )
    expect_false(is.na(res$eye_sim[1]))
    expect_true(is.na(res$perm_sim[1]))
    expect_equal(res$n_perm, c(0L, 0L))
    expect_true(is.na(res$eye_sim[2]))
  }
})

test_that("multiscale permutation baselines are averaged by sigma (#48)", {
  a  <- ms_density(sigma_density(c(1, 5, 2, 8, 3, 3), 20), sigma_density(c(2, 2, 3, 1, 0, 4), 10))
  b  <- ms_density(sigma_density(c(2, 1, 4, 3, 9, 0), 10), sigma_density(c(1, 3, 3, 2, 2, 5), 20))
  c3 <- ms_density(sigma_density(c(3, 1, 4, 1, 5, 9), 20), sigma_density(c(2, 6, 5, 3, 5, 8), 10))
  s  <- ms_density(sigma_density(c(2, 2, 3, 1, 0, 4), 20), sigma_density(c(3, 1, 4, 1, 5, 9), 10))
  ref <- tibble::tibble(k = 1:3, density = list(a, b, c3))
  src <- tibble::tibble(k = 1:3, density = list(s, a, a))

  res <- template_similarity(ref, src, match_on = "k", method = "pearson",
                             permutations = 5, multiscale_aggregation = "none")
  pb <- similarity(b, s, method = "pearson", multiscale_aggregation = "none")
  pc <- similarity(c3, s, method = "pearson", multiscale_aggregation = "none")
  expected <- c(sigma_20 = mean(c(pb[["sigma_20"]], pc[["sigma_20"]])),
                sigma_10 = mean(c(pb[["sigma_10"]], pc[["sigma_10"]])))

  expect_equal(res$perm_sim[[1]][names(expected)], expected)
})

test_that("cca_transform passes NA fit_by rows through unfitted (#49)", {
  set.seed(1)
  maps <- lapply(1:8, function(i) gen_density(1:3, 1:2, matrix(runif(6), 3, 2)))
  ref <- tibble::tibble(k = 1:8, g = rep(c("a", "b"), 4), density = maps)
  src <- ref
  src$g[1] <- NA

  expect_warning(cca_transform(ref, src, match_on = "k", comps = 2, fit_by = "g"),
                 "no CCA model for group\\(s\\) NA")
  res <- expect_warning(
    template_similarity(ref, src, match_on = "k", method = "pearson", permutations = 0,
                        similarity_transform = cca_transform,
                        similarity_transform_args = list(comps = 2, fit_by = "g")),
    "no CCA model"
  )
  expect_equal(nrow(res), 8L)
})

test_that("MultiMatch handles all-zero durations (#50)", {
  expect_equal(eyesim:::duration_diff_1d(c(0, 100, 50), c(0, 50, 50), cbind(1:3, 1:3)),
               c(0, 0.5, 0))

  skip_if_not(has_emd_backend(), "no EMD backend installed")
  fg1 <- fixation_group(x = c(100, 300, 500, 200), y = c(100, 400, 200, 300),
                        onset = c(0, 100, 200, 300), duration = c(0, 0, 0, 0))
  fg2 <- fixation_group(x = c(120, 280, 520, 210), y = c(110, 390, 220, 310),
                        onset = c(0, 100, 200, 300), duration = c(0, 0, 0, 0))
  res <- similarity(scanpath(fg1), scanpath(fg2), screensize = c(1000, 1000))
  expect_equal(res[["mm_duration"]], 1)
  # Positions differ, so equal-weight EMD must not report a perfect match.
  expect_lt(res[["mm_position_emd"]], 1)
})

test_that("similarity() for scanpaths rejects an unknown method (#51)", {
  fg1 <- fixation_group(x = c(100, 200, 300), y = c(100, 150, 200),
                        onset = c(0, 200, 400), duration = c(200, 200, 200))
  fg2 <- fixation_group(x = c(110, 210, 290), y = c(110, 140, 210),
                        onset = c(0, 200, 400), duration = c(200, 200, 200))
  expect_error(similarity(scanpath(fg1), scanpath(fg2), method = "not_a_method",
                          screensize = c(1000, 1000)),
               "should be")
})
