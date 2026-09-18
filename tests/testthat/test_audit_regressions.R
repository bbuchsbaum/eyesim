# Regression tests for defects found while auditing eyesim against the eyes4s
# parity baseline. Each block names the audit item it covers.

audit_fg <- function() {
  fixation_group(x = c(10, 30, 80), y = c(10, 20, 40),
                 onset = c(0, 100, 200), duration = c(1, 2, 4))
}

audit_density <- function(fg, ...) {
  suppressMessages(eye_density(fg, sigma = 10, xbounds = c(0, 100), ybounds = c(0, 50),
                               outdim = c(5, 3), ...))
}

# Audit item 1 ---------------------------------------------------------------
test_that("eye_density honours explicit fixation weights", {
  fg <- audit_fg()
  d0 <- audit_density(fg)
  dw <- audit_density(fg, weights = c(4, 2, 1))
  expect_false(isTRUE(all.equal(d0$z, dw$z)))

  # Explicit weights equal to the durations reproduce duration weighting.
  expect_equal(audit_density(fg, weights = c(1, 2, 4))$z,
               audit_density(fg, duration_weighted = TRUE)$z)

  # Explicit weights take precedence over duration weighting.
  expect_equal(audit_density(fg, weights = c(4, 2, 1), duration_weighted = TRUE)$z, dw$z)

  # A zero weight removes a fixation: the map equals one built without it.
  expect_equal(audit_density(fg, weights = c(1, 1, 0))$z,
               audit_density(fg[1:2, ])$z)

  # Weights are aligned with the rows that survive the window filter.
  expect_equal(audit_density(fg, weights = c(4, 2, 1), window = c(50, 300))$z,
               audit_density(fg[2:3, ], weights = c(2, 1))$z)

  expect_error(audit_density(fg, weights = c(1, 2)), "one value per fixation")
  expect_error(audit_density(fg, weights = c(1, -1, 2)), "non-negative")
})

# Audit item 2 ---------------------------------------------------------------
test_that("eye_density forwards extra arguments to ks::kde", {
  fg <- audit_fg()
  dens <- audit_density(fg, binned = FALSE)
  expect_s3_class(dens, "eye_density")

  ref <- ks::kde(cbind(fg$x, fg$y), H = diag(c(100, 100)), gridsize = c(5, 3),
                 xmin = c(0, 0), xmax = c(100, 50), binned = FALSE,
                 compute.cont = FALSE)$estimate
  expect_equal(dens$z, ref / sum(ref), tolerance = 1e-6)

  ms <- eye_density(fg, sigma = c(5, 10), xbounds = c(0, 100), ybounds = c(0, 50),
                    outdim = c(5, 3), binned = FALSE)
  expect_s3_class(ms, "eye_density_multiscale")
  expect_equal(ms[[2]]$z, dens$z)

  expect_error(audit_density(fg, not_a_kde_argument = 1), "Unsupported ks::kde")
  expect_error(audit_density(fg, gridsize = c(2, 2)), "Unsupported ks::kde")
  expect_error(audit_density(fg, kde_pkg = "MASS", binned = FALSE),
               "not supported when kde_pkg")
})

# Audit item 3 ---------------------------------------------------------------
test_that("weighted MASS densities are computed on non-square grids", {
  fg <- audit_fg()
  mass_density <- function(fg, ...) {
    suppressMessages(eye_density(fg, sigma = 40, xbounds = c(0, 100), ybounds = c(0, 50),
                                 outdim = c(5, 3), kde_pkg = "MASS", ...))
  }

  dw <- expect_no_warning(mass_density(fg, duration_weighted = TRUE))
  expect_s3_class(dw, "eye_density")
  expect_equal(dim(dw$z), c(5L, 3L))
  expect_false(isTRUE(all.equal(dw$z, mass_density(fg)$z)))

  # With uniform weights, the weighted kernel reproduces MASS::kde2d exactly.
  uw <- kde2d_weighted(fg$x, fg$y, h = 40, n = c(5, 3), lims = c(0, 100, 0, 50),
                       w = c(2, 2, 2))
  ref <- MASS::kde2d(fg$x, fg$y, h = 40, n = c(5, 3), lims = c(0, 100, 0, 50))
  expect_equal(uw$z, ref$z)

  # A zero weight removes a fixation from the weighted MASS map.
  expect_equal(mass_density(fg, weights = c(1, 1, 0))$z, mass_density(fg[1:2, ])$z)
})

# Audit item 5 ---------------------------------------------------------------
test_that("density maps do not depend on options(digits)", {
  fg <- audit_fg()
  old <- options(digits = 7)
  on.exit(options(old), add = TRUE)
  d7 <- audit_density(fg)$z
  options(digits = 17)
  d17 <- audit_density(fg)$z
  options(digits = 3)
  d3 <- audit_density(fg)$z
  expect_identical(d17, d7)
  expect_identical(d3, d7)

  # The geometric warp used by the affine and contract transforms is rounded
  # with the same fixed precision.
  dens <- audit_density(fg)
  A <- matrix(c(1.1, 0.05, 0, 0.9), 2)
  options(digits = 7)
  w7 <- warp_density_object(dens, A = A, t = c(1, -2))$z
  options(digits = 17)
  w17 <- warp_density_object(dens, A = A, t = c(1, -2))$z
  expect_identical(w17, w7)
})

# Audit items 6 and 8 --------------------------------------------------------
audit_unit_map <- function(i) {
  z <- rep(0, 4)
  z[i] <- 1
  gen_density(x = 1:2, y = 1:2, z = matrix(z, 2))
}

audit_perm_tables <- function(keys) {
  ref <- tibble::tibble(key = c("a", "b", "c"), participant = "p",
                        density = lapply(1:3, audit_unit_map))
  src <- tibble::tibble(key = keys, participant = "p",
                        density = lapply(match(keys, ref$key), audit_unit_map))
  list(ref = ref, src = src)
}

test_that("the true match is excluded before permutation candidates are sampled", {
  tabs <- audit_perm_tables(c("a", "b", "c"))
  for (method in c("cosine", "pearson")) {
    for (seed in 1:15) {
      set.seed(seed)
      res <- suppressMessages(template_similarity(
        tabs$ref, tabs$src, "key", permute_on = "participant",
        method = method, permutations = 2
      ))
      # Two non-matching candidates per row, so permutations = 2 is exhaustive.
      expect_equal(res$n_perm, rep(2L, 3), info = paste(method, seed))

      set.seed(seed)
      res1 <- suppressMessages(template_similarity(
        tabs$ref, tabs$src, "key", permute_on = "participant",
        method = method, permutations = 1
      ))
      expect_equal(res1$n_perm, rep(1L, 3), info = paste(method, seed))
      # The unit maps are orthogonal: a control is never the true match.
      control <- if (method == "cosine") 0 else -1 / 3
      expect_equal(res1$perm_sim, rep(control, 3), info = paste(method, seed))
    }
  }
})

test_that("duplicated focal keys never enter their own permutation baseline", {
  tabs <- audit_perm_tables(c("a", "a", "b", "c"))
  for (method in c("cosine", "pearson")) {
    res <- suppressMessages(template_similarity(
      tabs$ref, tabs$src, "key", permute_on = "participant",
      method = method, permutations = 100
    ))
    a_rows <- res$key == "a"
    control <- if (method == "cosine") 0 else -1 / 3
    expect_equal(res$perm_sim[a_rows], rep(control, 2), info = method)
    expect_equal(res$n_perm[a_rows], rep(2L, 2), info = method)
  }
})

# Audit item 7 (documentation) -----------------------------------------------
# The documented behaviour: sampling uses the session RNG, so set.seed() makes
# the baseline reproducible and different seeds can select different controls.
test_that("permutation baselines follow the session RNG", {
  mk <- function(i) {
    z <- rep(0, 4)
    z[i] <- 1
    z[(i %% 4) + 1] <- 0.5
    gen_density(x = 1:2, y = 1:2, z = matrix(z, 2))
  }
  ref <- tibble::tibble(key = letters[1:4], participant = "p", density = lapply(1:4, mk))
  run <- function(method, seed) {
    set.seed(seed)
    suppressMessages(template_similarity(ref, ref, "key", permute_on = "participant",
                                         method = method, permutations = 1))$perm_sim
  }
  for (method in c("cosine", "pearson")) {
    expect_identical(run(method, 1), run(method, 1), info = method)
    draws <- lapply(1:10, function(s) run(method, s))
    expect_gt(length(unique(draws)), 1L)
  }
})

# Audit item 24 --------------------------------------------------------------
test_that("template_similarity_cv leaves the caller's RNG stream untouched", {
  mk <- function(i) gen_density(x = 1:2, y = 1:2, z = matrix(c(i, 7 - i, (i %% 3) + 1, 2), 2))
  ref <- tibble::tibble(key = letters[1:6], participant = "p", density = lapply(1:6, mk))
  cv <- function(...) {
    suppressMessages(template_similarity_cv(ref, ref, "key", permute_on = "participant",
                                            method = "pearson", n_folds = 2, seed = 1, ...))
  }

  set.seed(42)
  expected <- runif(3)
  set.seed(42)
  res1 <- cv(permutations = 1)
  expect_identical(runif(3), expected)

  # Results depend on `seed` only, not on the caller's RNG state.
  set.seed(7)
  res2 <- cv(permutations = 1)
  expect_identical(res2$perm_sim, res1$perm_sim)

  # Without a prior .Random.seed, none is left behind.
  if (exists(".Random.seed", envir = globalenv())) {
    saved <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
    rm(".Random.seed", envir = globalenv())
  }
  cv(permutations = 0)
  expect_false(exists(".Random.seed", envir = globalenv()))

  # Controls come from the held-out fold only: 3 keys per fold, so 2 controls.
  full <- cv(permutations = 100)
  expect_equal(full$n_perm, rep(2L, 6))
})

# Audit item 9 ---------------------------------------------------------------
test_that("fisherz gives the same clamped value for every pair of identical maps", {
  z_max <- atanh(1 - .Machine$double.eps)
  constant <- rep(0.25, 4)
  varying <- c(0.1, 0.2, 0.3, 0.4)
  expect_equal(similarity(constant, constant, method = "fisherz"), z_max)
  expect_equal(similarity(varying, varying, method = "fisherz"), z_max)

  flat <- gen_density(x = 1:2, y = 1:2, z = matrix(constant, 2))
  expect_equal(similarity(flat, flat, method = "fisherz"), z_max)

  # Correlation-scale methods still report r = 1 for identical constant maps.
  expect_equal(similarity(constant, constant, method = "pearson"), 1)
  expect_equal(similarity(constant, constant, method = "spearman"), 1)
})

# Audit item 10 --------------------------------------------------------------
test_that("density maps on different lattices are refused", {
  z <- matrix(c(1, 2, 3, 4), 2)
  a <- gen_density(x = 1:2, y = 1:2, z = z)
  b <- gen_density(x = c(10, 20), y = c(10, 20), z = z)
  for (method in c("pearson", "cosine", "fisherz", "l1")) {
    expect_error(similarity(a, b, method = method), "different lattices", info = method)
  }
  expect_equal(similarity(a, gen_density(x = 1:2, y = 1:2, z = z), method = "pearson"), 1)
  expect_error(similarity(a, 1:5, method = "pearson"), "grid cells")

  # template_similarity refuses on both the fast cosine and the general path.
  ref <- tibble::tibble(key = c("k1", "k2"), density = list(a, a))
  src <- tibble::tibble(key = c("k1", "k2"), density = list(b, b))
  for (method in c("cosine", "pearson")) {
    expect_error(suppressMessages(template_similarity(ref, src, "key", method = method,
                                                      permutations = 0)),
                 "different lattices", info = method)
  }
})

# Audit item 11 --------------------------------------------------------------
test_that("similarity() and fixation_overlap() share the overlap threshold", {
  expect_identical(formals(eyesim:::similarity.fixation_group)$dthresh,
                   formals(fixation_overlap)$dthresh)

  # Fixations 50 units apart overlap under the documented default of 60.
  a <- fixation_group(x = c(0, 100), y = c(0, 0), onset = c(0, 100), duration = c(100, 100))
  b <- fixation_group(x = c(50, 150), y = c(0, 0), onset = c(0, 100), duration = c(100, 100))
  times <- c(0, 50, 100)
  expect_equal(similarity(a, b, method = "overlap", time_samples = times), 1)
  expect_equal(similarity(a, b, method = "overlap", time_samples = times),
               fixation_overlap(a, b, time_samples = times)$perc)
})

# Audit item 12 --------------------------------------------------------------
test_that("the default overlap time grid spans both fixation groups", {
  a <- fixation_group(x = c(0, 0), y = c(0, 0), onset = c(0, 1000), duration = c(1000, 1000))
  b <- fixation_group(x = c(0, 500), y = c(0, 0), onset = c(0, 100), duration = c(100, 100))
  expect_equal(fixation_overlap(a, b)$perc, fixation_overlap(b, a)$perc)
  expect_equal(fixation_overlap(a, b)$perc,
               fixation_overlap(a, b, time_samples = seq(0, 1000, by = 20))$perc)
})

# Audit item 13 --------------------------------------------------------------
test_that("sample_fixations holds the last fixation on both paths", {
  fg <- fixation_group(x = c(0, 1), y = c(0, 1), onset = c(0, 100), duration = c(100, 100))
  times <- c(-10, 0, 50, 100, 150, 200)
  fast <- sample_fixations(fg, times)
  slow <- sample_fixations(fg, times, fast = FALSE)
  expect_equal(fast$x, c(NA, 0, 0, 1, 1, 1))
  expect_equal(fast$y, c(NA, 0, 0, 1, 1, 1))
  expect_equal(fast$x, slow$x)
  expect_equal(fast$y, slow$y)

  # A single fixation no longer fails on the fast path.
  one <- fixation_group(x = 5, y = 6, onset = 10, duration = 100)
  expect_equal(sample_fixations(one, c(0, 10, 500))$x, c(NA, 5, 5))
  expect_equal(sample_fixations(one, c(0, 10, 500), fast = FALSE)$x, c(NA, 5, 5))

  # Density sampling over time follows the same rule.
  tmpl <- gen_density(x = c(0, 1), y = c(0, 1), z = matrix(c(0.1, 0.2, 0.3, 0.4), 2))
  res <- sample_density_time(tibble::tibble(k = "a", density = list(tmpl)),
                             tibble::tibble(k = "a", fixgroup = list(fg)), "k",
                             times = c(0, 50, 100, 150, 200))
  expect_equal(res$sampled[[1]]$z, c(0.1, 0.1, 0.4, 0.4, 0.4))
})

# Audit item 14 --------------------------------------------------------------
test_that("rep_fixations counts replicates without floating-point loss", {
  fg <- fixation_group(x = c(1, 2, 3, 4), y = rep(1, 4), onset = c(0, 1, 2, 3),
                       duration = c(0.29, 0.57, 0.295, 0.001))
  reps <- rep_fixations(fg, 100)
  expect_equal(as.vector(table(reps$index)), c(29L, 57L, 29L, 1L))
  expect_equal(nrow(rep_fixations(fixation_group(x = 5.5, y = 1, onset = 0, duration = 0.29), 100)), 29L)
})

# Audit item 16 --------------------------------------------------------------
test_that("sample_density_time bins are half-open, including the last", {
  tmpl <- gen_density(x = c(0, 1), y = c(0, 1), z = matrix(c(0.1, 0.2, 0.3, 0.4), 2))
  # Fixation at (0, 0) until t = 200, then at (1, 1).
  fg <- fixation_group(x = c(0, 1), y = c(0, 1), onset = c(0, 200), duration = c(200, 100))
  res <- sample_density_time(tibble::tibble(k = "a", density = list(tmpl)),
                             tibble::tibble(k = "a", fixgroup = list(fg)), "k",
                             times = c(0, 100, 200), time_bins = c(0, 100, 200),
                             aggregate_fun = function(v, na.rm) length(v))
  # t = 0 falls in [0, 100), t = 100 in [100, 200), and t = 200 in no bin.
  expect_equal(c(res$bin_1, res$bin_2), c(1, 1))
})

# Audit item 17 --------------------------------------------------------------
test_that("mm_position_emd compares every fixation, including the last", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("emdist")
  p1 <- scanpath(fixation_group(x = c(10, 40, 80), y = c(10, 20, 40),
                                onset = c(0, 100, 250), duration = c(80, 120, 100)))
  p2 <- scanpath(fixation_group(x = c(10, 40, 5), y = c(10, 20, 45),
                                onset = c(0, 100, 250), duration = c(80, 120, 100)))
  mm <- multi_match(p1, p2, screensize = c(100, 50))
  expect_lt(mm[["mm_position_emd"]], 1)

  emd <- emdist::emdw(cbind(p1$x, p1$y), p1$duration, cbind(p2$x, p2$y), p2$duration)
  expect_equal(mm[["mm_position_emd"]], 1 - emd / sqrt(100^2 + 50^2))
  expect_equal(multi_match(p1, p1, screensize = c(100, 50))[["mm_position_emd"]], 1)
})

# Audit item 18 --------------------------------------------------------------
test_that("template_multireg defaults to lm", {
  m <- tibble::tibble(
    response = list(list(z = matrix(c(1, 2, 3, 4, 5, 7), 3))),
    a = list(list(z = matrix(c(1, 0, 1, 2, 1, 1), 3))),
    b = list(list(z = matrix(c(0, 1, 1, 1, 2, 3), 3)))
  )
  default <- template_multireg(m, "response", c("a", "b"))
  explicit <- template_multireg(m, "response", c("a", "b"), method = "lm")
  expect_equal(default$multireg, explicit$multireg)
  expect_error(template_multireg(m, "response", c("a", "b"), method = "ols"), "should be one of")
})

# Audit item 19 --------------------------------------------------------------
test_that("template_regression refuses a duplicated baseline key", {
  mk <- function(v) gen_density(x = 1:2, y = 1:2, z = matrix(v, 2))
  ref <- tibble::tibble(key = c("a", "b"), density = list(mk(c(1, 2, 3, 4)), mk(c(4, 3, 2, 1))))
  src <- tibble::tibble(key = c("a", "b"), base = "x",
                        density = list(mk(c(1, 2, 3, 5)), mk(c(4, 3, 2, 2))))
  dup <- tibble::tibble(base = c("x", "x"), density = list(mk(c(1, 1, 2, 2)), mk(c(2, 2, 1, 1))))
  expect_error(template_regression(ref, src, "key", dup, "base"),
               "more than one row for base = x")

  # An unused duplicate does not block the rows that have a unique baseline.
  ok <- tibble::tibble(base = c("x", "y", "y"),
                       density = list(mk(c(1, 1, 2, 3)), mk(c(2, 2, 1, 1)), mk(c(2, 2, 1, 1))))
  res <- template_regression(ref, src, "key", ok, "base")
  expect_equal(nrow(res), 2L)
  expect_true(all(is.finite(res$beta_source)))
})
