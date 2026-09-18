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
