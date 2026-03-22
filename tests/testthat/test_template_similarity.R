library(testthat)
library(dplyr)

context("density_by")

test_that("template_similarity produces perfect similarity for identical patterns", {
  options(future.rng.onMisuse = "ignore")
  g1 <- tibble(fixgroup=lapply(1:10, function(i) {
    x <- runif(10)
    y <- runif(10)
    onset <- seq(1,length.out=length(x), by=50)
    duration <- rep(1,length(x))
    fixgroup <- fixation_group(x,y,onset,duration)
  }), image=1:10)

  g2 <- tibble(fixgroup=lapply(1:10, function(i) {
    x <- runif(10)
    y <- runif(10)
    onset <- seq(1,length.out=length(x), by=50)
    duration <- rep(1,length(x))
    fixgroup <- fixation_group(x,y,onset,duration)
  }), image=1:10)

  dens <- density_by(g1, "image", xbounds=c(0,1), ybounds=c(0,1))
  dens2 <- density_by(g2, "image", xbounds=c(0,1), ybounds=c(0,1))
  tsim <- template_similarity(dens, dens2, match_on="image", 
        method="spearman", permutations=3)
  expect_true(max(tsim$eye_sim) <=1)
  expect_true(min(tsim$eye_sim) >=-1)

})

test_that("template_similarity works for permute_on", {
  g1 <- tibble(fixgroup=lapply(1:100, function(i) {
    x <- runif(10)
    y <- runif(10)
    onset <- seq(1,length.out=length(x), by=50)
    duration <- rep(1,length(x))
    fixgroup <- fixation_group(x,y,onset,duration)
  }), image=1:100, subject=rep(1:10, each=10))

  g2 <- g1
  dens <- density_by(g1, "image", keep_vars="subject", xbounds=c(0,1), ybounds=c(0,1), duration_weighted=TRUE)
  dens2 <- density_by(g2, "image", keep_vars="subject", xbounds=c(0,1), ybounds=c(0,1), duration_weighted = TRUE)
  tsim <- template_similarity(dens, dens2, match_on="image", method="pearson", permute_on="subject",
                              permutations=6)
  expect_true(all(tsim$eye_sim >.99))

})

test_that("template_similarity cosine permutations match manual baseline computation", {
  make_density <- function(vals) {
    structure(
      list(z = matrix(vals, nrow = 2, ncol = 2, byrow = TRUE), x = 1:2, y = 1:2, sigma = 50),
      class = c("density", "eye_density")
    )
  }

  ref_tab <- tibble(
    image = 1:6,
    subject = rep(1:2, each = 3),
    density = list(
      make_density(c(1, 2, 3, 4)),
      make_density(c(2, 3, 4, 5)),
      make_density(c(3, 4, 5, 6)),
      make_density(c(6, 5, 4, 3)),
      make_density(c(5, 4, 3, 2)),
      make_density(c(4, 3, 2, 1))
    )
  )
  source_tab <- tibble(
    row_id = 1:6,
    image = 1:6,
    subject = rep(1:2, each = 3),
    density = list(
      make_density(c(1.1, 2.1, 3.1, 4.1)),
      make_density(c(2.1, 3.1, 4.1, 5.1)),
      make_density(c(3.1, 4.1, 5.1, 6.1)),
      make_density(c(6.1, 5.1, 4.1, 3.1)),
      make_density(c(5.1, 4.1, 3.1, 2.1)),
      make_density(c(4.1, 3.1, 2.1, 1.1))
    )
  )

  manual_perm <- function(ref_tab, source_tab, permutations) {
    matchind <- match(source_tab$image, ref_tab$image)
    match_split <- split(matchind, source_tab$subject)

    purrr::map_dfr(seq_len(nrow(source_tab)), function(i) {
      d1 <- ref_tab$density[[matchind[[i]]]]
      d2 <- source_tab$density[[i]]
      eye_sim <- similarity(d1, d2, method = "cosine")

      candidates <- match_split[[as.character(source_tab$subject[[i]])]]
      if (permutations < length(candidates)) {
        candidates <- sample(candidates, permutations)
      }
      match_pos <- match(matchind[[i]], candidates)
      if (!is.na(match_pos)) {
        candidates <- candidates[-match_pos]
      }

      perm_vals <- vapply(candidates, function(j) {
        similarity(ref_tab$density[[j]], d2, method = "cosine")
      }, numeric(1))
      perm_sim <- if (length(perm_vals) == 0L) NA_real_ else mean(perm_vals, na.rm = TRUE)

      tibble(
        row_id = source_tab$row_id[[i]],
        eye_sim = eye_sim,
        perm_sim = perm_sim,
        eye_sim_diff = eye_sim - perm_sim
      )
    })
  }

  set.seed(42)
  fast <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "image",
    permute_on = "subject",
    refvar = "density",
    sourcevar = "density",
    method = "cosine",
    permutations = 2
  )

  set.seed(42)
  manual <- manual_perm(ref_tab, source_tab, permutations = 2)

  fast <- fast[order(fast$row_id), c("row_id", "eye_sim", "perm_sim", "eye_sim_diff")]
  manual <- manual[order(manual$row_id), , drop = FALSE]

  expect_equal(fast$row_id, manual$row_id)
  expect_equal(fast$eye_sim, manual$eye_sim, tolerance = 1e-10)
  expect_equal(fast$perm_sim, manual$perm_sim, tolerance = 1e-10)
  expect_equal(fast$eye_sim_diff, manual$eye_sim_diff, tolerance = 1e-10)
})

test_that("template_similarity cosine uses exhaustive within-stratum baseline when permutations exceed candidates", {
  make_density <- function(vals) {
    structure(
      list(z = matrix(vals, nrow = 2, ncol = 2, byrow = TRUE), x = 1:2, y = 1:2, sigma = 50),
      class = c("density", "eye_density")
    )
  }

  ref_tab <- tibble(
    image = 1:6,
    subject = rep(1:2, each = 3),
    density = list(
      make_density(c(1, 0, 0, 1)),
      make_density(c(0, 1, 1, 0)),
      make_density(c(1, 1, 0, 0)),
      make_density(c(0, 0, 1, 1)),
      make_density(c(1, 0, 1, 0)),
      make_density(c(0, 1, 0, 1))
    )
  )
  source_tab <- tibble(
    row_id = 1:6,
    image = 1:6,
    subject = rep(1:2, each = 3),
    density = list(
      make_density(c(0.9, 0.1, 0.1, 0.9)),
      make_density(c(0.1, 0.9, 0.9, 0.1)),
      make_density(c(0.9, 0.9, 0.1, 0.1)),
      make_density(c(0.1, 0.1, 0.9, 0.9)),
      make_density(c(0.9, 0.1, 0.9, 0.1)),
      make_density(c(0.1, 0.9, 0.1, 0.9))
    )
  )

  manual <- purrr::map_dfr(seq_len(nrow(source_tab)), function(i) {
    d2 <- source_tab$density[[i]]
    match_idx <- match(source_tab$image[[i]], ref_tab$image)
    candidates <- ref_tab %>%
      filter(subject == source_tab$subject[[i]], image != source_tab$image[[i]])

    perm_vals <- vapply(candidates$density, function(d1) similarity(d1, d2, method = "cosine"), numeric(1))
    eye_sim <- similarity(ref_tab$density[[match_idx]], d2, method = "cosine")
    perm_sim <- mean(perm_vals)

    tibble(
      row_id = source_tab$row_id[[i]],
      eye_sim = eye_sim,
      perm_sim = perm_sim,
      eye_sim_diff = eye_sim - perm_sim
    )
  })

  res <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "image",
    permute_on = "subject",
    refvar = "density",
    sourcevar = "density",
    method = "cosine",
    permutations = 99
  )

  res <- res[order(res$row_id), c("row_id", "eye_sim", "perm_sim", "eye_sim_diff")]
  manual <- manual[order(manual$row_id), , drop = FALSE]

  expect_equal(res$row_id, manual$row_id)
  expect_equal(res$eye_sim, manual$eye_sim, tolerance = 1e-10)
  expect_equal(res$perm_sim, manual$perm_sim, tolerance = 1e-10)
  expect_equal(res$eye_sim_diff, manual$eye_sim_diff, tolerance = 1e-10)
})

test_that("template_similarity cosine exhaustive permutation results are invariant to row order", {
  make_density <- function(vals) {
    structure(
      list(z = matrix(vals, nrow = 2, ncol = 2, byrow = TRUE), x = 1:2, y = 1:2, sigma = 50),
      class = c("density", "eye_density")
    )
  }

  ref_tab <- tibble(
    image = 1:6,
    subject = rep(1:2, each = 3),
    density = list(
      make_density(c(1, 2, 3, 4)),
      make_density(c(2, 3, 4, 5)),
      make_density(c(3, 4, 5, 6)),
      make_density(c(6, 5, 4, 3)),
      make_density(c(5, 4, 3, 2)),
      make_density(c(4, 3, 2, 1))
    )
  )
  source_tab <- tibble(
    row_id = 1:6,
    image = 1:6,
    subject = rep(1:2, each = 3),
    density = list(
      make_density(c(1.1, 2.1, 3.1, 4.1)),
      make_density(c(2.1, 3.1, 4.1, 5.1)),
      make_density(c(3.1, 4.1, 5.1, 6.1)),
      make_density(c(6.1, 5.1, 4.1, 3.1)),
      make_density(c(5.1, 4.1, 3.1, 2.1)),
      make_density(c(4.1, 3.1, 2.1, 1.1))
    )
  )
  shuffled_ref <- ref_tab[sample(seq_len(nrow(ref_tab))), , drop = FALSE]
  shuffled_source <- source_tab[sample(seq_len(nrow(source_tab))), , drop = FALSE]

  set.seed(7)
  ordered <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "image",
    permute_on = "subject",
    refvar = "density",
    sourcevar = "density",
    method = "cosine",
    permutations = 99
  )

  set.seed(7)
  shuffled <- template_similarity(
    shuffled_ref,
    shuffled_source,
    match_on = "image",
    permute_on = "subject",
    refvar = "density",
    sourcevar = "density",
    method = "cosine",
    permutations = 99
  )

  ordered <- ordered[order(ordered$row_id), c("row_id", "eye_sim", "perm_sim", "eye_sim_diff")]
  shuffled <- shuffled[order(shuffled$row_id), c("row_id", "eye_sim", "perm_sim", "eye_sim_diff")]

  expect_equal(ordered$row_id, shuffled$row_id)
  expect_equal(ordered$eye_sim, shuffled$eye_sim, tolerance = 1e-10)
  expect_equal(ordered$perm_sim, shuffled$perm_sim, tolerance = 1e-10)
  expect_equal(ordered$eye_sim_diff, shuffled$eye_sim_diff, tolerance = 1e-10)
})

test_that("template_similarity cosine preserves degenerate zero-vector behavior", {
  make_density <- function(vals) {
    structure(
      list(z = matrix(vals, nrow = 2, ncol = 2, byrow = TRUE), x = 1:2, y = 1:2, sigma = 50),
      class = c("density", "eye_density")
    )
  }

  ref_tab <- tibble(
    image = 1:4,
    density = list(
      make_density(c(0, 0, 0, 0)),
      make_density(c(1, 0, 0, 1)),
      make_density(c(0, 1, 1, 0)),
      make_density(c(1, 1, 1, 1))
    )
  )
  source_tab <- tibble(
    row_id = 1:4,
    image = 1:4,
    density = list(
      make_density(c(0, 0, 0, 0)),
      make_density(c(1, 0, 0, 1)),
      make_density(c(0, 1, 1, 0)),
      make_density(c(0, 0, 0, 0))
    )
  )

  res <- template_similarity(
    ref_tab,
    source_tab,
    match_on = "image",
    refvar = "density",
    sourcevar = "density",
    method = "cosine",
    permutations = 0
  )

  manual <- purrr::map_dbl(seq_len(nrow(source_tab)), function(i) {
    similarity(ref_tab$density[[i]], source_tab$density[[i]], method = "cosine")
  })

  expect_equal(res$eye_sim, manual, tolerance = 1e-10)
  expect_equal(res$eye_sim[[1]], 1)
  expect_equal(res$eye_sim[[4]], 0)
})

test_that("compute density with variable name other than 'fixgroup'", {
  g1 <- tibble(fg=lapply(1:100, function(i) {
    x <- runif(10)
    y <- runif(10)
    onset <- seq(1,length.out=length(x), by=50)
    duration <- rep(1,length(x))
    fixgroup <- fixation_group(x,y,onset,duration)
  }), image=1:100, subject=rep(1:10, each=10))

  dens <- density_by(g1, "image", keep_vars="subject", xbounds=c(0,1), ybounds=c(0,1), fixvar="fg")
  expect_true(!is.null(dens$fg))

})
