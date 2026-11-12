test_that("references equal independent membership averages on incomplete designs", {
  set.seed(53)
  for (iteration in 1:12) {
    tab <- reference_fixture()
    if (iteration > 1) tab <- tab[-sample.int(nrow(tab), 2), ]
    if (iteration > 6) tab$density[sample.int(nrow(tab), 1)] <- list(NULL)
    for (space in c("density", "log")) {
      refs <- suppressWarnings(reference_maps(tab, space = space, floor = .01))
      for (i in seq_len(nrow(tab))) for (ref in refs$refs) {
        expected <- reference_oracle(tab, i, ref, space)
        actual <- reference_map(refs, i, ref)
        if (is.null(expected)) expect_null(actual) else {
          expect_equal(actual$z, expected, tolerance = 1e-12)
          points <- expand.grid(x = 0:2, y = 0:1)
          points$subject <- tab$subject[i]
          points$item <- tab$item[i]
          expect_equal(sample_density(refs, points)[[ref]], as.vector(expected), tolerance = 1e-12)
        }
      }
    }
  }
})

test_that("leave-out references cannot see the excluded source", {
  tab <- reference_fixture()
  before <- reference_maps(tab, space = "density")
  tab$density[[1]]$z <- matrix(c(100, rep(0, 5)), 3)
  after <- reference_maps(tab, space = "density")
  for (ref in setdiff(before$refs, "self_item")) {
    expect_equal(reference_map(before, 1, ref)$z, reference_map(after, 1, ref)$z, tolerance = 1e-12)
  }
  expect_false(isTRUE(all.equal(reference_map(before, 1)$z, reference_map(after, 1)$z)))
})

test_that("reference domains, counts and missingness are explicit", {
  tab <- reference_fixture()
  tab$density[1] <- list(NULL)
  refs <- expect_warning(reference_maps(tab, floor = .01), "min_n")
  expect_equal(refs$counts[1, ], data.frame(self_item = 0L, self_other_items = 2L,
                                          others_item = 2L, others_other_items = 2L))
  expect_null(reference_map(refs, 1))
  map <- reference_map(refs, 2)
  expect_identical(map$domain, "log")
  expect_false(inherits(map, "density"))
  fix <- data.frame(x = 0, y = 0, onset = 0)
  expect_equal(sample_density(map, fix)$z, map$z[1, 1])
  expect_error(sample_density(map, fix, transform = "log", floor = .01), "already prepared")
  expect_error(reference_maps(reference_fixture()), "explicit positive")
  expect_error(reference_maps(rbind(tab, tab[1, ]), floor = .01), "unique")
  tab$density[[2]]$x <- c(1, 2, 3)
  expect_error(reference_maps(tab, floor = .01), "different lattices")
})

test_that("log before averaging is not log after averaging", {
  tab <- reference_fixture()
  density <- reference_maps(tab, space = "density")
  logrefs <- reference_maps(tab, floor = .01)
  log_after <- reference_maps(tab, space = "density", transform = "log", floor = .01)
  expect_equal(reference_map(log_after, 1, "self_other_items")$z,
               log(reference_map(density, 1, "self_other_items")$z + .01))
  expect_error(reference_maps(tab, transform = "log", floor = .01), "already log")
  expect_false(isTRUE(all.equal(reference_map(logrefs, 1, "self_other_items")$z,
                               log(reference_map(density, 1, "self_other_items")$z + .01))))
})

test_that("keys are matched without delimiter collisions and row-order dependence", {
  tab <- reference_fixture()
  tab$subject <- rep(c("s::1", "s", "s:"), 3)
  tab$item <- rep(c("a", ":a", "::a"), each = 3)
  ref <- reference_maps(tab, floor = .01)
  fix <- tab[c("subject", "item")]
  fix$x <- 0; fix$y <- 1
  scores <- sample_density(ref, fix)
  reverse <- reference_maps(tab[9:1, ], floor = .01)
  expect_equal(sample_density(reverse, fix), scores, tolerance = 1e-12)
  expect_equal(nrow(sample_density(ref, fix, format = "long")), 36)
  fix$item[1] <- "missing"
  expect_error(sample_density(ref, fix), "no matching")
})
