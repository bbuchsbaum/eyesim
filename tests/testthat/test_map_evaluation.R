test_that("lookup preparation and multiple-map adapters agree", {
  maps <- reference_fixture()$density[1:2]
  names(maps) <- c("one", "two")
  fix <- data.frame(x = c(0, 1, 2), y = c(0, 1, 0), onset = 1:3)
  wide <- sample_density(maps, fix, normalize = "sum", transform = "log", floor = .01, format = "wide")
  long <- sample_density(maps, fix, normalize = "sum", transform = "log", floor = .01)
  for (nm in names(maps)) {
    z <- maps[[nm]]$z / sum(maps[[nm]]$z)
    expected <- log(z[cbind(c(1, 2, 3), c(1, 2, 1))] + .01)
    expect_equal(wide[[nm]], expected)
    expect_equal(long$z[long$reference == nm], expected)
  }
  tab <- tibble::tibble(reference = names(maps), density = unname(maps))
  expect_equal(sample_density(tab, fix, normalize = "sum", transform = "log", floor = .01), long)
  expect_error(sample_density(maps[[1]], fix, transform = "log"), "floor")
  expect_error(sample_density(maps[[1]], fix, normalize = "zscore", transform = "log", floor = .01), "zscore")
  expect_error(sample_density(unname(maps), fix), "named list")
})

test_that("time sampling and template adapter share log preparation", {
  tab <- reference_fixture()
  tab$key <- seq_len(nrow(tab))
  fg <- fixation_group(c(0, 1), c(0, 1), onset = c(0, 100), duration = c(100, 100))
  source <- data.frame(key = 1:9)
  source$fixgroup <- rep(list(fg), 9)
  set.seed(91)
  fast <- sample_density_time(tab, source, "key", times = c(0, 100), permutations = 2,
                              normalize = "sum", transform = "log", floor = .01)
  for (i in 1:9) expect_equal(fast$sampled[[i]], sample_density(tab$density[[i]], fg, c(0,100),
                                                             normalize = "sum", transform = "log", floor = .01))
  source$template <- tab$density
  expect_equal(template_sample(source, "template", time = c(0,100), normalize = "sum",
                                transform = "log", floor = .01)$sample_out, fast$sampled)
  testthat::local_mocked_bindings(density_lookup_maps = function(dens_list, ...) {
    list(x = NULL, y = NULL, z = vector("list", length(dens_list)))
  })
  set.seed(91)
  slow <- sample_density_time(tab, source[c("key", "fixgroup")], "key", times = c(0,100),
                              permutations = 2, normalize = "sum", transform = "log", floor = .01)
  expect_equal(slow, fast)
})

test_that("candidate sampling handles singleton integer pools and preserves legacy order", {
  set.seed(77)
  expected <- sample(c(7L, 2L, 9L), 2)
  expected_rng <- .Random.seed
  set.seed(77)
  expect_equal(draw_candidates(c(7L, 2L, 9L), 2), expected)
  expect_identical(.Random.seed, expected_rng)
  expect_equal(draw_candidates(7L, 20, replace = TRUE, cap = FALSE), rep(7L, 20))
  expect_equal(reference_candidates(c(1L, 2L, 1L, 3L, 2L), 1L), c(2L, 3L))
})

test_that("legacy S3 sampling extensions receive the original default arguments", {
  registry <- get(".__S3MethodsTable__.", envir = asNamespace("eyesim"))
  assign("sample_density.strict_test", function(x, fix, times = NULL, normalize = "none") {
    data.frame(z = rep(x$value, length(times)), time = times)
  }, envir = registry)
  on.exit(rm("sample_density.strict_test", envir = registry))
  refs <- tibble::tibble(key = 1:2, density = list(structure(list(value = 3), class = "strict_test"),
                                                structure(list(value = 7), class = "strict_test")))
  fg <- fixation_group(c(0,1),c(0,1),onset=c(0,100),duration=c(100,100))
  src <- tibble::tibble(key = 1:2, fixgroup = list(fg,fg))
  out <- expect_no_warning(sample_density_time(refs,src,"key",times=c(0,100),permutations=1))
  expect_equal(out$sampled[[1]]$z,c(3,3))
  expect_equal(out$perm_sampled[[1]]$z,c(7,7))
})
