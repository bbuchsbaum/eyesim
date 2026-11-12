make_issue_map <- function(z = 1:4) {
  d <- gen_density(1:2, 1:2, matrix(z, 2))
  class(d) <- c("eye_density", class(d))
  d
}

test_that("ambiguous template keys warn and compound keys match participants (#31)", {
  a <- make_issue_map()
  b <- make_issue_map(4:1)
  tab <- tibble::tibble(participant = c("s1", "s2"), image = "i1",
                        key = c("s1_i1", "s2_i1"), density = list(a, b))
  expect_warning(template_similarity(tab, tab, match_on = "image", method = "pearson",
                                     permutations = 0), "not unique")
  out <- expect_no_warning(template_similarity(tab, tab, match_on = "key", method = "pearson",
                                               permutations = 0))
  expect_equal(out$eye_sim, c(1, 1))
})

test_that("splitting participants and matching images selects the intended repeats (#28)", {
  a <- make_issue_map()
  b <- make_issue_map(4:1)
  tab <- tibble::tibble(participant = rep(c("s1", "s2"), each = 4),
                        image = rep(c("a", "b"), 4),
                        phase = rep(rep(c("enc", "ret"), each = 2), 2),
                        density = list(a, b, a, b, b, a, b, a))
  out <- dplyr::bind_rows(lapply(split(tab, tab$participant), function(x) {
    repetitive_similarity(x, condition_var = "image", method = "pearson")
  }))
  expect_equal(out$repsim, rep(1, 8))
  expect_equal(out$othersim, rep(-1, 8))
})

test_that("density arithmetic distinguishes contrasts even when all cells are positive (#32)", {
  a <- make_issue_map(2:5)
  b <- make_issue_map(rep(1, 4))
  a$sigma <- b$sigma <- 10
  for (contrast in list(a - b, a / b)) {
    expect_s3_class(contrast, "eye_contrast")
    expect_false(inherits(contrast, "density"))
    expect_equal(contrast$sigma, 10)
    expect_error(fixation_entropy(contrast), "No fixation_entropy")
    expect_s3_class(plot(contrast), "ggplot")
  }
  avg <- a + b
  expect_s3_class(avg, "eye_density")
  expect_equal(avg$z, (a$z + b$z) / 2)
  expect_true(is.finite(fixation_entropy(avg)))
  expect_equal(avg$sigma, 10)
  expect_equal((a / b)$z, log(a$z / b$z))
  b$sigma <- 20
  expect_error(a + b, "bandwidth")
  b <- make_issue_map()
  b$x <- c(1, 3)
  expect_error(a - b, "grid")
})

test_that("density estimates record the actual kernel and refuse implicit fallback (#30)", {
  fg <- fixation_group(x = c(10, 30, 60), y = c(20, 40, 30),
                        onset = c(0, 100, 200), duration = c(100, 150, 200))
  for (backend in c("ks", "MASS")) {
    d <- eye_density(fg, sigma = 10, xbounds = c(0, 100), ybounds = c(0, 100),
                     outdim = c(10, 10), kde_pkg = backend)
    expect_equal(d$estimator, backend)
    expect_equal(d$kernel_sd, if (backend == "ks") 10 else 2.5)
  }
  # Exercise the missing-ks branch without altering the installed dependency.
  f <- eyesim:::.compute_single_eye_density
  env <- new.env(parent = environment(f))
  env$requireNamespace <- function(...) FALSE
  environment(f) <- env
  expect_error(f(fg, 10, c(0, 100), c(0, 100), c(10, 10),
                 TRUE, TRUE, fg$duration, as.matrix(fg[c("x", "y")]), "ks"),
               "no automatic")
})
