density_delta_source <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    gaze_weave_test_inst_path("validation", "gaze-weave-density-delta.R"),
    envir = environment
  )
  environment
}

density_delta_numeric_mass <- function(court, template, h, step = 2) {
  screen <- c(width = 800, height = 600)
  gx <- seq(step / 2, screen[["width"]] - step / 2, by = step)
  gy <- seq(step / 2, screen[["height"]] - step / 2, by = step)
  grid <- expand.grid(x = gx, y = gy)
  density <- court$density_delta_template_density(
    template, grid$x, grid$y, h, screen
  )
  sum(density) * step^2
}

test_that("screen-normalised kernels integrate to one on the screen", {
  court <- density_delta_source()
  for (centre in list(c(0, 0), c(400, 3), c(795, 598), c(400, 300))) {
    template <- list(x = centre[[1L]], y = centre[[2L]], w = 1)
    for (h in c(20, 80)) {
      expect_equal(density_delta_numeric_mass(court, template, h), 1,
                   tolerance = 2e-3)
    }
  }
  # An edge kernel is raised relative to the plane-normalised Gaussian by
  # exactly the reciprocal of its on-screen mass.
  mass <- court$density_delta_kernel_mass(0, 300, 40, c(width = 800, height = 600))
  expect_equal(mass, 0.5, tolerance = 1e-9)
  edge <- court$density_delta_template_density(
    list(x = 0, y = 300, w = 1), 0, 300, 40, c(width = 800, height = 600)
  )
  expect_equal(drop(edge), 2 / (2 * pi * 40^2), tolerance = 1e-9)
})

test_that("episodes are normalised before averaging", {
  court <- density_delta_source()
  long <- data.frame(x = c(100, 200), y = c(100, 100), duration = c(900, 900))
  short <- data.frame(x = 600, y = 500, duration = 100)
  template <- court$density_delta_episode_template(list(long, short))
  expect_equal(sum(template$w), 1)
  expect_equal(template$w, c(0.25, 0.25, 0.5))
  pooled <- court$density_delta_pool_templates(list(template, list(
    x = 1, y = 1, w = 7
  )))
  expect_equal(sum(pooled$w), 1)
  expect_equal(pooled$w[[4L]], 0.5)
  mixed <- density_delta_numeric_mass(court, template, 30)
  expect_equal(mixed, 1, tolerance = 2e-3)
})

test_that("mixture EM recovers known population weights", {
  court <- density_delta_source()
  set.seed(11)
  n <- 4000
  truth <- c(0.6, 0.3, 0.1)
  component <- sample.int(3L, n, replace = TRUE, prob = truth)
  x <- ifelse(component == 1L, stats::rnorm(n, -2), ifelse(
    component == 2L, stats::rnorm(n, 2), stats::runif(n, -10, 10)
  ))
  F <- cbind(stats::dnorm(x, -2), stats::dnorm(x, 2), 1 / 20)
  fit <- court$density_delta_em(F, rep(1 / n, n), u = 1 / 20, floor = 0,
                                max_iter = 5000L, tol = 1e-12)
  expect_true(fit$converged)
  expect_equal(fit$weights, truth, tolerance = 0.05)
})

test_that("delta is near zero when gaze has no own contribution", {
  court <- density_delta_source()
  config <- court$density_delta_config()
  sim <- court$density_delta_simulate(
    participants = 12L, items = 24L, strength = 0, overlap = "high", seed = 3L
  )
  folds <- court$density_delta_synthetic_folds(sim$trials, 3L)
  own <- court$density_delta_crossfit(
    court$density_delta_synthetic_design(sim, config), folds
  )
  pseudo <- court$density_delta_crossfit(
    court$density_delta_synthetic_design(sim, config, "pseudo"), folds
  )
  expect_setequal(own$scores$trial, seq_len(nrow(sim$trials)))
  expect_lt(abs(mean(own$scores$delta)), config$null_delta_tolerance)
  contrast <- court$density_delta_crossed_test(
    own$scores, own$scores$delta - pseudo$scores$delta, "own_minus_pseudo",
    draws = 199L, seed = 1L
  )
  expect_lt(abs(contrast$estimate), config$null_delta_tolerance)
  expect_false(contrast$reject_one_sided)
})

test_that("an injected own contribution yields a positive held-out delta", {
  court <- density_delta_source()
  config <- court$density_delta_config()
  sim <- court$density_delta_simulate(
    participants = 12L, items = 24L, strength = 0.3, overlap = "low", seed = 4L
  )
  folds <- court$density_delta_synthetic_folds(sim$trials, 4L)
  own <- court$density_delta_crossfit(
    court$density_delta_synthetic_design(sim, config), folds
  )
  test <- court$density_delta_crossed_test(
    own$scores, own$scores$delta, "delta", draws = 199L, seed = 1L
  )
  expect_gt(test$estimate, 0.05)
  expect_true(test$reject_one_sided)
  expect_true(all(own$fits$w1_own > 0.1))
})

test_that("fold fits use training rows only", {
  court <- density_delta_source()
  config <- court$density_delta_config()
  sim <- court$density_delta_simulate(
    participants = 12L, items = 24L, strength = 0.2, seed = 5L
  )
  folds <- court$density_delta_synthetic_folds(sim$trials, 5L)
  fold <- folds[[1L]]
  held_out <- which(sim$trials$participant %in% fold$eval_participants &
                      sim$trials$item %in% fold$eval_items)
  train <- which(!sim$trials$participant %in% fold$eval_participants &
                   !sim$trials$item %in% fold$eval_items)
  design <- court$density_delta_synthetic_design(sim, config)
  fit <- court$density_delta_fit(design, train)

  perturbed <- sim
  set.seed(99)
  for (t in held_out) {
    n <- nrow(perturbed$y[[t]])
    perturbed$y[[t]]$x <- stats::runif(n, 0, 800)
    perturbed$y[[t]]$y <- stats::runif(n, 0, 600)
  }
  perturbed_design <- court$density_delta_synthetic_design(perturbed, config)
  perturbed_fit <- court$density_delta_fit(perturbed_design, train)

  expect_identical(perturbed_fit, fit)
  before <- court$density_delta_score(design, fit, held_out)
  after <- court$density_delta_score(perturbed_design, perturbed_fit, held_out)
  expect_false(isTRUE(all.equal(before$delta, after$delta)))
})

test_that("crossed inference leaves the caller's RNG untouched", {
  court <- density_delta_source()
  tab <- expand.grid(participant = letters[1:6], item = 1:8)
  value <- stats::rnorm(nrow(tab))
  set.seed(42)
  expected <- stats::runif(3)
  set.seed(42)
  court$density_delta_crossed_test(tab, value, "x", draws = 99L, seed = 7L)
  expect_identical(stats::runif(3), expected)
})
