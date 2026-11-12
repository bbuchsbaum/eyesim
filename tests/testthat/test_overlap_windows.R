test_that("overlap density equals manually clipped intervals and weights", {
  fg <- fixation_group(c(0, 1, 2, 3), c(1, 2, 0, 1), onset = c(-20, 20, 80, 100),
                       duration = c(40, 30, 50, 20))
  clipped <- fg[1:3, ]
  clipped$duration <- c(20, 30, 20)
  for (backend in c("ks", "MASS")) {
    calc <- function(x, ...) eye_density(x, sigma = 1, xbounds = c(0,3), ybounds = c(0,3),
                                         outdim = c(5,4), kde_pkg = backend, ...)
    expect_equal(calc(fg, window = c(0,100), window_weight = "overlap", duration_weighted = TRUE)$z,
                 calc(clipped, duration_weighted = TRUE)$z)
    expect_equal(calc(fg, window = c(0,100), window_weight = "overlap")$z,
                 calc(clipped, duration_weighted = TRUE)$z)
    unit <- expect_warning(calc(fg, window = c(0,100), window_weight = "overlap",
                                duration_weighted = FALSE), "unit weight")
    expect_equal(unit$z, calc(clipped)$z)
    expect_equal(calc(fg, window = c(0,100), window_weight = "overlap", weights = c(2,3,5,7))$z,
                 calc(clipped, weights = c(1,3,2))$z)
    expect_equal(calc(fg, window = c(0,100), duration_weighted = TRUE)$z,
                 calc(fg[2:3, ], duration_weighted = TRUE)$z)
  }
  expect_error(eye_density(fg, window = c(0,Inf), window_weight = "overlap"), "finite")
  expect_equal(eye_density(fg, window = c(-Inf,Inf))$z, eye_density(fg)$z)
  fg$duration[1] <- NA_real_
  expect_error(eye_density(fg, window = c(0,100), window_weight = "overlap"), "finite")
})

test_that("density_by forwards overlap weighting and zero overlap returns NULL", {
  fg <- fixation_group(c(0,1,2), c(0,1,0), onset = c(-20,20,80), duration = c(40,30,50))
  tab <- tibble::tibble(subject = "s", fixgroup = list(fg))
  args <- list(sigma = 1, xbounds = c(0,3), ybounds = c(0,3), outdim = c(5,4),
               window = c(0,100), window_weight = "overlap", duration_weighted = TRUE)
  direct <- do.call(eye_density, c(list(x = fg), args))
  grouped <- do.call(density_by, c(list(x = tab, groups = "subject"), args))
  expect_equal(grouped$density[[1]]$z, direct$z)
  expect_null(expect_warning(eye_density(fg, window = c(200,300), window_weight = "overlap"), "No fixations"))
})
