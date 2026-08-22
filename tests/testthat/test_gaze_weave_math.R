test_that("gaze measures normalize duration and time with declared edge policies", {
  fg <- fixation_group(
    x = c(2, 0, 1, 3),
    y = c(0, 0, 1, 1),
    duration = c(2, 1, 0, 1),
    onset = c(2, 0, 1, 4)
  )
  measure <- eyesim:::as_gaze_measure(fg, gaze_local_order(0.25))

  expect_equal(sum(measure$mass), 1)
  expect_equal(nrow(measure$coords), 3L)
  expect_true(all(diff(measure$fixations$onset) > 0))
  expect_true(all(measure$time >= 0 & measure$time <= 1))
  expect_equal(diag(measure$relation), rep(0, 3))
  expect_true(all(measure$relation[lower.tri(measure$relation)] == 0))
})

test_that("gaze measures reject ambiguous and invalid duration inputs", {
  expect_error(
    eyesim:::as_gaze_measure(
      fixation_group(c(0, 1), c(0, 1), c(1, 1), c(0, 0)),
      gaze_local_order()
    ),
    "unique"
  )
  expect_error(
    eyesim:::as_gaze_measure(
      fixation_group(c(0, 1), c(0, 1), c(0, 0), c(0, 1)),
      gaze_local_order()
    ),
    "positive"
  )
})

test_that("uniform temporal dilation leaves local chronology unchanged", {
  coords <- rbind(c(0, 0), c(1, 1), c(2, 0))
  original <- make_gaze_fixations(coords, c(1, 2, 1), c(0, 1, 3))
  dilated <- make_gaze_fixations(coords, c(3, 6, 3), c(0, 3, 9))
  chronology <- gaze_local_order(0.2)

  m1 <- eyesim:::as_gaze_measure(original, chronology)
  m2 <- eyesim:::as_gaze_measure(dilated, chronology)

  expect_equal(m1$mass, m2$mass, tolerance = 1e-14)
  expect_equal(m1$time, m2$time, tolerance = 1e-14)
  expect_equal(m1$relation, m2$relation, tolerance = 1e-14)
})


test_that("Gaussian mixture cost is stable, zero at identity, and increasing", {
  spatial <- gaze_gaussian_mixture(c(0.5, 1.5), c(0.25, 0.75))
  ref <- matrix(c(0, 0), nrow = 1)
  source <- rbind(c(0, 0), c(1, 0), c(1e5, 0))
  cost <- eyesim:::gaze_spatial_cost(ref, source, spatial)

  expect_equal(cost[[1]], 0, tolerance = 1e-14)
  expect_gt(cost[[2]], cost[[1]])
  expect_gt(cost[[3]], cost[[2]])
  expect_true(all(is.finite(cost)))
})
