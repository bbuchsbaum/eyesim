
context("centering")
options(future.rng.onMisue = "ignore")

test_that("can center a fixation_group", {
  x <- runif(10)
  y <- runif(10)
  onset <- seq(1,length.out=length(x), by=50)
  duration <- rep(1,length(x))
  fixgroup <- fixation_group(x,y,onset,duration)
  cfix <- center(fixgroup)
  expect_equivalent(mean(cfix$x), 0)
  expect_equivalent(mean(cfix$y), 0)
})
