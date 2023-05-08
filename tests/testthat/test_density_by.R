library(testthat)
library(dplyr)

context("density_by")
options(future.rng.onMisue = "ignore")
test_that("density_by produces perfect similarity for identical patterns", {
  g1 <- tibble(fixgroup=lapply(1:100, function(i) {
    x <- runif(10)
    y <- runif(10)
    onset <- seq(1,length.out=length(x), by=50)
    duration <- rep(1,length(x))
    fixgroup <- fixation_group(x,y,onset,duration)
  }), image=1:100)

  dens <- density_by(g1, "image", xbounds=c(0,1), ybounds=c(0,1))
  dens2 <- density_by(g1, "image", xbounds=c(0,1), ybounds=c(0,1))
  tsim <- template_similarity(dens, dens2, match_on="image", method="spearman", permutations=30)
  expect_equal(tsim$eye_sim, rep(1,nrow(tsim)))

})


