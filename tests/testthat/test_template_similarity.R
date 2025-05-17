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


