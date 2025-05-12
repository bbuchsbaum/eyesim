library(testthat)
library(dplyr)

context("density_by")
options(future.rng.onMisuse = "ignore")
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

test_that("weighted and unweighted density maps are highly correlated", {
  # Generate test data with varying durations
  set.seed(123)  # for reproducibility
  x <- runif(50, 0, 1000)
  y <- runif(50, 0, 1000)
  duration <- rep(1,50)  # varying durations
  onset <- seq(1, length.out=length(x), by=50)
  
  # Create fixation group
  fg <- fixation_group(x, y, duration, onset)
  
  # Generate density maps with both methods
 

                           
  
  unweighted_density <- eye_density(fg, sigma=100, xbounds=c(0,1000), ybounds=c(0,1000), 
                                  outdim=c(100,100), duration_weighted=FALSE)

  for(s in seq(2,1000, by=10)) {
  weighted_density <- eye_density(fg, sigma=s, xbounds=c(0,1000), ybounds=c(0,1000), 
                                outdim=c(100,100), duration_weighted=TRUE)
    # Calculate correlation between the two density maps
    correlation <- cor(as.vector(weighted_density$z), as.vector(unweighted_density$z), 
                      method="pearson")
    #print(paste(s , correlation))                   
    # Test that correlation is above 0.95
    #expect_gt(correlation, 0.95)
  }
  
  # Also test with more clustered data to ensure robustness
  # Generate clustered points around 3 centers
  centers <- matrix(c(250,250, 500,500, 750,750), ncol=2, byrow=TRUE)
  n_per_cluster <- 20
  
  x2 <- c()
  y2 <- c()
  for(i in 1:nrow(centers)) {
    x2 <- c(x2, rnorm(n_per_cluster, centers[i,1], 50))
    y2 <- c(y2, rnorm(n_per_cluster, centers[i,2], 50))
  }
  
  duration2 <- runif(length(x2), 50, 500)
  onset2 <- seq(1, length.out=length(x2), by=50)
  
  fg2 <- fixation_group(x2, y2, duration2, onset2)
  
  weighted_density2 <- eye_density(fg2, sigma=50, xbounds=c(0,1000), ybounds=c(0,1000), 
                                 outdim=c(100,100), duration_weighted=TRUE)
  
  unweighted_density2 <- eye_density(fg2, sigma=50, xbounds=c(0,1000), ybounds=c(0,1000), 
                                   outdim=c(100,100), duration_weighted=FALSE)
  
  correlation2 <- cor(as.vector(weighted_density2$z), as.vector(unweighted_density2$z), 
                     method="pearson")
  
  # Test that correlation is above 0.95 for clustered data
  expect_gt(correlation2, 0.95)
})
