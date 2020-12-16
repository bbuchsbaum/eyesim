
context("centering")
test_that("can center an eye_frame", {
  g1 <- tibble(fixgroup=lapply(1:100, function(i) {
    x <- runif(10)
    y <- runif(10)
    onset <- seq(1,length.out=length(x), by=50)
    duration <- rep(1,length(x))
    fixgroup <- fixation_group(x,y,onset,duration)
  }), image=1:100, subject=rep(1:10, each=10))


  g2 <- g1 %>% mutate(cen_fixgroup=list(center(fixgroup[[1]])))

})
