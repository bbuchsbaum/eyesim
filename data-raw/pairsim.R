

ij1 <- read.table("data-raw/pair_confuseability.txt", header=TRUE) %>% arrange(ImageNumber)

## TODO add "window" arg to template_similarity and fixation_similarity
## what window maximizes correlation between memory and perception?
## what window maximizes correlation between split half perception/perception?
out <- wynn_study_image %>% ungroup() %>% nest_by(ImageNumber) %>% mutate(sim={
  similarity(data$fixgroup[[1]], data$fixgroup[[2]], tweight=0)
})

out <- out %>% arrange(as.integer(as.character(ImageNumber)))


outt1 <- wynn_study_image %>% ungroup() %>% nest_by(ImageNumber) %>% mutate(sim={
  similarity(data$fixgroup[[1]], data$fixgroup[[2]], tweight=1)
}) %>% arrange(as.integer(as.character(ImageNumber)))

outt1 <- outt1 %>% arrange(as.integer(as.character(ImageNumber)))

compute_tsim <- function(tweight=1, window=c(0,1000)) {
  outt1 <- wynn_study_image %>% ungroup() %>% nest_by(ImageNumber) %>% mutate(dsim={
    similarity(data$fixgroup[[1]], data$fixgroup[[2]], tweight=tweight, window=window)
  }) %>% arrange(as.integer(as.character(ImageNumber)))
}

compute_sim <- function(metric="pearson", window=c(0,1000), sigma=60) {

  outd1 <- wynn_study_image %>% ungroup() %>% nest_by(ImageNumber) %>% mutate(dsim={
    d1 <- eye_density(data$fixgroup[[1]], sigma=sigma, xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), window=window)
    d2 <- eye_density(data$fixgroup[[2]], sigma=sigma, xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), window=window)
    similarity.density(d1,d2, metric)
  }) %>% arrange(as.integer(as.character(ImageNumber)))
}

cvals <- unlist(lapply(seq(0, 2500, by=200), function(i) {
  print(i)
  s <- compute_sim(metric="spearman", window=c(i, i+300))
  cor(s$dsim, ij1$Accuracy)
}))

cvals2 <- unlist(lapply(seq(10,200, by=10), function(i) {
  print(i)
  s <- compute_sim(metric="spearman", sigma=i, window=c(300,2000))
  cor(s$dsim, ij1$Accuracy)
}))

cvals3 <- unlist(lapply(seq(0, 2500, by=200), function(i) {
  print(i)
  s <- compute_tsim(tweight=0, window=c(i, i+1000))
  cor(s$dsim, ij1$Accuracy)
}))

cvals4 <- do.call(rbind, lapply(seq(0, 10,by=1), function(i) {
  print(i)
  s <- compute_tsim(tweight=i, window=c(100,2500))
  cc=cor(s$dsim, ij1$Accuracy)
  cc
}))



