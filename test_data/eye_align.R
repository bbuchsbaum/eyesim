library(dplyr)
library(tibble)
library(tidyr)
library(imager)
library(mgcv)
library(ggplot2)
library(memoise)


devtools::load_all()

exclude_subs <- c(28, 32, 109)

## load testdelay fixation data
pcdelay <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/delay_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% exclude_subs ) & Subject < 300) %>% droplevels()

delay_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                      groupvar=c("Subject", "ImageVersion"), data=pcdelay,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber"))

pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% exclude_subs) & Subject < 300) %>% droplevels() %>% mutate(ImageNumber=as.integer(as.character(ImageNumber)))


odim <- as.integer(c(80,60)/8)
## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Subject", "ImageVersion"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))


## construct heatmaps for the study phase, averaged within subjects
study_dens <- density_by(study_tab, groups=c("Subject", "ImageVersion"),
                         xbounds=c(0,800), ybounds=c(0,600), outdim=odim,
                         duration_weighted=TRUE, sigma=40, result_name="study_density")


## construct heatmaps for the study phase, averaged within subjects
delay_dens <- density_by(delay_tab, groups=c("Subject", "ImageVersion"),
                         xbounds=c(0,800), ybounds=c(0,600), outdim=odim,
                         duration_weighted=TRUE, sigma=40, result_name="delay_density")

study_grouped <- study_dens %>% group_by(Subject) %>% dplyr::group_split()
delay_grouped <- delay_dens %>% group_by(Subject) %>% dplyr::group_split()

study_blocks <- study_grouped %>% purrr::map(function(g) {
  out <- do.call(rbind, g$study_density %>% purrr::map(~ as.vector(.$z)))
  row.names(out) <- g$ImageVersion
  out
})

delay_blocks <- delay_grouped %>% purrr::map(function(g) {
  out <- do.call(rbind, g$delay_density %>% purrr::map(~ as.vector(.$z)))
  row.names(out) <- g$ImageVersion
  out
})

#labels <- c(study_dens$ImageVersion, delay_dens$ImageVersion)
labels <- unlist(lapply(study_blocks, row.names))
Sl <- neighborweights::label_matrix2(labels, labels)
Sl <- Sl/RSpectra::eigs(Sl,1)$values[1]
#X <- Matrix::bdiag(c(study_blocks, delay_blocks))
X <- Matrix::bdiag(study_blocks)

lds_to_mat <- function(lvec) {
  n <- length(study_blocks) + length(delay_blocks)
  aa <- array(lvec, c(prod(odim), n))
  m1 <- rowMeans(aa[,1:length(study_blocks)])
  #m2 <- rowMeans(aa[,(length(study_blocks)+1):(n)])

  #list(study=matrix(m1, odim[1], odim[2]), delay=matrix(m2, odim[1], odim[2]))
  list(study=matrix(m1, odim[1], odim[2]))
}

gres <- neuroca::pca(X, ncomp=8, preproc=neuroca::pass())


