library(dplyr)
library(tibble)
library(tidyr)
library(energy)


pcstudy <- as_tibble(read.csv("~/Dropbox/New_pc_behavioural_data/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()

pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject", "Block"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))
study_tab <- study_tab %>% mutate(dummy="global")

subject_dens <- density_by(study_tab, groups=c("dummy"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                           duration_weighted=TRUE, sigma=60)

im22_dens <- study_tab %>% filter(ImageVersion == "96 B") %>% density_by(., groups=c("ImageSet"),
                                                                         xbounds=c(0,800), ybounds=c(0,600),
                                                                         outdim=c(80,60),
                                                                         duration_weighted=TRUE, sigma=60)

im22_dens_test <- test_tab %>% filter(ImageNumber == "96") %>% density_by(., groups=c("ImageSet"),
                                                                         xbounds=c(0,800), ybounds=c(0,600),
                                                                         outdim=c(80,60),
                                                                         duration_weighted=TRUE, sigma=60)


## load test data
pctest <- as_tibble(read.csv("~/Dropbox/New_pc_behavioural_data/testdelay_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()


## create table for each test trial
test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                      groupvar=c("Image", "Subject"), data=pctest,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageNumber"))


binned_density <- function(min_onset, max_onset, iversion) {

  print(min_onset)
  ## create table for each test trial
  pctest_binned <- pctest %>% filter(FixOffset >= min_onset & FixOffset < max_onset & ImageVersion == iversion)


  test_tab_binned <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                               groupvar=c("Image", "Subject"), data=pctest_binned,
                               clip_bounds=c(112, (112+800), 684, 84),
                               vars=c("ImageVersion", "Saliency", "Accuracy",
                                      "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber"))


  test_dens <- density_by(test_tab_binned, groups=c("ImageNumber"),xbounds=c(0,800), ybounds=c(0,600),
                          outdim=c(80,60), duration_weighted=TRUE, sigma=50)
  test_dens %>% mutate(time=(min_onset + max_onset)/2)

}

res <- lapply(seq(0, 3500, by=25), function(i) {
  bd <- binned_density(i, i+200, "96 B")
})

library(gifski)
td <- tempdir()
for (i in 1:length(res)) {
  png(filename=paste0(td, "/im_", i, ".png"))
  image(res[[i]]$density[[1]]$z^(1/4))
  dev.off()
}

fnames <- paste0(td, "/im_", 1:length(res), ".png")
gifski(fnames, "animation.gif", delay=.05, loop=FALSE)
