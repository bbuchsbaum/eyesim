library(dplyr)
library(tibble)
library(tidyr)
library(imager)
library(mgcv)
library(ggplot2)
library(memoise)

## load study data
pcstudy <- as_tibble(read.csv("~/Dropbox/New_pc_behavioural_data/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()

pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))

## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject", "Block"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))

get_mask <- function(im) {
  fname <- paste0("~/Dropbox/Jordana_experiments/Jordana_saliency_study/images/Mat_", gsub("jpeg", "jpeg.rds", im))

  mask <- if (file.exists(fname)) {
    message("getting fname", fname)
    m <- readRDS(fname)
    m <- t(apply(m, 1, rev))
    m/sum(m)
  } else {
    mask <- matrix(1, 80,60)
    mask/sum(mask)
  }

  mask
}

#get_mask <- memoise(get_mask)


## load test data
pctest <- as_tibble(read.csv("~/Dropbox/New_pc_behavioural_data/testdelay_fixations.csv")) %>%
  mutate(fix_onset=FixOffset) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()


## create table for each test trial
test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                      groupvar=c("Image", "Subject"), data=pctest,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageNumber", "ImageRepetition"))
test_tab <- test_tab %>% filter(Subject < 300)
study_tab <- study_tab %>% filter(Subject < 300)


mask_tab <- tibble(Image=levels(test_tab$Image)) %>% rowwise() %>% mutate(mask=list(get_mask(Image))) %>% rowwise() %>% do({

  idx <- which(.$mask > 0, arr.ind=TRUE)
  tibble(Image=.$Image, fixgroup=list(fixation_group(idx[,1]*10, idx[,2]*10, duration=1, onset=0)))
})


## construct heatmaps for the study phase, averaged within subjects
study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                         duration_weighted=TRUE, sigma=60)






