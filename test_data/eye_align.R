library(dplyr)
library(tibble)
library(tidyr)
library(imager)
library(mgcv)
library(ggplot2)
library(memoise)

exclude_subs <- c(28, 32, 109)

## load testdelay fixation data
pcdelay <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/delay_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% exclude_subs ) & Subject < 300) %>% droplevels()

delay_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                      groupvar=c("Subject", "ImageVersion"), data=pctest,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber", "testdelayOnset"))

pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% exclude_subs) & Subject < 300) %>% droplevels()

pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))


## construct heatmaps for the study phase, averaged within subjects
study_dens <- density_by(study_tab, groups=c("ImageVersion", "Subject"),
                         xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                         duration_weighted=TRUE, sigma=80, result_name="study_density")


## construct heatmaps for the study phase, averaged within subjects
delay_dens <- density_by(delay_tab, groups=c("ImageVersion", "Subject"),
                         xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                         duration_weighted=TRUE, sigma=80, result_name="study_density")

