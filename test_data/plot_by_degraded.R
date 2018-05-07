library(dplyr)
library(tibble)
library(tidyr)



## load study data
pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()

pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))

## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))


## load test data
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/delay_fixations.csv")) %>%
  mutate(fix_onset=FixStartTime - DelayOnset) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()

pctest2 <- pctest %>% group_by(Subject, Trial) %>% slice(1)
## create table for each test trial
test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="fix_onset",
                      groupvar=c("Image", "Saliency"), data=pctest,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageNumber", "DelayOnset"))

m <- test_tab %>% group_by(Subject) %>% rowwise() %>% do( {
  sversion <- as.character(filter(pcstudy, Subject == .$Subject & ImageNumber == .$ImageNumber)$ImageSet[1])
  Match <- if (sversion == as.character(.$ImageSet)) "match"  else "mismatch"
  data.frame(Match=Match)
})

test_tab$Match <- m$Match


doplot <- function(inum, saliency=20, version="A", alpha=c(.2,.5)) {
  dfx <- test_tab %>% filter(ImageVersion == paste(inum, version) & Saliency==saliency)
  iname <- paste0("~/Dropbox/Jordana_experiments/Jordana_saliency_study/images/", paste0(inum, "_", version, "_1_", saliency, ".jpeg"))
  fg <- dfx$fixgroup[[1]]
  plot(fg, type="density", bg_image=iname, alpha=alpha)
}

## construct heatmaps for the study phase, averaged over subjects
study_dens <- density_by(study_tab, groups="Image", xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)


