library(dplyr)
library(tibble)
library(tidyr)

## load study data
start_time <- 4000
pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_study_input.csv")) %>% filter(Image != ".") %>% droplevels()
pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


## load test data
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delaytest_input.csv")) %>%
  filter(Image != "." & FixStartTime < start_time) %>% droplevels()

## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject", "ImageRepetition"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))

## create table for each test trial
test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject"), data=pctest,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion", "Saliency", "Accuracy",
                              "ImageSet", "Block", "Duration", "ImageNumber"))

m <- test_tab %>% group_by(Subject) %>% rowwise() %>% do( {
  sversion <- as.character(filter(pcstudy, Subject == .$Subject & ImageNumber == .$ImageNumber)$ImageSet[1])
  Match <- if (sversion == as.character(.$ImageSet)) "match"  else "mismatch"
  data.frame(Match=Match)
})

test_tab$Match <- m$Match

## construct heatmaps for the study phase, averaged within subjects
study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)
study_rep_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject", "ImageRepetition"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)

avg_rep_sim <- aggregate(eye_sim_diff ~ Subject + ImageNumber, data=study_rep_sim, FUN=mean)
avg_rep_sim$Image_Subj <- paste0(avg_rep_sim$ImageNumber, "_", avg_rep_sim$Subject)
test_tab$Image_Subj <- paste0(test_tab$ImageNumber, "_", test_tab$Subject)
avg_rep_sim <- inner_join(avg_rep_sim, test_tab %>% dplyr::select(-fixgroup), by="Image_Subj")
avg_rep_agg <- aggregate(eye_sim_diff ~ Accuracy  + Match + Duration, FUN=mean, data=avg_rep_sim)

ggplot(aes(Duration, eye_sim_diff, colour=factor(Accuracy)), data=avg_rep_agg) + facet_wrap(~ Match) + geom_line()
