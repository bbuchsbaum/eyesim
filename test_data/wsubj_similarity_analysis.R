library(dplyr)
library(tibble)
library(tidyr)

## load study data
start_time <- 2000
pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_study_input.csv")) %>% filter(Image != ".") %>% droplevels()
pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


## load test data
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delay_input.csv")) %>%
  filter(Image != "." & FixStartTime > Duration & FixStartTime < start_time) %>% droplevels()

## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion", "ImageRepetition",
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

## construct heatmaps for the study phase, averaged over subjects
study_dens_all <- density_by(study_tab, groups=c("ImageVersion"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)

## compute grand mean density
study_dens_avg <- Reduce("+", lapply(study_dens$density, function(x) x$z))/length(study_dens)
study_dens_avg <- study_dens_avg/sum(study_dens_avg)

dens_avg_subj <- study_dens %>% group_by(Subject) %>% do({
  dens <- Reduce("+", lapply(.$density, function(x) x$z))/nrow(.)
  tibble(Subject=.$Subject[1], dens=list(dens))
})

test_dens <- density_by(test_tab, groups=c("ImageNumber", "Subject"),xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)
test_dens_all <- density_by(test_tab, groups=c("ImageVersion", "Subject"),xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)

test_dens$Image_Subj <- paste0(test_dens$Subject, "_", test_dens$ImageNumber)
study_dens$Image_Subj <- paste0(study_dens$Subject, "_", study_dens$ImageNumber)
test_tab$Image_Subj <- paste0(test_tab$Subject, "_", test_tab$ImageNumber)

## compute similarity between each trial and study image derived from group average
test_sim <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=10)
test_sim_all <- template_similarity(study_dens_all, test_dens_all, "ImageVersion", permutations=10)

test_sim_all$Image_Subj_Version <- paste0(test_sim_all$Subject, "_", test_sim_all$ImageVersion)
test_tab$Image_Subj_Version <- paste0(test_tab$Subject, "_", test_tab$ImageVersion)

## join eye_sim with test_tab
test_sim <- inner_join(test_sim, test_tab, by="Image_Subj")

## compute means over variables
test_sim_agg <- aggregate(eye_sim_diff ~ Saliency + Duration, data=test_sim, FUN=mean)

# aggregate over subjects
#subj_acc <- aggregate(Accuracy ~ Subject.x + Saliency + Duration, data=test_sim, FUN=mean)
#subj_sim <- aggregate(eye_sim_diff ~ Subject.x + Saliency + Duration, data=test_sim, FUN=mean)
#subj_acc$sim <- subj_sim$eye_sim_diff

test_sim_all <- inner_join(test_sim_all, test_tab, by="Image_Subj_Version")
test_sim_all_agg <- aggregate(eye_sim_diff ~ Saliency + Duration + Accuracy + Match, data=test_sim_all, FUN=mean)
test_sim_all_agg <- aggregate(eye_sim_diff ~ Saliency + Duration, data=test_sim_all, FUN=mean)

# aggregate over subjects
#subj_acc_all <- aggregate(Accuracy ~ Subject.x + Saliency + Duration, data=subset(test_sim_all), FUN=mean)
#subj_sim_all <- aggregate(eye_sim_diff ~ Subject.x + Saliency + Duration, data=subset(test_sim_all), FUN=mean)
#subj_acc_all$sim <- subj_sim_all$eye_sim_diff

library(ggplot2)

ggplot(aes(Saliency, eye_sim_diff, colour=factor(Duration)), data=test_sim_agg) + geom_line()
ggplot(aes(Saliency, eye_sim_diff, linetype=factor(Accuracy), colour=factor(Duration)), data=test_sim_all_agg) +  facet_wrap( ~ Match) + geom_line()
ggplot(aes(Saliency, eye_sim_diff, linetype=factor(Accuracy), colour=factor(Duration)), data=test_sim_agg) +  facet_wrap( ~ Match) + geom_line()

