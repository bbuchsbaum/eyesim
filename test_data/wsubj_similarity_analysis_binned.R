library(dplyr)
library(tibble)
library(tidyr)

exclude_subs <- c(28, 32, 109)

prepare_study <- function() {
  pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% exclude_subs) & Subject < 300) %>% droplevels()

  pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


  ## create table for each study trial (Subject/Image)
  study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject", "Block"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))


  subject_dens <- density_by(study_tab, groups=c("Subject"), xbounds=c(0,800), ybounds=c(0,600),
                           outdim=c(80,60), duration_weighted=TRUE, sigma=80)


  ## construct heatmaps for the study phase, averaged within subjects
  study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                         duration_weighted=TRUE, sigma=80)


  study_dens_subj_avg <- density_by(study_tab, groups=c("Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                                  duration_weighted=TRUE, sigma=80)

  ## construct heatmaps for the study phase, averaged over subjects
  study_dens_all <- density_by(study_tab, groups=c("ImageVersion"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)

  ## compute grand mean density
  study_dens_avg <- Reduce("+", lapply(study_dens$density, function(x) x$z))/length(study_dens)
  study_dens_avg <- study_dens_avg/sum(study_dens_avg)
  study_dens$Image_Subj <- paste0(study_dens$Subject, "_", study_dens$ImageNumber)

  dens_avg_subj <- study_dens %>% group_by(Subject) %>% do({
    dens <- Reduce("+", lapply(.$density, function(x) x$z))/nrow(.)
    tibble(Subject=.$Subject[1], dens=list(dens))
  })

  list(
    study_dens=study_dens,
    study_dens_subj_avg=study_dens_subj_avg,
    study_dens_all=study_dens_all,
    study_dens_avg=study_dens_avg,
    study_dens_avg_subj=dens_avg_subj
  )

}

stud <- prepare_study()




## load testdelay fixation data
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/testdelay_fixations.csv")) %>%
  #mutate(fix_onset=FixStartTime - testdelayOnset) %>%
  filter(Image != "." & !(Subject %in% exclude_subs) & Subject < 300) %>% droplevels()

binned_similarity <- function(min_onset, max_onset) {
  ## create table for each test trial
  pctest_binned <- pctest %>% filter(FixOffset >= min_onset & FixOffset < max_onset)



  test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                        groupvar=c("Image", "Subject"), data=pctest_binned,
                        clip_bounds=c(112, (112+800), 684, 84),
                        vars=c("ImageVersion", "Saliency", "Accuracy",
                               "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber", "testdelayOnset"))

  test_tab$Image_Subj <- paste0(test_tab$Subject, "_", test_tab$ImageNumber)


  test_tab$Match <- test_tab$ImageRepetition

  test_dens <- density_by(test_tab, groups=c("ImageNumber", "Subject"),xbounds=c(0,800), ybounds=c(0,600),
                          outdim=c(80,60), duration_weighted=TRUE, sigma=80)
  test_dens$Image_Subj <- paste0(test_dens$Subject, "_", test_dens$ImageNumber)

  ## compute similarity between each trial and study image derived from group average
  test_sim <- template_similarity(stud$study_dens, test_dens, "Image_Subj", permutations=50) %>% mutate(min_onset=min_onset, max_onset=max_onset)

  test_sim <- inner_join(test_sim, test_tab, by="Image_Subj")

}

library(purrr)
bin_onsets <- seq(0, 3500, by=250)
test_sim_binned <- bin_onsets %>% map(~ binned_similarity(., . + 250)) %>% map_df(bind_rows) %>% select(-fixgroup.x, -fixgroup.y, -density)
write.table(test_sim_binned, "test_sim_binned.txt", row.names=FALSE)



library(ggplot2)
binned_saliency <- test_sim_binned %>% group_by(min_onset, Saliency) %>% summarize(eye_sim_diff=mean(eye_sim_diff))

qplot(min_onset, eye_sim_diff, colour=factor(Saliency), data=binned_saliency, geom=c("point", "line"))

binned_saliency_oldnew <- test_sim_binned %>% group_by(min_onset, Match, Accuracy) %>% summarize(eye_sim_diff=mean(eye_sim_diff),
                                                                                   perm_sim=mean(perm_sim), eye_sim=mean(eye_sim))

qplot(min_onset, eye_sim_diff, colour=factor(Accuracy), facets= . ~ Match,
      data=subset(binned_saliency_oldnew), geom=c("point", "line"))

binned_saliency_acc_oldnew <- test_sim_binned %>% group_by(min_onset, Match, Saliency, Accuracy) %>% summarize(eye_sim_diff=mean(eye_sim_diff),
                                                                                                 perm_sim=mean(perm_sim), eye_sim=mean(eye_sim))

qplot(min_onset, eye_sim_diff, colour=factor(Accuracy), facets= Saliency ~ Match,
      data=subset(binned_saliency_acc_oldnew), geom=c("point", "line"))

