library(dplyr)
library(tibble)
library(tidyr)

## load study data
start_time <- 2500
pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_study_input.csv")) %>% filter(Image != ".") %>% droplevels()
pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


## load test data
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delay_input.csv")) %>%
  filter(Image != "." & FixStartTime < start_time) %>% droplevels()

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

## construct heatmaps for the study phase, averaged over subjects
study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)
study_dens_all <- density_by(study_tab, groups=c("ImageVersion"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)

study_dens_avg <- Reduce("+", lapply(study_dens$density, function(x) x$z))/length(study_dens)
study_dens_avg <- study_dens_avg/sum(study_dens_avg)

test_dens <- density_by(test_tab, groups=c("ImageNumber", "Subject"),xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)
test_dens_all <- density_by(test_tab, groups=c("ImageVersion", "Subject"),xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)

test_dens$Image_Subj <- paste0(test_dens$Subject, "_", test_dens$ImageNumber)
study_dens$Image_Subj <- paste0(study_dens$Subject, "_", study_dens$ImageNumber)
test_tab$Image_Subj <- paste0(test_tab$Subject, "_", test_tab$ImageNumber)

## compute similarity between each trial and study image derived from group average
#test_sim <- template_similarity(study_dens, test_dens, "Image_Subj", reference_density = study_dens_avg)
#test_sim_all <- template_similarity(study_dens_all, test_dens_all, "ImageVersion", reference_density = study_dens_avg)
test_sim <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=5)
test_sim_all <- template_similarity(study_dens_all, test_dens_all, "ImageVersion", permutations=5)


test_sim_all$Image_Subj_Version <- paste0(test_sim_all$Subject, "_", test_sim_all$ImageVersion)
test_tab$Image_Subj_Version <- paste0(test_tab$Subject, "_", test_tab$ImageVersion)

## add back in per-trial experimental variables of interest
#test_sim$saliency <- test_dens$Saliency
#test_sim$duration <- test_dens$Duration
#test_sim$accuracy <- test_dens$Accuracy
#test_sim$match <- test_dens$Match

test_sim <- inner_join(test_sim, test_tab, by="Image_Subj")
test_sim_agg <- aggregate(eye_sim ~ Saliency + Duration + Accuracy + Subject.x, data=test_sim, FUN=mean)
subj_acc <- aggregate(Accuracy ~ Subject.x, data=test_sim, FUN=mean)
subj_sim <- aggregate(eye_sim ~ Subject.x, data=test_sim, FUN=mean)


test_sim_all <- inner_join(test_sim_all, test_tab, by="Image_Subj_Version")
test_sim_all_agg <- aggregate(eye_sim ~ Saliency + Duration + Accuracy + Match, data=test_sim_all, FUN=mean)
test_sim_all_agg <- aggregate(eye_sim ~ Saliency + Duration, data=test_sim_all, FUN=mean)
subj_acc_all <- aggregate(Accuracy ~ Subject.x, data=subset(test_sim_all), FUN=mean)
subj_sim_all <- aggregate(eye_sim ~ Subject.x, data=subset(test_sim_all), FUN=mean)


library(ggplot2)
ggplot(aes(Saliency, eye_sim, colour=factor(Duration)), data=test_sim_all_agg) + geom_line()
ggplot(aes(Saliency, eye_sim, linetype=factor(Accuracy), colour=factor(Duration)), data=test_sim_all_agg) +  facet_wrap( ~ Match) + geom_line()
ggplot(aes(Saliency, eye_sim, linetype=factor(Accuracy), colour=factor(Duration)), data=test_sim_agg) +  facet_wrap( ~ Match) + geom_line()

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) + facet_wrap( ~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("vis", "vis_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))


ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) + facet_wrap( ~ Saliency) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("vis", "vis_other"))) + facet_wrap( ~ Saliency) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("vis", "vis_other"))) + facet_wrap( ~ Saliency) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))


ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("novis", "novis_other"))) + facet_wrap( ~ Saliency) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("novis", "novis_other") & Duration == 250)) + facet_wrap( ~ factor(Accuracy)) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("novis", "novis_other") & Saliency < 40)) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("novis", "novis_other"))) + facet_wrap( ~ Duration) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("novis", "novis_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))


ggplot(aes(fixgroup.onset, totvis, colour=factor(Accuracy)), data=sal_out) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))


ggplot(aes(fixgroup.onset, vis, colour=factor(Saliency)), data=sal_out) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, vis, colour=factor(Duration), linetype=Match), data=sal_out) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=7, fx=TRUE))

ggplot(aes(fixgroup.onset, vis, linetype=factor(Match)), data=sal_out) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))
