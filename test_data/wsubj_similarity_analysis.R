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


subject_dens <- density_by(study_tab, groups=c("Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=60, angular=TRUE)


## load test data
pctest <- as_tibble(read.csv("~/Dropbox/New_pc_behavioural_data/testdelay_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()


## create table for each test trial
test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                       groupvar=c("Image", "Subject"), data=pctest,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion", "Saliency", "Accuracy",
                              "ImageSet", "Trial", "Duration", "ImageNumber"))

m <- test_tab %>% group_by(Subject) %>% rowwise() %>% do( {
  sversion <- as.character(filter(pcstudy, Subject == .$Subject & ImageNumber == .$ImageNumber)$ImageSet[1])
  Match <- if (sversion == as.character(.$ImageSet)) "match"  else "mismatch"
  data.frame(Match=Match)
})

test_tab$Image_Subj <- paste0(test_tab$Subject, "_", test_tab$ImageNumber)
test_tab$Image_Subj_Version <- paste0(test_tab$Subject, "_", test_tab$ImageVersion)
test_tab$Match <- m$Match


outdim <- c(16,12)
## construct heatmaps for the study phase, averaged within subjects
study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=outdim,
                         duration_weighted=TRUE, sigma=60)

study_dens_ang <- density_by(study_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=outdim,
                         duration_weighted=TRUE, sigma=60, angular=TRUE, angle_bins=12)



## construct heatmaps for the study phase, averaged over subjects
study_dens_all <- density_by(study_tab, groups=c("ImageVersion"), xbounds=c(0,800), ybounds=c(0,600), outdim=outdim, duration_weighted=TRUE, sigma=60)

## compute grand mean density
study_dens_avg <- Reduce("+", lapply(study_dens$density, function(x) x$z))/length(study_dens)
study_dens_avg <- study_dens_avg/sum(study_dens_avg)

dens_avg_subj <- study_dens %>% group_by(Subject) %>% do({
  dens <- Reduce("+", lapply(.$density, function(x) x$z))/nrow(.)
  tibble(Subject=.$Subject[1], dens=list(dens))
})

test_dens <- density_by(test_tab, groups=c("ImageNumber", "Subject"),xbounds=c(0,800), ybounds=c(0,600), outdim=outdim, duration_weighted=TRUE, sigma=60)
test_dens_ang <- density_by(test_tab, groups=c("ImageNumber", "Subject"),xbounds=c(0,800), ybounds=c(0,600), outdim=outdim, duration_weighted=TRUE, sigma=60,
                            angular=TRUE, angle_bins=12)


test_dens_all <- density_by(test_tab, groups=c("ImageVersion", "Subject"),xbounds=c(0,800), ybounds=c(0,600), outdim=outdim, duration_weighted=TRUE, sigma=60)

test_dens$Image_Subj <- paste0(test_dens$Subject, "_", test_dens$ImageNumber)
study_dens$Image_Subj <- paste0(study_dens$Subject, "_", study_dens$ImageNumber)

test_dens_ang$Image_Subj <- paste0(test_dens_ang$Subject, "_", test_dens_ang$ImageNumber)
study_dens_ang$Image_Subj <- paste0(study_dens$Subject, "_", study_dens_ang$ImageNumber)



binned_similarity <- function(min_onset, max_onset, method="cosine") {
  print(min_onset)
  ## create table for each test trial
  pctest_binned <- pctest %>% filter(FixOffset >= min_onset & FixOffset < max_onset)


  test_tab_binned <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                        groupvar=c("Image", "Subject"), data=pctest_binned,
                        clip_bounds=c(112, (112+800), 684, 84),
                        vars=c("ImageVersion", "Saliency", "Accuracy",
                               "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber", "testdelayOnset"))

  test_tab_binned$Image_Subj <- paste0(test_tab_binned$Subject, "_", test_tab_binned$ImageNumber)


  test_tab_binned$Match <- test_tab_binned$ImageRepetition

  test_dens <- density_by(test_tab_binned, groups=c("ImageNumber", "Subject"),xbounds=c(0,800), ybounds=c(0,600),
                          outdim=outdim, duration_weighted=TRUE, sigma=50)
  test_dens$Image_Subj <- paste0(test_dens$Subject, "_", test_dens$ImageNumber)

  ## compute similarity between each trial and study image derived from group average
  test_sim <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=10, method=method) %>% mutate(min_onset=min_onset, max_onset=max_onset)

  test_sim <- inner_join(test_sim, test_tab_binned, by="Image_Subj") %>% mutate(method=method)

}

library(purrr)
library(ggplot2)
bin_onsets <- seq(0, 3500, by=250)
sal_binned_cos <- bin_onsets %>% map(~ binned_similarity(., . + 250, method="jaccard")) %>% map_df(bind_rows) #%>% select(-fixgroup.x, -fixgroup.y, -density)
sal_binned_cos$sim <- sal_binned_cos$eye_sim_diff


sal_sum1 <- sal_binned_cos %>% filter(Subject.x < 300) %>% group_by(Accuracy,  ImageRepetition, min_onset) %>% summarize(eye_sim_diff=mean(eye_sim_diff), eye_sim=mean(eye_sim))
qplot(min_onset, eye_sim_diff, colour=factor(Accuracy), data=sal_sum1, facets = ~ ImageRepetition) + geom_line()


test_reg <- template_regression(study_dens, test_dens, "Image_Subj", subject_dens, baseline_key="Subject", method="rank")
test_reg$Image_Subj <- paste0(test_reg$Subject, "_", test_reg$ImageNumber)







cor_with_accuracy <- function(simtab, type) {
  #test_sim <- inner_join(simtab, test_tab, by="Image_Subj")
  test_sim <- subset(simtab, Subject.x < 300)

  res <- test_sim %>% group_by(Subject.x) %>% do({

    hits <- sum(.$Accuracy[.$Match == "old"])
    cr <- sum(.$Accuracy[.$Match == "new"])

    misses <- sum(.$Match == "old") - hits
    fa <- sum(.$Match == "new") - cr
    res <- dprime(hits, fa, misses, cr)

    sim_match <- mean(.$sim[.$Match == "old"])
    sim_mismatch <- mean(.$sim[.$Match == "new"])

    acc_match <- mean(.$Accuracy[.$Match == "old"])
    acc_mismatch <- mean(.$Accuracy[.$Match == "new"])
    data.frame(Subject=.$Subject.x[1], dprime=res$dprime, aprime=res$aprime, beta=res$beta, bppd=res$bppd, c=res$c,
               sim_match=sim_match, sim_mismatch=sim_mismatch, acc_match=acc_match, acc_mismatch=acc_mismatch)
  })

  data.frame(cormatch=cor(res$aprime, res$sim_match, use="complete.obs"),
             cormismatch=cor(res$aprime, res$sim_mismatch, use="complete.obs"),
             corhits=cor(res$sim_match, res$acc_match,use="complete.obs"),
             corcr=cor(res$sim_mismatch, res$acc_mismatch,use="complete.obs"), type=type)

}


## compute similarity between each trial and study image derived from group average
test_sim_spear <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=10, method="spearman")
test_sim_spear %>% mutate(sim = eye_sim_diff) %>% cor_with_accuracy(., type="spearman")
test_sim_spear %>% mutate(sim = eye_sim) %>% cor_with_accuracy(., type="spearman")

test_sim_cos <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=10, method="cosine")
test_sim_cos %>% mutate(sim = eye_sim_diff) %>% cor_with_accuracy(., type="cosine")
test_sim_cos %>% mutate(sim = eye_sim) %>% cor_with_accuracy(., type="cosine")

test_sim_jacc <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=10, method="jaccard")
test_sim_jacc %>% mutate(sim = eye_sim_diff) %>% cor_with_accuracy(., type="jaccard")
test_sim_jacc %>% mutate(sim = eye_sim) %>% cor_with_accuracy(., type="jaccard")

test_sim_l1 <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=10, method="l1")
test_sim_l1 %>% mutate(sim = eye_sim_diff) %>% cor_with_accuracy(., type="l1")
test_sim_l1 %>% mutate(sim = eye_sim) %>% cor_with_accuracy(., type="l1")

test_sim_dcov <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=10, method="dcov")
test_sim_dcov %>% mutate(sim = eye_sim_diff) %>% cor_with_accuracy(., type="dcov")
test_sim_dcov %>% mutate(sim = eye_sim) %>% cor_with_accuracy(., type="dcov")

test_sim_cos_ang <- template_similarity(study_dens_ang, test_dens_ang, "Image_Subj", permutations=10, method="cosine")
test_sim_cos_ang %>% mutate(sim = eye_sim_diff) %>% cor_with_accuracy(., type="cosine")
test_sim_cos_ang %>% mutate(sim = eye_sim) %>% cor_with_accuracy(., type="cosine")



test_sim_angle <- template_similarity(study_dens, test_dens, "Image_Subj", permutations=10, method="angular")

library(lme4)
test_sim <- inner_join(test_sim_cos_ang, test_tab, by="Image_Subj")
lme.1 <- glmer(Accuracy ~ eye_sim_diff*Match + (1 | Subject.x), family="binomial", data=subset(test_sim, Subject.x > 300))
summary(lme.1)

## join eye_sim with test_tab
test_sim <- inner_join(test_sim, test_tab, by="Image_Subj")

test_reg <- inner_join(test_reg, test_tab, by="Image_Subj")

## compute means over variables
test_sim_agg <- aggregate(eye_sim_diff ~ Saliency + Duration, data=test_sim, FUN=mean)
test_sim_agg1 <- aggregate(eye_sim_diff ~ Duration + Accuracy + Match, data=test_sim, FUN=mean)
test_sim_agg2 <- aggregate(eye_sim_diff ~ Saliency + Accuracy + Match, data=test_sim, FUN=mean)

test_reg_agg <- aggregate(beta_source ~ Saliency + Duration, data=test_reg, FUN=mean)
test_reg_agg1 <- aggregate(beta_source ~ Duration + Accuracy + Match, data=test_reg, FUN=mean)
test_reg_agg2 <- aggregate(beta_source ~ Saliency + Accuracy + Match, data=test_reg, FUN=mean)



# aggregate over subjects
subj_acc <- aggregate(Accuracy ~ Subject.x + Match, data=test_sim, FUN=mean)
subj_sim <- aggregate(eye_sim_diff ~ Subject.x + Match, data=test_sim, FUN=mean)
subj_acc$sim <- subj_sim$eye_sim_diff

subj_acc <- aggregate(Accuracy ~ Subject.x + Match, data=test_reg, FUN=mean)
subj_sim <- aggregate(beta_source ~ Subject.x + Match, data=test_reg, FUN=mean)
subj_acc$sim <- subj_sim$beta_source

subj_acc <- aggregate(Accuracy.x ~ Subject.x + Match.x, data=test_sim_all, FUN=mean)
subj_sim <- aggregate(eye_sim_diff ~ Subject.x + Match.x, data=test_sim_all, FUN=mean)
subj_acc$sim <- subj_sim$eye_sim_diff



test_sim_all <- inner_join(test_sim_all, test_tab, by="Image_Subj_Version")
test_sim_all_agg <- aggregate(eye_sim_diff ~ Saliency.x + Duration.x, data=test_sim_all, FUN=mean)
test_sim_all_agg1 <- aggregate(eye_sim_diff ~ Duration.x + Accuracy.x + Match.x, data=test_sim_all, FUN=mean)
test_sim_all_agg2 <- aggregate(eye_sim_diff ~ Saliency.x + Accuracy.x + Match.x, data=test_sim_all, FUN=mean)
test_sim_all_agg3 <- aggregate(eye_sim_diff ~ Duration.x + Accuracy.x + Match.x, data=test_sim_all, FUN=mean)
test_sim_all_agg4 <- aggregate(eye_sim_diff ~ Duration.x, data=test_sim_all, FUN=mean)
test_sim_all_agg5 <- aggregate(eye_sim_diff ~ Saliency.x + Match.x, data=test_sim_all, FUN=mean)
# aggregate over subjects
subj_acc_all <- aggregate(Accuracy ~ Subject.x + Match, data=subset(test_sim_all), FUN=mean)
subj_sim_all <- aggregate(eye_sim_diff ~ Subject.x + Match, data=subset(test_sim_all), FUN=mean)
subj_acc_all$sim <- subj_sim_all$eye_sim_diff

library(ggplot2)

ggplot(aes(Saliency, eye_sim_diff, colour=factor(Duration)), data=test_sim_agg) + geom_line()
ggplot(aes(Saliency, eye_sim_diff, colour=factor(Accuracy)), data=test_sim_agg2) + geom_line() + facet_wrap(~ Duration)

ggplot(aes(Saliency, eye_sim_diff, linetype=factor(Accuracy), colour=factor(Duration)), data=test_sim_all_agg) +  facet_wrap( ~ Match) + geom_line()
ggplot(aes(Saliency, eye_sim_diff, linetype=factor(Accuracy), colour=factor(Duration)), data=test_sim_agg) +  facet_wrap( ~ Match) + geom_line()

