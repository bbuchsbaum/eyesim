library(dplyr)
library(tibble)
library(tidyr)



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



## construct heatmaps for the study phase, averaged within subjects
study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                         duration_weighted=TRUE, sigma=60)


#study_dens_subj_avg <- density_by(study_tab, groups=c("Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
#                                  duration_weighted=TRUE, sigma=80)


study_dens_avg <- Reduce("+", lapply(study_dens$density, function(x) x$z))/length(study_dens)
study_dens_avg <- study_dens_avg/sum(study_dens_avg)
study_dens_sqw <- study_dens_avg^(1/2)
study_dens_sqw <- study_dens_sqw/sum(study_dens_sqw)

#sigma <- .1
#weights <- exp(-study_dens_avg^2/(2 * sigma^2))

test_tab <- test_tab %>% filter(Subject < 300)

saliency <- study_dens %>% rowwise() %>% do({

  zdens <- .$density$z/sum(.$density$z)
  zsqw <-  .$density$z^(1/2)
  zsqw <- zsqw / sum(zsqw)

  zrank <- rank(.$density$z)
  zrank <- zrank/sum(zrank)

  gg <- expand.grid(x=1:80, y=1:60)
  tibble(Subject=.$Subject, Image=.$Image, zdens=list(zdens), zrank=list(matrix(zrank, 80,60)), zsqw=list(zsqw))
})

#write.table(saliency, "~/Dropbox/Jordana_experiments/Jordana_saliency_study/saliency_grid.txt", row.names=FALSE)

library(imager)
maskset <- lapply(levels(pctest$Image), function(im) {
  fname <- paste0("~/Dropbox/Jordana_experiments/Jordana_saliency_study/images/Mat_", gsub("jpeg", "jpeg.rds", im))
  if (file.exists(fname)) {
    print("got it")
    readRDS(fname)
  } else {
    print(paste("no ", im))
    NULL
  }
})

names(maskset) <- levels(pctest$Image)

binned_visibility <- function(min_onset, max_onset) {
  print(min_onset)
  pctest_binned <- pctest %>% filter(FixOffset >= min_onset & FixOffset < max_onset)
  test_tab_binned <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                               groupvar=c("Image", "Subject"), data=pctest_binned,
                               clip_bounds=c(112, (112+800), 684, 84),
                               vars=c("ImageVersion", "Saliency", "Accuracy",
                                      "ImageSet", "Trial", "Duration", "ImageNumber", "ImageRepetition"))

  #binned_visibility <- function(min, max) {
  sal_out <- test_tab_binned %>% rowwise() %>% do({

    #print(as.character(.$Image))

    if (.$Saliency == 100) {
      im <- paste0(strsplit(as.character(.$Image), "_")[[1]][1:2], collapse="_")
      im <- paste0(im, "_1")
    } else {
      im <- paste0(strsplit(as.character(.$Image), "_")[[1]][1:3], collapse="_")
    }

    im <- paste0(im, ".jpeg")
    sal <- saliency$zdens[[which(saliency$Image == .$ImageNumber & saliency$Subject == .$Subject)]]

    #browser()
    fix <- .$fixgroup
    fm <- round(cbind(fix$x, fix$y)/10)
    fm[,1] <- ifelse(fm[,1] < 1, 1, fm[,1])
    fm[,2] <- ifelse(fm[,2] < 1, 1, fm[,2])

    if (.$Saliency == 100) {
      mvals <- rep(1, nrow(fm))
    } else {
      mask <- maskset[[as.character(.$Image)]]
      mask <- t(apply(mask, 1, rev))
      mvals <- mask[fm]
    }

    #browser()

    ## the salience of the visible items
    vis <- ifelse(mvals, sal[fm], NA)

    ## the salience of the invisible items
    novis <- ifelse(mvals == 0, sal[fm], NA)

    pvis <- ifelse(mvals, 1, 0)
    pnovis <- ifelse(mvals == 0, 1, 0)

    ret <- data.frame(pvis=sum(pvis), pnovis=sum(pnovis), pcount=length(mvals), min_onset=min_onset, sactual=sum(mask)/length(mask))
    as_tibble(cbind(ret, .))
  }) %>% ungroup()

}


library(purrr)
library(ggplot2)

bin_onsets <- seq(0, 3500, by=250)
sal_binned <- bin_onsets %>% map(~ binned_visibility(., . + 250)) %>% map_df(bind_rows) #%>% select(-fixgroup.x, -fixgroup.y, -density)

sal_sum1 <- sal_binned %>% group_by(Saliency, ImageRepetition, min_onset) %>% dplyr::summarize(pvis= (sum(pvis)/sum(pcount)),
                                                                                               pnovis=sum(pnovis)/sum(pcount), count=sum(pcount), sactual=median(sactual))
sal_sum1$oddsrat <- (sal_sum1$pvis/(1-sal_sum1$pvis)) /(sal_sum1$sactual/(1-sal_sum1$sactual))

sal_out <- gather(sal_out, key=measure, value=sim, vis, novis, totvis, pvis, pnovis)

#sal_out %>% group_by(Saliency, Duration) %>% summarize(vis=mean(vis), novis=mean(novis), totvis=mean(totvis))
library(mgcv)
library(ggplot2)
gam.1 <- gam(totvis ~ s(fixgroup.onset), data=sal_out)

ggplot(aes(fixgroup.onset, sim, linetype=measure), data=subset(sal_out, measure %in% c("vis", "novis")))  +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE)) + facet_wrap(~ Match)

ggplot(aes(fixgroup.onset, sim, colour=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis") )) + facet_wrap(Saliency ~ Match, ncol=5) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=8, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis") )) + facet_wrap(~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=8, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=factor(Duration)), data=subset(sal_out, measure %in% c("totvis") )) + facet_wrap(~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=8, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=factor(Saliency)), data=subset(sal_out, measure %in% c("totvis") )) + facet_wrap(~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=8, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis") & Saliency < 40 & Duration < 700)) + facet_wrap(~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=6, fx=TRUE))




ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("totvis"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis"))) + facet_wrap( ~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=factor(Saliency)), data=subset(sal_out, measure %in% c("totvis"))) +  facet_wrap(~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))


ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  facet_wrap(~ Match, nrow=2) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  facet_wrap(Match ~ Duration, nrow=2) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))


ggplot(aes(fixgroup.onset, sim, colour=factor(Duration), linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("vis"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE)) + facet_wrap(~ Match)


ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis"))) +
  facet_wrap( ~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("novis", "novis_other") & Duration == 250)) + facet_wrap( ~ factor(Accuracy)) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("novis", "novis_other") & Saliency < 40)) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Match)), data=subset(sal_out, measure %in% c("novis"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE)) + facet_wrap(Duration ~ Saliency, nrow=3)

ggplot(aes(fixgroup.onset, sim, colour=measure), data=subset(sal_out, measure %in% c("novis", "novis_other",
                                                                                                     "totvis", "totvis_other",
                                                                                                     "vis", "vis_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE)) + facet_wrap(Match ~ Duration, nrow=2)


ggplot(aes(fixgroup.onset, totvis, colour=factor(Accuracy)), data=sal_out) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))


ggplot(aes(fixgroup.onset, vis, colour=factor(Saliency)), data=sal_out) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, vis, colour=factor(Duration), linetype=Match), data=sal_out) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=7, fx=TRUE))

ggplot(aes(fixgroup.onset, vis, linetype=factor(Match)), data=sal_out) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))
