library(dplyr)
library(tibble)
library(tidyr)
library(imager)
library(mgcv)
library(ggplot2)

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
test_tab <- test_tab %>% filter(Subject < 300)
study_tab <- study_tab %>% filter(Subject < 300)

## construct heatmaps for the study phase, averaged within subjects
study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                         duration_weighted=TRUE, sigma=60)


## get the masks for each image
maskset <- lapply(levels(pctest$Image), function(im) {
  fname <- paste0("~/Dropbox/Jordana_experiments/Jordana_saliency_study/images/Mat_", gsub("jpeg", "jpeg.rds", im))
  if (file.exists(fname)) {
    readRDS(fname)
  } else {
    print(paste("no ", im))
    NULL
  }
})

names(maskset) <- levels(pctest$Image)

saliency <- study_dens %>% rowwise() %>% do({

  zdens <- .$density$z/sum(.$density$z)
  zsqw <-  .$density$z^(1/2)
  zsqw <- zsqw / sum(zsqw)

  zrank <- rank(.$density$z)
  zrank <- zrank/sum(zrank)

  gg <- expand.grid(x=1:80, y=1:60)
  tibble(Subject=.$Subject, Image=.$Image, zdens=list(zdens), zrank=list(matrix(zrank, 80,60)), zsqw=list(zsqw))
})




sal_out <- test_tab %>% rowwise() %>% do({


  print(as.character(.$Image))

  if (.$Saliency == 100) {
    im <- paste0(strsplit(as.character(.$Image), "_")[[1]][1:2], collapse="_")
    im <- paste0(im, "_1")
  } else {
    im <- paste0(strsplit(as.character(.$Image), "_")[[1]][1:3], collapse="_")
  }

  other_version <- if (as.character(.$ImageSet[1]) == "A") "B" else "A"

  im <- paste0(im, ".jpeg")
  im_other <- gsub(as.character(.$ImageSet[1]), other_version, im)

  sal <- saliency$zdens[[which(saliency$Image == .$ImageNumber & saliency$Subject == .$Subject)]]

  fix <- .$fixgroup
  fm <- round(cbind(fix$x, fix$y)/10)
  fm[,1] <- ifelse(fm[,1] < 1, 1, fm[,1])
  fm[,2] <- ifelse(fm[,2] < 1, 1, fm[,2])

  if (.$Saliency == 100) {
    mvals <- rep(1, nrow(fm))
  } else {
    mask <- maskset[[as.character(.$Image)]]
    #mask <- t(apply(mask, 1, rev))
    mvals <- mask[fm]

    #mask <- maskset[[as.character(.$Image)]]
    #mvals <- mask[fm]
  }

  ## the salience of the visible items
  vis <- ifelse(mvals, sal[fm], NA)

  ## the salience of the invisible items
  novis <- ifelse(mvals == 0, sal[fm], NA)

  ## the total salience
  tot <- sal[fm]

  ret <- data.frame(vis=vis, novis=novis, totvis=tot)
  as_tibble(cbind(ret, .))
}) %>% ungroup()

sal_out <- gather(sal_out, key=measure, value=sim, vis, novis, totvis)

#sal_out %>% group_by(Saliency, Duration) %>% summarize(vis=mean(vis), novis=mean(novis), totvis=mean(totvis))

gam.1 <- gam(sim ~ s(fixgroup.onset), data=subset(sal_out, measure=="novis"))

ggplot(aes(fixgroup.onset, sim, colour=ImageRepetition), data=subset(sal_out, measure %in% c("novis") & Accuracy == 1)) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=ImageRepetition), data=subset(sal_out, measure %in% c("vis") & Accuracy == 1)) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=12, fx=TRUE))



ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=factor(Accuracy)), data=subset(sal_out, measure %in% c("totvis", "tot_other"))) + facet_wrap( ~ Match) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("vis", "vis_other"))) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=10, fx=TRUE))

