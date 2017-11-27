library(dplyr)
library(tibble)
library(tidyr)

## load study data
pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_study_input.csv")) %>% filter(Image != ".") %>% droplevels()
pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


## load test data
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delay_input.csv")) %>% filter(Image != ".") %>% droplevels()


#pctest$StudyVersion <- as.character(pcstudy$ImageSet[match_ind])
#pctest$Match <- ifelse(pctest$StudyVersion == pctest$ImageSet, "match", "mismatch")

study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion", "ImageRepetition",
                              "ImageSet", "Block", "Image", "ImageNumber"))

test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject"), data=pctest,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion", "Saliency", "Accuracy",
                              "ImageSet", "Block", "Duration", "ImageNumber"))

tmp <- test_tab %>% group_by(Subject) %>% rowwise() %>% do( {
  sversion <- as.character(filter(pcstudy, Subject == .$Subject & ImageNumber == .$ImageNumber)$ImageSet[1])
  Match <- if (sversion == as.character(.$ImageSet)) {
    "match"
  } else {
    "mismatch"
  }
  data.frame(Match=Match)
})

test_tab$Match <- tmp$Match

## construct heatmaps for the study phase, averaged over subjects
study_dens <- density_by(study_tab, groups="Image", xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)

study_dens_avg <- Reduce("+", lapply(study_dens$density, function(x) x$z))/length(study_dens)
study_dens_avg <- study_dens_avg/sum(study_dens_avg)
study_dens_sqw <- study_dens_avg^(1/2)
study_dens_sqw <- study_dens_sqw/sum(study_dens_sqw)

#sigma <- .1
#weights <- exp(-study_dens_avg^2/(2 * sigma^2))

saliency <- study_dens %>% rowwise() %>% do({
  zcubed <- .$density$z^(1/3)
  zcubed <- zcubed/sum(zcubed)

  zsqw <-  .$density$z^(1/2)
  zsqw <- zsqw / sum(zsqw)

  zrank <- rank(.$density$z)
  zrank <- zrank/sum(zrank)

  gg <- expand.grid(x=1:80, y=1:60)
  tibble(Image=.$Image, zcuberoot = list(zcubed), zrank=list(matrix(zrank, 80,60)), zsqw=list(zsqw))
})

#write.table(saliency, "~/Dropbox/Jordana_experiments/Jordana_saliency_study/saliency_grid.txt", row.names=FALSE)

library(imager)
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

  sal <- saliency$zsqw[[which(as.character(saliency$Image) == im)]]
  sal_other <- saliency$zsqw[[which(as.character(saliency$Image) == im_other)]]

  fix <- .$fixgroup
  fm <- round(cbind(fix$x, fix$y)/10)
  fm[,1] <- ifelse(fm[,1] < 1, 1, fm[,1])
  fm[,2] <- ifelse(fm[,2] < 1, 1, fm[,2])

  if (.$Saliency == 100) {
    mvals <- rep(1, nrow(fm))
  } else {
    mask <- maskset[[as.character(.$Image)]]
    mvals <- mask[fm]
  }


  vis <- ifelse(mvals, sal[fm], NA)
  novis <- ifelse(mvals == 0, sal[fm], NA)
  tot <- sal[fm]

  vis_other <- ifelse(mvals, sal_other[fm], NA)
  novis_other <- ifelse(mvals==0, sal_other[fm], NA)
  tot_other <- sal_other[fm]

  bvis <- study_dens_sqw[fm]

  pvis <- ifelse(mvals, 1, 0)
  pnovis <- ifelse(mvals == 0, 1, 0)

  ret <- data.frame(vis=vis-bvis, novis=novis-bvis, totvis=tot-bvis, vis_other=vis_other-bvis, novis_other=novis_other-bvis, tot_other=tot_other-bvis, pvis=pvis, pnovis=pnovis)
  as_tibble(cbind(ret, .))
}) %>% ungroup()

sal_out <- gather(sal_out, key=measure, value=sim, vis, novis, totvis, vis_other, novis_other, tot_other)

#sal_out %>% group_by(Saliency, Duration) %>% summarize(vis=mean(vis), novis=mean(novis), totvis=mean(totvis))
library(mgcv)
library(ggplot2)
gam.1 <- gam(totvis ~ s(fixgroup.onset), data=sal_out)

ggplot(aes(fixgroup.onset, sim, colour=measure, linetype=Match), data=subset(sal_out, measure %in% c("totvis", "tot_other") )) + facet_wrap( ~ Duration) +
  geom_smooth(se=FALSE, method=gam, formula = y ~ s(x, k=7, fx=TRUE))

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
