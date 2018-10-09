library(dplyr)
library(tibble)
library(tidyr)
library(ppcor)  ### install.packages("ppcor")
library(MASS)

## load study data
pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()



## load test data to get accuracy values (only taking first fixation row of each trial -- see: 'slice(1)')
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/delay_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>%
  group_by(Subject, Trial) %>% slice(1) %>% droplevels() %>% ungroup()

pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject", "Block"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))

study_tab_all <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("ImageVersion"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageSet", "Block", "Image", "ImageNumber"))


## compute average fixation map for each subject (to be used as baseline) when computing encoding similarity
subject_dens <- density_by(study_tab, groups=c("Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                           duration_weighted=TRUE, sigma=80)

# dfx <- study_tab %>% filter(Subject %in% c(115,116,117,118,119,120))
# dfxdens <- density_by(study_tab, groups=c("Subject", "ImageVersion"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
#                       duration_weighted=TRUE, sigma=80)
# imat <- do.call(rbind, purrr::map(dfxdens$density, ~ as.vector(.$z)))


## compute encoding sim for all subjects/images.
## loop over subject/imagenumber cominations
encoding_sim <- study_tab %>% group_by(Subject, ImageNumber) %>% do({
  cmb <- combn(length(.$fixgroup),2)

  ## compute density map for each trial separately
  dens <- lapply(.$fixgroup, function(fg) {
    eye_density(fg, sigma=80,xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE)
  })

  ## get the subject average density map for the current subject
  sd <- subject_dens %>% filter(Subject == .$Subject[1])

  ## loop over all (6) pairwise combinations
  savg <- apply(cmb, 2, function(ind) {

    ## density map of pair 1
    x1 <- dens[[ind[1]]]

    ## density map of pair 2
    x2 <- dens[[ind[2]]]


    csim <- as.vector(x1$z) %*% as.vector(x2$z)
    avg <- as.vector(sd$density[[1]]$z)
    csim2 <-  ((as.vector(x1$z) %*% avg) + (as.vector(x2$z) %*% avg))/2

    #df1 <- data.frame(x1=as.vector(x1$z), x2=as.vector(x2$z), avg=as.vector(sd$density[[1]]$z))

    ## multiple regression where the average density map ('avg') is a covariate
    #res <- lm(x1 ~ avg + x2, data=df1)
    #coef(res)[2:3]
    #browser()
    #res <- pcor(df1, "spearman")
    #res <- pcor(df1)
    #res$estimate[3:2,1]
    c(csim=csim, csim_avg=csim2, csim_diff = csim-csim2)
  })



  fix <- sapply(.$fixgroup, nrow)
  ## total number of fixations over the four repetitions
  total_fix <- sum(fix)

  m <- do.call(rbind, .$fixgroup)
  ## dispersion of fixations
  sd_fix <- sd(m$x) * sd(m$y)

  ## slope of fixatons over trials
  slope <- coef(lm(fix ~ seq(1,length(fix))))[2]
  slope_cor <- cor(fix, seq(1,length(fix)), method="spearman", use="complete.obs")

  rmns <- rowMeans(savg)
  data.frame(Subject=.$Subject[1], ImageNumber=.$ImageNumber[1], sim_wavg=rmns[2],
             sim_within=rmns[1], sim_wdiff=rmns[3], total_fix=total_fix,
             slope=slope, slope_cor=slope_cor, sd_fix=sd_fix)
})



## code match trials
m <- pctest %>% group_by(Subject) %>% rowwise() %>% do( {
  sversion <- as.character(filter(pcstudy, Subject == .$Subject & ImageNumber == .$ImageNumber)$ImageSet[1])
  Match <- if (sversion == as.character(.$ImageSet)) "match"  else "mismatch"
  data.frame(Match=Match)
})

pctest$Match <- m$Match

## join encoding sim and test data
pctest2 <- inner_join(pctest, encoding_sim, by=c("Subject", "ImageNumber"))
pctest2$sid <- factor(pctest2$Subject)


## experiment with within-subject lmers... after adding ImageNumber as random effect, sig. effects are gone.
lme.1 <- lmer(Accuracy ~ sim_wavg  + Saliency + Duration + (1 | sid), data=subset(pctest2, Match=="mismatch"))
lme.2 <- lmer(Accuracy ~ sim_within + sim_wavg + Saliency + Duration + (1 | sid), data=subset(pctest2, Match=="mismatch"))
lme.3 <- lmer(Accuracy ~ sim_within + sim_wavg + Saliency + Duration + (1 | sid) + (1 | ImageNumber), data=subset(pctest2, Match=="mismatch"))

lme.4 <- lmer(Accuracy ~ sim_wdiff + Saliency + Duration + (1 | sid) + (1 | ImageNumber), data=subset(pctest2, Match=="mismatch"))


## for between subject analysis, average accuracy over subject
pctestacc <- pctest %>% group_by(Subject) %>% summarize(acc=mean(Accuracy))

## compute subject averages for encoding_sim
pcstudysim <- encoding_sim %>% group_by(Subject) %>% summarize(sim_wavg=mean(sim_wavg), sim_within=mean(sim_within), sim_wdiff=mean(sim_wdiff),
                                                               sd_fix=mean(sd_fix),
                                                               total_fix=mean(total_fix), slope=mean(slope), slope_cor=mean(slope_cor, na.rm=TRUE))
## join accuracy and encoding_sim averages
dfx <- inner_join(pctestacc, pcstudysim, by="Subject")

## linear model
## sim_wavg = the average similarity with the subject's own mean fixation map.
## sim_within = average similarity among repetitions after partialing out the average_map
lm.1 <- lm(acc ~ sim_wavg, data=dfx, subset=Subject < 300 & Subject != 16)
lm.2 <- lm(acc ~ sim_within, data=dfx, subset=Subject < 300 & Subject != 16)
lm.3 <- lm(acc ~ sim_wavg + sim_within, data=dfx, subset=Subject < 300 & Subject != 16)
lm.4 <- lm(acc ~ sim_wdiff, data=dfx, subset=Subject < 300 & Subject != 16)


lm.4 <- lm(acc ~ sd_fix, data=dfx, subset=Subject < 300 & Subject != 16)
lm.5 <- lm(acc ~ total_fix, data=dfx, subset=Subject < 300 & Subject != 16)
lm.6 <- lm(acc ~ slope, data=dfx, subset=Subject < 300 & Subject != 16)


## plot
ggplot(aes(sim_within, acc), data=dfx) + geom_point() + stat_smooth(method=rlm)

ggplot(aes(sim_wavg, acc), data=dfx) + geom_point() + stat_smooth(method=rlm)

