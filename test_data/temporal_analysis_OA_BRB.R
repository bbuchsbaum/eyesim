library(imager)
library(dplyr)
library(tibble)
library(tidyr)
library(mgcv)
library(ggplot2)
library(memoise)
library(eyesim)
library(wesanderson)
#need to run dplyr after imager (or any other package that contains plyr) because we want to make sure we are using
#"mutate" from the dplyr package and not the plyr package

exclude_subs <- c(28, 32, 109, 10, 103,323,333,329,300,325)

## this function retrieves the "input mask" associated with an given image.
## the mask contains a 1 of the cell is visible, otherwise it contains a 0.
get_mask <- function(im, transpose=TRUE, rev=TRUE) {
  print(im)
  fname <- paste0("~/Dropbox/Jordana_experiments/Jordana_saliency_study/images/Mat_", gsub("jpeg", "jpeg.rds", im))

  mask <- if (file.exists(fname)) {

    m <- readRDS(fname)

    if (rev) {
      m <- apply(m, 1, rev)
    }

    if (transpose) {
      m <- t(m)
    }
    m/sum(m)
  } else {
    mask <- matrix(1, 80,60)
    mask/sum(mask)
  }

  mask
}


## load testdelay fixation data
pctest <- as_tibble(read.csv("test_data/testdelay_fixations.csv")) %>%
  #mutate(fix_onset=FixStartTime - testdelayOnset) %>%
  filter(Image != "." & !(Subject %in% exclude_subs)& Subject >200) %>% droplevels()

pctest <- subset(pctest, pctest$FixDuration >80)
#pctest <- subset(pctest, pctest$Subject < 300)

test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                      groupvar=c("Image"), data=pctest,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber", "testdelayOnset"))

prepare_study <- function() {
  pcstudy <- as_tibble(read.csv("test_data/study_fixations_all.csv")) %>%
  filter(Image != "." & !(Subject %in% exclude_subs) & Subject >200) %>% droplevels()

  pcstudy <- subset(pcstudy, pcstudy$FixDuration >80)

  pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


  ## create table for each study trial (Subject/Image)
  study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject", "Block"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))


  ## construct heatmaps for the study phase, averaged within subjects
  study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject"),
                         xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                         duration_weighted=TRUE, sigma=80, result_name="study_density")

  ## construct heatmaps for the average fixation map for each subject
  study_dens_subj_avg <- density_by(study_tab, groups=c("Subject"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                                  duration_weighted=TRUE, sigma=80, result_name="subject_density")

  ## construct heatmaps for the image version, averaged over subjects
  study_dens_all <- density_by(study_tab, groups=c("ImageVersion"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                               duration_weighted=TRUE, sigma=80, result_name="image_density")

  ## compute grand mean density
  study_dens_avg <- Reduce("+", lapply(study_dens$study_density, function(x) x$z))/length(study_dens)
  study_dens_avg <- study_dens_avg/sum(study_dens_avg)
  study_dens$Image_Subj <- paste0(study_dens$Subject, "_", study_dens$ImageNumber)

  dens_avg_subj <- study_dens %>% group_by(Subject) %>% do({
    dens <- Reduce("+", lapply(.$density, function(x) x$z))/nrow(.)
    tibble(Subject=.$Subject[1], dens=list(dens))
  })

  ## load mask and create pseudo-fixations over the visible squares in each mask image.
  ## in other words for every pixel in the mask that is visible add a fixation at that location.
  ## This can then be used as a control model that assumes that subjects uniformly fixate the visible cells.
  mask_tab <- tibble(Image=sort(unique(test_tab$Image))) %>% rowwise() %>% dplyr::mutate(mask=list(get_mask(Image))) %>% rowwise() %>% do({
    #browser()
    idx <- which(.$mask > 0, arr.ind=TRUE)
    tibble(Image=.$Image, fixgroup=list(fixation_group(idx[,1]*10, idx[,2]*10, duration=rep(1, nrow(idx)), onset=rep(0, nrow(idx)))))
  })

  ## construct density maps for each mask
  ## given the pseudo-fixations create a density map for each mask. This effectively smooths the mask.
  mask_dens <- density_by(mask_tab, groups=c("Image"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                          duration_weighted=TRUE, sigma=80, result_name="mask_density")

  d1 <- mask_dens$mask_density[[1]]
  cbias <- gen_density(x=d1$x, y=d1$y, z=study_dens_avg)

  mask_dens <- mask_dens %>% rowwise() %>% mutate(centerbias=list(cbias))
                                                  #mod_mask_dens = list(gen_density(x=density$x,
                                                  #                        y=density$y,
                                                  #                        z= density$z * sqrt(study_dens_avg))))
    #mutate(mod_mask_dens = list(list(x=mod_mask_dens$x,
    #                                 y=mod_mask_dens$y,
    #                                 z=mod_mask_dens$z/sum(mod_mask_dens$z))))


  list(
    study_dens=study_dens,
    study_dens_subj_avg=study_dens_subj_avg,
    study_dens_all=study_dens_all,
    study_dens_avg=study_dens_avg,
    study_dens_avg_subj=dens_avg_subj,
    mask_dens=mask_dens
  )

}

stud <- prepare_study()

## construct test_tab
test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                      groupvar=c("Image", "Subject"), data=filter(pctest, Subject<400),
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber", "testdelayOnset"))

## create Subject_ImageNumber variable
test_tab$Image_Subj <- paste0(test_tab$Subject, "_", test_tab$ImageNumber)
test_tab$Match <- test_tab$ImageRepetition

## density map for grouped by Subject and Image Number
test_dens <- density_by(test_tab, groups=c("ImageNumber", "Subject"), xbounds=c(0,800), ybounds=c(0,600),
                        outdim=c(80,60), duration_weighted=TRUE, sigma=80, result_name="test_density",
                        keep_vars=c("Image", "Saliency", "Duration", "ImageRepetition", "Accuracy", "ImageVersion"))


test_dens$Image_Subj <- paste0(test_dens$Subject, "_", test_dens$ImageNumber)

test_dens <- left_join(test_dens, stud$study_dens, by="Image_Subj", suffix=c(".test", ".stud")) #%>% rename(Image=Image.x)
test_dens <- left_join(test_dens, stud$mask_dens, by="Image")
test_dens <- left_join(test_dens, stud$study_dens_all, by="ImageVersion")


### ???? what is this, it fails for BRB, commenting out
## test_dens <- subset(test_dens, test_dens(select=-c()))

library(purrr)
## create a mask desnity map that is weighted by the grouped derived image salience.
# image density = heatmap for the image version, averaged over subjects
mod_density <- test_dens %>% pmap(function(mask_density, image_density, ...) {
  #browser()
  xx <- mask_density$z * image_density$z
  xx <- xx/sum(xx)
  gen_density(x=mask_density$x, y=mask_density$y, z=xx)
})


test_dens <- test_dens %>% add_column(mod_density=mod_density) #%>%
  #add_column(nomod_density=nomod_density) %>%
  #add_column(mod_density2=mod_density2) %>%
  #add_column(cell_density=cell_density)

test_dens2 <- subset(test_dens, select=c("mod_density","Image_Subj"))

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################


## create a mask density map that is weighted by the square root of the centerbias.
cell_density <- test_dens %>% pmap(function(mask_density, centerbias, ...) {
  xx <- mask_density$z * centerbias$z^(1/2)
  xx <- xx/sum(xx)
  gen_density(x=mask_density$x, y=mask_density$y, z=xx)
})

## compute density maps for the non-visible area.
nomod_density <- test_dens %>% pmap(function(mask_density, study_density, Saliency, ...) {
  q=quantile(mask_density$z, max((1-Saliency/100) - .04,.0001))
  d <- mask_density$z
  d <- ifelse(d > q, 0, 1-d)
  xx <- d * study_density$z
  xx <- xx/sum(xx)
  gen_density(x=mask_density$x, y=mask_density$y, z=xx)
})

test_dens <- test_dens %>% ## (BRB removed this line) add_column(mod_density=mod_density) %>%
  add_column(nomod_density=nomod_density) %>%
  #add_column(mod_density2=mod_density2) %>%
  add_column(cell_density=cell_density)

library(purrr)

## template sample extracts the density for any arbitrary time point in a trial. It simply extract the value of the density map for the fixation at time t.
## we repeat this for each of 5 density models.
# at each time-point (from 0 to 3800), we extract the value at the location fixated at that point in time (fixations are contained in test_dens$fixgroup.test).
# The value we extract is one the models defined below.

## this is different than correlating two heat haps

## 1. the mask density modulated by the center bias.
res <- template_sample(test_dens, "mod_density", fixgroup="fixgroup.test", time=seq(0,3800,by=50), outcol="mod_sam")
## 2. the density of the non-visible area
res <- template_sample(res, "nomod_density", fixgroup="fixgroup.test", time=seq(0,3800,by=50), outcol="nomod_sam")
## 3. a simple center bias model
res <- template_sample(res, "centerbias", fixgroup="fixgroup.test", time=seq(0,3800,by=50), outcol="center_sam")
## 4. the study density (standard model)
res <- template_sample(res, "study_density", fixgroup="fixgroup.test", time=seq(0,3800,by=50), outcol="study_sam")
## 5. the image density averaged over subjects
res <- template_sample(res, "image_density", fixgroup="fixgroup.test", time=seq(0,3800,by=50), outcol="image_sam")
## 6. mask density multiplied by centerbias
res <- template_sample(res, "cell_density", fixgroup="fixgroup.test", time=seq(0,3800,by=50), outcol="cell_sam")

res <- res %>% select(Subject.test, Image, Saliency, Duration, ImageRepetition, Accuracy, mod_sam, nomod_sam, center_sam, study_sam,image_sam, cell_sam) %>%
  rowwise() %>% do( {
    tibble(Image=.$Image, Saliency=.$Saliency, Duration=.$Duration, ImageRepetition=.$ImageRepetition,
           Subject=.$Subject.test,
           Accuracy=.$Accuracy, time=.$mod_sam$time,
           mod_sam=.$mod_sam$z,
           nomod_sam=.$nomod_sam$z,
           center_sam=.$center_sam$z,
           cell_sam=.$cell_sam$z,
           study_sam=.$study_sam$z,
           image_sam=.$image_sam$z)
  })

setwd('/Users/jwynn/Dropbox/New_pc_behavioural_data')
write.csv(res, "res_OA.csv")
res <- read.csv("res.csv")

## compute mean summaries of the density over 4 variables
bsum <- res %>% group_by(Duration, Subject,ImageRepetition, Accuracy, time) %>% summarize(
  mod_sam=mean(mod_sam, na.rm=TRUE),
  nomod_sam=mean(nomod_sam, na.rm=TRUE),
  cell_sam=mean(cell_sam, na.rm=TRUE),
  center_sam=mean(center_sam, na.rm=TRUE),
  image_sam=mean(image_sam, na.rm=TRUE),
  study_sam=mean(study_sam, na.rm=TRUE))

#colnames(bsum)[which(names(bsum) == "study_sam")] <- "Subject sim study"
#colnames(bsum)[which(names(bsum) == "image_sam")] <- "Image sim study"
#colnames(bsum)[which(names(bsum) == "mod_sam")] <- "Image sim test"

#change names of measures

bsum <- gather(bsum, measure, value, -Subject, -Duration, -ImageRepetition, -Accuracy, -time)
bsum <- bsum %>% mutate(Acc=ifelse(Accuracy==0, "Incorrect", "Correct"))

#bsum_avg <- bsum %>% filter(Saliency <= 101) %>%
#               group_by(time, Acc, ImageRepetition,Duration,measure) %>%
#               summarize(value=mean(value, na.rm=TRUE))

bsum2 <- subset(bsum, bsum$time > 0)
bsum2 <- subset(bsum2, bsum2$time < 3340)
### median latency for the first test saccade is 286 so may want to remove bottom quartile or so

bsum2 <- subset(bsum2, bsum2$measure== "study_sam" |
                  bsum2$measure== "image_sam"|
                  bsum2$measure== "mod_sam")

bsum_250 <- subset(bsum2, bsum2$Duration ==250)
bsum_250$time <- bsum_250$time-250
bsum_500 <- subset(bsum2, bsum2$Duration ==500)
bsum_500$time <- bsum_500$time-500
bsum_750 <- subset(bsum2, bsum2$Duration ==750)
bsum_750$time <- bsum_750$time-750

bsum_new_times <- rbind(bsum_250,bsum_500,bsum_750)

#reorder
bsum_500$ImageRepetition <- factor(bsum_500$ImageRepetition, levels=c("old","new"), labels=c("Old", "Lure"))
bsum2$ImageRepetition <- factor(bsum2$ImageRepetition, levels=c("old","new"), labels=c("Old", "Lure"))
bsum_new_times$ImageRepetition <- factor(bsum_new_times$ImageRepetition, levels=c("old","new"), labels=c("Old", "Lure"))

bsum_new_times2 <- aggregate(value~Subject+ImageRepetition+Acc+time+measure, data=bsum_new_times,FUN=mean)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$time < 3050)
bsum_new_times2$Age <- ifelse(bsum_new_times2$Subject > 200, "OA","YA")
#bsum_500$time <- as.factor(bsum_500$time)
#bsum_new_times2$value <- range01(bsum_new_times2$value)
#range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#remove YA subs who got <50% correct (also <2.5sd from mean)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 28)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 32)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 109)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 13)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 10)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 103)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 323)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 333)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 329)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 300)
bsum_new_times2 <- subset(bsum_new_times2, bsum_new_times2$Subject != 325)

bsum_new_times2$measure <- factor(bsum_new_times2$measure, levels=c("study_sam","image_sam","mod_sam"), labels=c("Gaze","Image", "Probe"))

plot <- ggplot(bsum_new_times2, aes(x=time, y=value,
  colour=measure, linetype=Acc)) +
  geom_smooth(method=gam, se=F, formula=y ~ s(x, k=15), size=1.5) +
  geom_vline(xintercept = 0, color="grey40",size=.5) +
  #geom_vline(xintercept = 500, color="purple",size=20, alpha=.2) +
  #geom_text(aes(x=500, y=0.0011, label = "Delay Onset"), colour= "grey40", angle=0, vjust=1, text=element_text(size=10)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_x_continuous (breaks = c(0, 1000, 2000, 3000),
                     labels = paste0(c("Delay Onset","1000","2000","3000"))) +
  #scale_x_continuous(breaks=c(500, 1500, 2500, 3500)) +
  scale_color_manual(values=wes_palette("Darjeeling1"))+
  theme(axis.text.x = element_text(size=12, color="gray30", angle=-35, vjust=-.05)) +
  theme(axis.text.y = element_text(size=12, color="gray30")) +
  theme(axis.title.x = element_text(size=16, color="black")) +
  theme(axis.title.y = element_text(size=16, color="black")) +
  theme(legend.title=element_text(size=16, color="black")) +
  theme(legend.text=element_text(size=16, color="gray30")) +
  theme(legend.position = "right") +
  theme(legend.title=element_blank()) +
  theme(legend.text=element_text(size=12, color="gray30")) +
  theme(axis.ticks=element_blank()) +
  theme(axis.line.x=element_blank()) +
  theme(axis.line.y=element_blank()) +
  xlab("Time") +
  ylab("Density")

plot + facet_grid(Age ~ ImageRepetition) + guides(linetype=FALSE) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16))
  #theme(strip.background = element_rect(fill="white"))


### ANOVA
library(afex)
library(lsmeans)

bsum_anova <- subset(bsum_new_times2, bsum_new_times2$time>-.99)

bsum_anova$groupA <- ifelse(bsum_anova$time >-.99 &
                           bsum_anova$time < 999,1,
                           ifelse(bsum_anova$time > 999 &
                           bsum_anova$time < 1999,2,
                           ifelse(bsum_anova$time >1999 &
                           bsum_anova$time < 3001,3,0)))

bsum_anova$groupB <- ifelse(bsum_anova$time >-.99 &
                             bsum_anova$time < 499,1,
                           ifelse(bsum_anova$time > 499 &
                                    bsum_anova$time < 999,2,
                                  ifelse(bsum_anova$time >999 &
                                           bsum_anova$time < 1499,3,
                                         ifelse(bsum_anova$time >1499 &
                                                  bsum_anova$time < 1999,4,
                                                ifelse(bsum_anova$time >1999 &
                                                         bsum_anova$time < 2499,5,
                                                       ifelse(bsum_anova$time >2499 &
                                                                bsum_anova$time < 3050,6,0))))))


#bsum_anova2 <- aggregate(value~Subject+ImageRepetition+group+measure+Age, data=bsum_anova, FUN=mean)
#bsum_anova2$value <- bsum_anova2$V1

gp_A1 <- subset(bsum_anova, bsum_anova$groupA==1)
gp_A1 <- aggregate(value~Subject+measure+ImageRepetition+Age, data=gp_A1, FUN=mean)
gp_A2 <- subset(bsum_anova, bsum_anova$groupA==2)
gp_A2 <- aggregate(value~Subject+measure+ImageRepetition+Age, data=gp_A2, FUN=mean)
gp_A3 <- subset(bsum_anova, bsum_anova$groupA==3)
gp_A3 <- aggregate(value~Subject+measure+ImageRepetition+Age, data=gp_A3, FUN=mean)

gp_B1 <- subset(bsum_anova, bsum_anova$groupB==1)
gp_B1 <- aggregate(value~Subject+measure+ImageRepetition+Age+Acc, data=gp_B1, FUN=mean)
gp_B2 <- subset(bsum_anova, bsum_anova$groupB==2)
gp_B2 <- aggregate(value~Subject+measure+ImageRepetition+Age+Acc, data=gp_B2, FUN=mean)
gp_B3 <- subset(bsum_anova, bsum_anova$groupB==3)
gp_B3 <- aggregate(value~Subject+measure+ImageRepetition+Age+Acc, data=gp_B3, FUN=mean)
gp_B4 <- subset(bsum_anova, bsum_anova$groupB==4)
gp_B4 <- aggregate(value~Subject+measure+ImageRepetition+Age+Acc, data=gp_B4, FUN=mean)
gp_B5 <- subset(bsum_anova, bsum_anova$groupB==5)
gp_B5 <- aggregate(value~Subject+measure+ImageRepetition+Age+Acc, data=gp_B5, FUN=mean)
gp_B6 <- subset(bsum_anova, bsum_anova$groupB==6)
gp_B6 <- aggregate(value~Subject+measure+ImageRepetition+Age+Acc, data=gp_B6, FUN=mean)

#we want to compare gp_B1 image sam to gp_A2 or gp_A3 mod sam
#get median image sam for gp_B1 and median mod sam for gp_A2/A3 and do median split
#compare percorr for people who are in the top group for each or for both

gp_1_OA <- subset(gp_A1, gp_A1$Age=="OA")
gp_2_OA <- subset(gp_A2, gp_A2$Age=="OA")
gp_3_OA <- subset(gp_A3, gp_A3$Age=="OA")
gp_1_YA <- subset(gp_A1, gp_A1$Age=="YA")
gp_2_YA <- subset(gp_A2, gp_A2$Age=="YA")
gp_3_YA <- subset(gp_A3, gp_A3$Age=="YA")


gp_1_OA_lure_image <- subset(gp_1_OA, gp_1_OA$ImageRepetition=="Lure" &
                               gp_1_OA$measure=="image_sam")
gp_1_OA_lure_image$image_sam <- gp_1_OA_lure_image$value

gp_2_OA_lure_mod <- subset(gp_3_OA, gp_3_OA$ImageRepetition=="Lure" &
                               gp_3_OA$measure=="mod_sam" &
                               gp_3_OA$Acc=="Correct")
gp_2_OA_lure_mod$mod_sam <- gp_2_OA_lure_mod$value

lure_autocor <- merge(gp_1_OA_lure_image, gp_2_OA_lure_mod, by=c("Subject"))
### everything seems to be positively correlated with everything

Mixed.aov.1<-aov_car(value ~ measure*Acc*ImageRepetition*Age + Error(Subject/(measure*Acc*ImageRepetition)), data=gp_A3,
                     anova_table = list(es = "pes",correction = "none"))

Mixed.aov.1
summary(Mixed.aov.1)

##Main effects
emmeans(Mixed.aov.1, ~measure)

##Interactions
Mixed_Fitted_measureAcc<-emmeans(Mixed.aov.1, ~Acc|measure)
Mixed_Fitted_measureAcc

##pairwise comparison with bonferroni correction.
pairs(Mixed_Fitted_measureAcc, adjust = "bon")

#########################################################
qplot(time, value, colour=measure, linetype=factor(Accuracy),
      data=subset(bsum_500, measure %in% c("Within-subject study-delay similarity", "Across-subject study-delay similarity", "Across-subject test-delay similarity")
      & time > 50),
      facets = ~ ImageRepetition, geom=c("smooth"), se=FALSE, method=gam, formula=y ~ s(x, k=15))

# qplot(time, value, colour=measure, linetype=factor(Accuracy),
#       data=subset(bsum, measure %in% c("Within-subject study-delay similarity", "Across-subject study-delay similarity", "Across-subject test-delay similarity")
#       & time > 50),
#       facets = Duration ~ ImageRepetition, geom=c("line"))

bsum_avg_old <- subset(bsum_avg, bsum_avg$ImageRepetition=="old")
bsum_avg_new <- subset(bsum_avg, bsum_avg$ImageRepetition=="new")

# qplot(time, value, colour=measure, linetype=Acc,
#       data=subset(bsum, measure %in% c("Within-subject study-delay similarity", "Across-subject study-delay similarity", "Across-subject test-delay similarity")
#       & time > 50),
#       facets = .  ~ Duration, geom=c("line")) + theme_cowplot()
#

## do the same, but also add "Subject"
bsumS <- res %>% group_by(Saliency, ImageRepetition, Accuracy, time, Subject) %>% summarize(
  mod_sam=mean(mod_sam, na.rm=TRUE),
  mod_sam2=mean(mod_sam2, na.rm=TRUE),
  nomod_sam=mean(nomod_sam, na.rm=TRUE),
  center_sam=mean(center_sam, na.rm=TRUE),
  image_sam=mean(image_sam, na.rm=TRUE),
  study_sam=mean(study_sam, na.rm=TRUE))

bsumS <- gather(bsumS, measure, value, -Saliency, -ImageRepetition, -Accuracy, -time)

library(tidyr)

