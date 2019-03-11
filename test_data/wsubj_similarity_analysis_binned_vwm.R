library(dplyr)
library(tibble)
library(tidyr)
library(imager)
library(mgcv)
library(ggplot2)
library(memoise)


exclude_subs <- c(28, 32, 109)



## this function retrieves the "input mask" associated with an given image.
## the mask contains a 1 of the cell is visible, otherwise it contains a 0.
get_mask <- function(im, transpose=TRUE, rev=TRUE) {
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
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/testdelay_fixations.csv")) %>%
  #mutate(fix_onset=FixStartTime - testdelayOnset) %>%
  filter(Image != "." & !(Subject %in% exclude_subs)) %>% droplevels()


test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                      groupvar=c("Image"), data=pctest,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber", "testdelayOnset"))



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
  mask_tab <- tibble(Image=levels(test_tab$Image)) %>% rowwise() %>% mutate(mask=list(get_mask(Image))) %>% rowwise() %>% do({
    idx <- which(.$mask > 0, arr.ind=TRUE)
    tibble(Image=.$Image, fixgroup=list(fixation_group(idx[,1]*10, idx[,2]*10, duration=1, onset=0)))
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
                      groupvar=c("Image", "Subject"), data=filter(pctest, Subject<300),
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
test_dens <- test_dens %>% rename(Subject=Subject.test)
test_dens <- left_join(test_dens, stud$mask_dens, by="Image")
test_dens <- left_join(test_dens, stud$study_dens_all, by="ImageVersion")
test_dens <- left_join(test_dens, stud$study_dens_subj_avg, by="Subject")

library(purrr)
## create a mask desnity map that is weighted by the grouped derived image salience.
mod_density <- test_dens %>% pmap(function(mask_density, image_density, ...) {
  xx <- mask_density$z * image_density$z
  xx <- xx/sum(xx)
  gen_density(x=mask_density$x, y=mask_density$y, z=xx)
})



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

test_dens <- test_dens %>% add_column(mod_density=mod_density) %>%
  add_column(nomod_density=nomod_density) %>%
  add_column(mod_density2=mod_density2) %>%
  add_column(cell_density=cell_density)

library(purrr)

## template sample extracts the density for any arbitrary time point in a trial. It simply extract the value of the density map for the fixation at time t.
## we repeat this for each of 5 density models.

## 1. the mask density modulated by the center bias.
res <- template_sample(test_dens, "mod_density", fixgroup="fixgroup.test", time=seq(0,3500,by=50), outcol="mod_sam")
## 2. the density of the non-visible area
res <- template_sample(res, "nomod_density", fixgroup="fixgroup.test", time=seq(0,3500,by=50), outcol="nomod_sam")
## 3. a simple center bias model
res <- template_sample(res, "centerbias", fixgroup="fixgroup.test", time=seq(0,3500,by=50), outcol="center_sam")
## 4. the study density (standard model)
res <- template_sample(res, "study_density", fixgroup="fixgroup.test", time=seq(0,3500,by=50), outcol="study_sam")
## 5. the image density averaged over subjects
res <- template_sample(res, "image_density", fixgroup="fixgroup.test", time=seq(0,3500,by=50), outcol="image_sam")
## 6. mask density multiplied by centerbias
res <- template_sample(res, "cell_density", fixgroup="fixgroup.test", time=seq(0,3500,by=50), outcol="cell_sam")
## 7. subject center bias
res <- template_sample(res, "subject_density", fixgroup="fixgroup.test", time=seq(0,3500,by=50), outcol="subject_sam")

res <- res %>% select(Subject, Image, Saliency, Duration, ImageRepetition, Accuracy, mod_sam, nomod_sam, center_sam, subject_sam, study_sam,image_sam, cell_sam) %>%
  rowwise() %>% do( {
    tibble(Image=.$Image, Saliency=.$Saliency, Duration=.$Duration, ImageRepetition=.$ImageRepetition,
           Subject=.$Subject,
           Accuracy=.$Accuracy, time=.$mod_sam$time,
           mod_sam=.$mod_sam$z,
           nomod_sam=.$nomod_sam$z,
           center_sam=.$center_sam$z,
           cell_sam=.$cell_sam$z,
           study_sam=.$study_sam$z,
           subject_sam=.$subject_sam$z,
           image_sam=.$image_sam$z)
  })

##lm.1 <- lmer(nomod_sam ~ Accuracy*ImageRepetition + (1 | Subject) + (0 + Accuracy + ImageRepetition | Subject), data=subset(res, time > 800  & Saliency <= 60))

## compute mean summaries of the density over 4 variables
bsum <- res %>% group_by(Saliency, ImageRepetition, Accuracy, time) %>% summarize(
  mod_sam=mean(mod_sam, na.rm=TRUE),
  nomod_sam=mean(nomod_sam, na.rm=TRUE),
  cell_sam=mean(cell_sam, na.rm=TRUE),
  center_sam=mean(center_sam, na.rm=TRUE),
  subject_sam=mean(subject_sam, na.rm=TRUE),
  image_sam=mean(image_sam, na.rm=TRUE),
  study_sam=mean(study_sam, na.rm=TRUE))

library(cowplot)

bsum <- bsum %>% mutate(Acc=ifelse(Accuracy==0, "incorrect", "correct"))
bsum <- gather(bsum, measure, value, -Saliency, -ImageRepetition, -Acc, -time)
bsum_avg <- bsum %>% filter(Saliency <= 101) %>%
               group_by(time, Acc, ImageRepetition,measure) %>%
               summarize(value=mean(value, na.rm=TRUE))

qplot(time, value, colour=measure, linetype=factor(Accuracy),
      data=subset(bsum, measure %in% c("study_sam", "mod_sam", "cell_sam") & time > 50),
      facets = Saliency ~ ImageRepetition, geom=c("smooth"), se=FALSE, method=gam, formula=y ~ s(x, k=15))

qplot(time, value, colour=measure, linetype=factor(Accuracy),
      data=subset(bsum, measure %in% c("study_sam", "mod_sam", "cell_sam") & time > 50),
      facets = Saliency ~ ImageRepetition, geom=c("line"))

qplot(time, value, colour=measure, linetype=Acc,
      data=subset(bsum_avg, measure %in% c("study_sam", "mod_sam", "subject_sam") & time > 50 & time < 3001),
      facets = .  ~ ImageRepetition, geom=c("line")) + theme_cowplot()


## do the same, but also add "Subject"
bsumS <- res %>% group_by(Saliency, ImageRepetition, Accuracy, time, Subject) %>% summarize(
  mod_sam=mean(mod_sam, na.rm=TRUE),
  nomod_sam=mean(nomod_sam, na.rm=TRUE),
  center_sam=mean(center_sam, na.rm=TRUE),
  image_sam=mean(image_sam, na.rm=TRUE),
  study_sam=mean(study_sam, na.rm=TRUE))

bsumS <- gather(bsumS, measure, value, -Saliency, -ImageRepetition, -Accuracy, -time)

library(tidyr)

