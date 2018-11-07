library(dplyr)
library(tibble)
library(tidyr)
library(imager)
library(mgcv)
library(ggplot2)
library(memoise)


exclude_subs <- c(28, 32, 109)


get_mask <- function(im, transpose=TRUE, rev=TRUE) {
  fname <- paste0("~/Dropbox/Jordana_experiments/Jordana_saliency_study/images/Mat_", gsub("jpeg", "jpeg.rds", im))

  mask <- if (file.exists(fname)) {
    message("getting fname", fname)
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


  subject_dens <- density_by(study_tab, groups=c("Subject"), xbounds=c(0,800), ybounds=c(0,600),
                           outdim=c(80,60), duration_weighted=TRUE, sigma=80)


  ## construct heatmaps for the study phase, averaged within subjects
  study_dens <- density_by(study_tab, groups=c("ImageNumber", "Subject", "Image"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
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


  mask_tab <- tibble(Image=levels(test_tab$Image)) %>% rowwise() %>% mutate(mask=list(get_mask(Image))) %>% rowwise() %>% do({
    idx <- which(.$mask > 0, arr.ind=TRUE)
    tibble(Image=.$Image, fixgroup=list(fixation_group(idx[,1]*10, idx[,2]*10, duration=1, onset=0)))
  })

  mask_dens <- density_by(mask_tab, groups=c("Image"), xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60),
                          duration_weighted=TRUE, sigma=160)


  d1 <- mask_dens$density[[1]]
  cbias <- list(x=d1$x, y=d1$y, z=study_dens_avg)
  mask_dens <- mask_dens %>% rowwise %>% mutate(centerbias=list(cbias),
                                                mod_mask_dens = list(list(x=density$x, y=density$y, z= density$z * sqrt(study_dens_avg)))) %>%
    rename(mask_density=density)


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






res <- do.call(rbind, lapply(as.character(test_dens_20$Image), function(im) {
  #mask1 <- get_mask(im)
  mask1 <- filter(mask_dens, Image ==im)$density[[1]]$z

  #mask2 <- get_mask(im, transpose=FALSE, rev=FALSE)
  #mask3 <- get_mask(im, transpose=TRUE, rev=FALSE)
  #mask4 <- get_mask(im, transpose=FALSE, rev=TRUE)
  id <- which(as.character(test_dens_20$Image) == im)
  c1 <- cor(as.vector(mask1), as.vector(test_dens_20$density[[id]]$z), method="spearman")
  #c2 <- cor(as.vector(mask2), as.vector(test_dens_20$density[[id]]$z), method="spearman")
  #c3 <- cor(as.vector(mask3), as.vector(test_dens_20$density[[id]]$z), method="spearman")
  #c4 <- cor(as.vector(mask4), as.vector(test_dens_20$density[[id]]$z),method="spearman")
  data.frame(im=im, c1=c1)
}))



binned_similarity <- function(min_onset, max_onset, type=c("mask", "fixations"), method="cosine") {
  type <- match.arg(type)
  ## create table for each test trial
  pctest_binned <- pctest %>% filter(FixOffset >= min_onset & FixOffset < max_onset)

  test_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                        groupvar=c("Image", "Subject"), data=pctest_binned,
                        clip_bounds=c(112, (112+800), 684, 84),
                        vars=c("ImageVersion", "Saliency", "Accuracy",
                               "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber", "testdelayOnset"))

  test_tab$Image_Subj <- paste0(test_tab$Subject, "_", test_tab$ImageNumber)
  test_tab$Match <- test_tab$ImageRepetition

  test_dens <- density_by(test_tab, groups=c("ImageNumber", "Subject", "Image"),xbounds=c(0,800), ybounds=c(0,600),
                          outdim=c(80,60), duration_weighted=TRUE, sigma=80) %>% rename(test_density=density)

  test_dens$Image_Subj <- paste0(test_dens$Subject, "_", test_dens$ImageNumber)

  test_dens <- inner_join(test_dens, stud$study_dens, by="Image_Subj") %>% rename(Image=Image.x)
  test_dens <- inner_join(test_dens, stud$mask_dens, by="Image")

  #test_dens$Image <- test_tab$Image

  ## compute similarity between each trial and study image derived from group average
  if (type == "fixations") {
    test_sim <- template_similarity(stud$study_dens, test_dens, "Image_Subj", permutations=50, method=method) %>%
      mutate(min_onset=min_onset, max_onset=max_onset)
  } else {
    test_sim <- template_similarity(stud$mask_dens, test_dens, "Image", permutations=10, method=method) %>%
      mutate(min_onset=min_onset, max_onset=max_onset)
  }


  test_sim <- inner_join(test_sim, test_tab, by="Image_Subj")

}

library(purrr)
interval=250
bin_onsets <- seq(0, 3500, by=interval)
test_sim_binned <- bin_onsets %>% map(~ binned_similarity(., . + interval, type="mask", method="jaccard")) %>%
  map_df(bind_rows) #%>%
  #select(-fixgroup.x, -fixgroup.y, -density)


write.table(test_sim_binned, "test_sim_binned.txt", row.names=FALSE)



library(ggplot2)
binned_saliency <- test_sim_binned %>% group_by(min_onset, ImageRepetition,Duration, Saliency) %>% summarize(eye_sim=mean(eye_sim),
                                                                                                    eye_sim_diff=mean(eye_sim_diff))

qplot(min_onset, eye_sim, colour=factor(Saliency), data=binned_saliency, geom=c("point", "line"))

binned_saliency_oldnew <- test_sim_binned %>% group_by(min_onset, Match, Accuracy) %>% summarize(eye_sim_diff=mean(eye_sim_diff),
                                                                                   perm_sim=mean(perm_sim), eye_sim=mean(eye_sim))

qplot(min_onset, eye_sim_diff, colour=factor(Accuracy), facets= . ~ Match,
      data=subset(binned_saliency_oldnew), geom=c("point", "line"))

binned_saliency_acc_oldnew <- test_sim_binned %>% group_by(min_onset, Match, Saliency, Accuracy) %>% summarize(eye_sim_diff=mean(eye_sim_diff),
                                                                                                 perm_sim=mean(perm_sim), eye_sim=mean(eye_sim))

qplot(min_onset, eye_sim_diff, colour=factor(Accuracy), facets= Saliency ~ Match,
      data=subset(binned_saliency_acc_oldnew), geom=c("point", "line"))

