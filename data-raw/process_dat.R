library(dplyr)
library(eyesim)

df1 <- read.table("data-raw/study_fix_report_all.xls", header=TRUE) %>% filter(imagestudy != "UNDEFINEDnull") %>%
  mutate(
  Subject=readr::parse_number(as.character(RECORDING_SESSION_LABEL)),
  ImageNumber=imagenumberstudy,
  Trial=TRIAL_INDEX,
  Block=blockstudy,
  Version=sapply(strsplit(as.character(imagestudy), "_"), "[[", 2),
  ImageVersion=paste0(imagenumberstudy, "_", sapply(strsplit(as.character(imagestudy), "_"), "[[", 2))
)

## create table for each test trial
wynn_study <- eye_table("CURRENT_FIX_X", "CURRENT_FIX_Y", duration="CURRENT_FIX_DURATION", onset="CURRENT_FIX_START",
                      groupvar=c("ImageNumber", "Subject"), data=df1,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion",
                             "Version", "Trial", "ImageNumber", "Block"))


wynn_study <- wynn_study %>% mutate(ImageNumberS = paste0(Subject, "_", ImageNumber))


wynn_study_image <- eye_table("CURRENT_FIX_X", "CURRENT_FIX_Y", duration="CURRENT_FIX_DURATION", onset="CURRENT_FIX_START",
                        groupvar=c("ImageVersion"), data=df1,
                        clip_bounds=c(112, (112+800), 684, 84),
                        vars=c("ImageVersion",
                               "Version", "Trial", "ImageNumber", "Block"))






usethis::use_data(wynn_study, overwrite=TRUE)
usethis::use_data(wynn_study_image, overwrite=TRUE)





df2 <- read.table("data-raw/testdelay_fix_report_all.xls", header=TRUE) %>% mutate(
  Subject=readr::parse_number(as.character(RECORDING_SESSION_LABEL)),
  Saliency=saliencytest,
  Trial=TRIAL_INDEX,
  ImageNumber=imagenumbertest,
  Probetype=repetitiontest,
  Duration=durationtest,
  Accuracy=Accuracy,
  Version=sapply(strsplit(as.character(imagetest), "_"), "[[", 2),
  ImageVersion=paste0(imagenumbertest, "_", sapply(strsplit(as.character(imagetest), "_"), "[[", 2))
)



## create table for each test trial
wynn_test <- eye_table("CURRENT_FIX_X", "CURRENT_FIX_Y", duration="CURRENT_FIX_DURATION", onset="CURRENT_FIX_START",
                      groupvar=c("ImageNumber", "Subject"), data=df2,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "Version", "Trial", "Duration", "ImageNumber"))

wynn_test <- wynn_test %>% mutate(ImageNumberS = paste0(Subject, "_", ImageNumber))


wynn_test_image <- eye_table("CURRENT_FIX_X", "CURRENT_FIX_Y", duration="CURRENT_FIX_DURATION", onset="CURRENT_FIX_START",
                       groupvar=c("ImageVersion"), data=df2,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion", "Saliency", "Accuracy",
                              "Version", "Trial", "Duration", "ImageNumber"))



mids <- match(wynn_test$ImageNumberS, wynn_study$ImageNumberS)
mids2 <-  match(wynn_study$ImageNumberS, wynn_test$ImageNumberS)

ptype <- ifelse(wynn_test$ImageVersion == wynn_study$ImageVersion[mids], "old", "lure")
wynn_test$probe_type=ptype

usethis::use_data(wynn_test, overwrite=TRUE)



sinkdist <- function(x1, x2, xdenom=1000, ydenom=1000, tdenom=3000, tweight=.1, lambda=120) {
  xy1 <- cbind(x1$x/xdenom, x1$y/ydenom)
  xy2 <- cbind(x2$x/xdenom, x2$y/ydenom)

  spd <- proxy::dist(xy1,xy2)
  td <- proxy::dist(x1$onset/tdenom, x2$onset/tdenom)
  d <- spd + tweight*td

  T4transport::sinkhornD(d,wx=x1$duration, wy=x2$duration, lambda=lambda)$distance
  #T4transport::wassersteinD(d,wx=x1$duration, wy=x2$duration)$distance
}

idx <- which(wynn_study$ImageVersion == "1_A")
N=2000
tweight=0

ret <- do.call(rbind, lapply(c(0,.01,.1,.2, .5, .75, 1), function(tweight) {
  do.call(rbind, lapply(idx, function(i) {
    print(i)
    d <- sapply(1:N, function(j) sinkdist(wynn_study$fixgroup[[i]],
                                     wynn_study$fixgroup[[j]],
                                     tweight=tweight, lambda=50))

    keep <- rep(TRUE, N)
    keep[i] <- FALSE
    data.frame(d=d, ImageVersion=wynn_study$ImageVersion[1:N], keep=keep, tweight=tweight)
  }))
}))


