library(dplyr)
library(tibble)
library(tidyr)
library(energy)

exclude_subs <- c(28, 32, 109, 10, 103)

pcstudy <- as_tibble(read.csv("~/Dropbox/New_pc_behavioural_data/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109, 10, 103))) %>% droplevels()


pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))


## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))




plot(study_tab$fixgroup[[29]], bg_image = "~/Dropbox/Jordana_experiments/Jordana_saliency_study/images/111_A_1.jpeg",
     type="contour", show_points=FALSE, bandwidth=c(100,100), bins=50,xlim=c(0,800), ylim=c(0,600), alpha=.12)
