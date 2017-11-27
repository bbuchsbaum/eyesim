library(tibble)
## read fixation report tables (does not included creation of study density map)
pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_study_input.csv"))
pcdelay <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delay_input.csv"))
pcdelay_test <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delaytest_input.csv"))
pcdelay_test$onset <- rep(0, nrow(pcdelay_test))
#pcdelay_nofirst <- pcdelay %>% group_by(Subject, ImageNumber, ImageSet) %>% dplyr::filter(row_number() != 1)
pcdelay_nofirst <- pcdelay_test %>% group_by(Subject, ImageNumber, ImageSet)

## create eye fixation table
delay_nof_tab <- eye_table(x="FixX", y="FixY", duration="FixDuration", onset="onset",
                                  vars=c("Saliency", "Duration", "Accuracy"),
                                  groupvar=c("Subject", "ImageNumber", "ImageSet"), data=pcdelay_test)

## create density map for each "group" using a sigma of 100
delay_nof_density <- density_by(delay_nof_tab, groups=c("Subject", "ImageNumber", "ImageSet"), sigma=100)

## add unique identifier for group
delay_nof_density$ImageVersion <- paste0(delay_nof_density$ImageNumber, "_", delay_nof_density$ImageSet)

## compute similarity between each trial and study image derived from group average
delay_nof_sim <- template_similarity(study_density, delay_nof_density, "ImageVersion")

## add back in per-trial experimental variables of interest
delay_nof_sim$saliency <- delay_nof_tab$Saliency
delay_nof_sim$duration <- delay_nof_tab$Duration
delay_nof_sim$accuracy <- delay_nof_tab$Accuracy

## summarize data over variables and plot
sim_means <- delay_nof_sim %>% group_by(saliency, duration) %>% dplyr::summarize(eye_sim=mean(eye_sim))
qplot(saliency, eye_sim, colour=factor(duration), data=sim_means, geom=c("point", "line"))

## by accuracy
sim_means_acc <- delay_nof_sim %>% group_by(saliency, accuracy) %>% dplyr::summarize(eye_sim=mean(eye_sim))
qplot(saliency, eye_sim, colour=factor(accuracy), data=sim_means_acc, geom=c("point", "line"))

## correct only
sim_means_correct <- delay_nof_sim %>% filter(accuracy==1) %>% group_by(saliency, duration) %>% dplyr::summarize(eye_sim=mean(eye_sim))
qplot(saliency, eye_sim, colour=factor(duration), data=sim_means_correct, geom=c("point", "line"))


## 1
pcdelay_test <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delaytest_input.csv"))

## 2
delay_test_tab <- create_eye_table(x="FixX", y="FixY", duration="FixDuration", onset="onset", groupvar=c("Subject", "ImageNumber", "ImageSet"), data=pcdelay)

## 3 ...
