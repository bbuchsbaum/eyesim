library(dplyr)
library(tibble)
library(tidyr)
library(ppcor)
library(MASS)
library(mgcv)

## load test data to get accuracy values (only taking first fixation row of each trial -- see: 'slice(1)')
pctest <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/delay_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()

## divide the Fixations into 10 time bins, then compute the standard deviation of the fixations in x and y dimensions and multiply them together.
pctest_sd <- pctest %>% mutate(time_bin=ntile(FixOffset, n=15)) %>% group_by(Subject, Duration, time_bin) %>% dplyr::summarise(sdxy=sd(FixX)*sd(FixY), tmid=median(FixOffset))

## compute the mean xy standard deviation over subjects.
pctest_sd_dur <- pctest_sd %>% group_by(time_bin, Duration) %>% dplyr::summarise(time=median(tmid), mean_sdxy=mean(sdxy, na.rm=TRUE))

## compute the mean xy standard deviation over subjects, average over Duration
pctest_sd_all <- pctest_sd %>% group_by(time_bin) %>% dplyr::summarise(time=median(tmid), mean_sdxy=mean(sdxy, na.rm=TRUE))


## plot the curve, split by duration
qplot(time, mean_sdxy, colour=factor(Duration), data=pctest_sd_dur) + geom_line()

## plot the grand average curve
qplot(time, mean_sdxy, data=pctest_sd_all) + geom_line()
