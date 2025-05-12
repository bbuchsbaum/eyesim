
# Import YA and OA all ----------------------------------------------------

setwd("~/OneDrive - University of Toronto/PhD/AJ-fMRI/AJ_DataAnalysis/design2")

library(broom)
library(tidyverse)
library(Rmisc)
library(ggpubr)
library(lme4)
library(gridExtra)
library(reshape2)
library(nnls)
library(tweenr)
library(remotes)
library(Matrix)
library(eyesim)
library(gganimate)
#install.packages("patchwork")
library(patchwork)
library(dplyr)
library(readxl)
library(emmeans)

YA_merged_df_eyesim <- read.csv("test_data/YA_merged_df_eyesim5.csv")
OA_merged_df_eyesim <- read.csv("test_data/OA_merged_df_eyesim3.csv")


#YA first clean df

# Remove bad trials

YA_merged_df_eyesim_clean <- YA_merged_df_eyesim[YA_merged_df_eyesim$bad != 1, ]
str(YA_merged_df_eyesim_clean)

YA_merged_df_eyesim_clean <- YA_merged_df_eyesim_clean %>% dplyr::rename(ons=onset)

#filter for images that repeat 4 times (rep 1 - 4) only

YA_filtered_df <- YA_merged_df_eyesim_clean %>%
  group_by(SID,image) %>%
  filter(all(c(1, 2, 3, 4) %in% rep_num)) %>%
  ungroup()

#create eye_table. Note for some subjects (e.g., 1049) repetition for a given image is dropped due to removing bad trials

YA_etab_filtered_original <- eyesim::eye_table("CURRENT_FIX_X", "CURRENT_FIX_Y", "CURRENT_FIX_DURATION", "CURRENT_FIX_START",
                                               groupvar=c("SID", "image","rep_num","run"),
                                               vars=c("ons", "sex", "trial", "set", "miniblock", "TRIAL_START_TIME", "Age_Reaction_Time"),
                                               data=YA_filtered_df,
                                               clip_bounds=c(272, 751, 144,623),
                                               relative_coords=FALSE)

YA_eyedens2 <- eyesim::density_by(YA_etab_filtered_original, groups = c("SID","image", "rep_num", "run"), sigma = 100, xbounds = c(272, 751), ybounds = c(144,623)) # sigma smooths density maps

# Repetitive Similarity  ---------------------------------------------

# Filter for rep 1 vs rep 2 -Younger Adults -------------------------------------------

set.seed(1234)

YA_eyedens2_rep <- YA_eyedens2 %>% mutate(subj_image = factor(paste(SID, image, sep = "_")))

YA_rep1_dens <- YA_eyedens2_rep %>% filter(rep_num == "1") #Filters YA_eyedens to create density data for repetition numbers 1 and 2 (YA_rep1_dens and YA_rep2_dens).

YA_rep2_dens <- YA_eyedens2_rep %>% filter(rep_num == "2")


# Identify mismatches between the two datasets. I had to do this because I would get an error message when running template_similarity:
#Error message:template_similarity: similarity metric is fisherz
#Warning message: In run_similarity_analysis(ref_tab, source_tab, match_on, permutations, did not find matching template map for all source maps. Removing non-matching elements.

common_subj_images <- intersect(YA_rep1_dens$subj_image, YA_rep2_dens$subj_image) #identify mis-match issue
setdiff(YA_rep1_dens$subj_image, YA_rep2_dens$subj_image) # Unique to rep1
setdiff(YA_rep2_dens$subj_image, YA_rep1_dens$subj_image) # Unique to rep2

YA_rep1_dens <- YA_rep1_dens[YA_rep1_dens$subj_image %in% common_subj_images, ]
YA_rep2_dens <- YA_rep2_dens[YA_rep2_dens$subj_image %in% common_subj_images, ]
setdiff(YA_rep1_dens$subj_image, YA_rep2_dens$subj_image) # remove sub-1049 mis-match, had flagged this person previously
setdiff(YA_rep2_dens$subj_image, YA_rep1_dens$subj_image) # remove sub-1049 mis-match

#run template_similarity for 1-2
YA_simres_rep1_2 <- template_similarity(YA_rep1_dens, YA_rep2_dens, match_on="subj_image", permute_on = "SID", method="fisherz", permutations=50)

YA_rep1_dens <- YA_rep1_dens %>% mutate(subj_image_run = factor(paste(SID, image, run, sep = "_")))
YA_rep2_dens <- YA_rep2_dens %>% mutate(subj_image_run = factor(paste(SID, image, run, sep = "_")))

YA_rep1_dens <- YA_rep1_dens %>% mutate(subj_run = factor(paste(SID, run, sep = "_")))
YA_rep2_dens <- YA_rep2_dens %>% mutate(subj_run = factor(paste(SID, run, sep = "_")))

YA_simres_rep1_2_test <- template_similarity(YA_rep1_dens, YA_rep2_dens, match_on="subj_image", permute_on = "subj_run", method ="fisherz", permutations= 50)

# Using aggregate to compute the mean of 'eye_sim', 'perm_sim', and 'eye_sim_diff' for each 'SID' and 'image'
mean_eye_sim1_2_image <- aggregate(eye_sim ~ SID + image, data = YA_simres_rep1_2, FUN = mean)
mean_perm_sim1_2_image <- aggregate(perm_sim ~ SID + image, data = YA_simres_rep1_2, FUN = mean)
mean_eye_sim_diff1_2_image <- aggregate(eye_sim_diff ~ SID +image, data = YA_simres_rep1_2, FUN = mean) #corrected eye similarity is the final repetitive similarity score to go with as it controls for across participant variability

#now aggregate across image
mean_eye_sim1_2 <- aggregate(eye_sim ~ SID, data = YA_simres_rep1_2, FUN = mean)
mean_perm_sim1_2 <- aggregate(perm_sim ~ SID , data = YA_simres_rep1_2, FUN = mean)
mean_eye_sim_diff1_2 <- aggregate(eye_sim_diff ~ SID, data = YA_simres_rep1_2, FUN = mean) # use this (eye_sim_diff as final rep sim value as per Jordana's paper)

# t-test
YA_t.test1_rep1v2 <- t.test(mean_eye_sim1_2$eye_sim) #sig
YA_t.test2_rep1v2 <- t.test(mean_perm_sim1_2$perm_sim) #sig
YA_t.test3_rep1v2 <- t.test(mean_eye_sim_diff1_2$eye_sim_diff) #sig

#par(mfrow=c(1,3))
quartz(width = 7, height = 7)  # For macOS
par(mfrow=c(1,3))

YA_p1_rep1v2 <- hist(mean_eye_sim1_2$eye_sim, main="YA raw eye movement similarity (rep 1 vs rep 2)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram shows the distribution of similarity scores for the matched pairs, where each pair represents the same image during rep 1 and rep 4 phases.

YA_p2_rep1v2 <- hist(mean_perm_sim1_2$perm_sim, main="YA image-permuted eye movement similarity (rep 1 vs rep 2)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

YA_p3_rep1v2 <- hist(mean_eye_sim_diff1_2$eye_sim_diff, main="YA corrected eye movement similarity (rep 1 vs rep 2)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram illustrates the corrected similarity scores, which account for the non-specific patterns by subtracting the permutated similarity scores from the raw similarity scores.



# Filter for rep 1 vs rep 3 -YA -------------------------------------------

set.seed(1234)

YA_rep1_dens <- YA_eyedens2_rep %>% filter(rep_num == "1") #Filters YA_eyedens to create density data for repetition numbers 1 and 2 (YA_rep1_dens and YA_rep2_dens).

YA_rep3_dens <- YA_eyedens2_rep %>% filter(rep_num == "3")

# Identify mismatches between the two datasets.
common_subj_images <- intersect(YA_rep1_dens$subj_image, YA_rep3_dens$subj_image)
setdiff(YA_rep1_dens$subj_image, YA_rep3_dens$subj_image) # Unique to rep1
setdiff(YA_rep3_dens$subj_image, YA_rep1_dens$subj_image) # Unique to rep2

YA_rep1_dens <- YA_rep1_dens[YA_rep1_dens$subj_image %in% common_subj_images, ] #remove mis-match
YA_rep3_dens <- YA_rep3_dens[YA_rep3_dens$subj_image %in% common_subj_images, ]

YA_simres_rep1_3 <- template_similarity(YA_rep1_dens, YA_rep3_dens, match_on="subj_image", permute_on = "SID", method="fisherz", permutations=50)

# Using aggregate to compute the mean of 'eye_sim', 'perm_sim', and 'eye_sim_diff' for each 'SID' and 'image'
mean_eye_sim1_3_image <- aggregate(eye_sim ~ SID + image, data = YA_simres_rep1_3, FUN = mean)
mean_perm_sim1_3_image <- aggregate(perm_sim ~ SID + image, data = YA_simres_rep1_3, FUN = mean)
mean_eye_sim_diff1_3_image <- aggregate(eye_sim_diff ~ SID + image, data = YA_simres_rep1_3, FUN = mean)

#aggregate across image for t-test
mean_eye_sim1_3 <- aggregate(eye_sim ~ SID , data = YA_simres_rep1_3, FUN = mean)
mean_perm_sim1_3 <- aggregate(perm_sim ~ SID, data = YA_simres_rep1_3, FUN = mean)
mean_eye_sim_diff1_3 <- aggregate(eye_sim_diff ~ SID, data = YA_simres_rep1_3, FUN = mean)

YA_t.test1_rep1v3<- t.test(mean_eye_sim1_3$eye_sim) #sig
YA_t.test2_rep1v3 <- t.test(mean_perm_sim1_3$perm_sim) #sig
YA_t.test3_rep1v3 <- t.test(mean_eye_sim_diff1_3$eye_sim_diff) #sig

#par(mfrow=c(1,3))
quartz(width = 7, height = 7)  # For macOS
par(mfrow=c(1,3))

YA_p1_rep1v3 <- hist(mean_eye_sim1_3$eye_sim, main="YA raw eye movement similarity (rep 1 vs rep 3)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram shows the distribution of similarity scores for the matched pairs, where each pair represents the same image during rep 1 and rep 4 phases.

YA_p2_rep1v3 <- hist(mean_perm_sim1_3$perm_sim, main="YA image-permuted eye movement similarity (rep 1 vs rep 3)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

YA_p3_rep1v3 <- hist(mean_eye_sim_diff1_3$eye_sim_diff, main="YA corrected eye movement similarity (rep 1 vs rep 3)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram illustrates the corrected similarity scores, which account for the non-specific patterns by subtracting the permutated similarity scores from the raw similarity scores.


# Filter for rep 1 vs rep 4 - YA ------------------------------------------

set.seed(1234)

YA_rep1_dens <- YA_eyedens2_rep %>% filter(rep_num == "1") #Filters YA_eyedens to create density data for repetition numbers 1 and 4 (YA_rep1_dens and YA_rep4_dens).

YA_rep4_dens <- YA_eyedens2_rep %>% filter(rep_num == "4")

# Identify mismatches between the two datasets.

common_subj_images <- intersect(YA_rep1_dens$subj_image, YA_rep4_dens$subj_image)
setdiff(YA_rep1_dens$subj_image, YA_rep4_dens$subj_image) # Unique to rep1
setdiff(YA_rep4_dens$subj_image, YA_rep1_dens$subj_image) # Unique to rep4

YA_rep1_dens <- YA_rep1_dens[YA_rep1_dens$subj_image %in% common_subj_images, ] #remove mis-matches
YA_rep4_dens <- YA_rep4_dens[YA_rep4_dens$subj_image %in% common_subj_images, ]


YA_simres_rep1_4 <- template_similarity(YA_rep1_dens, YA_rep4_dens, match_on="subj_image", permute_on ="SID", method="fisherz", permutations=50)


# Using aggregate to compute the mean of 'eye_sim', 'perm_sim', and 'eye_sim_diff' for each 'SID'and 'image'
mean_eye_sim_1_4_image <- aggregate(eye_sim ~ SID + image, data = YA_simres_rep1_4, FUN = mean)
mean_perm_sim_1_4_image <- aggregate(perm_sim ~ SID + image, data = YA_simres_rep1_4, FUN = mean)
mean_eye_sim_diff1_4_image <- aggregate(eye_sim_diff ~ SID + image, data = YA_simres_rep1_4, FUN = mean)

# Using aggregate to compute the mean of 'eye_sim', 'perm_sim', and 'eye_sim_diff' for each 'SID'
mean_eye_sim_1_4 <- aggregate(eye_sim ~ SID , data = YA_simres_rep1_4, FUN = mean)
mean_perm_sim_1_4 <- aggregate(perm_sim ~ SID , data = YA_simres_rep1_4, FUN = mean)
mean_eye_sim_diff1_4 <- aggregate(eye_sim_diff ~ SID , data = YA_simres_rep1_4, FUN = mean)


YA_t.test1_rep1v4 <- t.test(mean_eye_sim_1_4$eye_sim)  #significantly different than zero
YA_t.test2_rep1v4 <- t.test(mean_perm_sim_1_4$perm_sim)  #significantly different than zero
YA_t.test3_rep1v4 <- t.test(mean_eye_sim_diff1_4$eye_sim_diff) #significantly different than zero



#par(mfrow=c(1,3))
quartz(width = 7, height = 7)  # For macOS
par(mfrow=c(1,3))

YA_p1_rep1v4 <- hist(mean_eye_sim_1_4$eye_sim, main="YA raw eye movement similarity (rep 1 vs rep 4)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram shows the distribution of similarity scores for the matched pairs, where each pair represents the same image during rep 1 and rep 4 phases.

YA_p2_rep1v4 <- hist(mean_perm_sim_1_4$perm_sim, main="YA image-permuted eye movement similarity (rep 1 vs rep 4)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

YA_p3_rep1v4 <- hist(mean_eye_sim_diff1_4$eye_sim_diff, main="YA corrected eye movement similarity (rep 1 vs rep 4)", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram illustrates the corrected similarity scores, which account for the non-specific patterns by subtracting the permutated similarity scores from the raw similarity scores.

#combine all repetition contrasts (controled repetitive similarity scores)
# Add a column to indicate the repetition
mean_eye_sim_diff1_2$rep <- "corrected_rep1_2"
mean_eye_sim_diff1_3$rep <- "corrected_rep1_3"
mean_eye_sim_diff1_4$rep <- "corrected_rep1_4"

# Combine all dataframes into one
YA_corrected_rep_combined_df <- rbind(mean_eye_sim_diff1_2, mean_eye_sim_diff1_3, mean_eye_sim_diff1_4)
count(YA_corrected_rep_combined_df$SID)



## Idiosyncratic - YA ------------------------------------------------------

# Idiosyncratic - rep 1 to 1 - YA ----------------------------------------------
set.seed(1234)

#First you need to create a "subj_image" factor.
#This has a level for each combination of subject and image.

YA_eyedens2_idio <- YA_eyedens2 %>% mutate(subj_image = factor(paste(SID, image, sep = "_")))

#Then, for each repnum (using 1 as an example):

YA_rep1_dens <- YA_eyedens2_idio %>% filter(rep_num == 1)

YA_idiosim1 <- template_similarity(YA_rep1_dens, YA_rep1_dens, match_on="subj_image", permute_on="SID", method="fisherz", permutations=50)

#Now, the values you want are in "perm_sim", e.g.

YA_idiosim1$perm_sim

#This is because we are matching on the identical image ("subj_image")
#since we provide YA_rep1_dens as both reference and source.
#Since we "permute_on" SID, then these permutation comparisons are actually what you want.
#They are the average similarity between each image and all other images of the same repnum except itself.

YA_mean_perm_sim_idio1 <- aggregate(perm_sim ~ SID + image, data = YA_idiosim1, FUN = mean)
YA_mean_perm_sim_idio1_noimage <- aggregate(perm_sim ~ SID, data = YA_idiosim1, FUN = mean)

#t-test
YA_idiosim1_ttest <- t.test(YA_mean_perm_sim_idio1_noimage$perm_sim) #sig
YA_idiosim1_ttest

#histogram
YA_idiosim1_p1<- hist(YA_mean_perm_sim_idio1_noimage$perm_sim, main="YA idiosyncratic similarity - rep 1-1", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

YA_idiosim1_p1

#across subject control, remove permute_on = SID

#options(future.rng.onMisuse = "ignore") #if i dont put this there is an error message when running the next line (template_similarity). Is it enough to just rempve permute_on = SID?

YA_idiosim1_across <- template_similarity(YA_rep1_dens, YA_rep1_dens, match_on="subj_image", method="fisherz", permutations=50)
YA_idiosim1_across$perm_sim

YA_mean_perm_sim_idio1_across <- aggregate(perm_sim ~ SID + image, data = YA_idiosim1_across, FUN = mean)
YA_mean_perm_sim_idio1_across_noimage <- aggregate(perm_sim ~ SID, data = YA_idiosim1_across, FUN = mean)
YA_idiosim1_ttest_across <- t.test(YA_mean_perm_sim_idio1_across_noimage$perm_sim) #sig
YA_idiosim1_ttest_across

YA_idiosim1_p1_across <- hist(YA_mean_perm_sim_idio1_across_noimage$perm_sim, main="YA idiosyncratic similarity (across control)- rep 1", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

#subtract across perm_sim from within perm_sim control to obtain final perm_sim value

#rename perm_sim columns to merge dfs
YA_mean_perm_sim_idio1_across <- YA_mean_perm_sim_idio1_across %>% dplyr::rename(perm_sim_across=perm_sim)

YA_mean_perm_sim_idio1 <- YA_mean_perm_sim_idio1 %>% dplyr::rename(perm_sim_within=perm_sim)

#merge perm_sim within and across
YA_mean_perm_sim_idio1_across_within <- merge(YA_mean_perm_sim_idio1_across,YA_mean_perm_sim_idio1)

#subtract perm_sim across from perm_sim within
YA_mean_perm_sim_idio1_across_within$perm_sim_final <- YA_mean_perm_sim_idio1_across_within$perm_sim_within - YA_mean_perm_sim_idio1_across_within$perm_sim_across

#collapse across image
YA_mean_perm_sim_idio1_final_collapsed <- aggregate(perm_sim_final ~ SID, data = YA_mean_perm_sim_idio1_across_within, FUN = mean)

# Idiosyncratic - rep 2 to 2 - YA ----------------------------------------------


#Then, for each repnum (using 2 as an example): within control

YA_rep2_dens <- YA_eyedens2_idio %>% filter(rep_num == 2)
YA_idiosim2 <- template_similarity(YA_rep2_dens, YA_rep2_dens, match_on="subj_image", permute_on="SID", method="fisherz", permutations=50)

#Now, the values you want are in "perm_sim", e.g.

YA_idiosim2$perm_sim

YA_mean_perm_sim_idio2 <- aggregate(perm_sim ~ SID + image, data = YA_idiosim2, FUN = mean)
YA_mean_perm_sim_idio2_noimage <- aggregate(perm_sim ~ SID, data = YA_idiosim2, FUN = mean)
YA_idiosim2_ttest <- t.test(YA_mean_perm_sim_idio2_noimage$perm_sim) #sig
YA_idiosim2_ttest

YA_idiosim2_p2<- hist(YA_mean_perm_sim_idio2_noimage$perm_sim, main="YA idiosyncratic similarity - rep 2", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

#across subject control, remove permute_on = SID

#options(future.rng.onMisuse = "ignore") #if i dont put this there is an error message when running the next line


YA_idiosim2_across <- template_similarity(YA_rep2_dens, YA_rep2_dens, match_on="subj_image", method="fisherz", permutations=50)
YA_idiosim2_across$perm_sim

YA_mean_perm_sim_idio2_across <- aggregate(perm_sim ~ SID + image, data = YA_idiosim2_across, FUN = mean)
YA_mean_perm_sim_idio2_across_noimage <- aggregate(perm_sim ~ SID, data = YA_idiosim2_across, FUN = mean)
YA_idiosim2_ttest_across <- t.test(YA_mean_perm_sim_idio2_across_noimage$perm_sim) #sig
YA_idiosim2_ttest_across

YA_idiosim2_p1_across <- hist(YA_mean_perm_sim_idio2_across_noimage$perm_sim, main="YA idiosyncratic similarity (across control)- rep 2", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

#subtract across from within control
#rename perm_sim columns to merge dfs
YA_mean_perm_sim_idio2_across <- YA_mean_perm_sim_idio2_across %>% dplyr::rename(perm_sim_across=perm_sim)

YA_mean_perm_sim_idio2 <- YA_mean_perm_sim_idio2 %>% dplyr::rename(perm_sim_within=perm_sim)

YA_mean_perm_sim_idio2_across_within <- merge(YA_mean_perm_sim_idio2_across,YA_mean_perm_sim_idio2)

YA_mean_perm_sim_idio2_across_within$perm_sim_final <- YA_mean_perm_sim_idio2_across_within$perm_sim_within - YA_mean_perm_sim_idio2_across_within$perm_sim_across

YA_mean_perm_sim_idio2_final_collapsed <- aggregate(perm_sim_final ~ SID, data = YA_mean_perm_sim_idio2_across_within, FUN = mean)


# Idiosyncratic - rep 3 - YA ----------------------------------------------


#Then, for each repnum (using 3 as an example):

YA_rep3_dens <- YA_eyedens2_idio %>% filter(rep_num == 3)
YA_idiosim3 <- template_similarity(YA_rep3_dens, YA_rep3_dens, match_on="subj_image", permute_on="SID", method="fisherz", permutations=50)

#Now, the values you want are in "perm_sim", e.g.

YA_idiosim3$perm_sim

YA_mean_perm_sim_idio3 <- aggregate(perm_sim ~ SID + image, data = YA_idiosim3, FUN = mean)
YA_mean_perm_sim_idio3_noimage <- aggregate(perm_sim ~ SID, data = YA_idiosim3, FUN = mean)
YA_idiosim3_ttest <- t.test(YA_mean_perm_sim_idio3_noimage$perm_sim) #sig
YA_idiosim3_ttest

YA_idiosim3_p3<- hist(YA_mean_perm_sim_idio3_noimage$perm_sim, main="YA idiosyncratic similarity - rep 3", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

#across subject control, remove permute_on = SID

#options(future.rng.onMisuse = "ignore") #if i dont put this there is an error message when running the next line


YA_idiosim3_across <- template_similarity(YA_rep3_dens, YA_rep3_dens, match_on="subj_image", method="fisherz", permutations=50)
YA_idiosim3_across$perm_sim

YA_mean_perm_sim_idio3_across <- aggregate(perm_sim ~ SID + image, data = YA_idiosim3_across, FUN = mean)
YA_mean_perm_sim_idio3_across_noimage <- aggregate(perm_sim ~ SID, data = YA_idiosim3_across, FUN = mean)
YA_idiosim3_ttest_across <- t.test(YA_mean_perm_sim_idio3_across_noimage$perm_sim) #sig
YA_idiosim3_ttest_across

YA_idiosim3_p1_across <- hist(YA_idiosim3_across$perm_sim, main="YA idiosyncratic similarity (across control)- rep 3", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

#subtract across from within control
#rename perm_sim columns to merge dfs
YA_mean_perm_sim_idio3_across <- YA_mean_perm_sim_idio3_across %>% dplyr::rename(perm_sim_across=perm_sim)

YA_mean_perm_sim_idio3 <- YA_mean_perm_sim_idio3 %>% dplyr::rename(perm_sim_within=perm_sim)

YA_mean_perm_sim_idio3_across_within <- merge(YA_mean_perm_sim_idio3_across,YA_mean_perm_sim_idio3)

YA_mean_perm_sim_idio3_across_within$perm_sim_final <- YA_mean_perm_sim_idio3_across_within$perm_sim_within - YA_mean_perm_sim_idio3_across_within$perm_sim_across

YA_mean_perm_sim_idio3_final_collapsed <- aggregate(perm_sim_final ~ SID, data = YA_mean_perm_sim_idio3_across_within, FUN = mean)


# Idiosyncratic - rep 4 - YA ----------------------------------------------


#Then, for each repnum (using 4 as an example):

YA_rep4_dens <- YA_eyedens2_idio %>% filter(rep_num == 4)
YA_idiosim4 <- template_similarity(YA_rep4_dens, YA_rep4_dens, match_on="subj_image", permute_on="SID", method="fisherz", permutations=50)

#Now, the values you want are in "perm_sim", e.g.

YA_idiosim4$perm_sim

YA_mean_perm_sim_idio4 <- aggregate(perm_sim ~ SID + image, data = YA_idiosim4, FUN = mean)
YA_mean_perm_sim_idio4_noimage <- aggregate(perm_sim ~ SID, data = YA_idiosim4, FUN = mean)
YA_idiosim4_ttest <- t.test(YA_mean_perm_sim_idio4_noimage$perm_sim) #sig
YA_idiosim4_ttest

YA_idiosim4_p4<- hist(YA_mean_perm_sim_idio4_noimage$perm_sim, main="YA idiosyncratic similarity - rep 4", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.


#across subject control, remove permute_on = SID

#options(future.rng.onMisuse = "ignore") #if i dont put this there is an error message when running the next line


YA_idiosim4_across <- template_similarity(YA_rep4_dens, YA_rep4_dens, match_on="subj_image", method="fisherz", permutations=50)
YA_idiosim4_across$perm_sim

YA_mean_perm_sim_idio4_across <- aggregate(perm_sim ~ SID + image, data = YA_idiosim4_across, FUN = mean)
YA_mean_perm_sim_idio4_across_noimage <- aggregate(perm_sim ~ SID, data = YA_idiosim4_across, FUN = mean)
YA_idiosim4_ttest_across <- t.test(YA_mean_perm_sim_idio4_across_noimage$perm_sim) #sig
YA_idiosim4_ttest_across

YA_idiosim4_p1_across <- hist(YA_mean_perm_sim_idio4_across_noimage$perm_sim, main="YA idiosyncratic similarity (across control)- rep 4", cex.main=0.8, cex.axis=0.7, cex.lab=0.7) #This histogram displays the distribution of similarity scores for the permutated pairs, which are created by randomly shuffling the fixation data.

#subtract across from within control for idiosyncratic
#rename perm_sim columns to merge dfs
YA_mean_perm_sim_idio4_across <- YA_mean_perm_sim_idio4_across %>% dplyr::rename(perm_sim_across=perm_sim)

YA_mean_perm_sim_idio4 <- YA_mean_perm_sim_idio4 %>% dplyr::rename(perm_sim_within=perm_sim)

YA_mean_perm_sim_idio4_across_within <- merge(YA_mean_perm_sim_idio4_across,YA_mean_perm_sim_idio4)

YA_mean_perm_sim_idio4_across_within$perm_sim_final <- YA_mean_perm_sim_idio4_across_within$perm_sim_within - YA_mean_perm_sim_idio4_across_within$perm_sim_across

YA_mean_perm_sim_idio4_final_collapsed <- aggregate(perm_sim_final ~ SID, data = YA_mean_perm_sim_idio4_across_within, FUN = mean)

#combine all repetition contrasts (controlled idiosyncratic similarity scores)

# Add a column to indicate the repetition
YA_mean_perm_sim_idio1_final_collapsed$rep <- "rep1"
YA_mean_perm_sim_idio2_final_collapsed$rep <- "rep2"
YA_mean_perm_sim_idio3_final_collapsed$rep <- "rep3"
YA_mean_perm_sim_idio4_final_collapsed$rep <- "rep4"

# Combine all dataframes into one
YA_combined_df_final <- rbind(YA_mean_perm_sim_idio1_final_collapsed, YA_mean_perm_sim_idio2_final_collapsed, YA_mean_perm_sim_idio3_final_collapsed, YA_mean_perm_sim_idio4_final_collapsed)
