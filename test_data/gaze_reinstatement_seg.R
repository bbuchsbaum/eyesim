library(eyesim)
library(patchwork)
library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)
library(rstudioapi)

#### set the working directory ####
set.seed(2023)
path = dirname(getSourceEditorContext()$path)
setwd(path)

df <- read_excel('df_fixation_seg.xlsx')

#### Overall gaze reinstatement ####
# the syntax of eyesim is the same for the following windowed analysis, so I only
# included detailed comments here for simplicity's sake 

# convert df to eye table format for eyesim
eyetab <- eye_table("CURRENT_FIX_X", "CURRENT_FIX_Y", "CURRENT_FIX_DURATION",
                    "CURRENT_FIX_START",
                    groupvar=c("RECORDING_SESSION_LABEL", "AGE", "CONDITION", "IMAGE"),
                    clip_bounds=c(0, 1920, 1080, 0),
                    data=df)

# create the template encoding density map at encoding for each trial, and create
# respective <Subject_Image> labels 
template_subj_enc_dens <- eyetab %>%
  filter(CONDITION == "DESC") %>%
  density_by(groups=c("IMAGE", "RECORDING_SESSION_LABEL"),
             sigma=200, xbounds=c(0, 1920), ybounds=c(0, 1080))
template_subj_enc_dens <- template_subj_enc_dens %>%
  mutate(Subject_Image = interaction(IMAGE, RECORDING_SESSION_LABEL))

# create recall/retention density map to compare to encoding maps, also create
# the same <Subject_Image> labels
ret_dens <- eyetab %>%
  filter(CONDITION == "REC") %>%
  density_by(groups=c("IMAGE", "RECORDING_SESSION_LABEL"),
             sigma=200, xbounds=c(0, 1920), ybounds=c(0, 1080))
ret_dens_subj <- ret_dens %>%
  mutate(Subject_Image = interaction(IMAGE, RECORDING_SESSION_LABEL))

# here, we are moving in order of the recall sessions (<ret_dens_subj>),
# and for each recall session, finding up to 50 description sessions (only 9 in 
# this case because we only have 10 images) from the same participant but not the
# original image to do the permutation. The same for all template_similarity
# call following this
simres_subj <- template_similarity(template_subj_enc_dens, ret_dens_subj,
                                   match_on="Subject_Image", method="fisherz",
                                   permute_on='RECORDING_SESSION_LABEL',permutations=50)
simres_subj <- simres_subj %>%
  mutate(AGE_CODE=substr(simres_subj$RECORDING_SESSION_LABEL,1,1)) %>%
  mutate(AGE=case_when(AGE_CODE=='0'~'YA', AGE_CODE=='1'~'OA'))

# t-test if the gaze reinstatement difference is greater than 0
simres_subj_YA <- simres_subj%>%filter(AGE=='YA')
simres_subj_OA <- simres_subj%>%filter(AGE=='OA')
par(mfrow=c(1,3))
hist(simres_subj_YA$eye_sim, main="YA raw eye movement similarity")
hist(simres_subj_YA$perm_sim, main="YA image-permuted eye movement similarity")
hist(simres_subj_YA$eye_sim_diff, main="YA corrected eye movement similarity")
t.test(simres_subj_YA$eye_sim_diff)
par(mfrow=c(1,3))
hist(simres_subj_OA$eye_sim, main="OA raw eye movement similarity")
hist(simres_subj_OA$perm_sim, main="OA image-permuted eye movement similarity")
hist(simres_subj_OA$eye_sim_diff, main="OA corrected eye movement similarity")
t.test(simres_subj_OA$eye_sim_diff)
t.test(simres_subj_YA$eye_sim_diff, simres_subj_OA$eye_sim_diff)
write_xlsx(simres_subj, "output/gaze_reinstatement.xlsx")
