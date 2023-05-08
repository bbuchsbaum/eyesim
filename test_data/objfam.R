
library(dplyr)
library(lme4)

## load raw fixation data
df1 <- read.table("test_data/SimilarityAnalysis.txt", header=TRUE, stringsAsFactors=FALSE)
df1$RepID <- factor(ifelse(as.character(df1$RepID) == "Test", "test", as.character(df1$RepID)))
## fill in missing duration information with constant
#df1$duration <- rep(1, nrow(df1))

## create an 'imagegroup' variable
df1$imagegroup <- substr(df1$image_filename, 1,5)


## create an 'eye_table' data.frame that groups fixations by Subject/imagegroup/RepID
etab_study <- eye_table("FixX", "FixY", duration="FixDur", onset="FixTime", groupvar=c("Subject", "imagegroup", "RepID"),
                  vars=c("Trial", "Accuracy", "SubsequentMorph", "test_confidence", "RepID"),
                  data=subset(df1, RepID != "test" & Subject != "ofy09alb"), clip_bounds=c(212,812, 84, 684))

## create an 'eye_table' data.frame that groups fixations by Subject/imagegroup/RepID
etab_test<- eye_table("FixX", "FixY", duration="FixDur", onset="FixTime", groupvar=c("Subject", "imagegroup", "RepID"),
                        vars=c("Trial", "Accuracy", "SubsequentMorph", "test_confidence", "RepID"),
                        data=subset(df1, RepID == "test" & Subject != "ofy09alb"), clip_bounds=c(212,812, 84, 684))

## compute fixation density maps, with bins formed by Subject/imagegroup/RepID
#edtemp <- density_by(etab_study, groups=c("imagegroup"), sigma=100, xbounds=c(0, 600), ybounds=c(0,600))
## compute fixation density maps, with bins formed by Subject/imagegroup/RepID

edens_study <- density_by(etab_study, groups=c("Subject", "imagegroup"), sigma=80, xbounds=c(0, 600), ybounds=c(0,600),
                            keep_vars=c("SubsequentMorph", "Accuracy", "test_confidence"))


edens_test <- density_by(etab_test, groups=c("Subject", "imagegroup"), sigma=80, xbounds=c(0,600), ybounds=c(0,600),
                         keep_vars=c("SubsequentMorph", "Accuracy", "test_confidence"))

edens_test <- edens_test %>% mutate(correct = ifelse(Accuracy == "CorrectRejection" | Accuracy == "Hit", 1, 0))
gacc <- edens_test %>% group_by(Subject) %>% summarize(avgcor=mean(correct),
                                                       dprime=sum(Accuracy == "Hit")/sum(Accuracy=="Hit" | Accuracy == "Miss") -
                                                         sum(Accuracy == "FalseAlarm")/sum(Accuracy=="FalseAlarm" | Accuracy == "CorrectRejection")

                                                       ) %>% arrange(avgcor)


edens_study$Subj_Image <- paste0(edens_study$Subject, "_", edens_study$imagegroup)
edens_test$Subj_Image <- paste0(edens_test$Subject, "_", edens_test$imagegroup)

tsim <- template_similarity(edens_study, edens_test, match_on="Subj_Image", permute_on="Subject",
                            method="fisherz", permutations=80)

lme.1 <- lmer(eye_sim_diff ~ correct + (1 | Subject), data=subset(tsim, SubsequentMorph == "Old"))
lme.2 <- lmer(eye_sim_diff ~ correct*SubsequentMorph + (1 | Subject), data=subset(tsim, SubsequentMorph != "Old"))
lme.3 <- lmer(eye_sim_diff ~ correct + (1 | Subject), data=subset(tsim, SubsequentMorph == "Near"))
lme.4 <- lmer(eye_sim_diff ~ correct + (1 | Subject), data=subset(tsim, SubsequentMorph == "Far"))
lme.5 <- lmer(eye_sim_diff ~ Accuracy + (1 | Subject), data=subset(tsim, SubsequentMorph == "Far"))
lme.6 <- lmer(eye_sim_diff ~ Accuracy + (1 | Subject), data=subset(tsim, SubsequentMorph == "Near"))
lme.7 <- lmer(eye_sim_diff ~ SubsequentMorph*Accuracy + (1 | Subject), data=subset(tsim, SubsequentMorph != "Old" & Accuracy != "Guess" &
                                                                     Accuracy != "NoResponse"))

lme.8 <- lmer(eye_sim_diff ~ test_confidence + (1 | Subject), data=subset(tsim, SubsequentMorph != "Old" &  Accuracy != "NoResponse"))
lme.9 <- lmer(eye_sim_diff ~ test_confidence + (1 | Subject), data=subset(tsim, SubsequentMorph == "Old" &  Accuracy != "NoResponse"))
lme.10 <- lmer(eye_sim_diff ~ test_confidence*SubsequentMorph + (1 | Subject), data=subset(tsim, SubsequentMorph != "Old" &
                                                                                             Accuracy != "NoResponse"))

tsim2 <- tsim %>% filter(Accuracy != "NoResponse") %>% group_by(Subject) %>% mutate(sconf = scale(test_confidence)) %>% ungroup()
lme.11 <- lmer(eye_sim_diff ~ test_confidence*SubsequentMorph + (1 | Subject), data=subset(tsim,
                                                                                             Accuracy != "NoResponse"))
library(brms)

brm.1 <- brm(eye_sim_diff ~ test_confidence*SubsequentMorph + (1 | Subject), data=subset(tsim,
                                                                                         Accuracy != "NoResponse"))

lme.12 <- lmer(eye_sim_diff ~ sconf*SubsequentMorph + (1 | Subject), data=subset(tsim2,
                                                                                           Accuracy != "NoResponse" & Accuracy != "Guess"))

tsim3 <- tsim2 %>% group_by(Subject, SubsequentMorph) %>% summarize(eye_sim=mean(eye_sim), eye_sim_diff=mean(eye_sim_diff))
tsim4 <- inner_join(tsim3, gacc, by="Subject")

lm.1 <- lm(dprime ~ eye_sim_diff*SubsequentMorph, data=tsim4)


baseline_tab <- density_by(etab_study, groups=c("Subject"), sigma=80, xbounds=c(0, 600), ybounds=c(0,600))
treg <- template_regression(edens_study, edens_test, match_on="Subj_Image", baseline_tab, baseline_key="Subject")
lme.13 <- lmer(beta_source~ test_confidence*SubsequentMorph + (1 | Subject), data=subset(treg,
                                                                                 Accuracy != "NoResponse" & Accuracy != "Guess"))

lme.14 <- lmer(test_confidence ~ beta_source*SubsequentMorph + (1 | Subject), data=subset(treg,
                                                                                      Accuracy != "NoResponse" & Accuracy != "Guess"))
treg <- treg %>% mutate(morph = ifelse(treg$SubsequentMorph == "Near", 1, 0))
lme.15 <- glmer(morph ~ factor(test_confidence)*beta_source + (1 | Subject), family="binomial", data=subset(treg, SubsequentMorph != "Old" &
                                                                                                      Accuracy != "NoResponse" & Accuracy != "Guess"))

tmp <- tsim %>% filter(Accuracy != "NoResponse") %>% group_by(test_confidence, SubsequentMorph) %>% summarize(eye_sim=mean(eye_sim),
                                                                                                              eye_sim_diff=mean(eye_sim_diff))








## compute similarities among RepIDs during encoding
res <- edens %>% group_by(Subject, imagegroup) %>% do({
  if (nrow(.) < 4) {
    data.frame(senc=NA, sdev=NA)
  } else {

    s12 <- similarity(.$density[[1]], .$density[[2]])
    s13 <- similarity(.$density[[1]], .$density[[3]])
    s23 <- similarity(.$density[[2]], .$density[[3]])
    senc <- mean(c(s12,s13,s23))

    x <- c(.$fixgroup[[1]]$x, .$fixgroup[[2]]$x,.$fixgroup[[3]]$x)
    y <- c(.$fixgroup[[1]]$y, .$fixgroup[[2]]$y,.$fixgroup[[3]]$y)

    data.frame(senc=senc, sdev=sd(x)*sd(y))
  }
})

## compute fixation density maps, with bins formed by Subject/imagegroup/RepID
edens <- density_by(etab, groups=c("Subject", "imagegroup", "RepID"))




etab2 <- filter(etab, RepID == 1)
res2 <- inner_join(res, etab2, by=c("Subject", "imagegroup"))

res2 <- filter(res2, Accuracy != "NoResponse")
res2$Acc <- ifelse(res2$Accuracy %in% c("CorrectRejection", "Hit"), 1, 0)

library(lme4)

lme.1 <- glmer(Acc ~ test_confidence + scale(senc) + SubsequentMorph + (1 | Subject), data=res2, family="binomial")
lme.2 <- glmer(Acc ~ test_confidence + scale(senc)*SubsequentMorph + (1 | Subject), data=res2, family="binomial")

lm.3 <- glmer(Acc ~ scale(senc) + (1 | Subject), data=subset(res2, SubsequentMorph=="Old"), family="binomial")
lm.4 <- glmer(Acc ~ test_confidence + scale(senc) + (1 | Subject), data=subset(res2, SubsequentMorph=="Near"), family="binomial")
lm.5 <- glmer(Acc ~ test_confidence + scale(senc) + (1 | Subject), data=subset(res2, SubsequentMorph=="Far"), family="binomial")

drop1(lme.2, test="Chisq")


lme.1a <- glmer(Acc ~ test_confidence + scale(sdev) + SubsequentMorph + (1 | Subject), data=res2, family="binomial")
lme.2a <- glmer(Acc ~ scale(sdev)*SubsequentMorph + (1 | Subject), data=res2, family="binomial")


