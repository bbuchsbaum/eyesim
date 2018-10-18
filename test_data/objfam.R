
library(dplyr)

## load raw fixation data
df1 <- read.table("~/analysis/objfam/SimilarityAnalysis.txt", header=TRUE)
df1$RepID <- factor(ifelse(as.character(df1$RepID) == "Test", "test", as.character(df1$RepID)))
## fill in missing duration information with constant
df1$duration <- rep(1, nrow(df1))

## create an 'imagegroup' variable
df1$imagegroup <- substr(df1$image_filename, 1,5)


## create an 'eye_table' data.frame that groups fixations by Subject/imagegroup/RepID
etab_study <- eye_table("FixX", "FixY", "duration", onset="FixTime", groupvar=c("Subject", "imagegroup", "RepID"),
                  vars=c("Trial", "Accuracy", "SubsequentMorph", "test_confidence", "RepID"),
                  data=subset(df1, RepID != "test" & Subject != "ofy09alb"), clip_bounds=c(212,812, 84, 684))

## create an 'eye_table' data.frame that groups fixations by Subject/imagegroup/RepID
etab_test<- eye_table("FixX", "FixY", "duration", onset="FixTime", groupvar=c("Subject", "imagegroup", "RepID"),
                        vars=c("Trial", "Accuracy", "SubsequentMorph", "test_confidence", "RepID"),
                        data=subset(df1, RepID == "test" & Subject != "ofy09alb"), clip_bounds=c(212,812, 84, 684))

## compute fixation density maps, with bins formed by Subject/imagegroup/RepID
edens_study <- density_by(etab_study, groups=c("Subject", "imagegroup"))

## compute fixation density maps, with bins formed by Subject/imagegroup/RepID
edens_test <- density_by(etab_test, groups=c("Subject", "imagegroup"))

edens_study$Subj_Image <- paste0(edens_study$Subject, "_", edens_study$imagegroup)
edens_test$Subj_Image <- paste0(edens_test$Subject, "_", edens_test$imagegroup)

tsim <- template_similarity(edens_study, edens_test, match_on="Subj_Image", method="pearson", permutations=50)


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


