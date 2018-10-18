library(dplyr)
library(tibble)
library(tidyr)
library(ppcor)  ### install.packages("ppcor")
library(MASS)
library(neuroca)


## load study data
pcstudy <- as_tibble(read.csv("~/Dropbox/New_pc_behavioural_data/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()


## load test data to get accuracy values (only taking first fixation row of each trial -- see: 'slice(1)')
pctest <- as_tibble(read.csv("~/Dropbox/New_pc_behavioural_data/delay_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% c(28,32, 109))) %>% droplevels()

pcstudy$ImageNumber <- as.integer(as.character(pcstudy$ImageNumber))

## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Image", "Subject", "Block"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))

## create table for each study trial (Subject/Image)
delay_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                       groupvar=c("Image", "Subject"), data=pctest,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Image", "ImageNumber"))


fill_mat <- function(imat, Y, levs) {
  omat <- matrix(NA, length(levs), ncol(imat))
  ind <- match(Y, levs)
  omat[ind,] <- imat
  row.names(omat) <- levs
  omat
}

get_b <- function(sid) {
  dfx <- study_tab %>% filter(Subject %in% sid)
  dfxdens <- density_by(dfx, groups=c("ImageVersion"),
                      xbounds=c(0,800),
                      ybounds=c(0,600),
                      outdim=c(16,12),
                      duration_weighted=TRUE, sigma=50)

  dfx2 <- delay_tab %>% filter(Subject %in% sid)
  dfxdens2 <- density_by(dfx2, groups=c("ImageVersion"),
                       xbounds=c(0,800),
                       ybounds=c(0,600),
                       outdim=c(16,12),
                       duration_weighted=TRUE, sigma=50)



  imat <- do.call(rbind, purrr::map(dfxdens$density, ~ as.vector(.$z)))
  Y1 <- droplevels(dfxdens$ImageVersion)

  imat2 <- do.call(rbind, purrr::map(dfxdens2$density, ~ as.vector(.$z)))
  Y2 <- droplevels(dfxdens2$ImageVersion)

  X1 <- fill_mat(imat, Y1, levels(pcstudy$ImageVersion))
  X2 <- fill_mat(imat2, Y2, levels(pcstudy$ImageVersion))

  list(study=X1, test=X2)
}

get_avgs <- function(Xs) {
  res <- lapply(1:nrow(Xs[[1]]), function(i) {
    iset <- lapply(1:length(Xs), function(j) {
      Xs[[j]][i,]
    })

    isetavg <- colMeans(do.call(rbind, iset), na.rm = TRUE)
  })

  Xavg <- do.call(rbind, res)
  row.names(Xavg) <- row.names(Xs[[1]])
  Xavg
}

get_estimate <- function(Xs, Xavg) {
  Xout <- Xs
  for (i in 1:length(Xs)) {
    for (j in 1:nrow(Xs[[i]])) {
      if (any(is.na(Xs[[i]][j,]))) {
        Xout[[i]][j,] <- Xavg[j,]
      }
    }
  }
  Xout
}

Xl <- lapply(unique(pcstudy$Subject), function(i) {
  print(i)
  get_block(i)
})

Xstudy <- lapply(Xl, "[[", "study")
Xtest <- lapply(Xl, "[[", "test")
Xavgstudy <- get_avgs(Xstudy)
Xstart <- get_estimate(Xstudy, Xavgstudy)
Xstart <- neuroca::to_block_matrix(do.call(cbind, Xstart),
                                       block_lengths=rep(ncol(Xstart[[1]]), length(Xstart)))
#Xlstudy <- neuroca::to_block_matrix(do.call(cbind, Xstudy),
#                                    block_lengths=rep(ncol(Xstudy[[1]]), length(Xstudy)))

mfa.1 <- neuroca:::mfa(Xstart, 15, center=FALSE, normalization="MFA")

avg_loadings <- function(fit, i) {
  lds <- loadings(fit)[,i]
  mat <- matrix(lds, length(Xstudy), 16*12, byrow=TRUE)
  matrix(colMeans(mat),16,12)
}

res <- imputeMFA(as.data.frame(Xlstudy), group=rep(ncol(Xstudy[[1]]),length(Xstudy)), ncp=15, maxiter=25)
gids <- unlist(lapply(1:length(Xstudy), function(i) rep(i, ncol(Xstudy[[i]]))))

iset <- lapply(1:length(Xstudy), function(i) {
  matrix(unlist(res$completeObs[7,gids==i]), 16,12)
})

isetavg <- t(apply(Reduce("+", iset)/length(iset), 1, rev))
image(iset[[1]]^(1/2))
image(iset[[22]]^(1/2))
image(iset[[21]]^(1/2))
Xstudy[[22]][7,]

pres <- neuroca:::mmpca(imat, Y, center=FALSE, knn=10)
refscores <- scores(pres)
sdres <- sda(scores(pres), L=dfx$ImageVersion, diagonal=TRUE)

dfy <- delay_tab %>% filter(Subject %in% c(115))
dfydens <- density_by(dfy, groups=c("ImageVersion"),
                      xbounds=c(0,800),
                      ybounds=c(0,600),
                      outdim=c(40,30),
                      duration_weighted=TRUE, sigma=80)
dmat <- do.call(rbind, purrr::map(dfydens$density, ~ as.vector(.$z)))


## code match trials
m <- pctest %>% group_by(Subject) %>% rowwise() %>% do( {
  sversion <- as.character(filter(pcstudy, Subject == .$Subject & ImageNumber == .$ImageNumber)$ImageSet[1])
  Match <- if (sversion == as.character(.$ImageSet)) "match"  else "mismatch"
  data.frame(Match=Match)
})

pctest$Match <- m$Match

## join encoding sim and test data
pctest2 <- inner_join(pctest, encoding_sim, by=c("Subject", "ImageNumber"))
pctest2$sid <- factor(pctest2$Subject)


## experiment with within-subject lmers... after adding ImageNumber as random effect, sig. effects are gone.
lme.1 <- lmer(Accuracy ~ sim_wavg  + Saliency + Duration + (1 | sid), data=subset(pctest2, Match=="mismatch"))
lme.2 <- lmer(Accuracy ~ sim_within + sim_wavg + Saliency + Duration + (1 | sid), data=subset(pctest2, Match=="mismatch"))
lme.3 <- lmer(Accuracy ~ sim_within + sim_wavg + Saliency + Duration + (1 | sid) + (1 | ImageNumber), data=subset(pctest2, Match=="mismatch"))

lme.4 <- lmer(Accuracy ~ sim_wdiff + Saliency + Duration + (1 | sid) + (1 | ImageNumber), data=subset(pctest2, Match=="mismatch"))


## for between subject analysis, average accuracy over subject
pctestacc <- pctest %>% group_by(Subject) %>% summarize(acc=mean(Accuracy))

## compute subject averages for encoding_sim
pcstudysim <- encoding_sim %>% group_by(Subject) %>% summarize(sim_wavg=mean(sim_wavg), sim_within=mean(sim_within), sim_wdiff=mean(sim_wdiff),
                                                               sd_fix=mean(sd_fix),
                                                               total_fix=mean(total_fix), slope=mean(slope), slope_cor=mean(slope_cor, na.rm=TRUE))
## join accuracy and encoding_sim averages
dfx <- inner_join(pctestacc, pcstudysim, by="Subject")

## linear model
## sim_wavg = the average similarity with the subject's own mean fixation map.
## sim_within = average similarity among repetitions after partialing out the average_map
lm.1 <- lm(acc ~ sim_wavg, data=dfx, subset=Subject < 300 & Subject != 16)
lm.2 <- lm(acc ~ sim_within, data=dfx, subset=Subject < 300 & Subject != 16)
lm.3 <- lm(acc ~ sim_wavg + sim_within, data=dfx, subset=Subject < 300 & Subject != 16)
lm.4 <- lm(acc ~ sim_wdiff, data=dfx, subset=Subject < 300 & Subject != 16)


lm.4 <- lm(acc ~ sd_fix, data=dfx, subset=Subject < 300 & Subject != 16)
lm.5 <- lm(acc ~ total_fix, data=dfx, subset=Subject < 300 & Subject != 16)
lm.6 <- lm(acc ~ slope, data=dfx, subset=Subject < 300 & Subject != 16)


## plot
ggplot(aes(sim_within, acc), data=dfx) + geom_point() + stat_smooth(method=rlm)

ggplot(aes(sim_wavg, acc), data=dfx) + geom_point() + stat_smooth(method=rlm)

