library(dplyr)
library(tibble)
library(tidyr)
library(imager)
library(mgcv)
library(ggplot2)
library(memoise)
library(missMDA)

devtools::load_all()

#exclude_subs <- c(28, 32, 109)
exclude_subs <- c()

## load testdelay fixation data
pcdelay <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/delay_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% exclude_subs )) %>% droplevels()

delay_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixOffset",
                      groupvar=c("Subject", "ImageVersion"), data=pcdelay,
                      clip_bounds=c(112, (112+800), 684, 84),
                      vars=c("ImageVersion", "Saliency", "Accuracy",
                             "ImageSet", "Trial", "Duration", "ImageRepetition", "ImageNumber"))

pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/study_fixations.csv")) %>%
  filter(Image != "." & !(Subject %in% exclude_subs)) %>% droplevels() %>% mutate(ImageNumber=as.integer(as.character(ImageNumber)))


odim <- as.integer(c(40,30))
## create table for each study trial (Subject/Image)
study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                       groupvar=c("Subject", "ImageVersion"), data=pcstudy,
                       clip_bounds=c(112, (112+800), 684, 84),
                       vars=c("ImageVersion",
                              "ImageSet", "Block", "Image", "ImageNumber"))


## construct heatmaps for the study phase, averaged within subjects
study_dens <- density_by(study_tab, groups=c("Subject", "ImageVersion"),
                         xbounds=c(0,800), ybounds=c(0,600), outdim=odim,
                         duration_weighted=TRUE, sigma=60, result_name="study_density")


## construct heatmaps for the study phase, averaged within subjects
delay_dens <- density_by(delay_tab, groups=c("Subject", "ImageVersion"),
                         xbounds=c(0,800), ybounds=c(0,600), outdim=odim,
                         duration_weighted=TRUE, sigma=60, result_name="delay_density")

study_grouped <- study_dens %>% group_by(Subject) %>% dplyr::group_split()
delay_grouped <- delay_dens %>% group_by(Subject) %>% dplyr::group_split()

all_labels <- sort(unique(as.character(study_dens$ImageVersion)))

get_blocks <- function(grouped, dname) {
  blocks <- grouped %>% purrr::map(function(g) {

    #xfull <- matrix(NA, length(all_labels), prod(odim))
    lbs <- as.character(g$ImageVersion)
    #match_ind <- match(lbs, all_labels)
    #xfull[match_ind,] <- do.call(rbind, g[[dname]] %>% purrr::map(~ as.vector(.$z)))
    xfull <- do.call(rbind, g[[dname]] %>% purrr::map(~ as.vector(.$z)))
    #out <- do.call(rbind, g$study_density %>% purrr::map(~ as.vector(.$z)))
    #row.names(out) <- g$ImageVersion
    #row.names(xfull) <- all_labels
    row.names(xfull) <- lbs
    xfull
  })
  blocks
}

study_blocks <- get_blocks(study_grouped, "study_density")
delay_blocks <- get_blocks(delay_grouped, "delay_density")

library(neuroca)
Xcat_perc <- do.call(rbind, study_blocks)
Xcat_mem <- do.call(rbind, delay_blocks)

Y_perc <- factor(row.names(Xcat_perc))
Y_mem <- factor(row.names(Xcat_mem))
S_perc <- rep(1:length(study_blocks), sapply(study_blocks, nrow))
S_mem <- rep(1:length(delay_blocks), sapply(delay_blocks, nrow))

bres_perc <- bada(Y_perc, Xcat_perc, S=S_perc, ncomp=33, preproc=center)
bres_mem <- bada(Y_mem,   Xcat_mem,  S=S_mem, ncomp=33, preproc=center)

pids <- match(levels(Y_mem), levels(Y_perc))
Xr_perc <- bres_perc$Xr[pids,]
Xr_mem <- bres_mem$Xr
Xr_dual <- cbind(Xr_perc, Xr_mem)
bres_dual <- bada(factor(levels(Y_mem)), Xr_dual, ncomp=6, preproc=center)
ind_perc <- 1:prod(odim)
ind_mem <- (prod(odim)+1):ncol(Xr_dual)

crossval <- function(blocks, S, Y) {
  res <- do.call(rbind, lapply(unique(S), function(s) {
    print(s)
    Xtrain <- do.call(rbind, blocks[-s])
    Xtest <- blocks[[s]]
    Ytrain <- Y[S != s]
    Ytest <- Y[S == s]
    bres <- bada(Ytrain, Xtrain, S=S[S != s], ncomp=35, preproc=center)
    do.call(rbind, lapply(seq(5,35,by=5), function(nc) {
      p <- predict(bres, Xtest, ncomp=nc)
      data.frame(pred=p, actual=Ytest, S=s, ncomp=nc, correct=p == Ytest)
    }))
  }))
}

crossval_dual <- function() {
  res <- do.call(rbind, lapply(unique(S), function(s) {
    bres_perc <- bada(Y_perc[S_perc != s], Xcat_perc[S_perc != s,], S=S_perc[S_perc != s], ncomp=33, preproc=center)
    bres_mem <- bada(Y_mem[S_mem != s],   Xcat_mem[S_mem != s,],  S=S_mem[S_mem != s], ncomp=33, preproc=center)
    Xr_perc <- bres_perc$Xr[pids,]
    Xr_mem <- bres_mem$Xr
    Xr_dual <- cbind(Xr_perc, Xr_mem)
    bres_dual <- bada(factor(levels(Y_mem)), Xr_dual, ncomp=25, preproc=center)

    Xtest <- Xcat_mem[S_mem == s,]
    Ytest <- Y_mem[S == s]

    do.call(rbind, lapply(seq(1,25,by=3), function(nc) {
      p <- predict(bres_dual, Xtest, ncomp=nc, colind=ind_mem)
      data.frame(pred=p, actual=Ytest, S=s, ncomp=nc, correct=p == Ytest)
    }))
  }))
}


cval_perc <- crossval(study_blocks, S_perc, Y_perc)
cval_mem <- crossval(delay_blocks, S_mem, Y_mem)

rsum_mem = cval_mem %>% group_by(S, ncomp) %>% summarize(correct=mean(correct))

bres <- neuroca::bada(Y, Xcat, S, ncomp=10, preproc=standardize)

feature_sig <- function(sfac, Xcat, Y) {
  pvals <- lapply(1:ncol(Xcat), function(i) {
    x <- Xcat[,i]
    res <- kruskal.test(x ~ Y)
    res$statistic
  })
  matrix(unlist(pvals), odim[1], odim[2])
}

f_perc <- feature_sig(S_perc, Xcat_perc, Y_perc)
f_mem <- feature_sig(S_mem, Xcat_mem, Y_mem)

#labels <- c(study_dens$ImageVersion, delay_dens$ImageVersion)
labels <- unlist(lapply(study_blocks, row.names))
Sl <- neighborweights::label_matrix2(labels, labels)
Sl <- Sl/RSpectra::eigs(Sl,1)$values[1]
#X <- Matrix::bdiag(c(study_blocks, delay_blocks))
X <- Matrix::bdiag(study_blocks)

lds_to_mat <- function(lvec) {
  n <- length(study_blocks) + length(delay_blocks)
  aa <- array(lvec, c(prod(odim), n))
  m1 <- rowMeans(aa[,1:length(study_blocks)])
  #m2 <- rowMeans(aa[,(length(study_blocks)+1):(n)])

  #list(study=matrix(m1, odim[1], odim[2]), delay=matrix(m2, odim[1], odim[2]))
  list(study=matrix(m1, odim[1], odim[2]))
}

gres <- neuroca::pca(X, ncomp=8, preproc=neuroca::pass())


