library(dplyr)
library(tibble)
library(tidyr)
library(imager)
library(mgcv)
library(ggplot2)
library(memoise)
library(missMDA)
library(neuroca)
library(Rfast)

devtools::load_all()

exclude_subs <- c(28, 32, 109)
#exclude_subs <- c()

## 85a is kissing birds
get_image <- function(version) {
  version <- gsub(" ", "_", version)
  paste0("~/Dropbox/Jordana_experiments/Jordana_saliency_study/images/",version, "_1", ".jpeg")
}


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


odim <- as.integer(c(80,60))
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

get_blocks <- function(grouped, dname, full=FALSE) {
  blocks <- grouped %>% purrr::map(function(g) {
    if (!full) {
      lbs <- as.character(g$ImageVersion)
      xfull <- do.call(rbind, g[[dname]] %>% purrr::map(~ as.vector(.$z)))
      row.names(xfull) <- lbs
      xfull
    } else {
      xfull <- matrix(NA, length(all_labels), prod(odim))
      match_ind <- match(lbs, all_labels)
      xfull[match_ind,] <- do.call(rbind, g[[dname]] %>% purrr::map(~ as.vector(.$z)))
      out <- do.call(rbind, g$study_density %>% purrr::map(~ as.vector(.$z)))
      row.names(out) <- g$ImageVersion
      row.names(xfull) <- all_labels
      xfull
    }
  })
  blocks
}

sim_procrustes <- function(snum) {
  dstudy <- study_grouped[[snum]]
  dtest <- delay_grouped[[snum]]
  lbs_study <- as.character(dstudy$ImageVersion)
  lbs_test <- as.character(dtest$ImageVersion)
  mind <- match(lbs_test, lbs_study)
  lbs_study[mind[!is.na(mind)]]

  dstudy_mind <- dstudy[mind[!is.na(mind)],]
  dtest_mind <- dtest[mind[!is.na(mind)],]

  lapply(1:nrow(dstudy_mind), function(i) {
    dx1 <- dstudy_mind[-i,]
    dx2 <- dtest_mind[-i,]
    X1 <- do.call(rbind, purrr::map(dx1$study_density, function(x) as.vector(x$z)/sum(x$z)))
    X2 <- do.call(rbind, purrr::map(dx2$delay_density, function(x) as.vector(x$z)/sum(x$z)))
    Xcat <- block_matrix(list(X1,X2))
    mres <- mfa(Xcat, ncomp=30, normalization="MFA")

    xtest_mem <- as.vector(dtest_mind[i,]$delay_density[[1]]$z)
    xtest_perc <- as.vector(dstudy_mind[i,]$study_density[[1]]$z)

    p1 <- project(mres, xtest_perc, colind=ind_mem)
    p2 <- project(mres, xtest_mem, colind=ind_perc)
  })

study_blocks <- get_blocks(study_grouped, "study_density")
delay_blocks <- get_blocks(delay_grouped, "delay_density")

library(neuroca)
Xcat_perc <- do.call(rbind, study_blocks)
Xcat_mem <- do.call(rbind, delay_blocks)

Y_perc <- factor(row.names(Xcat_perc))
Y_mem <- factor(row.names(Xcat_mem))
S_perc <- rep(1:length(study_blocks), sapply(study_blocks, nrow))
S_mem <- rep(1:length(delay_blocks), sapply(delay_blocks, nrow))

f_perc <- Rfast::anovas(Xcat_perc, as.integer(Y_perc))[,1]
f_mem <- Rfast::anovas(Xcat_mem, as.integer(Y_mem))[,1]


vec_to_image <- function(vec, palette="buda") {
  cplane <- colorplane::IntensityColorPlane(vec, cols=scico(100, palette = palette))
  mat <- array(colorplane::map_colors(cplane)@clrs, odim)
  as.cimg(as.raster(t(mat)))
}

vec_to_mask <- function(vec, expo=2) {
  svec <- scale(vec)
  svec <- abs(svec)^expo
  svec <- (svec - min(svec)) / (max(svec) - min(svec))
  mat <- array(svec, odim)
  as.cimg(as.raster(t(mat)))
}

show_basis <- function(vec, pal="grayC") {
  img <- vec_to_image(vec, palette=pal)
}

overlay <- function(imversion, vec, expo=1.2) {
  img <- vec_to_image(vec)
  mask <- vec_to_mask(vec, expo=expo)
  fname <- get_image(imversion)
  bgimage <- renorm(imager::load.image(fname))
  img <- imager::resize(img, 800,600)
  mask <- imager::resize(mask, 800,600)
  g <- renorm(grayscale(mask), 0,1)

  #img.a <- imlist(img,g) %>% imappend("c")
  #bgimg.a <- imlist(bgimage,imfill(dim=dim(bgimage), val=1)) %>% imappend("c")
  g <- add.colour(g)

  out <- bgimage*(1-g) + rm.alpha(img) * g
  plot(out)
  #grid.draw(as.raster(bgimg.a))
}

bres_perc <- bada(Y_perc, Xcat_perc, S=S_perc, ncomp=33, preproc=center())
bres_mem <- bada(Y_mem,   Xcat_mem,  S=S_mem, ncomp=33, preproc=center())

pids <- match(levels(Y_mem), levels(Y_perc))
Xr_perc <- bres_perc$Xr[pids,]
Xr_mem <- bres_mem$Xr

bres_perc2 <- bada(factor(levels(Y_mem)), Xr_perc, ncomp=33, preproc=center())
bres_mem2 <- bada(factor(levels(Y_mem)),  Xr_mem,  ncomp=33, preproc=center())

#procres <- vegan::procrustes(neuroca::scores(bres_perc2), neuroca::scores(bres_mem2))
#show_basis(loadings(bres_mem2)[,1])
#tmp <- loadings(bres_mem2) %*% t(procres$rotation)

ind_perc <- 1:prod(odim)
ind_mem <- (prod(odim)+1):ncol(Xr_dual)

Xr_dual <- cbind(Xr_perc, Xr_mem)
bres_dual <- bada(factor(levels(Y_mem)), Xr_dual, ncomp=12, preproc=center())

crossval <- function(blocks, S, Y, weighted=TRUE) {
  res <- do.call(rbind, lapply(unique(S), function(s) {
    print(s)
    Xtrain <- do.call(rbind, blocks[-s])
    Xtest <- blocks[[s]]
    Ytrain <- Y[S != s]
    Ytest <- Y[S == s]

    if (!weighted) {
      bres <- bada(Ytrain, Xtrain, S=S[S != s], ncomp=35, preproc=center())
    } else {
      wts <- sqrt(log(Rfast::anovas(Xtrain, as.integer(Ytrain))[,1]))
      #wts <- 1/sqrt(matrixStats::colSds(Xtrain))
      #cscore <- sda::catscore(Xtrain, Ytrain)
      #wts=log(apply(cscore, 1, function(v) max(abs(v))))
      bres <- bada(Ytrain, Xtrain, S=S[S != s], ncomp=35, preproc=center(), A=Matrix::Diagonal(x=wts))
    }
    do.call(rbind, lapply(seq(5,35,by=5), function(nc) {
      p <- predict(bres, Xtest, ncomp=nc, type="prob")
      data.frame(pred=p, actual=Ytest, S=s, ncomp=nc, correct=p == Ytest, weighted=weighted)
    }))
  }))
}

crossval_dual <- function() {
  res <- do.call(rbind, lapply(unique(S_mem), function(s) {
    print(s)
    bres_perc <- bada(Y_perc[S_perc != s], Xcat_perc[S_perc != s,], S=S_perc[S_perc != s], ncomp=33, preproc=center())
    bres_mem <- bada(Y_mem[S_mem != s],   Xcat_mem[S_mem != s,],  S=S_mem[S_mem != s], ncomp=33, preproc=center())
    Xr_perc <- bres_perc$Xr[pids,]
    Xr_mem <- bres_mem$Xr
    Xr_dual_bl <- block_matrix(list(Xr_perc, Xr_mem))
    Xr_dual <- cbind(Xr_perc, Xr_mem)
    bres_dual <- bada(factor(levels(Y_mem)), Xr_dual, ncomp=25, preproc=center())
    mres_dual <- mfa(Xr_dual_bl, ncomp=25,preproc=center(), normalization="MFA")


    ## procrustes
    ## no transformation

    #pc1 <- pca(Xr_perc, ncomp=20)
    #pc2 <- pca(Xr_mem, ncomp=20)
    #procr <- vegan::procrustes(Xr_perc, Xr_mem)



    #pls_dual <- plsr(Xr_mem ~ Xr_perc, ncomp=5)

    Xtest <- Xcat_mem[S_mem == s,]
    Ytest <- Y_mem[S_mem == s]

    ret <- do.call(rbind, lapply(seq(1,25,by=3), function(nc) {
      #p <- predict(bres_dual, Xtest, ncomp=nc, colind=ind_mem)
      #p2 <- predict(bres_dual, Xtest, ncomp=nc, colind=ind_perc)
      cfier <- neuroca:::classifier.projector(mres_dual, labels=row.names(scores(mres_dual)), colind=ind_mem)
      cfier2 <- neuroca:::classifier.projector(mres_dual, labels=row.names(scores(mres_dual)), colind=ind_perc)
      p <- predict(cfier, Xtest, ncomp=nc, metric="euclidean")$class
      p2 <- predict(cfier2, Xtest, ncomp=nc, metric="euclidean")$class

      data.frame(pred_perc=p, prep_mem=p2, actual=Ytest, S=s, ncomp=nc, correct_mem=p == Ytest, correct_perc=p2 == Ytest)
    }))

    message("mem:", mean(ret$correct_mem))
    message("perc:", mean(ret$correct_perc))
    ret
  }))
}

crossval_procrustes <- function(nc=30) {
  res <- do.call(rbind, lapply(unique(S_mem), function(s) {
    print(s)
    bres_perc <- bada(Y_perc[S_perc != s], Xcat_perc[S_perc != s,], S=S_perc[S_perc != s], ncomp=33, preproc=center())
    bres_mem <- bada(Y_mem[S_mem != s],   Xcat_mem[S_mem != s,],  S=S_mem[S_mem != s], ncomp=33, preproc=center())
    Xr_perc <- bres_perc$Xr[pids,]
    Xr_mem <- bres_mem$Xr

    #pc1 <- pca(Xr_perc, ncomp=nc)
    #pc2 <- pca(Xr_mem, ncomp=nc)
    #procr <- vegan::procrustes(scores(pc1), scores(pc2))

    procr <- vegan::procrustes(Xr_perc, Xr_mem)

    Xtest <- Xcat_mem[S_mem == s,]
    Ytest <- Y_mem[S_mem == s]

    #Xrot <- procr$scale * project(pc2, Xtest, comp=1:nc) %*% procr$rotation
    Xrot <- procr$scale * Xtest %*% procr$rotation

    #pcos <- proxy::simil(scores(pc1), Xrot, method="cosine")
    pcos <- proxy::simil(Xr_perc, Xrot, method="cosine")
    pcos2 <- proxy::simil(Xr_perc, Xtest, method="cosine")
    pclass <- apply(pcos, 2, function(v) {
      row.names(Xr_perc)[which.max(v)]
    })

    pclass2 <- apply(pcos2, 2, function(v) {
      row.names(Xr_perc)[which.max(v)]
    })
    df1 <- data.frame(pred=pclass, actual=Ytest, S=s, ncomp=nc, correct= pclass == Ytest, correct2=pclass2 == Ytest)

    message("proc:", mean(df1$correct))
    message("nonproc:", mean(df1$correct2))
   df1
  }))
}

crossval_procrustes <- function(nc=30) {
  res <- do.call(rbind, lapply(unique(S_mem), function(s) {
    print(s)
    Xtest_mem <- Xcat_mem[S_mem == s,]
    Xtest_perc <- Xcat_perc[S_perc == s,]
    labs=sort(intersect(row.names(Xtest_perc), row.names(Xtest_mem)))
    id_perc=match(labs, row.names(Xtest_perc))
    id_mem=match(labs, row.names(Xtest_mem))

    Xtest_perc = Xtest_perc[id_perc,]
    Xtest_mem = Xtest_mem[id_mem,]

    pc1 <- pca(Xtest_perc, ncomp=nc)
    pc2 <- pca(Xtest_mem, ncomp=nc)

    ret <- FactoMineR::coeffRV(scores(pc1), scores(pc2))
    data.frame(s=s, rvstd=ret$rvstd, nc=nc, pval=ret$p.value)
  }))
}





cval_perc_wtd <- crossval(study_blocks, S_perc, Y_perc, wts=TRUE)
cval_perc <- crossval(study_blocks, S_perc, Y_perc, wts=FALSE)
cval_mem_wtd <- crossval(delay_blocks, S_mem, Y_mem, weighted=TRUE)
cval_mem <- crossval(delay_blocks, S_mem, Y_mem, weighted=FALSE)
cval_mem_dual <- crossval_dual()
cval_procrustes <- crossval_procrustes(nc=15)

rsum_perc_wtd = cval_perc_wtd %>% group_by(ncomp) %>% summarize(correct=mean(correct))
rsum_perc = cval_perc %>% group_by(ncomp) %>% summarize(correct=mean(correct))
rsum_mem_wtd = cval_mem_wtd %>% group_by(ncomp) %>% summarize(correct=mean(correct))
rsum_mem = cval_mem %>% group_by(ncomp) %>% summarize(correct=mean(correct))

rsum_mem_dual_perc <- cval_mem_dual %>% group_by(ncomp) %>% summarize(correct=mean(correct_perc))
rsum_mem_dual_mem <- cval_mem_dual %>% group_by(ncomp) %>% summarize(correct=mean(correct_mem))
