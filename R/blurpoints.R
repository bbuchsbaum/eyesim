



#' blurpoints
#'
#' @param coords
#' @param sigma
#' @param xbounds
#' @param ybounds
#' @param resolution
#' @param weights
#' @param normalize
#' @importFrom Matrix sparseMatrix
#' @importFrom rflann RadiusSearch
blurpoints <- function(coords, sigma=50, xbounds, ybounds, resolution=1, weights=rep(1, nrow(coords)), normalize=TRUE) {

  keep <- coords[,1] >= xbounds[1] & coords[,1] <= xbounds[2] & coords[,2] >= ybounds[1] & coords[,2] <=ybounds[2]

  coords <- coords[keep,]
  weights <- weights[keep]

  xc <- seq(xbounds[1]+resolution/2, xbounds[2] - resolution/2, by=resolution)
  yc <- seq(ybounds[1]+resolution/2, ybounds[2] - resolution/2, by=resolution)

  grid <- expand.grid(x=xc,
                      y=yc)

  xind <- seq_along(xc)
  yind <- seq_along(yc)

  igrid <- expand.grid(x=xind,
                       y=yind)

  ret <- lapply(1:nrow(coords), function(i) {
    ret <- rflann::RadiusSearch(coords[i,,drop=FALSE], grid, (sigma*3)^2, max_neighbour=100000)
    ind <- ret$indices[[1]]

    d <- sqrt(ret$distances[[1]])
    vals <- dnorm(d, sd=sigma) * weights[i]
    out <- sparseMatrix(x=vals, i=igrid[ind,1], j=igrid[ind,2], dims=c(length(xind), length(yind)))
  })

  out <- Reduce("+", ret)

  if (normalize) {
    out <- out/sum(out)
  }

  list(im=out, x=xc, y=yc)
}



pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_study_input.csv"))
pcdelay <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delay_input.csv"))

pcdelay_nofirst <- pcdelay %>% group_by(Subject, ImageNumber) %>% mutate(rnum=1:nrow(.)) %>% filter(rnum != 1)
delay_nof_tab <- create_eye_table(x="FixX", y="FixY", duration="FixDuration", onset="onset",
                                  vars=c("Saliency", "Duration", "Accuracy"),
                                  groupvar=c("Subject", "ImageNumber", "ImageSet"), data=pcdelay_nofirst)
delay_nof_density <- density_by(delay_nof_tab, groups=c("Subject", "ImageNumber", "ImageSet"), sigma=100)

delay_nof_density$ImageVersion <- paste0(delay_nof_density$ImageNumber, "_", delay_nof_density$ImageSet)
delay_nof_sim <- template_similarity(study_density, delay_nof_density, "ImageVersion")
delay_nof_sim$saliency <- delay_nof_tab$Saliency
delay_nof_sim$duration <- delay_nof_tab$Duration
delay_nof_sim$accuracy <- delay_nof_tab$Accuracy
sim_means <- delay_nof_sim %>% group_by(saliency, duration,accuracy) %>% dplyr::summarize(eye_sim=mean(eye_sim))


## 1
pcdelay_test <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_delaytest_input.csv"))

## 2
delay_test_tab <- create_eye_table(x="FixX", y="FixY", duration="FixDuration", onset="onset", groupvar=c("Subject", "ImageNumber", "ImageSet"), data=pcdelay)



