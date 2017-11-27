pcstudy <- as_tibble(read.csv("~/Dropbox/Jordana_experiments/Jordana_saliency_study/fix_report_study_input.csv"))

study_tab <- eye_table("FixX", "FixY", duration="FixDuration", onset="FixStartTime",
                  groupvar=c("Image"), data=df1,
                  clip_bounds=c(112, (112+800), 684, 84),
                  vars=c("ImageVersion", "ImageRepetition",
                         "ImageSet", "Block", "Image"))


study_dens <- density_by(study_tab, groups="Image", xbounds=c(0,800), ybounds=c(0,600), outdim=c(8,6), duration_weighted=TRUE, sigma=80)
#study_dens_100 <- density_by(study_tab, groups="Image", xbounds=c(0,800), ybounds=c(0,600), outdim=c(80,60), duration_weighted=TRUE, sigma=80)
study_dens_avg <- Reduce("+", lapply(study_dens$density, function(x) x$z))/length(study_dens)
study_dens_avg <- study_dens_avg/sum(study_dens_avg)
sigma <- .1
weights <- exp(-study_dens_avg^2/(2 * sigma^2))

saliency <- study_dens %>% rowwise() %>% do({
  zcubed <- .$density$z^(1/3)
  zcubed <- zcubed/sum(zcubed)

  zsqw <-  .$density$z^(1/2) * weights
  zsqw <- zsqw / sum(zsqw)

  zrank <- rank(.$density$z * weights)
  zrank <- zrank/sum(zrank)

  gg <- expand.grid(x=1:8, y=1:6)
  tibble(Image=.$Image, x=gg$x, y=gg$y, zcuberoot = as.vector(zcubed), zrank=as.vector(zrank), zwsqroot=as.vector(zsqw))
})

write.table(saliency, "~/Dropbox/Jordana_experiments/Jordana_saliency_study/saliency_grid.txt", row.names=FALSE)
