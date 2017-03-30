
heat_kernel <- function(x, sigma=1) {
  exp(-x/(2*sigma^2))
}


blurpoints <- function(coords, sigma=10, dim, resolution=1, weights=NULL) {
  grid <- expand.grid(x=seq(0+resolution/2, dim[1] - resolution/2, by=resolution),
                      y=seq(0+resolution/2, dim[2] - resolution/2, by=resolution)


  lapply(1:nrow(coords), function(i) {
    rflann::RadiusSearch(coords[i,,drop=FALSE], grid, (sigma*3)^2, max_neighbour=25)

  }

}
