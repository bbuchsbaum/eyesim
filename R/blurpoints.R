



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
#' @export
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



