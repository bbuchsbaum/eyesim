cart2pol <- function(x, y) {
  rho <- sqrt(x ^ 2 + y ^ 2)
  theta = atan2(y, x)
  cbind(rho, theta)
}

rad2deg <- function(rad) {(rad * 180) / (pi)}

calcangle <- function(x1, x2) {
  rad2deg(acos(sum(x1 * x2) / (norm(x1, "2") * norm(x2, "2"))))
}

#' @export
scanpath.fixation_group <- function(fg) {
  lenx <- diff(fg$x)
  leny <- diff(fg$y)

  polar <- cart2pol(lenx, leny)
  fg <- fg %>% mutate(lenx=c(NA,lenx), leny=c(NA,leny), rho=c(NA,polar[,1]), theta=c(NA,polar[,2]))
  class(fg) <- c("scanpath", class(fg))
  fg
}

