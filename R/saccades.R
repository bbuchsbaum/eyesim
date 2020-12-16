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
add_scanpath.data.frame <- function(x, outvar="scanpath", fixvar="fixgroup" ) {
  x %>% mutate(!!outvar := list(scanpath(.data[[fixvar]][[1]])))
}

#' @export
add_scanpath.eye_table <- function(x, outvar="scanpath", fixvar="fixgroup" ) {
  x %>% mutate(!!outvar := list(scanpath(.data[[fixvar]][[1]])))
}



#' @export
scanpath.fixation_group <- function(x) {

  lenx <- diff(x$x)
  leny <- diff(x$y)

  polar <- cart2pol(lenx, leny)
  x <- x %>%
    mutate(lenx=c(lenx,0), leny=c(leny,0), rho=c(polar[,1],0), theta=c(polar[,2], 0))
  class(x) <- c("scanpath", class(x))
  x
}



