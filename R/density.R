



#' @export
#' @import rlang
density_by <- function(x, groups, sigma=50, xbounds=c(0, 1000), ybounds=c(0, 1000), outdim=c(100,100),
                       duration_weighted=TRUE, window=NULL, keep_vars=NULL, result_name="density", ...) {

  ## TODO what happens if window produces fixations < 0?

  rname <- rlang::sym(result_name)
  vars <- c(groups, keep_vars)

  if (!missing(groups) && !is.null(groups) ) {
    ret <- x %>% group_by(.dots=groups) %>% do( {
      g <- do.call(rbind, .$fixgroup)
      cbind(.[1,vars],tibble(fixgroup=list(g)))
    }) %>% rowwise() %>% do( {
      d <- eye_density(.$fixgroup, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim,
                       duration_weighted=duration_weighted, window=window, origin=attr(x, "origin"), ...)
      cbind(as_tibble(.[vars]), tibble( fixgroup=list(.$fixgroup), !!rname := list(d)))
    })
  } else {
    #browser()
    fx <- do.call(rbind, x$fixgroup)
    d <- eye_density(fx, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim,
                     duration_weighted=duration_weighted, window=window,origin=attr(x, "origin"), ...)
    ret <- tibble(fixgroup=list(fx), !!rname := list(d))

  }

  ret

}



rank_trans <- scales::trans_new(name="rank",
                                transform=function(x) { rank(x) },
                                inverse=function(x) (length(x)+1) - rank(x))

cuberoot_trans <- scales::trans_new(name="rank",
                                    transform=function(x) { x^(1/3) },
                                    inverse=function(x) x^3)
