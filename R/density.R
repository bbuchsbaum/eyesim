



#' @export
density_by <- function(x, groups, sigma=50, xbounds=c(0, 1000), ybounds=c(0, 1000), outdim=c(100,100),
                       duration_weighted=TRUE, angular=FALSE, angle_bins=8, keep_vars=NULL, result_name="density", ...) {

  vars <- c(groups, keep_vars)
  ret <- x %>% group_by_(.dots=groups) %>% do( {
    g <- do.call(rbind, .$fixgroup)
    cbind(.[1,vars],tibble(fixgroup=list(g)))
  }) %>% rowwise() %>% do( {
    d <- eye_density(.$fixgroup, sigma, xbounds=xbounds, ybounds=ybounds, outdim=outdim,
                     duration_weighted=duration_weighted, angular=angular, angle_bins=angle_bins, origin=attr(x, "origin"), ...)
    cbind(as_tibble(.[vars]), tibble( fixgroup=list(.$fixgroup), density=list(d)))
  })

  ret

}



rank_trans <- scales::trans_new(name="rank",
                                transform=function(x) { rank(x) },
                                inverse=function(x) (length(x)+1) - rank(x))

cuberoot_trans <- scales::trans_new(name="rank",
                                    transform=function(x) { x^(1/3) },
                                    inverse=function(x) x^3)
