

estimate_scale <- function(x, y, lower=c(.1,.1), upper=c(10,10), window) {

  if (!is.null(window)) {
    y <- subset(x, onset >= window[1] & onset < window[2])
  }

  if (nrow(x) == 0 || nrow(y) == 0) {
    return(list(par=c(1,1)))
  }
  #if (!is.null(window_y)) {
  #  x <- subset(y, onset >= window_y[1] & onset < window_y[2])
  #}

  cx <- as.matrix(x[,1:2])
  cy <- as.matrix(y[,1:2])
  par <- c(1,1)
  f <- function(p) {
    #browser()
    newy <- cy %*% diag(c(p[1],p[2]))
    pracma::hausdorff_dist(newy,cx)
  }

  ret <- optim(par, f, lower=lower, upper=upper, method="L-BFGS")
}


#' @inheritParams template_similarity
match_scale <- function(ref_tab, source_tab, match_on,
                        refvar="fixgroup",sourcevar="fixgroup",
                        window,...) {
  if (!is.null(window) ) {
    assertthat::assert_that(window[2] > window[1])
  }

  matchind <- match(source_tab[[match_on]], ref_tab[[match_on]])

  source_tab <- source_tab %>% ungroup() %>% mutate(matchind=matchind)
  if (any(is.na(matchind))) {
    warning("did not find matching template map for all source maps. Removing non-matching elements.")
    source_tab <- source_tab %>% filter(!is.na(matchind))
    matchind <- matchind[!is.na(matchind)]
  }

  ret <- source_tab %>% purrr::pmap(function(...) {
    . <- list(...)
    d1 <- ref_tab[[refvar]][[.$matchind]]
    d2 <- .[[sourcevar]]
    res=estimate_scale(d1, d2, window=window)
    tibble(scale_x=res$par[1], scale_y=res$par[2])
  }) %>% bind_rows()

  source_tab %>% mutate(scale_x=ret$scale_x, scale_y=ret$scale_y)
}


#' @xport
rescale.fixation_group <- function(x, sx,sy) {
  x %>% mutate(x=x*sx, y=y*sy)
}




