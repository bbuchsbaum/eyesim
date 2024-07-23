
#' Run similarity analysis for fixation data
#'
#' This function compares the similarity between each fixation group in a \code{source_tab}
#' with a matching fixation group in \code{ref_tab} using the specified similarity metric.
#' Optionally, permutation tests can be performed for assessing the significance of similarity values.
#'
#' @param ref_tab A data frame containing the reference fixation groups.
#' @param source_tab A data frame containing the source fixation groups to be compared with the reference fixation groups.
#' @param match_on A column name in both \code{ref_tab} and \code{source_tab} used for matching the fixation groups.
#' @param permutations The number of permutations to perform for assessing the significance of similarity values (default: 0, no permutation tests).
#' @param permute_on An optional column name for limiting the matching indices in permutation tests (default: NULL).
#' @param method The similarity metric to use for comparing fixation groups (e.g., "sinkhorn", "overlap").
#' @param refvar A column name in \code{ref_tab} containing the reference fixation groups.
#' @param sourcevar A column name in \code{source_tab} containing the source fixation groups.
#' @param window An optional numeric vector specifying the temporal window for computing similarity (default: NULL).
#' @param ... Extra arguments passed to the similarity function.
#' @keywords internal
#' @importFrom rlang set_names
run_similarity_analysis <- function(ref_tab, source_tab, match_on, permutations, permute_on=NULL, method,
                                    refvar, sourcevar, window=NULL, ...) {
  args <- list(...)

  # Match indices between source and reference tables
  matchind <- match(source_tab[[match_on]], ref_tab[[match_on]])

  # Add the matched indices to the source table
  source_tab <- source_tab %>% ungroup() %>% mutate(matchind=matchind)

  # Remove rows with no matching template map
  if (any(is.na(matchind))) {
    warning("did not find matching template map for all source maps. Removing non-matching elements.")
    source_tab <- source_tab %>% filter(!is.na(matchind))
    matchind <- matchind[!is.na(matchind)]
  }

  # If permutation tests are performed, split match indices by permute_on
  if (!is.null(permute_on)) {
    assertthat::assert_that(permute_on %in% names(source_tab) && assertthat::assert_that(permute_on %in% names(ref_tab)))
    match_split <- split(matchind, source_tab[[permute_on]])
  }



  # Define a helper function to calculate similarity between d1 and d2 using the specified method and window
  do_sim <- function(d1,d2,method, window=NULL) {
    if (!is.null(window)) {
      p <- purrr::partial(similarity, d1, d2, method=method, window=window)
      sim <- unlist(do.call(p, args))
    } else {
      p <- purrr::partial(similarity, d1, d2, method=method)
      sim <- unlist(do.call(p, args))
    }
  }

  # Calculate similarities and permutation tests (if specified) for each row in the source table
  ret <- source_tab %>% furrr::future_pmap(function(...) {
    . <- list(...)
    d1 <- ref_tab[[refvar]][[.$matchind]]  # Reference data
    d2 <- .[[sourcevar]]                  # Source data

    # Calculate similarity between d1 and d2 using the specified method and window
    sim <- do_sim(d1,d2,method, window)

    # Perform permutation tests if the number of permutations is greater than 0
    if (permutations > 0) {

      # Limit matching indices to the permute variable if specified
      mind <- if (!is.null(permute_on)) {
        match_split[[as.character(.[[permute_on]])]]
      } else {
        matchind
      }

      # Randomly sample a subset of matching indices if the number of permutations is less than the length of mind
      if (permutations < length(mind)) {
        mind <- sample(mind, permutations)
      }

      # Remove the current element from the list of matching indices
      elnum <- match(.$matchind, mind)
      if (!is.na(elnum) && length(elnum) > 0) {
        mind <- mind[-elnum]
      }

      if (length(mind) == 0) {
        warning("no matching candidate indices for permutation test. Skipping.")
        return(tibble(eye_sim=NA, perm_sim=NA, eye_sim_diff=NA))
      }

      # Calculate permuted similarities for each remaining index in mind
      psim <- do.call(rbind, lapply(mind, function(i) {
        d1p <- ref_tab[[refvar]][[i]]
        do_sim(d1p, d2, method, window)
      }))

      #if (length(psim) == 0) {
      #  browser()
      #}

      # Calculate the mean permuted similarity and the difference between the observed and permuted similarities
      if (ncol(psim) > 1) {
        psim <- colMeans(psim)
        cnames <- c(names(sim), paste0("perm_", names(sim)), paste0(names(sim), "_diff"))
        c(unlist(sim), psim, unlist(sim)-psim) %>% rlang::set_names(cnames) %>% bind_rows()
      } else {
        # Create a tibble with the observed similarity, mean permuted similarity, and their difference
        tibble(eye_sim=sim, perm_sim=mean(psim), eye_sim_diff=sim - mean(psim))
      }
    } else {
      # If no permutation tests, return a tibble with the observed similarity
      if (length(sim) == 1) {
        tibble(eye_sim=sim)
      } else {
        sim %>% bind_rows()
      }
    }
  }, .options=furrr::furrr_options(seed = TRUE)) %>% bind_rows()  # Combine the results of each row in the source table into a single tibble

  # Bind the calculated similarity values to the source table and return the result
  source_tab %>% bind_cols(ret)
}



#' Fixation Similarity
#'
#' Compute the similarity between each fixation group in a \code{source_tab} and a matching fixation group in \code{ref_tab}.
#'
#' @param ref_tab The reference table containing the fixation groups to compare.
#' @param source_tab The source table containing the fixation groups to compare.
#' @param match_on The column name in both tables used to match fixation groups.
#' @param permutations The number of permutations to perform for permutation tests (default is 0, no permutations).
#' @param permute_on The column name on which to permute for permutation tests (default is NULL).
#' @param method The similarity metric to use; options are "sinkhorn" and "overlap" (default is "sinkhorn").
#' @param refvar The name of the column containing fixation groups in the reference table (default is "fixgroup").
#' @param sourcevar The name of the column containing fixation groups in the source table (default is "fixgroup").
#' @param window The temporal window over which to compute similarity (default is NULL).
#' @param ... Additional arguments to pass to the similarity metric function.
#'
#' @return A table containing the computed similarities between fixation groups.
#'
#' @examples
#' # Example usage of the fixation_similarity function
#' ref_table <- # reference table data
#' source_table <- # source table data
#' match_column <- # column name to match fixation groups
#' similarity_results <- fixation_similarity(ref_table, source_table, match_column)
#'
#' @export
fixation_similarity <- function(ref_tab, source_tab, match_on, permutations=0, permute_on=NULL,
                                method=c("sinkhorn", "overlap"),
                                refvar="fixgroup", sourcevar="fixgroup", window=NULL, ...) {
  if (!is.null(window) ) {
    assertthat::assert_that(window[2] > window[1])
  }
  message("fixation_similarity: similarity metric is ", method)

  method <- match.arg(method)
  run_similarity_analysis(ref_tab,source_tab, match_on, permutations, permute_on, method, refvar, sourcevar, window, ...)

}

#' @inheritParams template_similarity
#' @export
scanpath_similarity <- function(ref_tab, source_tab, match_on, permutations=0, permute_on=NULL,
                                method=c("multimatch"),
                                refvar="scanpath", sourcevar="scanpath", window=NULL, ...) {

  if (!is.null(window) ) {
    assertthat::assert_that(window[2] > window[1])
  }

  message("scan_similarity: similarity metric is ", method)
  run_similarity_analysis(ref_tab,source_tab, match_on, permutations, permute_on,
                          method, refvar, sourcevar, window, ...)

}



#' template_similarity
#'
#' Compute similarity between each density map in a \code{source_tab} with a matching ("template") density map in \code{ref_tab}.
#'
#' @param ref_tab A data frame or tibble containing reference density maps.
#' @param source_tab A data frame or tibble containing source density maps.
#' @param match_on A character string representing the variable used to match density maps between \code{ref_tab} and \code{source_tab}.
#' @param permute_on A character string representing the variable used to stratify permutations (default is NULL).
#' @param refvar A character string representing the name of the variable containing density maps in the reference table (default is "density").
#' @param sourcevar A character string representing the name of the variable containing density maps in the source table (default is "density").
#' @param method A character string specifying the similarity method to use. Possible values are "spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", and "dcov" (default is "spearman").
#' @param permutations A numeric value specifying the number of permutations for the baseline map (default is 10).
#' @param ... Extra arguments to pass to the `similarity` function.
#'
#'
#' @return A data frame or tibble containing the source table and additional columns with the similarity scores and permutation results.
#' @export
template_similarity <- function(ref_tab, source_tab, match_on, permute_on = NULL, refvar="density", sourcevar="density",
                                method=c("spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov"),
                                permutations=10, ...) {


  method <- match.arg(method)
  message("template_similarity: similarity metric is ", method)
  run_similarity_analysis(ref_tab,source_tab, match_on, permutations, permute_on, method, refvar, sourcevar,...)
}


#' Sample a smooth fixation density map with a set of discrete fixations.
#'
#' This function samples a smooth fixation density map represented by the object \code{x} with a set of discrete fixations provided in \code{fix}.
#'
#' @param x An object of class "density" representing the smooth fixation density map.
#' @param fix A data frame or tibble containing discrete fixations with columns "x", "y", and "onset".
#' @param times A vector of numeric values representing the time points at which the density map should be sampled (default is NULL).
#'
#' @details The function first checks if the \code{times} parameter is NULL. If so, it directly samples the density map using the coordinates of the fixations in the \code{fix} argument. If the \code{times} parameter is provided, the function first calls the \code{sample_fixations} function to generate a new fixation sequence with the specified time points, and then samples the density map using the coordinates of the new fixation sequence. The result is a data frame containing the sampled density values and the corresponding time points.
#'
#' @return A data frame with columns "z" and "time", where "z" contains the sampled density values and "time" contains the corresponding time points.
#' @export
sample_density.density <- function(x, fix, times=NULL) {
  if (is.null(times)) {
    cds <- as.matrix(cbind(fix$x, fix$y))
    ix <- sapply(cds[,1], function(p) which.min(abs(x$x - p)))
    iy <- sapply(cds[,2], function(p) which.min(abs(x$y - p)))
    data.frame(z=x$z[cbind(ix,iy)], time=fix$onset)
  } else {
    fg <- sample_fixations(fix, times)
    cds <- as.matrix(cbind(fg$x, fg$y))
    res <- lapply(1:nrow(cds), function(i) {
      i1 <- which.min(abs(x$x - cds[i,1]))
      i2 <- which.min(abs(x$y - cds[i,2]))
      if (length(i1) > 0) {
        x$z[i1,i2]
      } else {
        NA
      }
    })

    data.frame(z=unlist(res), time=times)
  }
}


#' This function creates a density object from the provided x, y, and z matrices. The density object is a list containing the x, y, and z values with a class attribute set to "density" and "list".
#'
#' @param x A numeric vector representing the x-axis values of the density map.
#' @param y A numeric vector representing the y-axis values of the density map.
#' @param z A matrix representing the density values at each (x, y) coordinate.
#'
#' @details The function first checks if the dimensions of the z matrix are equal to the length of the x and y vectors. If not, it throws an error. Then, it creates a list containing the x, y, and z values and sets the class attribute of the list to "density" and "list".
#'
#' @return A density object which is a list containing the x, y, and z values with a class attribute set to "density" and "list".
#' @export
gen_density <- function(x,y,z) {
  if (!all(dim(z) == c(length(x), length(y)))) {
    stop("length of x and y must equal nrow(z) and ncol(z)")
  }

  out <- list(
    x=x,
    y=y,
    z=z)

  class(out) <- c("density", "list")
  out
}

#' @export
get_density.eye_density <- function(x, ...) {
  x$z
}

#' Convert an eye_density object to a data.frame.
#'
#' This function converts an eye_density object into a data.frame with x, y, and z values.
#'
#' @param x An eye_density object to be converted into a data.frame.
#' @param ... Additional arguments passed to the method (currently not used).
#'
#' @details The function extracts the x and y values from the eye_density object, then creates a data.frame with all possible combinations of x and y using purrr::cross_df(). It then adds a new column 'z' to the data.frame with the density values from the eye_density object.
#'
#' @return A data.frame with columns x, y, and z representing the x-axis, y-axis, and density values, respectively.
#' @export
#' @importFrom purrr cross_df
#' @importFrom dplyr mutate
as.data.frame.eye_density <- function(x, ...) {
  z <- x$z
  kde_df <- x %>%
    .[c("x", "y")] %>%
    purrr::cross_df() %>%
    dplyr::mutate(z = as.vector(z))

  kde_df
}

#' @export
print.eye_density <- function(x,...) {
  cat("fixation density map", "\n")
  cat("xlim: ", range(x$x), "\n")
  cat("ylim: ", range(x$y), "\n")
  cat("z range: ", range(x$z), "\n")
}



#' Compute a density map for a fixation group.
#'
#' This function computes a density map for a given fixation group using kernel density estimation.
#'
#' @param x A fixation_group object.
#' @param sigma The standard deviation of the kernel. Default is 50.
#' @param xbounds The x-axis bounds. Default is the range of x values in the fixation group.
#' @param ybounds The y-axis bounds. Default is the range of y values in the fixation group.
#' @param outdim The output dimensions of the density map. Default is the difference between the xbounds and ybounds.
#' @param normalize Whether to normalize the output map. Default is TRUE.
#' @param duration_weighted Whether to weight the fixations by their duration. Default is FALSE.
#' @param window The temporal window over which to compute the density map. Default is NULL.
#' @param origin The origin of the coordinate system. Default is c(0,0).
#'
#' @details The function computes a density map for a given fixation group using the kde2d function from the MASS package. The density map is computed based on the x and y coordinates of the fixations, with optional weighting by their duration. The resulting density map can be normalized if desired.
#'
#' @return An object of class c("eye_density", "density", "list") containing the computed density map and other relevant information.
#' @export
#' @importFrom MASS kde2d
#' @family eye_density
eye_density.fixation_group <- function(x, sigma=50, xbounds=c(min(x$x), max(x$x)), ybounds=c(min(x$y), max(x$y)),
                                       outdim=c(xbounds[2] - xbounds[1], ybounds[2] - ybounds[1]),
                                       normalize=TRUE, duration_weighted=FALSE, window=NULL,  origin=c(0,0)) {

  if (!is.null(window)) {
    assertthat::assert_that(window[2] > window[1])
    assert_that("onset" %in% colnames(x))
    x <- filter(x, onset >= window[1] & onset < window[2])
    assertthat::assert_that(nrow(x) > 0)
  }


  wts <- if (duration_weighted) {
    x$duration
  } else {
    rep(1, nrow(x))
  }

  out <- if (duration_weighted || !is.null(window)) {
    xrep <- rep_fixations(x, 50)
    xrep <- x
    kde2d(xrep$x, xrep$y, h=sigma, n=outdim, lims=c(xbounds, ybounds))
  } else {
    kde2d(x$x, x$y, n=outdim, h=sigma, lims=c(xbounds, ybounds))
  }

  if (normalize) {
    out$z <- out$z/sum(out$z)
  }

  out$z <- zapsmall(out$z)

  out$fixgroup <- x

  class(out) <- c("eye_density", "density", "list")
  out

  #ret <- blurpoints(cbind(x$x, x$y), sigma=sigma, xbounds=xbounds, ybounds=ybounds, resolution=resolution, weights=wts, normalize=normalize)
}

Ops.eye_density <- function(e1,e2) {
    op = .Generic[[1]]
    switch(op,
           `-` = {
             delta <- e1$z - e2$z
             #delta <- (delta - min(delta))/(max(delta) - min(delta))
             structure(list(x=e1$x, y=e2$y, z=delta, fixgorup=rbind(e1$fixgroup, e2$fixgroup)),
                       class=c("eye_density_delta", "eye_density", "density", "list"))

           },
           `+` = {
             add = (e1$z + e2$z)/2
             structure(list(x=e1$x, y=e2$y, z=add, fixgorup=rbind(e1$fixgroup, e2$fixgroup)),
                       class=c("eye_density_add", "eye_density", "density", "list"))

           },
           `/` = {
             div = log(e1$z/e2$z)
             structure(list(x=e1$x, y=e2$y, z=div, fixgorup=rbind(e1$fixgroup, e2$fixgroup)),
                       class=c("eye_density_div", "eye_density", "density", "list"))

           },
           stop("undefined operation")
    )
}




#' @noRd
to_angle <- function(x, y) {
  r <- sqrt(x^2 + y^2)
  asin(x/r)
}


#' @noRd
sigmoid <- function (x, a = 1, b = 0)  {
  if (length(x) == 0)
    return(c())
  stopifnot(is.numeric(x), is.numeric(a), is.numeric(b))
  a <- a[1]
  b <- b[1]
  1/(1 + exp(-a * (x - b)))
}


#' Compute Similarity Between Scanpaths
#'
#' This function computes the similarity between two scanpaths using a specified method.
#'
#' @param x A data frame containing the first scanpath.
#' @param y A data frame containing the second scanpath.
#' @param method A character string specifying the method to compute the similarity (default is "multimatch").
#' @param window A numeric vector of length 2 specifying the time window to restrict the fixations in the input scanpaths (default is NULL, which considers all fixations).
#' @param screensize A numeric vector of length 2 specifying the dimensions of the screen (e.g., c(1000, 1000)). Required for the "multimatch" method.
#' @param ... Additional arguments passed to the similarity computation method.
#'
#' @return A numeric value representing the similarity between the two input scanpaths.
#'
#' @examples
#' # Example usage of the similarity.scanpath function
#' scanpath1 <- # first scanpath data
#' scanpath2 <- # second scanpath data
#' similarity_value <- similarity.scanpath(scanpath1, scanpath2, method = "multimatch", screensize = c(1000, 1000))
#'
#' @importFrom dplyr filter
#' @export
#' @family similarity
similarity.scanpath <- function(x, y, method=c("multimatch"),
                                      window=NULL,
                                      screensize=NULL,...) {

  if (!inherits(y, "scanpath")) {
    stop("`y` must be of type `scanpath`")
  }

  if (!is.null(window)) {
    #print(paste("window", window))
    x <- filter(x, onset >= window[1] & onset < window[2])
    y <- filter(y, onset >= window[1] & onset < window[2])
  }

  if (nrow(x) == 0) {
    warning("no observations in 'x'")
    return(NA)
  }

  if (nrow(y) == 0) {
    warning("no observations in 'y'")
    return(NA)
  }

  if (is.null(screensize)) {
    stop("method `multi_match` requires a `screensize` argument (e.g. c(1000,1000)")
  }

  multi_match(x,y,screensize)

}


#' @export
#' @family similarity
similarity.fixation_group <- function(x, y, method=c("sinkhorn", "overlap"),
                                      window=NULL,
                                xdenom=1000, ydenom=1000, tdenom=3000,
                                tweight=.8,  lambda=.1, dthresh=40,
                                time_samples=NULL, screensize=NULL,...) {
  method <- match.arg(method)

  if (!inherits(y, "fixation_group")) {
    stop("`y` must be of type `fixation_group`")
  }

  if (!is.null(window)) {
      #print(paste("window", window))
      x <- filter(x, onset >= window[1] & onset < window[2])
      y <- filter(y, onset >= window[1] & onset < window[2])
  }

  if (nrow(x) == 0) {
    warning("no observations in 'x'")
    return(NA)
  }

  if (nrow(y) == 0) {
    warning("no observations in 'y'")
    return(NA)
  }

  if (method == "sinkhorn") {
    xy1 <- cbind(x$x/xdenom, x$y/ydenom)
    xy2 <- cbind(y$x/xdenom, y$y/ydenom)

    xyt1 <- cbind(x$x/xdenom, x$y/ydenom, x$onset/tdenom * tweight)
    xyt2 <- cbind(y$x/xdenom, y$y/ydenom, y$onset/tdenom * tweight)

    d <- proxy::dist(xyt1, xyt2)

    #stw1 <- sigmoid(x$onset, a=a, b=b)
    #stw2 <- sigmoid(y$onset, a=a, b=b)
    xdur <- x$duration/sum(x$duration)
    ydur <- y$duration/sum(y$duration)

    d0 <- T4transport::sinkhornD(d,wx=xdur, wy=ydur, lambda=lambda)$distance
    1/(1+d0)
  } else if (method == "overlap") {
    if (is.null(time_samples)) {
      stop("method `overlap` requires a vector of `time_samples`")
    }
    fixation_overlap(x, y, dthresh=dthresh, time_samples=time_samples)
  }

}




#' @importFrom proxy simil
#' @export
similarity.density <- function(x, y, method=c("pearson", "spearman", "fisherz", "cosine",
                                              "l1", "jaccard", "dcov"), ...) {
  method=match.arg(method)

  if (inherits(y, "density")) {
    y <- y$z
  }

  compute_similarity(x$z, as.vector(y), method)

}

compute_similarity <- function(x,y, method=c("pearson", "spearman", "fisherz", "cosine", "l1", "jaccard", "dcov")) {
  method=match.arg(method)
  if (method=="pearson" || method == "spearman") {
    cor(as.vector(x), as.vector(y), method=method)
  } else if (method == "fisherz") {
    atanh(cor(as.vector(x), as.vector(y)))
  } else if (method == "cosine") {
    proxy::simil(as.vector(x), as.vector(y), method="cosine", by_rows=FALSE)[,]
  } else if (method == "l1") {
    z <- as.vector(x)
    y <- as.vector(y)
    x1 <- z/sum(z)
    x2 <- y/sum(y)
    1-(proxy::dist(x1, x2, method="Manhattan", by_rows=FALSE)[,])
  } else if (method == "jaccard") {
    proxy::simil(as.vector(x), as.vector(y), method="eJaccard", by_rows=FALSE)[,]
  } else if (method=="dcov") {
    z <- as.vector(x)
    y <- as.vector(y)
    x1 <- z/sum(z)
    x2 <- y/sum(y)
    energy::dcor(x1,x2)
  }

}



kde2d_weighted <- function (x, y, h, n = 25, lims = c(range(x), range(y)), w)
{
  nx <- length(x)
  if (length(y) != nx)
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  #browser()
  n <- rep(n, length.out = 2L)
  gx <- seq(lims[1], lims[2], length = n[1])
  gy <- seq(lims[3], lims[4], length = n[2])
  h <- if (missing(h))
    c(bandwidth.nrd(x), bandwidth.nrd(y))
  else rep(h, length.out = 2L)
  if (any(h <= 0))
    stop("bandwidths must be strictly positive")
  if (missing(w))
    w <- numeric(nx) + 1
  h <- h/4
  ax <- outer(gx, x, "-")/h[1]
  ay <- outer(gy, y, "-")/h[2]
  z <- (matrix(rep(w, n), nrow = n, ncol = nx, byrow = TRUE) *
          matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n,
                                                 nx))/(sum(w) * h[1] * h[2])
  return(list(x = gx, y = gy, z = z))
}
