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
#' @importFrom stats bw.nrd0
run_similarity_analysis <- function(ref_tab, source_tab, match_on, permutations, permute_on=NULL, method,
                                    refvar, sourcevar, window=NULL, multiscale_aggregation = "mean", ...) {
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

  # Calculate similarities and permutation tests (if specified) for each row in the source table
  ret <- source_tab %>% furrr::future_pmap(function(...) {
    . <- list(...)
    d1 <- ref_tab[[refvar]][[.$matchind]]  # Reference data
    d2 <- .[[sourcevar]]                  # Source data

    is_valid_density <- function(obj) {
      inherits(obj, c("density", "eye_density", "eye_density_multiscale")) && !is.null(obj)
    }

    if (!is_valid_density(d1) || !is_valid_density(d2)) {
      warning("Invalid or NULL density encountered in run_similarity_analysis(). Returning NA for this comparison.")
      return(tibble(eye_sim = NA_real_, perm_sim = NA_real_, eye_sim_diff = NA_real_))
    }

    # Define a helper function inline to calculate similarity, passing multiscale_aggregation
    calculate_sim <- function(d1, d2, method, window = NULL, ...) {
        args_list <- list(...)
        if (!is.null(window)) {
            p <- purrr::partial(similarity, d1, d2, method = method, window = window)
        } else {
            p <- purrr::partial(similarity, d1, d2, method = method)
        }
        # Pass multiscale_aggregation explicitly here
        do.call(p, c(args_list, list(multiscale_aggregation = multiscale_aggregation)))
    }

    sim <- calculate_sim(d1, d2, method, window, ... = args)

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
        calculate_sim(d1p, d2, method, window)
      }))

      #if (length(psim) == 0) {
      #  browser()
      #}

      # Calculate the mean permuted similarity and the difference between the observed and permuted similarities
      if (ncol(psim) > 1) {
        # Handle potential vector/list 'sim' when calculating means and differences
        # Calculate mean permuted similarity, preserving names
        perm_sim_mean <- colMeans(psim, na.rm = TRUE)

        # Ensure operations work correctly if sim is a vector/list
        # Need to handle element-wise subtraction potentially
        if (is.list(sim) || (is.vector(sim) && length(sim) > 1)) {
            # Assuming sim and perm_sim_mean have compatible structures/names for subtraction
            sim_vec <- unlist(sim) # Use unlist carefully or match by name if needed
            diff_val <- sim_vec - perm_sim_mean

            # Construct tibble ensuring list columns
            tibble::lst(
                eye_sim = list(sim_vec), # Wrap in list
                perm_sim = list(perm_sim_mean),
                eye_sim_diff = list(diff_val)
            ) %>% tibble::as_tibble()
        } else {
             # Original logic for scalar sim
             diff_val <- sim - perm_sim_mean
             tibble::tibble(eye_sim = sim, perm_sim = perm_sim_mean, eye_sim_diff = diff_val)
        }
      } else {
        # Create a tibble with the observed similarity, mean permuted similarity, and their difference
        tibble(eye_sim=sim, perm_sim=mean(psim), eye_sim_diff=sim - mean(psim))
      }
    } else {
      # If no permutation tests, return a tibble with the observed similarity
      if (length(sim) == 1) {
        tibble(eye_sim=sim)
      } else {
        # Ensure vector 'sim' is wrapped in a list to create a list column
        tibble(eye_sim = list(sim))
      }
    }
  }, .options=furrr::furrr_options(seed = TRUE)) %>% dplyr::bind_rows() # Combine the results of each row in the source table into a single tibble

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
#' @param multiscale_aggregation If the density maps are multiscale (i.e., `eye_density_multiscale` objects), this specifies how to aggregate similarities from different scales. Options: "mean" (default, returns the average similarity across scales), "none" (returns a list or vector of similarities, one per scale, within the result columns). See `similarity.eye_density_multiscale`.
#' @param ... Extra arguments to pass to the `similarity` function.
#'
#'
#' @return A data frame or tibble containing the source table and additional columns with the similarity scores and permutation results.
#' @export
template_similarity <- function(ref_tab, source_tab, match_on, permute_on = NULL, refvar="density", sourcevar="density",
                                method=c("spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov", "emd"),
                                permutations=10, multiscale_aggregation = "mean", ...) {


  method <- match.arg(method)
  message("template_similarity: similarity metric is ", method)
  run_similarity_analysis(ref_tab, source_tab, match_on, permutations, permute_on, method, refvar, sourcevar, multiscale_aggregation = multiscale_aggregation, ...)
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
sample_density.density <- function(x, fix, times = NULL) {
  nearest_index <- function(coord, grid) {
    ind <- round(approx(grid, seq_along(grid), coord, rule = 2)$y)
    ind[ind < 1L] <- 1L
    ind[ind > length(grid)] <- length(grid)
    ind
  }

  if (is.null(times)) {
    cds <- cbind(fix$x, fix$y)
    ix <- nearest_index(cds[, 1], x$x)
    iy <- nearest_index(cds[, 2], x$y)
    data.frame(z = x$z[cbind(ix, iy)], time = fix$onset)
  } else {
    fg <- sample_fixations(fix, times)
    cds <- cbind(fg$x, fg$y)
    ix <- nearest_index(cds[, 1], x$x)
    iy <- nearest_index(cds[, 2], x$y)
    data.frame(z = x$z[cbind(ix, iy)], time = times)
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
  cat("Sigma:", x$sigma, "\n")
  cat("xlim: ", range(x$x), "\n")
  cat("ylim: ", range(x$y), "\n")
  cat("z range: ", range(x$z), "\n")
}



#' Compute a density map for a fixation group.
#'
#' This function computes a density map for a given fixation group using kernel density estimation.
#'
#' @param x A fixation_group object.
#' @param sigma The standard deviation(s) of the kernel. Can be a single numeric value or a numeric vector. If a vector is provided, a multiscale density object (`eye_density_multiscale`) will be created. Default is 50.
#' @param xbounds The x-axis bounds. Default is the range of x values in the fixation group.
#' @param ybounds The y-axis bounds. Default is the range of y values in the fixation group.
#' @param outdim The output dimensions of the density map. Default is c(100, 100).
#' @param normalize Whether to normalize the output map. Default is TRUE.
#' @param duration_weighted Whether to weight the fixations by their duration. Default is FALSE.
#' @param window The temporal window over which to compute the density map. Default is NULL.
#' @param min_fixations Minimum number of fixations required to compute a density map.
#'   If fewer fixations are present after optional filtering, the function returns NULL.
#'   Default is 2.
#' @param origin The origin of the coordinate system. Default is c(0,0).
#'
#' @details The function computes a density map for a given fixation group using kernel density estimation. If `sigma` is a single value, it computes a standard density map. If `sigma` is a vector, it computes a density map for each value in `sigma` and returns them packaged as an `eye_density_multiscale` object, which is a list of individual `eye_density` objects.
#'
#' @return An object of class `eye_density` (inheriting from `density` and `list`) if
#'   `sigma` is a single value, or an object of class `eye_density_multiscale` (a
#'   list of `eye_density` objects) if `sigma` is a vector. Returns `NULL` if
#'   filtering by `window` leaves fewer than `min_fixations` fixations, or if
#'   density computation fails (e.g., due to zero weights).
#' @export
#' @importFrom ks kde
#' @importFrom dplyr filter
#' @importFrom assertthat assert_that
eye_density.fixation_group <- function(x, sigma = 50,
                                       xbounds = c(min(x$x), max(x$x)),
                                       ybounds = c(min(x$y), max(x$y)),
                                       outdim = c(100, 100),
                                       normalize = TRUE, duration_weighted = FALSE,
                                       window = NULL, min_fixations = 2,
                                       origin = c(0, 0),
                                       kde_pkg = "ks",
                                       ...) {

  # Filter by window if specified
  x_filtered <- x # Use a new variable for the potentially filtered data
  if (!is.null(window)) {
    assert_that(length(window) == 2,
                msg = "Window must be a numeric vector of length 2.")
    assert_that(window[2] > window[1],
                msg = "The second element of window must be greater than the first.")
    assert_that("onset" %in% colnames(x_filtered),
                msg = "The data frame must contain an 'onset' column.")

    x_filtered <- dplyr::filter(x_filtered, onset >= window[1] & onset < window[2])
    if (nrow(x_filtered) == 0) {
         warning("No fixations remain after applying the window filter. Returning NULL.")
         return(NULL) # Return NULL if no data left
    }
    # assert_that(nrow(x_filtered) > 0, # Check done above
    #            msg = "No fixations remain after applying the window filter.")
  }

  # Basic check for enough data points for KDE
  if (nrow(x_filtered) < min_fixations) {
      warning(paste0("Not enough fixations (need >= ", min_fixations,
                     ") to compute density. Returning NULL. Provided: ",
                     nrow(x_filtered)))
      return(NULL)
  }


  # Prepare data and weights (original logic, using x_filtered)
  data_matrix <- as.matrix(x_filtered[, c("x", "y")])
  current_weights <- if (duration_weighted) {
    assert_that("duration" %in% colnames(x_filtered),
                msg = "The data frame must contain a 'duration' column.")
    assert_that(is.numeric(x_filtered$duration),
                msg = "'duration' must be numeric.")
    # Allow zero duration? Maybe filter them out? Let's allow for now, ks::kde handles w=0.
    assert_that(all(x_filtered$duration >= 0),
                msg = "All durations must be non-negative.")
    x_filtered$duration
  } else {
    rep(1, nrow(x_filtered))
  }

  # Decide on weights processing based on kde_pkg
  processed_weights <- current_weights
  if (duration_weighted && !(requireNamespace("ks", quietly = TRUE) && kde_pkg == "ks")) {
      # If using custom kde2d_weighted, it might expect specific weight normalization.
      # The example implementation used sum(w) in denominator.
      # Let's ensure weights are positive sum if using this path.
      if (sum(processed_weights) <= 0) {
          warning("Sum of weights is zero or negative, cannot compute weighted density with non-ks method. Returning NULL.")
          return(NULL)
      }
      # Normalization like w/sum(w) * N might be needed depending on kde2d_weighted implementation.
      # Keeping raw weights for now, assuming kde2d_weighted handles it.
  } else if (duration_weighted && requireNamespace("ks", quietly = TRUE) && kde_pkg == "ks") {
      # ks::kde handles raw weights (counts, proportions, etc.)
      if (sum(processed_weights) <= 0) {
          warning("Sum of weights is zero or negative for ks::kde. Result might be zero density. Proceeding.")
      }
  }


  if (length(sigma) > 1) {
    # Multiscale request
    all_eye_densities <- lapply(sigma, function(s_val) {
      .compute_single_eye_density(
        x_data = x_filtered, sigma_val = s_val, xbounds = xbounds, ybounds = ybounds,
        outdim = outdim, normalize = normalize, duration_weighted = duration_weighted,
        weights = processed_weights, data_matrix = data_matrix, kde_pkg = kde_pkg, ...
      )
    })
    # Filter out NULLs if any density computation failed
    all_eye_densities <- Filter(Negate(is.null), all_eye_densities)

    if (length(all_eye_densities) == 0) {
      warning("All multiscale density computations failed.")
      return(NULL)
    }

    # Store original sigma vector and individual sigmas with each element
    # The individual sigmas are already stored by .compute_single_eye_density
    attr(all_eye_densities, "sigmas_vector") <- sapply(all_eye_densities, `[[`, "sigma") # Store successful sigmas
    class(all_eye_densities) <- c("eye_density_multiscale", "list")
    return(all_eye_densities)
  } else {
    # Single scale request (delegates to the helper too for consistency)
    return(.compute_single_eye_density(
      x_data = x_filtered, sigma_val = sigma, xbounds = xbounds, ybounds = ybounds,
      outdim = outdim, normalize = normalize, duration_weighted = duration_weighted,
      weights = processed_weights, data_matrix = data_matrix, kde_pkg = kde_pkg, ...
    ))
  }
}

# Internal function, not exported
.compute_single_eye_density <- function(x_data, sigma_val, xbounds, ybounds, outdim,
                                       normalize, duration_weighted, weights, data_matrix,
                                       kde_pkg = "ks") { # Added kde_pkg, default to ks if available, else MASS

  current_sigma <- sigma_val # Use the specific sigma for this scale
  gridsize <- outdim

  # Determine the weights to use based on duration_weighted flag
  final_weights <- if (duration_weighted) {
    # Use the pre-calculated (potentially duration-based) weights passed in
    weights
  } else {
    # For unweighted case, create uniform weights
    rep(1, nrow(data_matrix))
  }

  # Check if ks package should be used
  use_ks <- requireNamespace("ks", quietly = TRUE) && kde_pkg == "ks"

  if (use_ks) {
    # --- Use ks::kde ---
    # Bandwidth matrix H for ks::kde (isotropic Gaussian)
    H_mat <- diag(rep(current_sigma^2, 2))

    # Pre-scale weights for ks::kde to avoid warning (expects sum(w) == n)
    scaled_final_weights <- final_weights
    sum_w <- sum(final_weights)
    n_obs <- nrow(data_matrix)

    if (sum_w > .Machine$double.eps) {
      scale_factor <- n_obs / sum_w
      scaled_final_weights <- final_weights * scale_factor
    } # else: weights are zero/negative sum, ks::kde handles/warns

    # Check if weights became zero or negative after scaling (unlikely but possible)
    if (sum(scaled_final_weights) <= 0 && n_obs > 0) {
        warning("Sum of weights is zero or negative even after scaling for ks::kde. Result might be zero density. Sigma: ", current_sigma)
        # Proceed, ks::kde might handle this
    }


    kde_result <- tryCatch({
        ks::kde(x = data_matrix,
                H = H_mat,
                gridsize = gridsize,
                xmin = c(xbounds[1], ybounds[1]),
                xmax = c(xbounds[2], ybounds[2]),
                w = scaled_final_weights, # Pass scaled weights
                compute.cont = FALSE)
      }, error = function(e) {
         warning("ks::kde failed for sigma=", current_sigma, ". Error: ", e$message)
         NULL
      })

      if (is.null(kde_result)) return(NULL) # Return NULL if ks::kde failed

      eval_points_val <- kde_result$eval.points
      density_matrix_val <- kde_result$estimate

  } else {
    # --- Fallback to MASS or custom ---
    message("ks package not found or not selected. Using MASS::kde2d (or custom kde2d_weighted if applicable).")

    # Check if weights are non-uniform (relevant if duration_weighted was TRUE)
    is_weighted_fallback <- length(unique(final_weights)) > 1

    if (is_weighted_fallback && exists("kde2d_weighted", mode = "function")) {
        message("Using custom kde2d_weighted(). Ensure it handles weights appropriately.")
         kde_result <- tryCatch({
            kde2d_weighted(data_matrix[,1], data_matrix[,2], h = current_sigma, n = outdim, lims = c(xbounds, ybounds), w = final_weights) # Pass original determined weights
         }, error = function(e) {
            warning("kde2d_weighted failed for sigma=", current_sigma, ". Error: ", e$message)
            NULL
         })
         if (is.null(kde_result)) return(NULL)

    } else {
        if (is_weighted_fallback) {
            warning("Using MASS::kde2d without weights for a weighted request because ks/custom unavailable.")
        }
        # Use standard MASS::kde2d (ignores weights)
        kde_result <- tryCatch({
          MASS::kde2d(data_matrix[,1], data_matrix[,2], n = outdim, h = current_sigma, lims = c(xbounds, ybounds))
        }, error = function(e) {
           warning("MASS::kde2d failed for sigma=", current_sigma, ". Error: ", e$message)
           NULL
        })
         if (is.null(kde_result)) return(NULL)
    }

    eval_points_val <- list(kde_result$x, kde_result$y)
    density_matrix_val <- kde_result$z
  }

  # --- Normalization and Object Creation (common path) ---
  if (normalize) {
    sum_dens <- sum(density_matrix_val)
    if (sum_dens > .Machine$double.eps) {
        density_matrix_val <- density_matrix_val / sum_dens
    } else {
        warning("Sum of density matrix is near zero, cannot normalize. Sigma: ", current_sigma)
    }
  }
  density_matrix_val <- zapsmall(density_matrix_val)

  out_list <- list(x = eval_points_val[[1]],
                   y = eval_points_val[[2]],
                   z = density_matrix_val,
                   fixgroup = x_data,
                   sigma = current_sigma)
  class(out_list) <- c("eye_density", "density", "list")
  return(out_list)
}


#' @export
Ops.eye_density <- function(e1,e2) {
    stopifnot(all(e1$x == e2$x), all(e1$y == e2$y))
    op = .Generic[[1]]
    switch(op,
           `-` = {
             delta <- e1$z - e2$z
             #delta <- (delta - min(delta))/(max(delta) - min(delta))
             structure(list(x=e1$x, y=e1$y, z=delta, fixgroup=rbind(e1$fixgroup, e2$fixgroup)),
                       class=c("eye_density_delta", "eye_density", "density", "list"))

           },
           `+` = {
             add = (e1$z + e2$z)/2
             structure(list(x=e1$x, y=e1$y, z=add, fixgroup=rbind(e1$fixgroup, e2$fixgroup)),
                       class=c("eye_density_add", "eye_density", "density", "list"))

           },
           `/` = {
             div = log(e1$z/e2$z)
             structure(list(x=e1$x, y=e1$y, z=div, fixgroup=rbind(e1$fixgroup, e2$fixgroup)),
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
#' @param x A scanpath object containing the first scanpath.
#' @param y A scanpath object containing the second scanpath.
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

  if (!inherits(x, "scanpath")) {
    stop("`x` must be of type `scanpath`")
  }

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
similarity.density <- function(x, y,
                               method = c("pearson", "spearman", "fisherz",
                                           "cosine", "l1", "jaccard",
                                           "dcov", "emd"),
                               saliency_map = NULL, ...) {
  method <- match.arg(method)
  if (method == "emd") {
    compute_similarity(x, y, method = method, saliency_map = saliency_map)
  } else {
    if (inherits(y, "density")) {
      y <- y$z
    }
    compute_similarity(x$z, as.vector(y), method)
  }
}
compute_similarity <- function(x, y,
                              method = c("pearson", "spearman", "fisherz",
                                          "cosine", "l1", "jaccard", "dcov",
                                          "emd"),
                              saliency_map = NULL) {
  method <- match.arg(method)
  if (method == "emd") {
    if (!all(c("x", "y", "z") %in% names(x)) || !all(c("x", "y", "z") %in% names(y))) {
      stop("method 'emd' requires eye_density objects with x, y, z components")
    }
    coords <- as.matrix(expand.grid(x = x$x, y = x$y))
    wx <- as.vector(x$z)
    wy <- as.vector(y$z)
    if (!is.null(saliency_map)) {
      s_mat <- if (is.list(saliency_map) && !is.null(saliency_map$z)) saliency_map$z else saliency_map
      if (!all(dim(s_mat) == dim(x$z))) stop("saliency_map dimensions must match density maps")
      r1 <- x$z - s_mat
      r2 <- y$z - s_mat
      pos1 <- pmax(r1, 0); neg1 <- pmax(-r1, 0)
      pos2 <- pmax(r2, 0); neg2 <- pmax(-r2, 0)
      emd_pos <- emdw(coords, as.vector(pos1), coords, as.vector(pos2))
      emd_neg <- emdw(coords, as.vector(neg1), coords, as.vector(neg2))
      return(-(emd_pos + emd_neg))
    } else {
      emd_dist <- emdw(coords, wx, coords, wy)
      return(1/(1 + emd_dist))
    }
  }
  vx <- as.vector(x)
  vy <- as.vector(y)

  # Check for sufficient non-NA data and variance for correlation-based methods
  valid_x <- !is.na(vx)
  valid_y <- !is.na(vy)
  valid_common <- valid_x & valid_y

  if (sum(valid_common) < 2) {
    warning("Less than 2 common valid data points for similarity calculation.")
    return(NA_real_)
  }

  vx_common <- vx[valid_common]
  vy_common <- vy[valid_common]

  # Check variance for correlation methods
  var_x <- stats::var(vx_common)
  var_y <- stats::var(vy_common)

  if (method %in% c("pearson", "spearman", "fisherz", "dcov")) {
    is_zero_var_x <- is.na(var_x) || var_x < .Machine$double.eps
    is_zero_var_y <- is.na(var_y) || var_y < .Machine$double.eps

    if (is_zero_var_x || is_zero_var_y) {
      # Check if vectors are identical despite zero variance
      # Use a tolerance for floating point comparisons
      if (is_zero_var_x && is_zero_var_y && all(abs(vx_common - vy_common) < sqrt(.Machine$double.eps))) {
         return(1.0) # Perfect correlation for identical vectors
      } else {
         # If only one has zero variance, or they have zero variance but are different (e.g. one is all 0, other is all 1), correlation is undefined/NA
         warning(paste("Method", method, "requires variance in both inputs, or identical inputs if variance is zero. One or both have near-zero variance and are not identical."))
         return(NA_real_)
      }
    }
  }


  if (method=="pearson" || method == "spearman") {
    stats::cor(vx_common, vy_common, method=method)
  } else if (method == "fisherz") {
    cor_val <- stats::cor(vx_common, vy_common, method="pearson")
    # Ensure cor_val is within (-1, 1) for atanh
    cor_val <- max(min(cor_val, 1 - .Machine$double.eps), -1 + .Machine$double.eps)
    atanh(cor_val)
  } else if (method == "cosine") {
    # proxy::simil handles vectors directly
    proxy::simil(matrix(vx_common, ncol=1), matrix(vy_common, ncol=1), method="cosine", by_rows=FALSE)[1,1]
  } else if (method == "l1") {
    # Assumes non-negative densities for probability distribution interpretation
    if (any(vx_common < 0) || any(vy_common < 0)) {
        warning("L1 distance interpretation assumes non-negative densities.")
    }
    sum_x <- sum(vx_common)
    sum_y <- sum(vy_common)
    if (sum_x <= .Machine$double.eps || sum_y <= .Machine$double.eps) {
        warning("Cannot normalize for L1 distance; sum is too small.")
        return(NA_real_)
    }
    x1 <- vx_common / sum_x
    x2 <- vy_common / sum_y
    # L1 distance is sum(|x1_i - x2_i|). Similarity = 1 - 0.5 * L1_distance (range 0-1)
    l1_dist <- sum(abs(x1 - x2))
    1 - 0.5 * l1_dist # Scale to [0, 1] where 1 is identical
  } else if (method == "jaccard") {
     # proxy::simil eJaccard expects matrices where rows are features, cols are samples
     # Using vectors directly might not be standard interpretation.
     # Let's treat vectors as single samples.
     proxy::simil(matrix(vx_common, ncol=1), matrix(vy_common, ncol=1), method="eJaccard", by_rows=FALSE)[1,1]
  } else if (method=="dcov") {
     # energy::dcor requires matrices or vectors
     # Normalize first? The original code normalized. Let's follow that.
     sum_x <- sum(vx_common)
     sum_y <- sum(vy_common)
     if (sum_x <= .Machine$double.eps || sum_y <= .Machine$double.eps) {
        warning("Cannot normalize for dcov; sum is too small.")
        return(NA_real_)
     }
     x1 <- vx_common / sum_x
     x2 <- vy_common / sum_y
    energy::dcor(x1, x2)
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
    c(bw.nrd0(x), bw.nrd0(y))
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

# Add new S3 print method for multiscale
#' @export
print.eye_density_multiscale <- function(x, ...) {
  sigmas <- attr(x, "sigmas_vector")
  if (is.null(sigmas)) sigmas <- sapply(x, `[[`, "sigma") # Fallback
  cat("Multiscale fixation density map", "\n")
  cat("Number of scales:", length(x), "\n")
  cat("Sigmas:", paste(sigmas, collapse=", "), "\n")
  if (length(x) > 0) {
      cat("--- First scale details ---", "\n")
      print(x[[1]]) # Print details of the first scale
  }
}

#' @export
#' @family similarity
# New S3 method for eye_density_multiscale
similarity.eye_density_multiscale <- function(x, y, method = c("pearson", "spearman", "fisherz", "cosine", "l1", "jaccard", "dcov"),
                                              multiscale_aggregation = "mean", ...) {
  method <- match.arg(method)

  if (!inherits(y, "eye_density_multiscale")) {
    stop("`y` must also be of type `eye_density_multiscale`")
  }

  if (length(x) == 0 || length(y) == 0) {
    warning("One or both multiscale objects are empty.")
    return(NA_real_)
  }

  # Match scales based on sigma value for robustness
  sigmas_x <- sapply(x, `[[`, "sigma")
  sigmas_y <- sapply(y, `[[`, "sigma")
  common_sigmas <- intersect(sigmas_x, sigmas_y)

  if (length(common_sigmas) == 0) {
    warning("No common sigmas found between multiscale objects. Cannot compute similarity.")
    return(NA_real_)
  }

  if (length(common_sigmas) < length(sigmas_x) || length(common_sigmas) < length(sigmas_y)) {
     warning("Not all sigmas matched between multiscale objects. Comparing only common sigmas: ", paste(common_sigmas, collapse=", "))
  }

  # Create pairs of scales matching on sigma
  matched_pairs_x <- x[match(common_sigmas, sigmas_x)]
  matched_pairs_y <- y[match(common_sigmas, sigmas_y)]


  per_scale_similarities <- mapply(function(scale_x, scale_y) {
    # Each scale_x, scale_y is an 'eye_density' object
    tryCatch({
      similarity(scale_x, scale_y, method = method, ...) # Dispatch to similarity.density
    }, error = function(e) {
      warning("Error computing similarity for sigma=", scale_x$sigma, ": ", e$message)
      NA_real_
    })
  }, matched_pairs_x, matched_pairs_y, SIMPLIFY = TRUE) # SIMPLIFY = TRUE to get a vector if possible

  if (all(is.na(per_scale_similarities))) return(NA_real_)

  aggregation_method <- match.arg(tolower(multiscale_aggregation), c("mean", "none"))

  if (aggregation_method == "mean") {
    return(mean(per_scale_similarities, na.rm = TRUE))
  } else if (aggregation_method == "none") {
    # Return named vector if possible
    names(per_scale_similarities) <- paste0("sigma_", common_sigmas)
    return(per_scale_similarities) # Return vector of similarities
  }
  # Fallback should not be reached due to match.arg
}
