#' Repetitive Similarity Analysis for Density Maps
#'
#' Computes within-condition and between-condition similarity for density maps.
#' For each density map (trial), this function calculates its average similarity
#' to all other maps within the same condition (`repsim`) and its average similarity
#' to all maps from different conditions (`othersim`).
#'
#' @param tab A data frame or tibble containing the density maps and condition identifiers.
#' @param density_var A character string specifying the name of the column containing the density maps (must be of class "density" or compatible). All maps must share one lattice (same bounds and \code{outdim}); otherwise the call is refused. Default is "density".
#' @param condition_var A character string specifying the name of the column identifying the condition for each trial.
#' @param method A character string specifying the similarity method to use, passed to `similarity.density`.
#'        Possible values include "spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov". Default is "spearman".
#' @param pairwise A logical value indicating whether to return the raw pairwise similarity scores within the same condition for each trial. Default is FALSE. If TRUE, a list column named `pairwise_repsim` will be added.
#' @param multiscale_aggregation If densities are multiscale (i.e., `eye_density_multiscale` objects),
#'        this specifies how to aggregate similarities from different scales.
#'        Options include "mean" (default) or "none" (to get a list column of similarity vectors in `pairwise_repsim`). Note: `repsim` and `othersim` always report the mean similarity across scales.
#' @param ... Additional arguments passed to the `similarity.density` function.
#'
#' @return The input tibble `tab` augmented with the following columns:
#'   - `repsim`: The average similarity of the trial's density map to other maps within the same condition.
#'   - `othersim`: The average similarity of the trial's density map to maps from all other conditions.
#'   - `pairwise_repsim` (optional, if `pairwise = TRUE`): A list column containing vectors of similarity scores between the trial and each other trial within the same condition.
#'
#' @importFrom dplyr %>% group_by mutate row_number select bind_cols bind_rows
#' @importFrom purrr map map_dbl map_lgl
#' @importFrom rlang sym !! := env
#' @importFrom stats na.omit density rnorm var
#' @importFrom tibble as_tibble tibble
#' @export
#'
#' @examples
#' \donttest{
#'   # Generate a small synthetic dataset of density maps across two conditions.
#'   # Fixations in condition A cluster at (30, 30), those in B at (70, 70).
#'   set.seed(123)
#'   n_trials   <- 20
#'   conditions <- rep(c("A", "B"), each = n_trials / 2)
#'
#'   make_map <- function(condition) {
#'     centre <- if (condition == "A") 30 else 70
#'     fg <- fixation_group(x = rnorm(15, centre, 10), y = rnorm(15, centre, 10),
#'                          duration = rep(200, 15), onset = seq(0, by = 250, length.out = 15))
#'     # Common bounds and outdim put every map on the same lattice, which
#'     # similarity() requires.
#'     eye_density(fg, sigma = 10, xbounds = c(0, 100), ybounds = c(0, 100),
#'                 outdim = c(40, 40))
#'   }
#'
#'   my_data <- tibble::tibble(
#'     subject         = rep(1:4, length.out = n_trials),
#'     trial_condition = conditions,
#'     density_map     = lapply(conditions, make_map)
#'   )
#'
#'   # Compute within- and between-condition similarity.
#'   result <- repetitive_similarity(
#'     my_data,
#'     density_var   = "density_map",
#'     condition_var = "trial_condition",
#'     method        = "cosine"
#'   )
#'
#'   # Optionally, return the raw pairwise similarities.
#'   result_pairwise <- repetitive_similarity(
#'     my_data,
#'     density_var   = "density_map",
#'     condition_var = "trial_condition",
#'     method        = "cosine",
#'     pairwise      = TRUE
#'   )
#' }
repetitive_similarity <- function(tab,
                                  density_var = "density",
                                  condition_var,
                                  method = c("spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov", "emd"),
                                  pairwise = FALSE,
                                  multiscale_aggregation = "mean",
                                  ...) {

  method <- match.arg(method)
  density_sym <- rlang::sym(density_var)
  condition_sym <- rlang::sym(condition_var)

  if (!density_var %in% names(tab)) {
    stop(paste("Density variable", density_var, "not found in input table."))
  }
  if (!condition_var %in% names(tab)) {
    stop(paste("Condition variable", condition_var, "not found in input table."))
  }

  # Check the type of the first non-null density object
  first_valid_density_idx <- which(!sapply(tab[[density_var]], is.null))[1]
  if (is.na(first_valid_density_idx)) {
      warning("No valid density objects found in column: ", density_var)
  } else {
      first_density <- tab[[density_var]][[first_valid_density_idx]]
      if (!inherits(first_density, "density") && !inherits(first_density, "eye_density_multiscale")) {
         warning(paste("Column", density_var, "may not contain compatible \'density\' or \'eye_density_multiscale\' objects. Ensure it works with similarity(). Class found:", class(first_density)[1]))
      }
  }

  # Ensure we are not in a grouped/rowwise state (e.g., after `density_by()` which returns a rowwise tibble)
  # so that `row_number()` generates a unique sequential index for the whole table.
  tab <- tab %>%
    dplyr::ungroup() %>%   # drop any grouping structure
    dplyr::mutate(.row_id = dplyr::row_number())

  n_rows <- nrow(tab)

  # Pre-compute similarity matrix for all density pairs
  sim_matrix <- matrix(vector("list", n_rows * n_rows), nrow = n_rows, ncol = n_rows)

  if (n_rows > 1) {
    combn_idx <- utils::combn(n_rows, 2)
    # lapply() rather than purrr::map() so a refused pair surfaces as the plain
    # lattice error, not wrapped with a meaningless pair index.
    pair_vals <- lapply(seq_len(ncol(combn_idx)), function(k) {
      i <- combn_idx[1, k]; j <- combn_idx[2, k]
      d1 <- tab[[density_var]][[i]]
      d2 <- tab[[density_var]][[j]]
      if (is.null(d1) || is.null(d2)) {
        NA_real_
      } else {
        tryCatch({
          similarity(d1, d2, method = method,
                     multiscale_aggregation = multiscale_aggregation, ...)
        }, error = function(e) {
          # Maps on different lattices are an input error, not a numerical
          # failure of one pair: refuse the whole call.
          if (inherits(e, "eyesim_lattice_mismatch")) stop(e)
          warning("Error in similarity calculation for row ", i, " vs ", j, ": ", e$message)
          NA_real_
        })
      }
    })
    for (k in seq_len(ncol(combn_idx))) {
      i <- combn_idx[1, k]; j <- combn_idx[2, k]
      sim_matrix[[i, j]] <- pair_vals[[k]]
      sim_matrix[[j, i]] <- pair_vals[[k]]
    }
  }
  for (i in seq_len(n_rows)) sim_matrix[[i, i]] <- NA_real_

  calculate_mean_sim <- function(sim_list) {
    if (length(sim_list) == 0) return(NA_real_)
    comparison_means <- purrr::map_dbl(sim_list, function(sim_val) {
      if (all(is.na(sim_val))) return(NA_real_)
      mean(sim_val, na.rm = TRUE)
    })
    mean(comparison_means, na.rm = TRUE)
  }

  repsim_vec <- numeric(n_rows)
  othersim_vec <- numeric(n_rows)
  pairwise_list <- vector("list", n_rows)

  for (i in seq_len(n_rows)) {
    current_condition <- tab[[condition_var]][[i]]
    same_idx <- which(tab[[condition_var]] == current_condition & seq_len(n_rows) != i)
    other_idx <- which(tab[[condition_var]] != current_condition)

    repsim_vals <- sim_matrix[i, same_idx]
    othersim_vals <- sim_matrix[i, other_idx]

    repsim_vec[i] <- calculate_mean_sim(repsim_vals)
    othersim_vec[i] <- calculate_mean_sim(othersim_vals)

    if (pairwise) {
      processed_repsim <- purrr::map(repsim_vals, function(sim_val) {
        if (is.atomic(sim_val)) {
          stats::na.omit(sim_val)
        } else {
          sim_val
        }
      })
      pairwise_list[[i]] <- processed_repsim
    }
  }

  result_tab <- tibble::tibble(repsim = repsim_vec, othersim = othersim_vec)
  if (pairwise) result_tab$pairwise_repsim <- pairwise_list

  out_tab <- dplyr::bind_cols(tab, result_tab) %>%
    dplyr::select(-.row_id)

  return(out_tab)
}
