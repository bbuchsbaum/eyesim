#' Repetitive Similarity Analysis for Density Maps
#'
#' Computes within-condition and between-condition similarity for density maps.
#' For each density map (trial), this function calculates its average similarity
#' to all other maps within the same condition (`repsim`) and its average similarity
#' to all maps from different conditions (`othersim`).
#'
#' @param tab A data frame or tibble containing the density maps and condition identifiers.
#' @param density_var A character string specifying the name of the column containing the density maps (must be of class "density" or compatible). Default is "density".
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
#'   # Each "density_map" is created from normally-distributed random samples.
#'   set.seed(123)
#'   n_trials   <- 20
#'   conditions <- rep(c("A", "B"), each = n_trials / 2)
#'
#'   my_data <- tibble::tibble(
#'     subject         = rep(1:4, length.out = n_trials),
#'     trial_condition = conditions,
#'     density_map     = purrr::map(seq_len(n_trials), function(i) {
#'       x <- rnorm(100,
#'                  mean = ifelse(conditions[i] == "A", 0, 2),
#'                  sd   = 1)
#'       stats::density(x)
#'     })
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
                                  method = c("spearman", "pearson", "fisherz", "cosine", "l1", "jaccard", "dcov"),
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

  all_conditions <- unique(tab[[condition_var]])

  # Use map instead of map_dbl initially to handle potential vector returns
  results_list <- purrr::map(1:nrow(tab), function(i) {
    current_row <- tab[i, ]
    current_density <- current_row[[density_var]][[1]]
    current_condition <- current_row[[condition_var]]
    current_row_id <- current_row$.row_id

    # Indices for same condition (excluding self)
    same_cond_indices <- which(tab[[condition_var]] == current_condition & tab$.row_id != current_row_id)
    # Indices for other conditions
    other_cond_indices <- which(tab[[condition_var]] != current_condition)

    # Similarity function now returns a list of similarities (potentially vectors)
    sim_fun_list <- function(other_idx) {
       if (length(other_idx) == 0) return(list()) # Return empty list
       purrr::map(other_idx, function(j) { # Use map here
           other_density <- tab[[density_var]][[j]]
           if (is.null(current_density) || is.null(other_density)) return(NA_real_)
           # Pass multiscale_aggregation here
           tryCatch({
             similarity(current_density, other_density, method = method,
                        multiscale_aggregation = multiscale_aggregation, ...)
           }, error = function(e) {
             warning("Error in similarity calculation for row ", i, " vs ", j, ": ", e$message)
             NA_real_
           })
       })
    }

    repsim_values_list <- sim_fun_list(same_cond_indices)
    othersim_values_list <- sim_fun_list(other_cond_indices)

    # Calculate mean similarity across scales for repsim/othersim columns
    # Regardless of multiscale_aggregation, these columns report the mean.
    calculate_mean_sim <- function(sim_list) {
        if (length(sim_list) == 0) return(NA_real_)
        # Calculate mean for each comparison (element in sim_list)
        comparison_means <- purrr::map_dbl(sim_list, function(sim_val) {
            if (all(is.na(sim_val))) return(NA_real_)
            mean(sim_val, na.rm = TRUE)
        })
        # Calculate mean of those means
        mean(comparison_means, na.rm = TRUE)
    }

    repsim_mean <- calculate_mean_sim(repsim_values_list)
    othersim_mean <- calculate_mean_sim(othersim_values_list)

    # Create result list
    res_list <- list(
      repsim = repsim_mean,
      othersim = othersim_mean
    )

    if (pairwise) {
      # Store the potentially multi-value list; filter out NAs within each vector/scalar
      # This will result in a list column where each element is a list of similarity results
      # (each result being a scalar or vector depending on multiscale_aggregation)
       processed_repsim_list <- purrr::map(repsim_values_list, function(sim_val) {
             if (is.atomic(sim_val)) {
                  stats::na.omit(sim_val)
             } else {
                  sim_val # Keep structure if not atomic (e.g., already processed NAs)
             }
        })
       # Filter out comparisons that resulted in only NAs or were empty after na.omit
       # processed_repsim_list <- Filter(function(x) length(x) > 0, processed_repsim_list)

      # Wrap the processed list in an additional list so that downstream code/tests can
       # reliably access it with `[[1]][[1]]`, obtaining the list of pairwise comparisons.
       # This maintains backward-compatibility with existing expectations.
       res_list$pairwise_repsim <- list(processed_repsim_list) # Store list of similarity results (scalar or vector)
    }

    return(tibble::as_tibble(res_list))

  })

  # Combine results and remove temporary ID
  results_df <- dplyr::bind_rows(results_list)
  out_tab <- dplyr::bind_cols(tab, results_df) %>% # bind_cols handles list columns
             dplyr::select(-.row_id) # Remove temporary row id

  return(out_tab)
}
