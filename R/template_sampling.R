#' Sample density maps at fixation locations over time
#'
#' This function samples template density maps at source fixation locations across
#' specified time points. It is useful for analyzing how fixation patterns relate
#' to a reference density map over the course of a trial (e.g., "reinstatement
#' across time" analysis).
#'
#' @param template_tab A data frame or tibble containing template density maps
#'   (e.g., from an encoding phase).
#' @param source_tab A data frame or tibble containing fixation groups
#'   (e.g., from a retrieval phase).
#' @param match_on A character string specifying the column name used to match
#'   template density maps to source fixation groups.
#' @param times A numeric vector of time points (in ms) at which to sample the
#'   density. Default is \code{seq(0, 3000, by = 50)}.
#' @param time_bins An optional numeric vector specifying bin boundaries for
#'   aggregating samples. For example, \code{c(0, 1000, 2000, 3000)} creates
#'   3 bins: [0-1000), [1000-2000), [2000-3000). Every bin is half-open,
#'   including the last, so a time point equal to the final break (3000 here)
#'   is not assigned to any bin. Default is NULL (no binning).
#' @param template_var A character string specifying the name of the density
#'   column in \code{template_tab}. Default is "density".
#' @param source_var A character string specifying the name of the fixation
#'   group column in \code{source_tab}. Default is "fixgroup".
#' @param permutations An integer specifying the number of permutation samples
#'   for computing a baseline. Default is 0 (no permutations).
#' @param permute_on An optional character string specifying the column used to
#'   stratify permutations (e.g., "subject" to permute within subjects).
#' @inheritParams sample_density
#' @param aggregate_fun A function used to aggregate density values within time
#'   bins. Default is \code{mean}.
#' @param normalize A character string specifying how to normalize each density
#'   map before sampling. One of \code{"none"} (default), \code{"max"},
#'   \code{"sum"}, or \code{"zscore"}. See \code{\link{sample_density}} for details.
#'   Using \code{"zscore"} is recommended when comparing across conditions or
#'   groups that may differ in fixation count or density bandwidth.
#'
#' @details
#' The function matches each row in \code{source_tab} to a corresponding density

#' map in \code{template_tab} based on the \code{match_on} column. For each matched
#' pair, it samples the template density at the fixation coordinates interpolated
#' at each time point in \code{times}.
#'
#' This approach avoids the problem of having too few fixations to compute density
#' maps within short time windows. Instead, a single density map is created from
#' all fixations (e.g., during encoding), and the sampling is done over time
#' (e.g., during retrieval).
#'
#' If \code{time_bins} is provided, the sampled values are aggregated within each
#' bin using \code{aggregate_fun}. The result includes columns named \code{bin_1},
#' \code{bin_2}, etc.
#'
#' If \code{permutations > 0}, a baseline is computed by sampling from non-matching
#' density maps. As in \code{\link{template_similarity}}, the candidates are the
#' distinct templates matched by other source rows (within the \code{permute_on}
#' stratum, if given); every copy of the true match is excluded, each other
#' template counts once, and up to \code{permutations} of them are drawn without
#' replacement. The result includes \code{perm_sampled} (mean permuted trajectory)
#' and bin-specific permutation columns if binning is used.
#'
#' @return A tibble containing:
#' \itemize{
#'   \item All columns from \code{source_tab}
#'   \item \code{sampled}: A list column with data frames containing \code{z}
#'     (density values) and \code{time} columns for each row
#'   \item If \code{time_bins} is provided: columns \code{bin_1}, \code{bin_2}, etc.
#'     with aggregated values for each bin
#'   \item If \code{permutations > 0}: \code{perm_sampled} (list column with mean
#'     permuted trajectory) and \code{perm_bin_1}, \code{perm_bin_2}, etc.
#' }
#'
#' @examples
#' set.seed(1)
#' # Encoding fixations for three images, and one recall per image that
#' # revisits the encoded locations with some noise.
#' centres <- list(c(150, 150), c(400, 300), c(650, 150))
#' fixations <- function(centre, jitter) {
#'   fixation_group(x = centre[1] + rnorm(6, 0, jitter),
#'                  y = centre[2] + rnorm(6, 0, jitter),
#'                  onset = seq(0, 1000, by = 200), duration = rep(200, 6))
#' }
#' encoding <- tibble::tibble(image = 1:3,
#'                            fixgroup = lapply(centres, fixations, jitter = 40))
#' recall <- tibble::tibble(image = 1:3,
#'                          fixgroup = lapply(centres, fixations, jitter = 60))
#' templates <- density_by(encoding, groups = "image", sigma = 60,
#'                         xbounds = c(0, 800), ybounds = c(0, 450),
#'                         outdim = c(80, 45))
#' res <- sample_density_time(templates, recall, match_on = "image",
#'                            times = seq(0, 1000, by = 100),
#'                            time_bins = c(0, 500, 1000), permutations = 2)
#' res[, c("image", "bin_1", "bin_2", "perm_bin_1", "perm_bin_2")]
#'
#' @md
#' @export
#' @importFrom dplyr bind_cols mutate filter select
#' @importFrom purrr map map_dbl
#' @importFrom tibble tibble
#' @importFrom rlang sym
sample_density_time <- function(template_tab,
                                source_tab,
                                match_on,
                                times = seq(0, 3000, by = 50),
                                time_bins = NULL,
                                template_var = "density",
                                source_var = "fixgroup",
                                permutations = 0,
                                permute_on = NULL,
                                aggregate_fun = mean,
                                normalize = c("none", "max", "sum", "zscore"),
                                transform = c("none", "log"), floor = NULL) {

  normalize <- match.arg(normalize)
  transform <- match.arg(transform)
  if (transform == "log") check_log_floor(floor)


  # Validate inputs

  assertthat::assert_that(match_on %in% names(template_tab),
                          msg = paste("match_on column", match_on, "not found in template_tab"))
  assertthat::assert_that(match_on %in% names(source_tab),
                          msg = paste("match_on column", match_on, "not found in source_tab"))
  assertthat::assert_that(template_var %in% names(template_tab),
                          msg = paste("template_var column", template_var, "not found in template_tab"))
  assertthat::assert_that(source_var %in% names(source_tab),
                          msg = paste("source_var column", source_var, "not found in source_tab"))

  if (!is.null(time_bins)) {
    assertthat::assert_that(length(time_bins) >= 2,
                            msg = "time_bins must have at least 2 values to define bin boundaries")
    assertthat::assert_that(all(diff(time_bins) > 0),
                            msg = "time_bins must be monotonically increasing")
  }

  # Ungroup source_tab to avoid issues with grouped operations

  source_tab <- dplyr::ungroup(source_tab)

  # Match source to template
  matchind <- match_map_rows(source_tab, template_tab, match_on)

  # Handle unmatched rows
  if (any(is.na(matchind))) {
    n_missing <- sum(is.na(matchind))
    warning("Did not find matching template for ", n_missing, " source rows. Removing non-matching elements.")
    keep_rows <- !is.na(matchind)
    source_tab <- source_tab[keep_rows, , drop = FALSE]
    matchind <- matchind[keep_rows]
  }

  # If permutations requested, set up match index splits
  match_split <- NULL
  if (permutations > 0 && !is.null(permute_on)) {
    assertthat::assert_that(permute_on %in% names(source_tab),
                            msg = paste("permute_on column", permute_on, "not found in source_tab"))
    assertthat::assert_that(permute_on %in% names(template_tab),
                            msg = paste("permute_on column", permute_on, "not found in template_tab"))
    match_split <- split(matchind, source_tab[[permute_on]])
  }

  # Helper function to aggregate samples into bins
  aggregate_bins <- function(sampled_df, time_bins, aggregate_fun) {
    if (is.null(sampled_df) || nrow(sampled_df) == 0) {
      return(rep(NA_real_, length(time_bins) - 1))
    }
    # Every bin is half-open, [lower, upper), including the last one, so a time
    # equal to the final break falls outside all bins (label NA).
    bin_labels <- cut(sampled_df$time, breaks = time_bins, right = FALSE,
                      labels = FALSE)
    vapply(seq_len(length(time_bins) - 1), function(b) {
      vals <- sampled_df$z[which(bin_labels == b)]
      if (length(vals) == 0 || all(is.na(vals))) NA_real_ else aggregate_fun(vals, na.rm = TRUE)
    }, numeric(1))
  }

  # Process each row
  template_data <- template_tab[[template_var]]
  source_data <- source_tab[[source_var]]

  # Fast path: templates handled by sample_density.density() that share one
  # lattice are normalized once here, and each source row's sampled grid
  # indices are computed once, so the observed value and every permutation are
  # plain lookups. Any other template (or a row whose indices cannot be
  # computed cleanly) goes through sample_density() exactly as before.
  lookup <- density_lookup_maps(template_data, normalize, which = matchind,
                                enabled = !is.null(times), transform = transform, floor = floor)

  results <- lapply(seq_along(matchind), function(i) {
    template_dens <- template_data[[matchind[i]]]
    source_fix <- source_data[[i]]

    idx <- NULL
    if (!is.null(lookup$x) && !is.null(source_fix)) {
      # Warnings and errors are left to the per-call path to report as before.
      idx <- tryCatch(
        density_grid_index(sample_fixations(source_fix, times), lookup$x, lookup$y),
        warning = function(w) NULL,
        error = function(e) NULL
      )
    }

    # Check for valid inputs
    if (is.null(template_dens) || is.null(source_fix)) {
      sampled <- data.frame(z = rep(NA_real_, length(times)), time = times)
    } else if (!is.null(lookup$z[[matchind[i]]]) && !is.null(idx)) {
      sampled <- data.frame(z = lookup$z[[matchind[i]]][idx], time = times)
    } else {
      sampled <- tryCatch({
        sample_map_dispatch(template_dens, source_fix, times, normalize, transform, floor)
      }, error = function(e) {
        warning("Error sampling density for row ", i, ": ", e$message)
        data.frame(z = rep(NA_real_, length(times)), time = times)
      })
    }

    result <- list(sampled = sampled)

    # Compute bin aggregates if requested
    if (!is.null(time_bins)) {
      bin_vals <- aggregate_bins(sampled, time_bins, aggregate_fun)
      for (b in seq_along(bin_vals)) {
        result[[paste0("bin_", b)]] <- bin_vals[b]
      }
    }

    # Compute permutation baseline if requested
    if (permutations > 0) {
      # Get candidate indices for permutation
      if (!is.null(permute_on)) {
        perm_key <- as.character(source_tab[[permute_on]][i])
        mind <- match_split[[perm_key]]
      } else {
        mind <- matchind
      }

      # Remove every copy of the current match and keep each other template
      # once, as in template_similarity()
      current_match <- matchind[i]
      mind <- reference_candidates(mind, current_match)

      if (length(mind) == 0) {
        # No candidates for permutation
        perm_sampled <- data.frame(z = rep(NA_real_, length(times)), time = times)
        if (!is.null(time_bins)) {
          for (b in seq_len(length(time_bins) - 1)) {
            result[[paste0("perm_bin_", b)]] <- NA_real_
          }
        }
      } else {
        # Sample permutation candidates
        mind <- draw_candidates(mind, permutations)

        # Compute permuted samples
        perm_samples <- lapply(mind, function(j) {
          zmat <- lookup$z[[j]]
          if (!is.null(zmat) && !is.null(idx)) {
            return(zmat[idx])
          }
          perm_dens <- template_data[[j]]
          if (is.null(perm_dens) || is.null(source_fix)) {
            return(rep(NA_real_, length(times)))
          }
          tryCatch({
            sample_map_dispatch(perm_dens, source_fix, times, normalize, transform, floor)$z
          }, error = function(e) {
            rep(NA_real_, length(times))
          })
        })

        # Average across permutations
        perm_matrix <- do.call(rbind, perm_samples)
        perm_mean <- colMeans(perm_matrix, na.rm = TRUE)
        perm_sampled <- data.frame(z = perm_mean, time = times)

        result$perm_sampled <- perm_sampled

        # Compute permuted bin aggregates
        if (!is.null(time_bins)) {
          perm_bin_vals <- aggregate_bins(perm_sampled, time_bins, aggregate_fun)
          for (b in seq_along(perm_bin_vals)) {
            result[[paste0("perm_bin_", b)]] <- perm_bin_vals[b]
          }
        }
      }

      if (is.null(result$perm_sampled)) {
        result$perm_sampled <- data.frame(z = rep(NA_real_, length(times)), time = times)
      }
    }

    result
  })

  # Assemble output tibble
  out <- source_tab

  # Add sampled column (list column)
  out$sampled <- lapply(results, `[[`, "sampled")

  # Add bin columns if present
  if (!is.null(time_bins)) {
    n_bins <- length(time_bins) - 1
    for (b in seq_len(n_bins)) {
      col_name <- paste0("bin_", b)
      out[[col_name]] <- vapply(results, function(r) r[[col_name]] %||% NA_real_, numeric(1))
    }
  }

  # Add permutation columns if present
  if (permutations > 0) {
    out$perm_sampled <- lapply(results, `[[`, "perm_sampled")

    if (!is.null(time_bins)) {
      n_bins <- length(time_bins) - 1
      for (b in seq_len(n_bins)) {
        col_name <- paste0("perm_bin_", b)
        out[[col_name]] <- vapply(results, function(r) r[[col_name]] %||% NA_real_, numeric(1))
      }

      # Add difference columns (observed - permuted)
      for (b in seq_len(n_bins)) {
        obs_col <- paste0("bin_", b)
        perm_col <- paste0("perm_bin_", b)
        diff_col <- paste0("diff_bin_", b)
        out[[diff_col]] <- out[[obs_col]] - out[[perm_col]]
      }
    }
  }

  tibble::as_tibble(out)
}
