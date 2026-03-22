#' Latent-space transforms for template-based similarity
#'
#' These helpers convert densities to low-dimensional vectors before similarity
#' is computed. They are intended to be passed via the `similarity_transform`
#' argument of `template_similarity`.
#'
#' @param ref_tab,source_tab Data frames passed from `template_similarity`.
#' @param match_on Column used for matching rows between tables.
#' @param refvar,sourcevar Density column names to transform.
#' @param comps Number of latent components to retain (post-PCA).
#' @param center,scale. Logical flags passed to `stats::prcomp`.
#' @param shrink Ridge term for covariance regularization.
#' @param fit_by Optional character vector of column names used to fit CCA
#'   separately within strata shared by `ref_tab` and `source_tab`.
#' @param unique_match_only Logical; if `TRUE`, CCA fitting uses at most one
#'   matched pair per `match_on` key within each fitting stratum.
#' @param ... Additional arguments (currently unused; reserved for extension).
#' 
#' @return A list with `ref_tab`, `source_tab`, `refvar`, `sourcevar`, and `info`.
#' @export
latent_pca_transform <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                                 comps = 30, center = TRUE, scale. = FALSE, ...) {
  proj <- latent_pca_projection(ref_tab, source_tab, refvar, sourcevar, comps = comps, center = center, scale. = scale.)

  ref_tab[[refvar]] <- split_rows(proj$ref_scores)
  source_tab[[sourcevar]] <- split_rows(proj$source_scores)

  list(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = refvar,
    sourcevar = sourcevar,
    info = list(transform = "latent_pca", comps = proj$k, center = center, scale = scale.)
  )
}

#' @rdname latent_pca_transform
#' @export
coral_transform <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                            comps = 30, center = TRUE, scale. = FALSE, shrink = 1e-3, ...) {
  proj <- latent_pca_projection(ref_tab, source_tab, refvar, sourcevar, comps = comps, center = center, scale. = scale.)

  cov_ref <- stats::cov(proj$ref_scores) + diag(shrink, proj$k)
  cov_src <- stats::cov(proj$source_scores) + diag(shrink, proj$k)

  mat_sqrt <- function(m) {
    ev <- eigen(m)
    ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*% t(ev$vectors)
  }
  mat_inv_sqrt <- function(m) {
    ev <- eigen(m)
    ev$vectors %*% diag(1 / sqrt(pmax(ev$values, shrink))) %*% t(ev$vectors)
  }

  adapt <- mat_inv_sqrt(cov_src) %*% mat_sqrt(cov_ref)
  adapted_src <- t(adapt %*% t(proj$source_scores))

  ref_tab[[refvar]] <- split_rows(proj$ref_scores)
  source_tab[[sourcevar]] <- split_rows(adapted_src)

  list(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = refvar,
    sourcevar = sourcevar,
    info = list(transform = "coral", comps = proj$k, center = center, scale = scale., shrink = shrink)
  )
}

#' @rdname latent_pca_transform
#' @export
cca_transform <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                          comps = 10, center = TRUE, scale. = FALSE, shrink = 1e-3,
                          fit_by = NULL, unique_match_only = FALSE, ...) {
  proj <- latent_pca_projection(ref_tab, source_tab, refvar, sourcevar, comps = comps, center = center, scale. = scale.)
  k_target <- min(proj$k, comps)

  if (k_target < 1L) {
    warning("Not enough components for CCA; returning PCA-transformed scores.")
    ref_tab[[refvar]] <- split_rows(proj$ref_scores)
    source_tab[[sourcevar]] <- split_rows(proj$source_scores)
    return(list(
      ref_tab = ref_tab,
      source_tab = source_tab,
      refvar = refvar,
      sourcevar = sourcevar,
      info = list(transform = "cca", comps = 0, note = "insufficient components")
    ))
  }

  ref_latent <- proj$ref_scores[, seq_len(k_target), drop = FALSE]
  src_latent <- proj$source_scores[, seq_len(k_target), drop = FALSE]
  group_keys <- transform_group_keys(ref_tab, source_tab, fit_by)
  group_levels <- union(unique(group_keys$ref), unique(group_keys$source))
  group_info <- vector("list", length(group_levels))
  total_pairs <- 0L

  for (g in seq_along(group_levels)) {
    group_key <- group_levels[[g]]
    ref_rows <- which(group_keys$ref == group_key)
    src_rows <- which(group_keys$source == group_key)
    group_label <- if (identical(group_key, "__all__")) "all" else group_key

    if (length(ref_rows) == 0L || length(src_rows) == 0L) {
      group_info[[g]] <- list(group = group_label, matched_pairs = 0L, comps = 0L, note = "missing group rows")
      next
    }

    paired_idx <- matched_transform_indices(
      ref_tab[ref_rows, , drop = FALSE],
      source_tab[src_rows, , drop = FALSE],
      match_on = match_on,
      unique_match_only = unique_match_only
    )
    n_pairs <- length(paired_idx$source)
    total_pairs <- total_pairs + n_pairs

    if (n_pairs < 2L) {
      group_info[[g]] <- list(group = group_label, matched_pairs = n_pairs, comps = 0L, note = "insufficient matched observations")
      next
    }

    k_use <- min(k_target, n_pairs - 1L)
    ref_train <- proj$ref_scores[ref_rows[paired_idx$ref], seq_len(k_use), drop = FALSE]
    src_train <- proj$source_scores[src_rows[paired_idx$source], seq_len(k_use), drop = FALSE]

    ref_fit <- fit_scale_with_shrink(ref_train, shrink)
    src_fit <- fit_scale_with_shrink(src_train, shrink)
    cca_fit <- tryCatch(stats::cancor(ref_fit$scores, src_fit$scores), error = function(e) NULL)

    if (is.null(cca_fit)) {
      group_info[[g]] <- list(group = group_label, matched_pairs = n_pairs, comps = 0L, note = "cancor failed")
      next
    }

    k_group <- min(k_use, length(cca_fit$cor))
    if (k_group < 1L) {
      group_info[[g]] <- list(group = group_label, matched_pairs = n_pairs, comps = 0L, note = "no canonical components")
      next
    }

    ref_all <- apply_scale_with_shrink(proj$ref_scores[ref_rows, seq_len(k_group), drop = FALSE], ref_fit)
    src_all <- apply_scale_with_shrink(proj$source_scores[src_rows, seq_len(k_group), drop = FALSE], src_fit)

    ref_latent[ref_rows, seq_len(k_group)] <- ref_all %*% cca_fit$xcoef[, seq_len(k_group), drop = FALSE]
    src_latent[src_rows, seq_len(k_group)] <- src_all %*% cca_fit$ycoef[, seq_len(k_group), drop = FALSE]

    group_info[[g]] <- list(
      group = group_label,
      matched_pairs = n_pairs,
      comps = k_group,
      cor = cca_fit$cor[seq_len(k_group)],
      note = "ok"
    )
  }

  ref_tab[[refvar]] <- split_rows(ref_latent)
  source_tab[[sourcevar]] <- split_rows(src_latent)

  list(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = refvar,
    sourcevar = sourcevar,
    info = list(
      transform = "cca",
      comps = k_target,
      center = center,
      scale = scale.,
      shrink = shrink,
      fit_by = fit_by,
      unique_match_only = unique_match_only,
      matched_pairs = total_pairs,
      groups = group_info
    )
  )
}

# Internal helpers -------------------------------------------------------

vectorize_density <- function(obj) {
  if (inherits(obj, "eye_density_multiscale")) {
    lens <- vapply(obj, function(el) length(as.vector(el$z)), integer(1))
    if (length(unique(lens)) != 1) {
      stop("All scales in eye_density_multiscale must have the same grid dimensions for latent transforms.")
    }
    obj <- obj[[1]]
  }
  if (inherits(obj, c("density", "eye_density"))) {
    return(as.vector(obj$z))
  }
  if (is.numeric(obj) || is.matrix(obj)) {
    return(as.vector(obj))
  }
  stop("Unsupported object passed to latent transform; expected density or numeric.")
}

split_rows <- function(mat) {
  lapply(seq_len(nrow(mat)), function(i) mat[i, , drop = TRUE])
}

stack_rows <- function(vecs) {
  if (length(vecs) == 0) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }

  matrix(
    unlist(vecs, use.names = FALSE),
    nrow = length(vecs),
    byrow = TRUE
  )
}

matched_transform_indices <- function(ref_tab, source_tab, match_on, unique_match_only = FALSE) {
  if (!match_on %in% names(ref_tab) || !match_on %in% names(source_tab)) {
    stop("match_on column must exist in both ref_tab and source_tab for latent transforms.")
  }

  ref_match <- match(source_tab[[match_on]], ref_tab[[match_on]])
  src_idx <- which(!is.na(ref_match))
  ref_idx <- ref_match[src_idx]

  if (unique_match_only && length(src_idx) > 0L) {
    keep <- !duplicated(source_tab[[match_on]][src_idx])
    src_idx <- src_idx[keep]
    ref_idx <- ref_idx[keep]
  }

  list(
    ref = ref_idx,
    source = src_idx
  )
}

transform_group_keys <- function(ref_tab, source_tab, fit_by = NULL) {
  if (is.null(fit_by)) {
    return(list(
      ref = rep("__all__", nrow(ref_tab)),
      source = rep("__all__", nrow(source_tab))
    ))
  }

  missing_ref <- setdiff(fit_by, names(ref_tab))
  missing_src <- setdiff(fit_by, names(source_tab))
  if (length(missing_ref) > 0L || length(missing_src) > 0L) {
    stop("fit_by columns must exist in both ref_tab and source_tab for CCA transforms.")
  }

  make_key <- function(tab) {
    interaction(tab[fit_by], drop = TRUE, lex.order = TRUE, sep = "::")
  }

  list(
    ref = as.character(make_key(ref_tab)),
    source = as.character(make_key(source_tab))
  )
}

fit_scale_with_shrink <- function(mat, shrink) {
  center_vec <- colMeans(mat)
  centered <- sweep(mat, 2, center_vec, FUN = "-")
  vars <- if (nrow(centered) > 1L) {
    colSums(centered^2) / (nrow(centered) - 1L)
  } else {
    rep(0, ncol(centered))
  }
  scale_vec <- sqrt(pmax(vars, 0) + shrink)
  scale_vec[scale_vec == 0] <- 1

  list(
    scores = sweep(centered, 2, scale_vec, FUN = "/"),
    center = center_vec,
    scale = scale_vec
  )
}

apply_scale_with_shrink <- function(mat, fit) {
  centered <- sweep(mat, 2, fit$center, FUN = "-")
  sweep(centered, 2, fit$scale, FUN = "/")
}

latent_pca_projection <- function(ref_tab, source_tab, refvar, sourcevar, comps, center = TRUE, scale. = FALSE) {
  ref_vecs <- lapply(ref_tab[[refvar]], vectorize_density)
  source_vecs <- lapply(source_tab[[sourcevar]], vectorize_density)

   ref_len <- vapply(ref_vecs, length, integer(1))
   src_len <- vapply(source_vecs, length, integer(1))
   if (length(unique(c(ref_len, src_len))) != 1) {
     stop("All density vectors must share the same length for latent transforms.")
   }

  ref_mat <- stack_rows(ref_vecs)
  src_mat <- stack_rows(source_vecs)

  combined <- rbind(ref_mat, src_mat)
  rank_target <- min(comps, nrow(combined), ncol(combined))
  pc_args <- list(x = combined, center = center, scale. = scale.)
  if ("rank." %in% names(formals(stats:::prcomp.default))) {
    pc_args$rank. <- rank_target
  }
  pc <- do.call(stats::prcomp, pc_args)

  k <- min(comps, ncol(pc$x))
  scores <- pc$x[, seq_len(k), drop = FALSE]

  list(
    ref_scores = scores[seq_len(nrow(ref_mat)), , drop = FALSE],
    source_scores = scores[nrow(ref_mat) + seq_len(nrow(src_mat)), , drop = FALSE],
    k = k,
    basis = pc
  )
}
