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
                          comps = 10, center = TRUE, scale. = FALSE, shrink = 1e-3, ...) {
  proj <- latent_pca_projection(ref_tab, source_tab, refvar, sourcevar, comps = comps, center = center, scale. = scale.)
  k_use <- min(proj$k, comps, nrow(proj$ref_scores) - 1L, nrow(proj$source_scores) - 1L)

  if (k_use < 1) {
    warning("Not enough observations for CCA; returning PCA-transformed scores.")
    ref_tab[[refvar]] <- split_rows(proj$ref_scores)
    source_tab[[sourcevar]] <- split_rows(proj$source_scores)
    return(list(ref_tab = ref_tab, source_tab = source_tab, refvar = refvar, sourcevar = sourcevar,
                info = list(transform = "cca", comps = 0, note = "insufficient observations")))
  }

  scale_with_shrink <- function(mat) {
    vars <- apply(mat, 2, stats::var)
    sds <- sqrt(pmax(vars, 0) + shrink)
    scale(mat[, seq_len(k_use), drop = FALSE], center = TRUE, scale = sds)
  }

  X <- scale_with_shrink(proj$ref_scores)
  Y <- scale_with_shrink(proj$source_scores)

  cca_fit <- tryCatch(stats::cancor(X, Y), error = function(e) NULL)

  if (is.null(cca_fit)) {
    warning("CCA failed to converge; returning PCA-transformed scores.")
    ref_tab[[refvar]] <- split_rows(proj$ref_scores)
    source_tab[[sourcevar]] <- split_rows(proj$source_scores)
    return(list(ref_tab = ref_tab, source_tab = source_tab, refvar = refvar, sourcevar = sourcevar,
                info = list(transform = "cca", comps = 0, note = "cancor failed")))
  }

  k_use <- min(k_use, length(cca_fit$cor))
  ref_latent <- X %*% cca_fit$xcoef[, seq_len(k_use), drop = FALSE]
  src_latent <- Y %*% cca_fit$ycoef[, seq_len(k_use), drop = FALSE]

  ref_tab[[refvar]] <- split_rows(ref_latent)
  source_tab[[sourcevar]] <- split_rows(src_latent)

  list(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = refvar,
    sourcevar = sourcevar,
    info = list(transform = "cca", comps = k_use, center = center, scale = scale., shrink = shrink,
                cor = cca_fit$cor[seq_len(k_use)])
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

latent_pca_projection <- function(ref_tab, source_tab, refvar, sourcevar, comps, center = TRUE, scale. = FALSE) {
  ref_vecs <- lapply(ref_tab[[refvar]], vectorize_density)
  source_vecs <- lapply(source_tab[[sourcevar]], vectorize_density)

   ref_len <- vapply(ref_vecs, length, integer(1))
   src_len <- vapply(source_vecs, length, integer(1))
   if (length(unique(c(ref_len, src_len))) != 1) {
     stop("All density vectors must share the same length for latent transforms.")
   }

  ref_mat <- do.call(rbind, ref_vecs)
  src_mat <- do.call(rbind, source_vecs)

  combined <- rbind(ref_mat, src_mat)
  pc <- stats::prcomp(combined, center = center, scale. = scale.)

  k <- min(comps, ncol(pc$x))
  scores <- pc$x[, seq_len(k), drop = FALSE]

  list(
    ref_scores = scores[seq_len(nrow(ref_mat)), , drop = FALSE],
    source_scores = scores[nrow(ref_mat) + seq_len(nrow(src_mat)), , drop = FALSE],
    k = k,
    basis = pc
  )
}
