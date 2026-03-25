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
#' @param fit_by Optional character vector of column names used to fit CORAL,
#'   CCA, or geometric transforms separately within strata shared by `ref_tab`
#'   and `source_tab`.
#' @param unique_match_only Logical; if `TRUE`, supervised fitting uses at most
#'   one matched pair per `match_on` key within each fitting stratum.
#' @param ... Additional arguments (currently unused; reserved for extension).
#' 
#' @return A list with `ref_tab`, `source_tab`, `refvar`, `sourcevar`, and `info`.
#' @export
latent_pca_transform <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                                 comps = 30, center = TRUE, scale. = FALSE, ...) {
  model <- fit_similarity_transform_model(
    latent_pca_transform,
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = match_on,
    refvar = refvar,
    sourcevar = sourcevar,
    comps = comps,
    center = center,
    scale. = scale.,
    ...
  )
  apply_similarity_transform_model(model, ref_tab = ref_tab, source_tab = source_tab)
}

#' @rdname latent_pca_transform
#' @export
contract_transform <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                               shrink = 1e-6, fit_by = NULL, unique_match_only = FALSE, ...) {
  model <- fit_similarity_transform_model(
    contract_transform,
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = match_on,
    refvar = refvar,
    sourcevar = sourcevar,
    shrink = shrink,
    fit_by = fit_by,
    unique_match_only = unique_match_only,
    ...
  )
  apply_similarity_transform_model(model, ref_tab = ref_tab, source_tab = source_tab)
}

#' @rdname latent_pca_transform
#' @export
affine_transform <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                             shrink = 1e-6, fit_by = NULL, unique_match_only = FALSE, ...) {
  model <- fit_similarity_transform_model(
    affine_transform,
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = match_on,
    refvar = refvar,
    sourcevar = sourcevar,
    shrink = shrink,
    fit_by = fit_by,
    unique_match_only = unique_match_only,
    ...
  )
  apply_similarity_transform_model(model, ref_tab = ref_tab, source_tab = source_tab)
}

#' @rdname latent_pca_transform
#' @export
coral_transform <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                            comps = 30, center = TRUE, scale. = FALSE, shrink = 1e-3,
                            fit_by = NULL, ...) {
  model <- fit_similarity_transform_model(
    coral_transform,
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = match_on,
    refvar = refvar,
    sourcevar = sourcevar,
    comps = comps,
    center = center,
    scale. = scale.,
    shrink = shrink,
    fit_by = fit_by,
    ...
  )
  apply_similarity_transform_model(model, ref_tab = ref_tab, source_tab = source_tab)
}

#' @rdname latent_pca_transform
#' @export
cca_transform <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                          comps = 10, center = TRUE, scale. = FALSE, shrink = 1e-3,
                          fit_by = NULL, unique_match_only = FALSE, ...) {
  model <- fit_similarity_transform_model(
    cca_transform,
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = match_on,
    refvar = refvar,
    sourcevar = sourcevar,
    comps = comps,
    center = center,
    scale. = scale.,
    shrink = shrink,
    fit_by = fit_by,
    unique_match_only = unique_match_only,
    ...
  )
  apply_similarity_transform_model(model, ref_tab = ref_tab, source_tab = source_tab)
}

# Internal helpers -------------------------------------------------------

fit_similarity_transform_model <- function(similarity_transform, ref_tab, source_tab, match_on,
                                           refvar = "density", sourcevar = "density", ...) {
  transform_fun <- resolve_similarity_transform_function(similarity_transform)
  args <- list(...)

  if (identical(transform_fun, latent_pca_transform)) {
    return(do.call(fit_latent_pca_model, c(
      list(ref_tab = ref_tab, source_tab = source_tab, refvar = refvar, sourcevar = sourcevar),
      args
    )))
  }

  if (identical(transform_fun, contract_transform)) {
    return(do.call(fit_contract_model, c(
      list(ref_tab = ref_tab, source_tab = source_tab, match_on = match_on, refvar = refvar, sourcevar = sourcevar),
      args
    )))
  }

  if (identical(transform_fun, affine_transform)) {
    return(do.call(fit_affine_model, c(
      list(ref_tab = ref_tab, source_tab = source_tab, match_on = match_on, refvar = refvar, sourcevar = sourcevar),
      args
    )))
  }

  if (identical(transform_fun, coral_transform)) {
    return(do.call(fit_coral_model, c(
      list(ref_tab = ref_tab, source_tab = source_tab, refvar = refvar, sourcevar = sourcevar),
      args
    )))
  }

  if (identical(transform_fun, cca_transform)) {
    return(do.call(fit_cca_model, c(
      list(ref_tab = ref_tab, source_tab = source_tab, match_on = match_on, refvar = refvar, sourcevar = sourcevar),
      args
    )))
  }

  stop("Out-of-sample transform fitting currently supports latent_pca_transform, contract_transform, affine_transform, coral_transform, and cca_transform only.")
}

apply_similarity_transform_model <- function(model, ref_tab, source_tab) {
  if (is.null(model$transform)) {
    stop("Invalid transform model: missing transform name.")
  }

  if (identical(model$transform, "latent_pca")) {
    ref_scores <- project_transform_scores(model, ref_tab, model$refvar)
    src_scores <- project_transform_scores(model, source_tab, model$sourcevar)
    ref_tab[[model$refvar]] <- split_rows(ref_scores)
    source_tab[[model$sourcevar]] <- split_rows(src_scores)
  } else if (identical(model$transform, "contract") || identical(model$transform, "affine")) {
    transformed <- apply_geometric_transform_model(model, ref_tab, source_tab)
    ref_tab <- transformed$ref_tab
    source_tab <- transformed$source_tab
  } else if (identical(model$transform, "coral")) {
    ref_scores <- project_transform_scores(model$pca_model, ref_tab, model$refvar)
    src_scores <- project_transform_scores(model$pca_model, source_tab, model$sourcevar)

    adapted_src <- src_scores
    if (nrow(src_scores) > 0L) {
      group_keys <- transform_group_keys(ref_tab, source_tab, model$fit_by)
      group_levels <- unique(group_keys$source)

      for (group_key in group_levels) {
        group_label <- if (identical(group_key, "__all__")) "all" else group_key
        group_model <- model$group_models[[group_label]]
        if (is.null(group_model)) {
          next
        }

        src_rows <- which(group_keys$source == group_key)
        adapted_src[src_rows, seq_len(group_model$k)] <- t(group_model$adapt %*% t(src_scores[src_rows, seq_len(group_model$k), drop = FALSE]))
      }
    }

    ref_tab[[model$refvar]] <- split_rows(ref_scores)
    source_tab[[model$sourcevar]] <- split_rows(adapted_src)
  } else if (identical(model$transform, "cca")) {
    transformed <- apply_cca_model(model, ref_tab, source_tab)
    ref_tab <- transformed$ref_tab
    source_tab <- transformed$source_tab
  } else {
    stop("Unsupported transform model: ", model$transform)
  }

  list(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = model$refvar,
    sourcevar = model$sourcevar,
    info = model$info
  )
}

fit_contract_model <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                               shrink = 1e-6, fit_by = NULL, unique_match_only = FALSE, ...) {
  fit_geometric_density_model(
    transform = "contract",
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = match_on,
    refvar = refvar,
    sourcevar = sourcevar,
    shrink = shrink,
    fit_by = fit_by,
    unique_match_only = unique_match_only
  )
}

fit_affine_model <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                             shrink = 1e-6, fit_by = NULL, unique_match_only = FALSE, ...) {
  fit_geometric_density_model(
    transform = "affine",
    ref_tab = ref_tab,
    source_tab = source_tab,
    match_on = match_on,
    refvar = refvar,
    sourcevar = sourcevar,
    shrink = shrink,
    fit_by = fit_by,
    unique_match_only = unique_match_only
  )
}

fit_geometric_density_model <- function(transform, ref_tab, source_tab, match_on, refvar, sourcevar,
                                        shrink = 1e-6, fit_by = NULL, unique_match_only = FALSE) {
  group_keys <- transform_group_keys(ref_tab, source_tab, fit_by)
  group_levels <- union(unique(group_keys$ref), unique(group_keys$source))
  group_models <- list()
  group_info <- vector("list", length(group_levels))

  for (g in seq_along(group_levels)) {
    group_key <- group_levels[[g]]
    ref_rows <- which(group_keys$ref == group_key)
    src_rows <- which(group_keys$source == group_key)
    group_label <- if (identical(group_key, "__all__")) "all" else group_key

    if (length(ref_rows) == 0L || length(src_rows) == 0L) {
      group_info[[g]] <- list(group = group_label, matched_pairs = 0L, note = "missing group rows")
      next
    }

    paired_idx <- matched_transform_indices(
      ref_tab[ref_rows, , drop = FALSE],
      source_tab[src_rows, , drop = FALSE],
      match_on = match_on,
      unique_match_only = unique_match_only
    )
    n_pairs <- length(paired_idx$source)

    if (n_pairs < 1L) {
      group_info[[g]] <- list(group = group_label, matched_pairs = 0L, note = "no matched observations")
      next
    }

    ref_objs <- ref_tab[[refvar]][ref_rows[paired_idx$ref]]
    src_objs <- source_tab[[sourcevar]][src_rows[paired_idx$source]]
    ref_mom <- aggregate_density_moments(ref_objs)
    src_mom <- aggregate_density_moments(src_objs)

    if (identical(transform, "contract")) {
      src_radius <- sum(diag(src_mom$cov))
      ref_radius <- sum(diag(ref_mom$cov))
      scale <- sqrt((ref_radius + shrink) / (src_radius + shrink))
      A <- diag(scale, 2L)
      extra <- list(scale = scale)
    } else {
      A <- mat_sqrt(src_mom$cov * 0 + ref_mom$cov + diag(shrink, 2L)) %*%
        mat_inv_sqrt(src_mom$cov + diag(shrink, 2L), shrink)
      extra <- list()
    }

    t_vec <- as.numeric(ref_mom$mean - A %*% src_mom$mean)
    group_models[[group_label]] <- c(
      list(A = A, t = t_vec, matched_pairs = n_pairs),
      extra
    )
    group_info[[g]] <- c(
      list(group = group_label, matched_pairs = n_pairs, note = "ok"),
      extra
    )
  }

  structure(
    list(
      transform = transform,
      group_models = group_models,
      fit_by = fit_by,
      refvar = refvar,
      sourcevar = sourcevar,
      info = list(
        transform = transform,
        shrink = shrink,
        fit_by = fit_by,
        unique_match_only = unique_match_only,
        groups = group_info
      )
    ),
    class = c("similarity_transform_model", "list")
  )
}

apply_geometric_transform_model <- function(model, ref_tab, source_tab) {
  group_keys <- transform_group_keys(ref_tab, source_tab, model$fit_by)
  source_tab[[model$sourcevar]] <- lapply(seq_len(nrow(source_tab)), function(i) {
    group_key <- group_keys$source[[i]]
    group_label <- if (identical(group_key, "__all__")) "all" else group_key
    group_model <- model$group_models[[group_label]]
    if (is.null(group_model)) {
      return(source_tab[[model$sourcevar]][[i]])
    }

    warp_density_object(source_tab[[model$sourcevar]][[i]], A = group_model$A, t = group_model$t)
  })

  list(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = model$refvar,
    sourcevar = model$sourcevar,
    info = model$info
  )
}

resolve_similarity_transform_function <- function(similarity_transform) {
  if (is.function(similarity_transform)) {
    return(similarity_transform)
  }

  if (is.character(similarity_transform) && length(similarity_transform) == 1L) {
    return(get(similarity_transform, mode = "function", inherits = TRUE))
  }

  stop("similarity_transform must be a function or a single function name.")
}

fit_latent_pca_model <- function(ref_tab, source_tab, refvar = "density", sourcevar = "density",
                                 comps = 30, center = TRUE, scale. = FALSE, ...) {
  ref_mat <- vectorize_density_tab(ref_tab, refvar)
  src_mat <- vectorize_density_tab(source_tab, sourcevar, expected_len = ncol(ref_mat))

  combined <- rbind(ref_mat, src_mat)
  if (nrow(combined) == 0L || ncol(combined) == 0L) {
    stop("Cannot fit latent transform on empty density matrices.")
  }

  rank_target <- min(comps, nrow(combined), ncol(combined))
  pc_args <- list(x = combined, center = center, scale. = scale.)
  prcomp_default <- utils::getS3method("prcomp", "default", optional = TRUE)
  if (!is.null(prcomp_default) && "rank." %in% names(formals(prcomp_default))) {
    pc_args$rank. <- rank_target
  }
  pc <- do.call(stats::prcomp, pc_args)
  k <- min(comps, ncol(pc$rotation))

  structure(
    list(
      transform = "latent_pca",
      basis = pc,
      k = k,
      refvar = refvar,
      sourcevar = sourcevar,
      input_dim = ncol(ref_mat),
      info = list(transform = "latent_pca", comps = k, center = center, scale = scale.)
    ),
    class = c("similarity_transform_model", "list")
  )
}

fit_coral_model <- function(ref_tab, source_tab, refvar = "density", sourcevar = "density",
                            comps = 30, center = TRUE, scale. = FALSE, shrink = 1e-3,
                            fit_by = NULL, ...) {
  pca_model <- fit_latent_pca_model(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = refvar,
    sourcevar = sourcevar,
    comps = comps,
    center = center,
    scale. = scale.
  )

  ref_scores <- project_transform_scores(pca_model, ref_tab, refvar)
  src_scores <- project_transform_scores(pca_model, source_tab, sourcevar)
  group_keys <- transform_group_keys(ref_tab, source_tab, fit_by)
  group_levels <- union(unique(group_keys$ref), unique(group_keys$source))
  group_models <- list()
  group_info <- vector("list", length(group_levels))

  for (g in seq_along(group_levels)) {
    group_key <- group_levels[[g]]
    ref_rows <- which(group_keys$ref == group_key)
    src_rows <- which(group_keys$source == group_key)
    group_label <- if (identical(group_key, "__all__")) "all" else group_key

    if (length(ref_rows) == 0L || length(src_rows) == 0L) {
      group_info[[g]] <- list(group = group_label, ref_n = length(ref_rows), source_n = length(src_rows), note = "missing group rows")
      next
    }

    cov_ref <- cov_with_shrink(ref_scores[ref_rows, , drop = FALSE], shrink)
    cov_src <- cov_with_shrink(src_scores[src_rows, , drop = FALSE], shrink)
    group_models[[group_label]] <- list(
      k = pca_model$k,
      adapt = mat_inv_sqrt(cov_src, shrink) %*% mat_sqrt(cov_ref)
    )
    group_info[[g]] <- list(
      group = group_label,
      ref_n = length(ref_rows),
      source_n = length(src_rows),
      note = "ok"
    )
  }

  structure(
    list(
      transform = "coral",
      pca_model = pca_model,
      group_models = group_models,
      k = pca_model$k,
      fit_by = fit_by,
      refvar = refvar,
      sourcevar = sourcevar,
      info = list(
        transform = "coral",
        comps = pca_model$k,
        center = center,
        scale = scale.,
        shrink = shrink,
        fit_by = fit_by,
        groups = group_info
      )
    ),
    class = c("similarity_transform_model", "list")
  )
}

fit_cca_model <- function(ref_tab, source_tab, match_on, refvar = "density", sourcevar = "density",
                          comps = 10, center = TRUE, scale. = FALSE, shrink = 1e-3,
                          fit_by = NULL, unique_match_only = FALSE, ...) {
  pca_model <- fit_latent_pca_model(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = refvar,
    sourcevar = sourcevar,
    comps = comps,
    center = center,
    scale. = scale.
  )
  k_target <- pca_model$k

  if (k_target < 1L) {
    return(structure(
      list(
        transform = "cca",
        pca_model = pca_model,
        group_models = list(),
        k = 0L,
        fit_by = fit_by,
        refvar = refvar,
        sourcevar = sourcevar,
        info = list(transform = "cca", comps = 0, note = "insufficient components")
      ),
      class = c("similarity_transform_model", "list")
    ))
  }

  ref_scores <- project_transform_scores(pca_model, ref_tab, refvar)
  src_scores <- project_transform_scores(pca_model, source_tab, sourcevar)
  group_keys <- transform_group_keys(ref_tab, source_tab, fit_by)
  group_levels <- union(unique(group_keys$ref), unique(group_keys$source))
  group_models <- list()
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
    ref_train <- ref_scores[ref_rows[paired_idx$ref], seq_len(k_use), drop = FALSE]
    src_train <- src_scores[src_rows[paired_idx$source], seq_len(k_use), drop = FALSE]
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

    group_models[[group_label]] <- list(
      k = k_group,
      ref_fit = ref_fit,
      src_fit = src_fit,
      xcoef = cca_fit$xcoef[, seq_len(k_group), drop = FALSE],
      ycoef = cca_fit$ycoef[, seq_len(k_group), drop = FALSE]
    )
    group_info[[g]] <- list(
      group = group_label,
      matched_pairs = n_pairs,
      comps = k_group,
      cor = cca_fit$cor[seq_len(k_group)],
      note = "ok"
    )
  }

  structure(
    list(
      transform = "cca",
      pca_model = pca_model,
      group_models = group_models,
      k = k_target,
      fit_by = fit_by,
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
    ),
    class = c("similarity_transform_model", "list")
  )
}

apply_cca_model <- function(model, ref_tab, source_tab) {
  ref_scores <- project_transform_scores(model$pca_model, ref_tab, model$refvar)
  src_scores <- project_transform_scores(model$pca_model, source_tab, model$sourcevar)
  ref_latent <- ref_scores
  src_latent <- src_scores
  group_keys <- transform_group_keys(ref_tab, source_tab, model$fit_by)
  group_levels <- union(unique(group_keys$ref), unique(group_keys$source))

  for (group_key in group_levels) {
    group_label <- if (identical(group_key, "__all__")) "all" else group_key
    group_model <- model$group_models[[group_label]]
    if (is.null(group_model)) {
      next
    }

    ref_rows <- which(group_keys$ref == group_key)
    src_rows <- which(group_keys$source == group_key)

    if (length(ref_rows) > 0L) {
      ref_all <- apply_scale_with_shrink(ref_scores[ref_rows, seq_len(group_model$k), drop = FALSE], group_model$ref_fit)
      ref_latent[ref_rows, seq_len(group_model$k)] <- ref_all %*% group_model$xcoef
    }

    if (length(src_rows) > 0L) {
      src_all <- apply_scale_with_shrink(src_scores[src_rows, seq_len(group_model$k), drop = FALSE], group_model$src_fit)
      src_latent[src_rows, seq_len(group_model$k)] <- src_all %*% group_model$ycoef
    }
  }

  ref_tab[[model$refvar]] <- split_rows(ref_latent)
  source_tab[[model$sourcevar]] <- split_rows(src_latent)

  list(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = model$refvar,
    sourcevar = model$sourcevar,
    info = model$info
  )
}

project_transform_scores <- function(model, tab, var) {
  if (nrow(tab) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = model$k))
  }

  mat <- vectorize_density_tab(tab, var, expected_len = model$input_dim)
  apply_prcomp_projection(mat, model$basis, model$k)
}

apply_prcomp_projection <- function(mat, basis, k) {
  scaled <- scale(
    mat,
    center = if (is.null(basis$center)) FALSE else basis$center,
    scale = if (is.null(basis$scale)) FALSE else basis$scale
  )
  scaled %*% basis$rotation[, seq_len(k), drop = FALSE]
}

vectorize_density_tab <- function(tab, var, expected_len = NULL) {
  vecs <- lapply(tab[[var]], vectorize_density)
  if (length(vecs) == 0L) {
    ncol <- if (is.null(expected_len)) 0L else expected_len
    return(matrix(numeric(0), nrow = 0L, ncol = ncol))
  }

  lens <- vapply(vecs, length, integer(1))
  if (length(unique(lens)) != 1L) {
    stop("All density vectors must share the same length for latent transforms.")
  }
  if (!is.null(expected_len) && unique(lens) != expected_len) {
    stop("Density vectors must match the transform training dimensionality.")
  }

  stack_rows(vecs)
}

mat_sqrt <- function(m) {
  ev <- eigen(m)
  ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*% t(ev$vectors)
}

mat_inv_sqrt <- function(m, shrink) {
  ev <- eigen(m)
  ev$vectors %*% diag(1 / sqrt(pmax(ev$values, shrink))) %*% t(ev$vectors)
}

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
    stop("fit_by columns must exist in both ref_tab and source_tab for grouped latent transforms.")
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

cov_with_shrink <- function(mat, shrink) {
  k <- ncol(mat)
  if (k == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }

  centered <- sweep(mat, 2, colMeans(mat), FUN = "-")
  cov_mat <- if (nrow(centered) > 1L) {
    crossprod(centered) / (nrow(centered) - 1L)
  } else {
    matrix(0, nrow = k, ncol = k)
  }

  cov_mat + diag(shrink, k)
}

aggregate_density_moments <- function(objs) {
  moments <- lapply(objs, density_object_moments)
  means <- do.call(rbind, lapply(moments, `[[`, "mean"))
  grand_mean <- colMeans(means)
  grand_cov <- Reduce(`+`, lapply(moments, function(m) {
    centered_mean <- m$mean - grand_mean
    m$cov + tcrossprod(centered_mean)
  })) / length(moments)

  list(mean = grand_mean, cov = grand_cov)
}

density_object_moments <- function(obj) {
  dens <- primary_density_scale(obj)
  if (!inherits(dens, c("density", "eye_density"))) {
    stop("Geometric transforms require eye_density or density objects.")
  }

  coords <- as.matrix(expand.grid(x = dens$x, y = dens$y))
  weights <- as.numeric(dens$z)
  sum_w <- sum(weights)

  if (!is.finite(sum_w) || sum_w <= .Machine$double.eps) {
    mean_vec <- c(mean(dens$x), mean(dens$y))
    cov_mat <- matrix(0, nrow = 2L, ncol = 2L)
    return(list(mean = mean_vec, cov = cov_mat))
  }

  weights <- weights / sum_w
  mean_vec <- colSums(coords * weights)
  centered <- sweep(coords, 2, mean_vec, FUN = "-")
  cov_mat <- t(centered * weights) %*% centered

  list(mean = mean_vec, cov = cov_mat)
}

primary_density_scale <- function(obj) {
  if (inherits(obj, "eye_density_multiscale")) {
    if (length(obj) == 0L) {
      stop("Geometric transforms cannot use empty eye_density_multiscale objects.")
    }
    return(obj[[1]])
  }
  obj
}

warp_density_object <- function(obj, A, t) {
  if (inherits(obj, "eye_density_multiscale")) {
    warped <- lapply(obj, function(el) warp_density_object(el, A = A, t = t))
    attr(warped, "sigmas_vector") <- attr(obj, "sigmas_vector")
    class(warped) <- class(obj)
    return(warped)
  }

  dens <- primary_density_scale(obj)
  if (!inherits(dens, c("density", "eye_density"))) {
    stop("Geometric transforms require eye_density or density objects.")
  }

  inv_A <- solve(A)
  grid_xy <- expand.grid(x = dens$x, y = dens$y)
  source_xy <- t(inv_A %*% t(sweep(as.matrix(grid_xy), 2, t, FUN = "-")))
  z_new <- bilinear_sample_density(
    x_grid = dens$x,
    y_grid = dens$y,
    z = dens$z,
    xq = source_xy[, 1],
    yq = source_xy[, 2]
  )
  z_new <- matrix(z_new, nrow = length(dens$x), ncol = length(dens$y))
  sum_z <- sum(z_new)
  if (is.finite(sum_z) && sum_z > .Machine$double.eps) {
    z_new <- z_new / sum_z
  }

  dens$z <- zapsmall(z_new)
  dens
}

bilinear_sample_density <- function(x_grid, y_grid, z, xq, yq) {
  nx <- length(x_grid)
  ny <- length(y_grid)
  out <- rep(0, length(xq))

  valid <- xq >= min(x_grid) & xq <= max(x_grid) & yq >= min(y_grid) & yq <= max(y_grid)
  if (!any(valid)) {
    return(out)
  }

  xqv <- xq[valid]
  yqv <- yq[valid]
  ix <- pmin(findInterval(xqv, x_grid, left.open = FALSE), nx - 1L)
  iy <- pmin(findInterval(yqv, y_grid, left.open = FALSE), ny - 1L)
  ix[ix < 1L] <- 1L
  iy[iy < 1L] <- 1L

  x1 <- x_grid[ix]
  x2 <- x_grid[ix + 1L]
  y1 <- y_grid[iy]
  y2 <- y_grid[iy + 1L]
  wx <- ifelse(x2 > x1, (xqv - x1) / (x2 - x1), 0)
  wy <- ifelse(y2 > y1, (yqv - y1) / (y2 - y1), 0)

  z11 <- z[cbind(ix, iy)]
  z21 <- z[cbind(ix + 1L, iy)]
  z12 <- z[cbind(ix, iy + 1L)]
  z22 <- z[cbind(ix + 1L, iy + 1L)]

  out[valid] <- (1 - wx) * (1 - wy) * z11 +
    wx * (1 - wy) * z21 +
    (1 - wx) * wy * z12 +
    wx * wy * z22

  out
}

apply_scale_with_shrink <- function(mat, fit) {
  centered <- sweep(mat, 2, fit$center, FUN = "-")
  sweep(centered, 2, fit$scale, FUN = "/")
}

latent_pca_projection <- function(ref_tab, source_tab, refvar, sourcevar, comps, center = TRUE, scale. = FALSE) {
  model <- fit_latent_pca_model(
    ref_tab = ref_tab,
    source_tab = source_tab,
    refvar = refvar,
    sourcevar = sourcevar,
    comps = comps,
    center = center,
    scale. = scale.
  )
  ref_scores <- project_transform_scores(model, ref_tab, refvar)
  source_scores <- project_transform_scores(model, source_tab, sourcevar)

  list(ref_scores = ref_scores, source_scores = source_scores, k = model$k, basis = model$basis)
}
